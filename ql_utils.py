from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import pathlib
import logging
import os

def_pvar = ('lat', 'lon', 'cldmask', 'illum', 'solar_zenith_view_no1',
            'cth', 'cot', 'aot550', 'aer', 'qcflag', 'lsflag', 'niter', 'costjm')

def_svar = ('reflectance_in_channel_no_1', 'reflectance_in_channel_no_2',
            'reflectance_in_channel_no_3', 'brightness_temperature_in_channel_no_9')


class OutputLims:
    def __init__(self,
                 min_aod=0.,
                 max_aod=5.,
                 min_cth=0.,
                 max_cth=15.,
                 min_cer=0.,
                 max_cer=50.,
                 min_dcth=0,
                 max_dcth=10.,
                 min_cod=0,
                 max_cod=100.):

        # Minimum aerosol optical depth for plotting
        self.min_aod = min_aod
        # Maximum aerosol optical depth for plotting
        self.max_aod = max_aod
        # Minimum cloud top height for plotting
        self.min_cth = min_cth
        # Maximum cloud top height for plotting
        self.max_cth = max_cth
        # Minimum cloud effective radius for plotting
        self.min_cer = min_cer
        # Maximum cloud effective radius for plotting
        self.max_cer = max_cer
        # Minimum cloud height uncertainty for plotting
        self.min_dcth = min_dcth
        # Maximum cloud height uncertainty for plotting
        self.max_dcth = max_dcth
        # Minimum cloud optical depth for plotting
        self.min_cod = min_cod
        # Maximum cloud optical depth for plotting
        self.max_cod = max_cod


class QuickLookOpts:
    def __init__(self,
                 in_dtstr='202109101100',
                 indir='./',
                 outdir_top='./TEST/',
                 coast_dir=None,
                 clobber=False,
                 aerosol_qc=284,
                 dist2cloud=3.,
                 aerosol_landsea=False,
                 cesium=True,
                 flip_data=True,
                 auto_out=True,
                 out_img_pix=(1700, 597),
                 out_img_scl_cs=(358, 192),
                 out_img_ll=(-71.8154, 29.1062, 30.1846, 64.7465),
                 pvar=def_pvar,
                 svar=def_svar,
                 sza_thresh=70.,
                 perc_max=99.,
                 res_meth='bilin'):

        # Directory containing the ORAC pri + sec files
        self.indir = indir
        # Top level directory for the output files
        self.outdir_top = outdir_top

        # Directory containing coastline shapefiles
        self.coast_dir = coast_dir

        # Timeslot to search for, YYYYMMDDHHMM
        self.dtstr = in_dtstr
        self.dater = datetime.strptime(self.dtstr, "%Y%m%d%H%M")

        # Whether to overwrite files
        self.clobber = clobber

        # Aerosol quality control flag
        self.aerosol_qc = aerosol_qc

        # Distance to cloud for flagging
        self.dist2cloud = dist2cloud

        # Aerosol land sea mask flag
        self.aerosol_landsea = aerosol_landsea

        # Are we doing cesium or quicklooks
        self.cesium = cesium

        # Set output resolution in degrees for cesium
        self.out_img_scl_cs = out_img_scl_cs

        # Should we flip data? ORAC outputs transposed in both directions
        self.flip_data = flip_data

        # Shall we make subdirectories based on date
        self.auto_out = auto_out

        # Set image output size in pixels
        self.out_img_pix = out_img_pix

        # Set lat/lon limits of output image
        self.out_img_ll = out_img_ll

        # List of variables to be read from primary file
        self.pvar = pvar

        # List of variables to be read from secondary file
        self.svar = svar

        # Solar zenith threshold for merging false color + IR data
        self.sza_thresh = sza_thresh

        # Scale factor for reflectance normalisation
        self.perc_max = perc_max

        # Resampling method, 'bilin' or 'near' supported
        self.res_meth = res_meth


def load_orac(in_file, var_list, flipper):
    """Read ORAC variables from file into dict.
    Inputs:
     - in_file: String, filename to read.
     - var_list: List, set of variables to read.
    Outputs:
     - var_dict: Dictionary of variable, data, pairs.
     """
    try:
        fid = Dataset(in_file, 'r')
    except FileNotFoundError:
        print(f"Error, can't open input file! {in_file}")
        raise
    var_dict = {}
    for variable in var_list:
        try:
            data = np.array(fid[variable])
            if flipper:
                data = np.fliplr(np.flipud(data))
            # Remove fill value pixels
            data = np.where(data > -100, data, np.nan)
            var_dict[variable] = data
            logging.info(f' - Read {variable}')
        except IndexError:
            logging.info(f' - Variable {variable} not found in file {in_file}')
    fid.close()
    return var_dict


def set_output_dir(odir_top, indate):
    """Define and create output directory based on the date.
    Inputs:
     - outdir_top: String, top-level output directory.
     - indate: Datetime, current processing timeslot.
    Returns:
     - outdir: String, correct directory for saving output."""
    outdir = f'{odir_top}/{indate.strftime("%Y/%m/%d")}/'
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    return outdir


def set_output_files_ql(odir, pri_fname, offset=17, var_out_list=('CTH', 'COT', 'dCTH', 'FC')):
    """Define and create output filenames for quicklooks based on the date.
    Inputs:
     - odir: String, output directory.
     - pri_fname: String, the filename of the ORAC primary file.
     - offset: Int, filename position offset from 'SEVIRI_ORAC' that gives timestamp.
     - var_out_list: List of strings, variables to save.
    Returns:
     - out_fnames: Dictionary, output filenames for saving.
     - need_proc: Boolean, do we need to do processing. False if all files already present."""
    # Find base filename
    base_fname = os.path.basename(pri_fname)
    pos = base_fname.find('SEVIRI_ORAC_MSG')
    dtstr = base_fname[pos+offset:pos+offset+12]
    out_fnames = {}
    for var in var_out_list:
        out_fnames[var] = f'{odir}/{dtstr}_{var}.png'
    return out_fnames


def set_output_files_cs(odir, pri_fname, var_out_list=('CTH', 'COT', 'FC', 'AOD')):
    """Define and create output filenames for cesium based on the date.
    Inputs:
     - odir: String, output directory.
     - pri_fname: String, the filename of the ORAC primary file.
     - var_out_list: List of strings, variables to save.
    Returns:
     - out_fnames: Dictionary, output filenames for saving.
     - need_proc: Boolean, do we need to do processing. False if all files already present."""

    # Find base filename
    base_fname = os.path.basename(pri_fname)
    pos = base_fname.find('.primary.nc')
    base_fname = base_fname[:pos]
    out_fnames = {}
    for var in var_out_list:
        out_fnames[var] = f'{odir}/{base_fname}_{var}.png'
    return out_fnames


def test_files(fdict):
    """Determine if all required output files are present.
    Inputs:
     - flist: Dictinary of filenames to test.
    Returns:
     - Boolean, true of all files are present."""
    if all([os.path.isfile(fdict[k]) for k in fdict.keys()]):
        return True
    else:
        return False


def apply_basic_qc(pri_data):
    """Apply basic QC to data and filter out bad pixels.
    Inputs
        -    pri_data: Dict, data from ORAC primary file.
    Returns
        -    pri_data: Dict, same data but with QC applied.
        """
    return None
