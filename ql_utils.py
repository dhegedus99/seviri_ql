from datetime import datetime
from netCDF4 import Dataset
from matplotlib import cm
import numpy as np
import pathlib
import logging
import os

def_pvar = ('lat', 'lon', 'cldmask', 'illum', 'solar_zenith_view_no1',
            'cth', 'cth_uncertainty', 'cot', 'aot550', 'aer', 'qcflag', 'lsflag', 'niter', 'costjm')

def_svar = ('reflectance_in_channel_no_1', 'reflectance_in_channel_no_2',
            'reflectance_in_channel_no_3', 'brightness_temperature_in_channel_no_9')
            
def_fvar = ('toa_swdn', 'toa_swup', 'toa_lwup', 'boa_swdn', 'boa_swup', 'boa_lwup', 'boa_lwdn')


def limdict(min_aod=0.01,
            max_aod=1.2,
            min_cth=0.02,
            max_cth=15.,
            min_cer=0.,
            max_cer=50.,
            min_dcth=0,
            max_dcth=10.,
            min_cod=0.1,
            max_cod=100.):

    lim_dict = {'AOD': [min_aod, max_aod],
                'CTH': [min_cth, max_cth],
                'CER': [min_cer, max_cer],
                'COT': [min_cod, max_cod],
                'dCTH': [min_dcth, max_dcth]}

    return lim_dict


class QuickLookOpts:
    def __init__(self,
                 in_dtstr='202109101100',
                 indir='/gws/pw/j07/rsgnceo/Data/seviri_msg3/nrt_processing/l2b/',
                 cache_dir='/gws/pw/j07/rsgnceo/Data/seviri_msg3/nrt_processing/cache_dir/',
                 outdir_top='./TEST/',
                 #outdir_top='/gws/pw/j07/rsgnceo/public/nrt/nrt_part_seviri_msg3/quick_look_cesium//',
                 coast_dir=None,
                 clobber=False,
                 aerosol_qc=284,
                 dist2cloud=3.,
                 aerosol_landsea=False,
                 cesium=True,
                 flip_data=True,
                 auto_out=True,
                 out_img_pix=(1420, 601),
                 out_img_scl_cs=(358, 192),
                 out_img_ll=(-55, 29, 30, 65),
                 pvar=def_pvar,
                 svar=def_svar,
                 fvar=def_fvar,
                 sza_thresh=70.,
                 perc_max=99.,
                 res_meth='nearest',
                 logscl=False,
                 keyticks=None,
                 keytitle='',
                 title_stub='',
                 title='',
                 outlims=limdict(),
                 platform='',
                 varname='',
                 cmap=cm.viridis,
                 add_coast=False,
                 fill_value=-1):
        # Directory containing the ORAC pri + sec files
        self.indir = indir
        # Top level directory for the output files
        self.outdir_top = outdir_top

        # Directory containing coastline shapefiles
        self.coast_dir = coast_dir
        
        # Directory to store resampling cache, speeds up subsequent runs
        self.cache_dir = cache_dir

        # Timeslot to search for, YYYYMMDDHHMM
        self.dtstr = in_dtstr
        self.dater = datetime.strptime(self.dtstr, "%Y%m%d%H%M")
        self.subdir = self.dater.strftime("%Y/%m/%d/")

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
        
        # List of variables to be read from flux file
        self.fvar = fvar

        # Solar zenith threshold for merging false color + IR data
        self.sza_thresh = sza_thresh

        # Scale factor for reflectance normalisation
        self.perc_max = perc_max

        # Resampling method, 'bilin' or 'near' supported
        self.res_meth = res_meth

        # Plotting scale, log if True or linear if False
        self.logscl = logscl

        # Tickmarks for colorbar, list of values
        self.keyticks = keyticks

        # Legend title string
        self.keytitle = keytitle

        # Plot title stub common across plots
        self.title_stub = title_stub

        # Plot title string
        self.title = title

        # Max min values for plotting. If not supplied, use default class.
        self.outlims = outlims

        # Platform name
        self.platform = platform

        # Name of variable being plotted
        self.varname = varname

        # Colormap for plotting
        self.cmap = cmap

        # Flag for whether we plot coastlines
        self.add_coast = add_coast

        # Set the fill value for filtering data
        self.fill_value = fill_value


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
            data = np.squeeze(np.array(fid[variable]))
            if flipper:
                data = np.fliplr(np.flipud(data))
            # Remove fill value pixels
            data = np.where(data > -100, data, 0)
            var_dict[variable] = data
            logging.info(f' - Read {variable}')
        except IndexError:
            logging.info(f' - Variable {variable} not found in file {in_file}')
    plat = fid.Platform
    fid.close()
    return var_dict, plat


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
    dtstr = base_fname[pos + offset:pos + offset + 12]
    out_fnames = {}
    for var in var_out_list:
        out_fnames[var] = f'{odir}/{dtstr}_{var}.png'
    return out_fnames
    
def set_output_files_flux_ql(odir, flx_fname, offset=17, var_out_list=('toa_swdn', 'toa_swup', 'toa_lwup', 'boa_swdn', 'boa_swup', 'boa_lwup', 'boa_lwdn')):
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
    base_fname = os.path.basename(flx_fname)
    pos = base_fname.find('SEVIRI_ORAC_MSG')
    dtstr = base_fname[pos + offset:pos + offset + 12]
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
