from datetime import datetime
from netCDF4 import Dataset
from matplotlib import cm
import numpy as np
import pathlib
import logging
import os

def_pvar = ('lat', 'lon', 'cldmask', 'illum', 'solar_zenith_view_no1', 'cer', 'phase',
            'cth', 'cth_uncertainty', 'cot', 'aot550', 'aer', 'qcflag', 'lsflag', 'niter', 'costjm')

def_svar = ('reflectance_in_channel_no_1', 'reflectance_in_channel_no_2',
            'reflectance_in_channel_no_3', 'brightness_temperature_in_channel_no_9')
            
def_fvar = ('toa_swdn', 'toa_swup', 'toa_lwup', 'boa_swdn', 'boa_swup', 'boa_lwup', 'boa_lwdn')


def limdict(min_aod=0.,
            max_aod=1.2,
            min_cth=0.02,
            max_cth=15.,
            min_cer=0.,
            max_cer=50.,
            min_dcth=0.,
            max_dcth=10.,
            min_cod=0.1,
            max_cod=100.,
            min_toa_swdn=0, 
            max_toa_swdn=1250, #check these
            min_toa_swup=0, 
            max_toa_swup=900,
            min_toa_lwup=100, 
            max_toa_lwup=350,
            min_boa_lwdn=200, 
            max_boa_lwdn=500,
            min_boa_lwup=200, #check these
            max_boa_lwup=800,
            min_boa_swup=0, 
            max_boa_swup=500,
            min_boa_swdn=0, 
            max_boa_swdn=1000):

    lim_dict = {'AOD': [min_aod, max_aod],
                'CTH': [min_cth, max_cth],
                'CER': [min_cer, max_cer],
                'COT': [min_cod, max_cod],
                'dCTH': [min_dcth, max_dcth],
                'toa_swup': [min_toa_swup, max_toa_swup],
                'toa_swdn': [min_toa_swdn, max_toa_swdn],
                'toa_lwup': [min_toa_lwup, max_toa_lwup],
                'boa_swup': [min_boa_swup, max_boa_swup],
                'boa_swdn': [min_boa_swdn, max_boa_swdn],
                'boa_lwup': [min_boa_lwup, max_boa_lwup],
                'boa_lwdn': [min_boa_lwdn, max_boa_lwdn]}

    return lim_dict


class QuickLookOpts:
    def __init__(self,
                 in_dtstr='202109101100',
                 use_aerosol = True,
                 indir='/gws/pw/j07/rsgnceo/from_j05/public/nrt/nrt_part_seviri_msg3_ext/data/lv2/',
                 #indir='/home/users/dhegedus/seviri_redo/Data/seviri_msg3/nrt_processing/l2b/',
                 #indir='/gws/pw/j07/rsgnceo/Data/seviri_msg3/nrt_processing/l2b/',
                 cache_dir='/gws/pw/j07/rsgnceo/Data/seviri_msg3/nrt_processing/cache_dir/',
                 outdir_top='./TEST/',
                 #outdir_top='/gws/pw/j07/rsgnceo/public/nrt/nrt_part_seviri_msg3/quick_look_cesium//',
                 coast_dir='/home/users/dhegedus/seviri_ql/coast_shp/', 
                 clobber=False,
                 aerosol_qc=284,
                 dist2cloud=3.,
                 aerosol_landsea=False,
                 cesium=True,
                 flip_data=True,
                 auto_out=True,
                 out_img_res=0.06,
                 #out_img_pix=(2840,2170),
                 out_img_scl_cs=(358, 192),
                 #out_img_ll=(-55, 0, 30, 65),
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
                 cmap_cld=cm.viridis,
                 cldcmappath='/home/users/dhegedus/seviri_ql/cube1_0-1.csv',
                 cmap_aer=cm.inferno,
                 add_coast=True,
                 cbar_path='/gws/pw/j07/rsg_share/public/rsgnceo/nrt/nrt_part_seviri_msg3/quick_look_cesium/colour_bars/',
                 fill_value=-1,
                 font='/home/users/dhegedus/seviri_ql/Verdana.ttf'):
        # Directory containing the ORAC pri + sec files
        self.indir = indir
        self.use_aerosol = use_aerosol
            
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
        
        # Set image output resolution in degrees
        self.out_img_res = out_img_res
        
        # Set image output size in pixels
        #self.out_img_pix = out_img_pix

        # Set lat/lon limits of output image
        #self.out_img_ll = out_img_ll

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
        self.cmap_cld = cmap_cld
        self.cldcmappath = cldcmappath
        self.cmap_aer = cmap_aer

        # Flag for whether we plot coastlines
        self.add_coast = add_coast
        
        # Location of colourbar for cesium
        self.cbar_path = cbar_path

        # Set the fill value for filtering data
        self.fill_value = fill_value

        self.font = font
        

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


def set_output_dir(odir_top, indate, cesium=True):
    """Define and create output directory based on the date.
    Inputs:
     - outdir_top: String, top-level output directory.
     - indate: Datetime, current processing timeslot.
    Returns:
     - outdir: String, correct directory for saving output."""
    if cesium:
        outdir = f'{odir_top}quick_look_cesium/{indate.strftime("%Y/%m/%d")}/'
    else:
        #outdir = f'/gws/pw/j07/rsgnceo/from_j05/public/nrt/nrt_part_seviri_msg3_ext/quick_look_hires/{indate.strftime("%Y/%m/%d")}/'
        outdir = f'{odir_top}quick_look_hires/{indate.strftime("%Y/%m/%d")}/'
        
    os.makedirs(outdir, 0o775, exist_ok=True)
    
    return outdir


def set_output_files_ql(odir, pri_fname, offset=17, var_out_list=('CTH', 'COT', 'dCTH', 'FC', 'AOD', 'CER', 'PHS', 'FC')):
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
    pos = base_fname.find('.primary.nc')
    base_fname = base_fname[:pos]
    out_fnames = {}
    for var in var_out_list:
        out_fnames[var] = f'{odir}/{base_fname}_{var}.png'
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
    pos = base_fname.find('.bugsrad.nc')
    base_fname = base_fname[:pos]
    out_fnames = {}
    for var in var_out_list:
        out_fnames[var] = f'{odir}/{base_fname}_{var}.png'
    return out_fnames

def set_output_files_flux_cs(odir, flx_fname, var_out_list=('toa_swdn', 'toa_swup', 'toa_lwup', 'boa_swdn', 'boa_swup', 'boa_lwup', 'boa_lwdn')):
    """Define and create output filenames for cesium based on the date.
    Inputs:
     - odir: String, output directory.
     - pri_fname: String, the filename of the ORAC primary file.
     - var_out_list: List of strings, variables to save.
    Returns:
     - out_fnames: Dictionary, output filenames for saving.
     - need_proc: Boolean, do we need to do processing. False if all files already present."""

    # Find base filename
    base_fname = os.path.basename(flx_fname)
    pos = base_fname.find('.bugsrad.nc')
    base_fname = base_fname[:pos]
    out_fnames = {}
    for var in var_out_list:
        out_fnames[var] = f'{odir}/{base_fname}_{var}.png'
    return out_fnames


def set_output_files_cs(odir, pri_fname, var_out_list=('CTH', 'COT', 'FC', 'AOD', 'CER', 'PHS')):
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
        out_fnames[var] = f'{odir}{base_fname}_{var}.png'
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
