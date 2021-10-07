#!/usr/bin/env python

import seviri_additional_cloud_tests as extra_tests
from pycoast import ContourWriterAGG
from datetime import datetime
import plotting as sev_plot
from netCDF4 import Dataset
from PIL import Image
from glob import glob
import numpy as np
import ql_utils
import warnings
import pathlib
import logging
import os

start_time = datetime.utcnow()

warnings.filterwarnings('ignore')

logfile = './NRT_QL.log'
logging.basicConfig(filename=logfile, level=logging.DEBUG)


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


def set_output_files(odir, pri_fname, cesium):
    """Define and create output directory based on the date.
    Inputs:
     - odir: String, output directory.
     - pri_fname: String, the filename of the ORAC primary file.
     - cesium: Bool, if true then set cesium filenames otherwise set QL.
    Returns:
     - out_fnames: Dictionary, output filenames for saving.
     - need_proc: Boolean, do we need to do processing. False if all files already present."""
    var_out_list = ['CTH', 'COD', 'FC', 'AOD']
    # Find base filename
    base_fname = os.path.basename(pri_fname)
    pos = base_fname.find('.primary.nc')
    base_fname = base_fname[:pos]
    out_fnames = {}
    for var in var_out_list:
        out_fnames[var] = f'{odir}/{base_fname}_{var}'
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


def main(opts):
    logging.info(f'Beginning processing for {opts.dater.strftime("%Y-%m-%d %H:%M")}')

    # Locate primary and secondary files
    inf_pri = glob(f'{opts.indir}/*{opts.dtstr}*.primary.nc')
    inf_sec = glob(f'{opts.indir}/*{opts.dtstr}*.secondary.nc')

    # Check that files exist, raise error if not
    if len(inf_pri) < 1:
        raise OSError(f"Error: Cannot find primary file for {opts.dtstr}!")
    else:
        inf_pri = inf_pri[0]
    if len(inf_sec) < 1:
        raise OSError(f"Error: Cannot find secondary file for {opts.dtstr}!")
    else:
        inf_sec = inf_sec[0]

    logging.info(f'SEVIRI_MK_NRT_QUICKLOOKS: directories')
    logging.info(f' - Input: {opts.indir}')
    logging.info(f' - Output: {opts.outdir_top}')

    # Set up directories
    logging.info(f'Setting output files and directories')
    outdir = set_output_dir(opts.outdir_top, opts.dater)
    # Set output filenames
    outfiles_ql = set_output_files(outdir, inf_pri, False)
    outfiles_cs = set_output_files(outdir, inf_pri, True)

    # Are files present? If so, no need to process
    all_done = test_files({**outfiles_ql, **outfiles_cs})
    if all_done:
        logging.info(f'All files are present. Halting processing.')
        return

    # Load the primary and secondary variables from ORAC output
    logging.info(f'Reading primary file: {inf_pri}')
    pri_data = load_orac(inf_pri, opts.pvar, opts.flip_data)
    logging.info(f'Reading secondary file: {inf_sec}')
    sec_data = load_orac(inf_sec, opts.svar, opts.flip_data)

    if opts.aerosol_qc is not None:
        logging.info(f'Applying additional aerosol quality filtering.')
        extra_tests.seviri_additional_cloud_tests(pri_data, opts.aerosol_qc, opts.dist2cloud)

    odata = sev_plot.retr_fc(pri_data, sec_data)
    res_data, res_area = sev_plot.resample_data(odata, pri_data, opts.out_img_pix, opts.out_img_ll)
    area_def = (res_area.proj4_string, res_area.area_extent)
    save_fc = sev_plot.make_alpha(res_data[:, :, ::-1])
    img = Image.fromarray(save_fc)

    cw = ContourWriterAGG(opts.coast_dir)
    cw.add_coastlines(img, area_def, resolution='l', level=4)
    cw.add_borders(img, area_def)
    img.save("out_res_border.png")

    return


main_opts = ql_utils.QuickLook_Opts()

main(main_opts)

end_time = datetime.utcnow()

print(end_time - start_time)
