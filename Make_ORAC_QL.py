#!/usr/bin/env python

import seviri_additional_cloud_tests as extra_tests
from datetime import datetime
import plotting as sev_plot
from glob import glob
import ql_utils
import warnings
import logging
import sys

start_time = datetime.utcnow()

warnings.filterwarnings('ignore')

logfile = './NRT_QL.log'
logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)


def main(opts):
    logging.info(f'Beginning processing for {opts.dater.strftime("%Y-%m-%d %H:%M")}')

    # Locate primary and secondary files
    inf_pri = glob(f'{opts.indir}/{opts.subdir}/*{opts.dtstr}*.primary.nc')
    inf_sec = glob(f'{opts.indir}/{opts.subdir}/*{opts.dtstr}*.secondary.nc')

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
    outdir = ql_utils.set_output_dir(opts.outdir_top, opts.dater)
    # Set output filenames
    outfiles_ql = ql_utils.set_output_files_ql(outdir, inf_pri)
    outfiles_cs = ql_utils.set_output_files_cs(outdir, inf_pri)

    # Are files present? If so, no need to process
    all_done = ql_utils.test_files({**outfiles_ql, **outfiles_cs})
    if all_done:
        logging.info(f'All files are present. Halting processing.')
        return

    # Load the primary and secondary variables from ORAC output
    logging.info(f'Reading primary file: {inf_pri}')
    pri_data, priplat = ql_utils.load_orac(inf_pri, opts.pvar, opts.flip_data)
    logging.info(f'Reading secondary file: {inf_sec}')
    sec_data, secplat = ql_utils.load_orac(inf_sec, opts.svar, opts.flip_data)

    if priplat != secplat:
        raise ValueError("Platforms do not match!", priplat, secplat)
    else:
        opts.platform = priplat
        opts.title_stub = f'{opts.platform}-SEVIRI {opts.dater.strftime("%Y %m %d %H:%M")}'

    if opts.aerosol_qc is not None:
        logging.info(f'Applying additional aerosol quality filtering.')
        pri_data = extra_tests.seviri_additional_cloud_tests(pri_data, opts.aerosol_qc, opts.dist2cloud)

    pri_data = extra_tests.apply_filters(pri_data, opts)

    # False color image
    # Get the output data scaled into byte range
    odata_fc = sev_plot.retr_fc(pri_data, sec_data, opts.perc_max, opts.sza_thresh)
    # Resample onto appropriate grid for cesium
    res_data_fc, res_area = sev_plot.resample_data(odata_fc, pri_data, opts)
    # Find area definition, needed for coastlines
    area_def = (res_area.proj4_string, res_area.area_extent)
    # Save the output to disk
    sev_plot.save_plot_fc(outfiles_cs['FC'], res_data_fc, opts, area_def)

    # Cloud top height image
    opts.logscl = False
    opts.varname = 'CTH'
    opts.title = opts.title_stub + opts.varname
    opts.keytitle = 'Cloud-top height (km)'
    opts.keyticks = ['0', '3', '6', '9', '12', '>15']
    res_data_cod, res_area = sev_plot.resample_data(pri_data['cth'], pri_data, opts)
    sev_plot.save_plot_cmap(outfiles_cs[opts.varname], res_data_cod, opts)

    # Cloud Optical depth image
    opts.logscl = True
    opts.varname = 'COT'
    opts.title = opts.title_stub + opts.varname
    opts.keytitle = 'Cloud optical depth (550 nm)'
    opts.keyticks = ['0.1', '1', '10', '>100']
    res_data_cod, res_area = sev_plot.resample_data(pri_data['cot'], pri_data, opts)
    sev_plot.save_plot_cmap(outfiles_cs[opts.varname], res_data_cod, opts)

    # Aerosol Optical depth image
    opts.logscl = False
    opts.varname = 'AOD'
    opts.title = opts.title_stub + opts.varname
    opts.keytitle = 'Aerosol optical depth (550 nm)'
    opts.keyticks = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '>1.2']
    res_data_cod, res_area = sev_plot.resample_data(pri_data['aot550'], pri_data, opts)
    sev_plot.save_plot_cmap(outfiles_cs[opts.varname], res_data_cod, opts)


if len(sys.argv) < 2:
    raise ValueError("You did not supply a datetime for processing. Please supply as command line argument in the "
                     "format YYYYmmddHHMM.")
in_dt = sys.argv[1]

main_opts = ql_utils.QuickLookOpts(coast_dir=None, res_meth='bilinear', in_dtstr=in_dt)
main(main_opts)
end_time = datetime.utcnow()
print(f'Time taken: {(end_time - start_time).total_seconds():5.3f} sec')
