#!/usr/bin/env python

import seviri_additional_cloud_tests as extra_tests
from datetime import datetime
import plotting as sev_plot
from glob import glob
import ql_utils
import warnings
import logging

start_time = datetime.utcnow()

warnings.filterwarnings('ignore')

logfile = './NRT_QL.log'
logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)


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
    pri_data = ql_utils.load_orac(inf_pri, opts.pvar, opts.flip_data)
    logging.info(f'Reading secondary file: {inf_sec}')
    sec_data = ql_utils.load_orac(inf_sec, opts.svar, opts.flip_data)

    print(pri_data.keys())

    if opts.aerosol_qc is not None:
        logging.info(f'Applying additional aerosol quality filtering.')
       # extra_tests.seviri_additional_cloud_tests(pri_data, opts.aerosol_qc, opts.dist2cloud)

    # Get the output data scaled into byte range
    odata_fc = sev_plot.retr_fc(pri_data, sec_data, opts.perc_max, opts.sza_thresh)

    # Resample onto appropriate grid for cesium
    res_data, res_area = sev_plot.resample_data(odata_fc, pri_data, opts)
    # Find area definition, needed for coastlines
    area_def = (res_area.proj4_string, res_area.area_extent)
    # Save the output to disk
    sev_plot.save_plot(outfiles_cs['FC'], res_data, opts, area_def)

    return


main_opts = ql_utils.QuickLookOpts(coast_dir='C:/Users/EUMETCAST#/Documents/coast/', res_meth='bilin')
#main_opts = ql_utils.QuickLookOpts(coast_dir=None, res_meth='bilin')

main(main_opts)

end_time = datetime.utcnow()

print(f'Time taken: {(end_time - start_time).total_seconds():5.3f} sec')
