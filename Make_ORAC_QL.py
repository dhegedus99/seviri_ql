#!/usr/bin/env python

import seviri_additional_cloud_tests as extra_tests
from pycoast import ContourWriterAGG
from datetime import datetime
import plotting as sev_plot
from PIL import Image
from glob import glob
import ql_utils
import warnings
import logging

start_time = datetime.utcnow()

warnings.filterwarnings('ignore')

logfile = './NRT_QL.log'
logging.basicConfig(filename=logfile, level=logging.DEBUG)


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
    outfiles_ql = ql_utils.set_output_files(outdir, inf_pri, False)
    outfiles_cs = ql_utils.set_output_files(outdir, inf_pri, True)

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
