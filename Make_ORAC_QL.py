#!/usr/bin/env python
import sys
sys.path.append('/home/users/dhegedus/seviri_ql')
import seviri_additional_cloud_tests as extra_tests
from datetime import datetime
import plotting as sev_plot
from getopt import getopt
import cProfile
from glob import glob
import ql_utils
import warnings
import logging
import sys
import numpy as np

start_time = datetime.utcnow()

warnings.filterwarnings('ignore')

logfile = './NRT_QL.log'
logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)


def main(opts):
    logging.info(f'Beginning processing for {opts.dater.strftime("%Y-%m-%d %H:%M")}')
    
    # Locate primary and secondary files
    if opts.use_aerosol:
            opts.indir+='withaer'
    else:
            opts.indir+='noaer'
    inf_pri = glob(f'{opts.indir}/{opts.subdir}*{opts.dtstr}*.primary.nc')
    inf_sec = glob(f'{opts.indir}/{opts.subdir}*{opts.dtstr}*.secondary.nc')
    inf_flx = glob(f'{opts.indir}/{opts.subdir}*{opts.dtstr}*.bugsrad.nc')[0]
    
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
    outdir_cs = ql_utils.set_output_dir(opts.outdir_top, opts.dater, cesium=True)
    outdir_ql = ql_utils.set_output_dir(opts.outdir_top, opts.dater, cesium=False)
    # Set output filenames
    outfiles_ql = ql_utils.set_output_files_ql(outdir_ql, inf_pri)
    outfiles_flx_cs = ql_utils.set_output_files_flux_cs(outdir_cs, inf_flx)
    outfiles_flx_ql = ql_utils.set_output_files_flux_ql(outdir_ql, inf_flx)
    outfiles_cs = ql_utils.set_output_files_cs(outdir_cs, inf_pri)

    # Are files present? If so, no need to process
    all_done = ql_utils.test_files({**outfiles_ql, **outfiles_flx_ql, **outfiles_cs})
    if all_done:
        logging.info(f'All files are present. Halting processing.')
        return

    # Load the primary and secondary variables from ORAC output #ADDED: Load flux variables from ORAC output
    logging.info(f'Reading primary file: {inf_pri}')
    pri_data, priplat = ql_utils.load_orac(inf_pri, opts.pvar, opts.flip_data)
    logging.info(f'Reading secondary file: {inf_sec}')
    sec_data, secplat = ql_utils.load_orac(inf_sec, opts.svar, opts.flip_data)
    logging.info(f'Reading flux file: {inf_flx}')
    flx_data, flxplat = ql_utils.load_orac(inf_flx, opts.fvar, opts.flip_data)
    
    if priplat != secplat:
        raise ValueError("Platforms do not match!", priplat, secplat)
    else:
        opts.platform = priplat
        opts.title_stub = f'{opts.platform}-SEVIRI {opts.dater.strftime("%Y %m %d %H:%M")}'

    if opts.aerosol_qc is not None:
        logging.info(f'Applying additional aerosol quality filtering.')
        pri_data = extra_tests.seviri_additional_cloud_tests(pri_data, opts.aerosol_qc, opts.dist2cloud)

    pri_data = extra_tests.apply_filters(pri_data, opts)
    print('Start plotting')
    # False color image
    # Get the output data scaled into byte range
    odata_fc = sev_plot.retr_fc(pri_data, sec_data, opts.perc_max, opts.sza_thresh)
    # Resample onto appropriate grid for cesium
    res_data_fc, res_area, area_ext = sev_plot.resample_data(odata_fc, pri_data, opts)
    # Find area definition, needed for coastlines
    area_def = (res_area.proj4_string, res_area.area_extent)
    # Save the output to disk
    sev_plot.save_plot_fc(outfiles_cs['FC'], res_data_fc, opts, area_def)
    #sev_plot.save_plot_fc(outfiles_cs['FC'], res_data_fc, opts, area_def, addcoast=True)
    
    # Phase image
    # Get the output data scaled into byte range
    odata_phs = sev_plot.retr_phs(pri_data)
    # Resample onto appropriate grid for cesium
    res_data_phs, res_area, area_ext = sev_plot.resample_data(odata_phs, pri_data, opts)
    # Save the output to disk
    im = sev_plot.save_plot_phs(outfiles_cs['PHS'], res_data_phs, opts, area_def)
    sev_plot.save_plot_phs_ql(outfiles_ql['PHS'], res_data_phs, opts, im, area_ext)

    # Cloud top height image
    opts.logscl = False
    opts.varname = 'CTH'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Cloud-top height (km)'
    opts.keyticks = ['0', '2.5', '5', '7.5', '10', '12.5', '>15']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(pri_data['cth'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_ql[opts.varname], res_data_cod, opts, im, area_ext)

    # Cloud Optical depth image
    opts.logscl = True
    opts.varname = 'COT'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Cloud optical depth (550 nm)'
    opts.keyticks = ['0.1', '1', '10', '>100']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(pri_data['cot'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_ql[opts.varname], res_data_cod, opts, im, area_ext)
    
    # Cloud effective radius image
    opts.logscl = False
    opts.varname = 'CER'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Cloud effective radius (micrometer)'
    opts.keyticks = ['0',  '25', '50']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(pri_data['cer'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_ql[opts.varname], res_data_cod, opts, im, area_ext)

    # Aerosol Optical depth image
    opts.logscl = False
    opts.varname = 'AOD'
    opts.aerosol_landsea = True
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Aerosol optical depth (550 nm)'
    opts.keyticks = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '>1.2']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(pri_data['aot550'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_cs[opts.varname], res_data_cod, opts, area_ext=area_ext)
    sev_plot.save_plot_cmap_ql(outfiles_ql[opts.varname], res_data_cod, opts, im, area_ext)

    # ToA upwelling SW radiation image
    opts.logscl = False
    opts.varname = 'toa_swup'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Top of Atmosphere upwelling SW radiation'
    opts.keyticks = ['0.0', '200', '400', '600', '800', '1000']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(flx_data['toa_swup'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_flx_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_flx_ql[opts.varname], res_data_cod, opts, im, area_ext)

    # ToA downwelling SW radiation image
    opts.logscl = False
    opts.varname = 'toa_swdn'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Top of Atmosphere downwelling SW radiation'
    opts.keyticks = ['0.0', '200', '400', '600', '800', '1000','1200', '1400']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(flx_data['toa_swdn'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_flx_ql[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_flx_ql[opts.varname], res_data_cod, opts, im, area_ext)
    
    # ToA upwelling LW radiation image
    opts.logscl = False
    opts.varname = 'toa_lwup'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Top of Atmosphere upwelling LW radiation'
    opts.keyticks = ['0.0', '100', '200', '300', '400', '500']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(flx_data['toa_lwup'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_flx_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_flx_ql[opts.varname], res_data_cod, opts, im, area_ext)
    
    # BoA upwelling SW radiation image
    opts.logscl = False
    opts.varname = 'boa_swup'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Bottom of Atmosphere upwelling SW radiation'
    opts.keyticks = ['0.0', '200', '400', '600', '800', '1000']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(flx_data['boa_swup'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_flx_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_flx_ql[opts.varname], res_data_cod, opts, im, area_ext)
    
    # BoA downwelling SW radiation image
    opts.logscl = False
    opts.varname = 'boa_swdn'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Bottom of Atmosphere downwelling SW radiation'
    opts.keyticks = ['0.0', '200', '400', '600', '800', '1000']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(flx_data['boa_swdn'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_flx_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_flx_ql[opts.varname], res_data_cod, opts, im, area_ext)
    
    # BoA upwelling LW radiation image
    opts.logscl = False
    opts.varname = 'boa_lwup'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Bottom of Atmosphere upwelling LW radiation'
    opts.keyticks = ['0.0', '200', '400', '600', '800', '1000']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(flx_data['boa_lwup'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_flx_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_flx_ql[opts.varname], res_data_cod, opts, im, area_ext)
    
    # BoA downwelling LW radiation image
    opts.logscl = False
    opts.varname = 'boa_lwdn'
    opts.aerosol_landsea = False
    opts.title = opts.title_stub + ' ' + opts.varname
    opts.keytitle = 'Bottom of Atmosphere downwelling LW radiation'
    opts.keyticks = ['0.0', '200', '400', '600', '800', '1000']
    res_data_cod, res_area, area_ext = sev_plot.resample_data(flx_data['boa_lwdn'], pri_data, opts)
    im = sev_plot.save_plot_cmap(outfiles_flx_cs[opts.varname], res_data_cod, opts)
    sev_plot.save_plot_cmap_ql(outfiles_flx_ql[opts.varname], res_data_cod, opts, im, area_ext)
    

if len(sys.argv) < 2:
    raise ValueError("You did not supply a datetime for processing. Please supply as command line argument in the "
                     "format YYYYmmddHHMM.")
in_dt = sys.argv[-1]
options, operands = getopt(sys.argv[1:], "", ["outdir_top=", 'use_aerosol='])
for o,v in options:
    if o == "--outdir_top":
        outdir = v
    elif o == "--use_aerosol":
        use_aerosol = v

main_opts = ql_utils.QuickLookOpts(res_meth='bilinear', in_dtstr=in_dt, outdir_top=outdir, use_aerosol=use_aerosol)
main(main_opts)
end_time = datetime.utcnow()
print(f'Time taken: {(end_time - start_time).total_seconds():5.3f} sec')
