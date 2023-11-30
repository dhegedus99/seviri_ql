#!/usr/bin/env python
import sys
sys.path.append('/home/users/dhegedus/seviri_ql')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import plotting as sev_plot
import ql_utils


def main(opts):
    opts.logscl=False
    opts.varname='AOD'
    opts.keyticks = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '>1.2']  
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'CTH'
    opts.keyticks = ['0', '2.5', '5', '7.5', '10', '12.5', '15']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'CER'
    opts.keyticks = ['0',  '10', '20', '30', '40', '50']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'COT'
    opts.keyticks = ['-1', '-0.5', '0.0', '0.5', '1.0', '1.5', '2.0']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'toa_swup'
    opts.keyticks = ['0', '200', '400', '600', '800', '1000']
    sev_plot.colorbar_plotting(opts)

    
    opts.logscl = False
    opts.varname = 'toa_lwup'
    opts.keyticks = ['0', '100', '200', '300', '400', '500']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'boa_swup'
    opts.keyticks = ['0', '200', '400', '600', '800', '1000']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'boa_swdn'
    opts.keyticks = ['0', '200', '400', '600', '800', '1000']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'boa_lwup'
    opts.keyticks = ['0', '200', '400', '600', '800', '1000']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'boa_lwdn'
    opts.keyticks = ['0', '200', '400', '600', '800', '1000']
    sev_plot.colorbar_plotting(opts)
    
    opts.logscl = False
    opts.varname = 'PHS'
    sev_plot.colorbar_phs_plotting(opts)
    

main_opts = ql_utils.QuickLookOpts(title_stub = 'NCEO-L2-CLOUD-AEROSOL-SEVIRI_ORAC_MSG4_')
main(main_opts)

