#from pycoast import ContourWriterAGG
from PIL import Image, ImageOps, ImageFont
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import aggdraw
import cartopy
import regionmask
    
def assign_cmap_cld(filepath):
    rgb_array = np.loadtxt(filepath, delimiter=',') # csv file with RGB values 0-1
    b3 = rgb_array[:,2] # value of blue at sample n
    b2 = rgb_array[:,2] # value of blue at sample n
    b1 = np.linspace(0, 1, len(b2)) # position of sample n - ranges from 0 to 1
    
    # Setting up columns for tuples
    g3 = rgb_array[:,1]
    g2 = rgb_array[:,1]
    g1 = np.linspace(0,1,len(g2))
    
    r3 = rgb_array[:,0]
    r2 = rgb_array[:,0]
    r1 = np.linspace(0,1,len(r2))
    
    # Creating tuples
    R = zip(r1,r2,r3)
    G = zip(g1,g2,g3)
    B = zip(b1,b2,b3)
    
    # Transposing
    RGB = zip(R,G,B)
    rgb = zip(*RGB)
    
    # Creating dictionary
    k = ['red', 'green', 'blue']
    cube1 = dict(zip(k,rgb))
    cmap = matplotlib.colors.LinearSegmentedColormap('cloudcmap', cube1)
    return cmap
    
def assign_cmap_aer():
    '''Create colormap for AOD - mimicking the Gareth's IDL colormap.
    Interpolates between the red, green and blue colors given.'''
    red=  2.55*np.array([30,  0, 100,100])
    green=2.55*np.array([0,  0, 80,  0])
    blue= 2.55*np.array([30, 80, 30,  0])
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('aercmap', np.array([red, green, blue]).T/255)
    return cmap

def adjust_lightness(color, amount=0.5):
    '''Lighten or darken an rgb color.
    color: can be rgb tuple, matplotlib color string, hex string.
    amount: darkens for <1, lightens for >1.'''
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])
    
def create_cmap_flx(cmap):
    '''Create colormap for fluxes - like cloud colormap, but darken the blue/purple side.'''  
    new_rows = cmap(np.arange(cmap.N))
    darken = np.linspace(0.5, 1, 128)
    for row in np.arange(0, 128):
        new_rows[row,:3] = adjust_lightness(cmap(row), darken[row])
    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('flxcmap', new_rows)
    return my_cmap
    
def colorbar_plotting(opts):
    fig, ax = plt.subplots(figsize=(1, 3))
    fig.patch.set_alpha(0)
    if opts.varname == 'AOD':
        cur_cmap = assign_cmap_aer()
    else:
        cur_cmap = assign_cmap_cld(opts.cldcmappath)
              
        if 'toa' in opts.varname or 'boa' in opts.varname:
            cur_cmap = create_cmap_flx(cur_cmap)
        cmaplist = [cur_cmap(i) for i in range(cur_cmap.N)]
        cur_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                      'discmap', cmaplist, (len(opts.keyticks)-1)*2) 
        
    if opts.logscl:
        rng_min = np.log10(opts.outlims[opts.varname][0])
        rng_max = np.log10(opts.outlims[opts.varname][1])
        cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=rng_min, vmax=rng_max), cmap=cur_cmap),
                 cax=ax, orientation='vertical', fraction=100, shrink=1.1, aspect=5)
        cbar.ax.set_yticks(np.arange(rng_min, rng_max+1, 1), color='white')
        cbar.ax.set_yticklabels(opts.keyticks, color='white', fontsize=10)
        cbar.outline.set_edgecolor('white')
        cbar.outline.set_linewidth(0.5)
        ax.yaxis.set_ticks_position('left')
        cbar.ax.tick_params(axis='y',colors='white')
    else:
        rng_min = opts.outlims[opts.varname][0]
        rng_max = opts.outlims[opts.varname][1]
        print(rng_min, rng_max)
        print((rng_max-rng_min)/(len(opts.keyticks)-1))
        print(np.arange(rng_min, rng_max+(rng_max-rng_min)/(len(opts.keyticks)-1), (rng_max-rng_min)/(len(opts.keyticks)-1)))
        cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=rng_min, vmax=rng_max), cmap=cur_cmap),
                 cax=ax, orientation='vertical', fraction=100, shrink=1.1, aspect=5)
        cbar.outline.set_edgecolor('white')
        cbar.outline.set_linewidth(0.5)
        cbar.add_lines(levels=np.arange(rng_min, rng_max+(rng_max-rng_min)/(len(opts.keyticks)-1), (rng_max-rng_min)/(len(opts.keyticks)-1)/2), 
                      colors=np.repeat('white', len(opts.keyticks)*2), 
                      linewidths=np.repeat(0.5,len(opts.keyticks)*2))
        ax.yaxis.set_ticks_position('left')
        cbar.ax.tick_params(axis='y',colors='white')
        cbar.ax.set_yticks(np.arange(rng_min, rng_max+(rng_max-rng_min)/(len(opts.keyticks)-1), (rng_max-rng_min)/(len(opts.keyticks)-1)), color='white')
        cbar.ax.set_yticklabels(opts.keyticks, color='white', fontsize=10)
    
    fig.tight_layout()
    fig.savefig(opts.cbar_path+opts.title_stub+opts.varname+'_colourbar.png', bbox_inches='tight', pad_inches = 0)   
    
def colorbar_phs_plotting(opts):
    fig, ax = plt.subplots(figsize=(1, 3))
    fig.patch.set_alpha(0)
    cur_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('phasecmap', [[100,0,0], [0,0,100]], 2)
    cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0, vmax=1), cmap=cur_cmap),
             cax=ax, orientation='vertical', fraction=100, shrink=1.1, aspect=5)
    cbar.outline.set_edgecolor('white')
    cbar.outline.set_linewidth(0.5)
    cbar.add_lines(levels=np.arange(0.5, 1, 0.5), 
                      colors=np.repeat('white',1), 
                      linewidths=np.repeat(0.5,1))
    ax.yaxis.set_ticks_position('left')
    cbar.ax.tick_params(axis='y',colors='white')
    cbar.ax.set_yticks(np.arange(0.25, 1.25, 0.5), color='white')
    cbar.ax.set_yticklabels( ['liquid', 'ice'])
    
    fig.tight_layout()
    fig.savefig(opts.cbar_path+opts.title_stub+opts.varname+'_colourbar.png', bbox_inches='tight', pad_inches = 0) 

def save_plot_fc(fname, data, opts, area_ext, area_def, fill_value=-999, addcoast=False):
    """Save plot-ready false color/phase class data to cesium file."""
    save_fc = make_alpha(data[:, :, ::-1])
    img = Image.fromarray(save_fc)
    if addcoast:
        img = Image.fromarray(save_fc)
        cw = ContourWriterAGG(opts.coast_dir)
        cw.add_coastlines(img, area_def, resolution='l', level=4)
        #cw.add_borders(img, area_def)
    
    #img.save(fname)
    fig, ax = plt.subplots(dpi=300)
    ax.imshow(img)
    ax.axis('off')
    fig.savefig(fname, bbox_inches='tight', transparent=True)
    return img
    
def save_plot_fc_ql(fname, data, opts, im, area_ext):
    """Save plot-ready phase class data to quicklook file."""
    import matplotlib.ticker as mticker
    import matplotlib.patches as mpatches
    
    # Set up plot
    fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    ax.coastlines(lw=0.3, color='r')
    ax.set_xticks(np.arange(round(area_ext[0].item(),-1),round(area_ext[1].item(), -1),20), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(round(area_ext[2].item(),-1),math.ceil(area_ext[3].item()/10.0)*10,10), crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_title(opts.title)
    
    # Overlay image of data
    ax.imshow(im, extent=area_ext)
    
    # Save quicklook figure
    fig.savefig(fname, bbox_inches='tight')

def save_plot_phs(fname, data, opts, area_def, addcoast=False):
    """Save plot-ready phase classification data to file."""
    save_fc = make_alpha(data[:, :, ::-1])
    save_fc[:,:,3] = np.where((save_fc[:, :,0]==0) & (save_fc[:, :,1]==0) & (save_fc[:, :,2]==0), 0, save_fc[:,:,3])
    #save_fc[:,:,3] = np.where((save_fc[:, :,0]==25) & (save_fc[:, :,1]==25) & (save_fc[:, :,2]==25), 0, save_fc[:,:,3])
    img = Image.fromarray(save_fc)
    if addcoast:
        cw = ContourWriterAGG(opts.coast_dir)
        cw.add_coastlines(img, area_def, resolution='l', level=4)
        #cw.add_borders(img, area_def)        
    #img.save(fname, optimize=True, quality=85)
    fig, ax = plt.subplots(dpi=300)
    ax.imshow(img)
    ax.axis('off')
    fig.savefig(fname, bbox_inches='tight', transparent=True)
    return img

def save_plot_phs_ql(fname, data, opts, im, area_ext):
    """Save plot-ready phase class data to quicklook file."""
    import matplotlib.ticker as mticker
    import matplotlib.patches as mpatches
    
    # Set up plot
    fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    cax = fig.add_axes([0.97, 0.175, 0.02, 0.64]) 
    ax.coastlines(lw=0.3)
    ax.set_xticks(np.arange(round(area_ext[0].item(),-1),round(area_ext[1].item(), -1),20), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(round(area_ext[2].item(),-1),math.ceil(area_ext[3].item()/10.0)*10,10), crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_title(opts.title)
    cur_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('phasecmap', [[100,0,0], [0,0,100]], 2)
    cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0, vmax=1), cmap=cur_cmap),
             cax=cax, orientation='vertical', label=opts.varname)
    cbar.ax.set_yticks(np.arange(0.25, 1.25, 0.5))
    cbar.ax.set_yticklabels( ['liquid', 'ice'])
    
    # Overlay image of data
    ax.imshow(im, extent=area_ext)
    
    # Save quicklook figure
    fig.savefig(fname, bbox_inches='tight')
        
        

def save_plot_cmap(fname, data, opts, area_ext=None, fill_value=-999, data_filt=None):
    """Save plot-ready AOD/cloud data to cesium file."""
    data_proc = np.copy(data)
    
    if data_filt is not None:
        data_proc = np.where(data_filt == 0, data_proc, fill_value)
    
    # Get the colormap
    if opts.varname == 'AOD':
        cur_cmap = assign_cmap_aer()
        #cur_cmap = opts.cmap_aer.copy()
    else:
        cur_cmap = assign_cmap_cld(opts.cldcmappath)
        #cur_cmap = opts.cmap_cld.copy()
        if 'toa' in opts.varname or 'boa' in opts.varname:
            cur_cmap = create_cmap_flx(cur_cmap)
            
    # Find the correct range limits for a given variable, in log scale if needed
    if opts.logscl:
        rng_min = np.log10(opts.outlims[opts.varname][0])
        rng_max = np.log10(opts.outlims[opts.varname][1])
        data_proc = np.log10(data_proc)
        data_proc = np.where(np.isfinite(data_proc), data_proc, fill_value)
    else:
        rng_min = opts.outlims[opts.varname][0]
        rng_max = opts.outlims[opts.varname][1]

    # Set data lims for plotting and init mask
    #mask = data_proc.copy()
    data_proc = np.where(data_proc < 0, fill_value, data_proc)
    data_proc = np.where((data_proc < rng_min) & (data_proc > 0), rng_min, data_proc)
    data_proc = np.where(data_proc > rng_max, rng_max, data_proc)
    mask = data_proc.copy()
    # Populate mask
    mask = np.where(mask == fill_value, 0, 255)
    
    # Normalise data
    data_proc = data_proc / rng_max
    
    
    # Make the image and save
    im = np.uint8(cur_cmap(data_proc) * 255)
    im[:, :, 3] = mask
    
    if opts.aerosol_landsea:
        land = regionmask.defined_regions.natural_earth_v5_0_0.land_10
        lat = np.linspace(area_ext[2].item(),area_ext[3].item(), data_proc.shape[0])
        lon = np.linspace(area_ext[0].item(),area_ext[1].item(), data_proc.shape[1])
        mask2 = land.mask(lon, lat)
        mask = np.where(np.flipud(mask2)==0, 0,mask)
        
    im[:, :, 3] = mask   
    img = Image.fromarray(im)
    #img.save(fname, optimize=True)
    
    fig, ax = plt.subplots(dpi=300)
    ax.imshow(img)
    ax.axis('off')
    fig.savefig(fname, bbox_inches='tight', transparent=True)
    return im
    
def save_plot_cmap_ql(fname, data, opts, im, area_ext):
    """Save plot-ready AOD/clouddata to quicklook file."""
    import matplotlib.ticker as mticker
    fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    cax = fig.add_axes([0.97, 0.175, 0.02, 0.64]) 
    ax.coastlines(lw=0.3)
    ax.set_xticks(np.arange(round(area_ext[0].item(),-1),round(area_ext[1].item(), -1),20), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(round(area_ext[2].item(),-1),math.ceil(area_ext[3].item()/10.0)*10,10), crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_title(opts.title) 
    
    # Get the colormap
    if opts.varname == 'AOD':
        cur_cmap = assign_cmap_aer()
        #cur_cmap = opts.cmap_aer.copy()
    else:
        cur_cmap = assign_cmap_cld(opts.cldcmappath)
        #cur_cmap = opts.cmap_cld.copy()
        if 'toa' in opts.varname or 'boa' in opts.varname:
            cur_cmap = create_cmap_flx(cur_cmap)
            
    # Find the correct range limits for a given variable, in log scale if needed
    if opts.logscl:
        rng_min = np.log10(opts.outlims[opts.varname][0])
        rng_max = np.log10(opts.outlims[opts.varname][1])
        cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=rng_min, vmax=rng_max), cmap=cur_cmap),
             cax=cax, orientation='vertical', label=opts.varname)
        cbar.ax.set_yticks(np.arange(rng_min, rng_max+1, 1))
        cbar.ax.set_yticklabels(opts.keyticks)
    else:
        rng_min = opts.outlims[opts.varname][0]
        rng_max = opts.outlims[opts.varname][1]
        cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=rng_min, vmax=rng_max), cmap=cur_cmap),
             cax=cax, orientation='vertical', label=opts.keytitle)
    
    # # Overlay image of data
    ax.imshow(im, extent=area_ext, vmin=rng_min, vmax=rng_max)
    
    # Save quicklook figure   
    fig.savefig(fname, bbox_inches='tight')



def retr_fc(pridata, secdata, perc_max=99, sza_thresh=70.):
    """Retrieve false colour image in original projection.
    Inputs:
        -   pridata: Dict, ORAC primary file data.
        -   secdata: Dict, ORAC secondary file data.
        -   perc_max: Float, maximum reflectance scaling value.
        -   sza_thresh: Float, solar zenith threshold for blending IR data
    Returns:
        -   img: 3d float array, false colour image from given ORAC data.
    """
    # Find bands and scale to percentile
    b1 = secdata['reflectance_in_channel_no_1']
    b2 = secdata['reflectance_in_channel_no_2']
    b3 = secdata['reflectance_in_channel_no_3']

    b1 = b1 / np.nanpercentile(b1, perc_max)
    b2 = b2 / np.nanpercentile(b2, perc_max)
    b3 = b3 / np.nanpercentile(b3, perc_max)

    img = np.dstack((b1, b2, b3))
    pts = (img < 0).nonzero()
    img[pts] = 0
    sza = pridata['solar_zenith_view_no1']
    if np.nanmax(sza > sza_thresh):
        img = np.where(np.isnan(img), 0, img)
        irbt = secdata['brightness_temperature_in_channel_no_9']
        pts = (irbt < 0).nonzero()
        irbt[pts] = np.nan
        irbt = np.nanpercentile(irbt, perc_max) - irbt
        irbt = irbt / np.nanpercentile(irbt, perc_max)
        irbt[pts] = 0.
        scale = (sza - sza_thresh) / 20.0
        scale = np.where(sza > 90, 1.0, scale)
        scale = np.where(sza < sza_thresh, 0.0, scale)
        img[:, :, 0] = img[:, :, 0]+scale * irbt
        img[:, :, 1] = img[:, :, 1]+scale * irbt
        img[:, :, 2] = img[:, :, 2]+scale * irbt
    img = np.where(img<0,0, img)
    img = np.round(np.where(img > 1, 1., img) * 255)
    img = img.astype(np.ubyte)

    return img


def retr_phs(pridata):
    """Retrieve phase in original projection.
    Inputs:
        -   pridata: Dict, ORAC primary file data.
    Returns:
        -   img: 3d float array, phase image from given ORAC data.
    """
    # Find bands and scale to percentile
    b1 = pridata['phase']
    b2 = pridata['phase']
    b3 = pridata['phase']

    b1 = np.where(b1==2, 1, 0) #Blue = ice
    b2 = np.where(b2==0, np.nan, 0) #Green = clear
    b3 = np.where(b3==1, 1, 0) #Red = liquid
    
    
    img = np.dstack((b1, b2, b3))
    pts = (img < 0).nonzero()
    img[pts] = 0
    
    img = np.round(np.where(img > 1, 1., img) * 255)
    img = img.astype(np.ubyte)
    return img

def resample_data(indata, pridata, opts, roi=50000, fill_value=-999):
    """Transform raw SEVIRI data into required output projection.
    Inputs:
        -   indata: 2d or 3d numpy array, ORAC/false colour data in native projection.
        -   pridata: Dict,  ORAC primary data.
        -   img_size: Tuple, x and y output image size
        -   img_bnds: Tuple, lat/lon boundaries for output (lon_0, lat_0, lon_1, lat_1)
        -   meth: String, resampling method. Only bilinear and nearest supported.
    Outputs:
        -   res_img: 2d/3d numpy array, resampled to desired projection
    """
    from pyresample import create_area_def, geometry, image
    from satpy import resample
    import xarray as xr
    lats = pridata['lat']
    lats = np.where(lats > -90, lats, 180.)
    lats = np.where(np.isfinite(lats), lats, 380.)
    lons = pridata['lon']
    lons = np.where(lons > -180, lons, 380.)
    lons = np.where(np.isfinite(lons), lons, 380.)
    lons = xr.DataArray(lons, dims=["y", "x"])
    lats = xr.DataArray(lats, dims=["y", "x"])
    
    lat_max = lats.max().values
    lat_min = lats.min().values
    lon_max = lons.max().values
    lon_min = lons.min().values
    print(lon_min, lon_max, lat_min, lat_max)
    pix_height = math.floor((lat_max-lat_min)/opts.out_img_res)
    pix_width = math.floor((lon_max-lon_min)/opts.out_img_res)
    
    indata_def = geometry.SwathDefinition(lats=lats, lons=lons)
    area_def = create_area_def('test_area',
                               {'proj': 'latlong', 'lon_0': 0},
                               area_extent=(lon_min, lat_min, lon_max, lat_max),
                               width=pix_width,
                               height=pix_height)
    if len(indata.shape) == 3:
        data_xr1 = xr.DataArray(indata[:,:,0], dims=["y", "x"])
        data_xr2 = xr.DataArray(indata[:,:,1], dims=["y", "x"])
        data_xr3 = xr.DataArray(indata[:,:,2], dims=["y", "x"])
        
                               
        res1 = resample.resample(indata_def,
                                 data_xr1,
                                 area_def,
                                 resampler=opts.res_meth,
                                 reduce_data=False,
                                 radius_of_influence=roi,
                                 fill_value=fill_value,
                                 cache_dir=opts.cache_dir)
        res2 = resample.resample(indata_def,
                                 data_xr2,
                                 area_def,
                                 resampler=opts.res_meth,
                                 reduce_data=False,
                                 radius_of_influence=roi,
                                 fill_value=fill_value,
                                 cache_dir=opts.cache_dir)
        res3 = resample.resample(indata_def,
                                 data_xr3,
                                 area_def,
                                 resampler=opts.res_meth,
                                 reduce_data=False,
                                 radius_of_influence=roi,
                                 fill_value=fill_value,
                                 cache_dir=opts.cache_dir)
        res = np.dstack((res1, res2, res3))
    else:
        data_xr = xr.DataArray(indata[:,:], dims=["y", "x"])
        
        res = resample.resample(indata_def,
                                data_xr,
                                area_def,
                                resampler=opts.res_meth,
                                reduce_data=False,
                                radius_of_influence=roi,
                                fill_value=fill_value,
                                cache_dir=opts.cache_dir)

    return res, area_def, (lon_min, lon_max, lat_min, lat_max)
'''
def resample_data(indata, pridata, opts, roi=50000):
    """Transform raw SEVIRI data into required output projection.
    Inputs:
        -   indata: 2d or 3d numpy array, ORAC/false colour data in native projection.
        -   pridata: Dict,  ORAC primary data.
        -   img_size: Tuple, x and y output image size
        -   img_bnds: Tuple, lat/lon boundaries for output (lon_0, lat_0, lon_1, lat_1)
        -   meth: String, resampling method. Only bilinear and nearest supported.
    Outputs:
        -   res_img: 2d/3d numpy array, resampled to desired projection
    """
    from pyresample import create_area_def, geometry, image

    lats = pridata['lat']
    lats = np.where(lats > -90, lats, 180.)
    lats = np.where(np.isfinite(lats), lats, 380.)
    lons = pridata['lon']
    lons = np.where(lons > -180, lons, 380.)
    lons = np.where(np.isfinite(lons), lons, 380.)

    indata_def = geometry.SwathDefinition(lats=lats, lons=lons)
    if opts.res_meth == 'near':
        swath_con = image.ImageContainerNearest(indata, indata_def, radius_of_influence=roi)
    elif opts.res_meth == 'bilin':
        swath_con = image.ImageContainerBilinear(indata, indata_def, radius_of_influence=roi)
    else:
        raise NotImplementedError('Only nearest (near) and bilinear (bilin) resampling are supported!')

    area_def = create_area_def('test_area',
                               {'proj': 'latlong', 'lon_0': 0},
                               area_extent=opts.out_img_ll,
                               width=opts.out_img_pix[0],
                               height=opts.out_img_pix[1],
                               units='degrees',)
    area_con = swath_con.resample(area_def)
    result = area_con.image_data

    return result, area_def
'''

def make_alpha(inarr, fill_value=-999):
    """Add an alpha channel to a numpy array.
    Inputs:
        -   inarr: Numpy array, the original data to mask.
        -   mask_val: Integer, the value that triggers masking.
    Outputs:
        -   outarr: Numpy array, input data with extra dimension for mask."""
    # Assume RGB if multiple bands
    if len(inarr.shape) > 2:
        mask = ~np.all(inarr == [fill_value, fill_value, fill_value], axis=-1)
    outarr = np.dstack((inarr, mask.astype(np.uint8) * 255.)).astype(np.uint8)
    return outarr
