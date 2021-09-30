import matplotlib.pyplot as plt
import numpy as np


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
        irbt = secdata['brightness_temperature_in_channel_no_9']
        pts = (irbt < 0).nonzero()
        irbt[pts] = np.nan

        irbt = np.nanpercentile(irbt, perc_max) - irbt
        irbt = irbt / np.nanpercentile(irbt, perc_max)
        irbt[pts] = 0.
        scale = (sza - sza_thresh) / 20.0
        scale = np.where(sza > 90, 1.0, scale)
        scale = np.where(sza < sza_thresh, 0.0, scale)
        img[:, :, 0] = img[:, :, 0] + scale * irbt
        img[:, :, 1] = img[:, :, 1] + scale * irbt
        img[:, :, 2] = img[:, :, 2] + scale * irbt

    print(np.nanmin(img), np.nanmean(img), np.nanmax(img))
    img = np.round(np.where(img > 1, 1., img) * 255)
    print(np.nanmin(img), np.nanmean(img), np.nanmax(img))
    img = img.astype(np.ubyte)
    print(np.nanmin(img), np.nanmean(img), np.nanmax(img))

    return img


def resample_data(indata, pridata, img_size, img_bnds, meth='bilin'):
    """Transform raw SEVIRI data into required output projection.
    Inputs:
        -   indata: 2d or 3d numpy array, ORAc/false colour data in native projection.
        -   pridata: Dict,  ORAC primary data.
        -   img_size: Tuple, x and y output image size
        -   img_bnds: Tuple, lat/lon boundaries for output (lon_0, lat_0, lon_1, lat_1)
        -   meth: String, determine resampling method. Currently only bilinear supported.
    Outputs:
        -   res_img: 2d/3d numpy array, resampled to desired projection
    """
    from pyresample import create_area_def, geometry, image

    indata_def = geometry.SwathDefinition(lats=pridata['lat'], lons=pridata['lon'])
    swath_con = image.ImageContainerNearest(indata, indata_def, radius_of_influence=25000)

    area_def = create_area_def('test_area',
                               {'proj': 'latlong', 'lon_0': 0},
                               area_extent=img_bnds,
                               width=img_size[0],
                               height=img_size[1],
                               units='degrees',)
    area_con = swath_con.resample(area_def)
    result = area_con.image_data
    return result, area_def


def make_alpha(inarr, mask_val=0):
    """Add an alpha channel to a numpy array.
    Inputs:
        -   inarr: Numpy array, the original data to mask.
        -   mask_val: Integer, the value that triggers masking.
    Outputs:
        -   outarr: Numpy array, input data with extra dimension for mask."""
    # Assume RGB if multiple bands
    if len(inarr.shape) > 2:
        mask = ~np.all(inarr == [mask_val, mask_val, mask_val], axis=-1)
    outarr = np.dstack((inarr, mask.astype(np.uint8) * 255.)).astype(np.uint8)
    return outarr
