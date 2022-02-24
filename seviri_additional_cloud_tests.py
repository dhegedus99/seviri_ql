from scipy.ndimage.filters import generic_filter
from scipy.signal import convolve2d
from numba import jit, prange
import numpy as np
import logging


@jit(nopython=True)
def _make_dist_arr(window_size):
    """Make an array of window_size x window_size containing distances to central pixel."""
    sizer = int((window_size-1)/2)
    dist_arr = np.zeros((int(window_size), int(window_size)))
    for x in range(0, int(window_size)):
        for y in range(0, int(window_size)):
            dist_arr[x, y] = np.sqrt((x - sizer)*(x - sizer) + (y - sizer)*(y - sizer))
    return dist_arr, sizer


def _bit_check(var, bit_n):
    """Check if a given bit is set.
    Inputs:
        -   var: Input array (int).
        -   bit_n: Nth bit to check.
    Returns:
        - Array showing whether bit 'n' is set.
    """
    return (var & (1 << bit_n)) > 0


def fmi_cloud_qc_1km(cld_data, aod_data, winsize=3, aot_thresh=0.1, tot_thresh=3):
    """Find cloud edges"""

    kernel = np.ones((winsize, winsize))
    totals = convolve2d(cld_data, kernel, mode='same', boundary='fill', fillvalue=0)

    edge = np.where((totals > tot_thresh) & (totals < 9), 1, 0)
    aod_data = np.copy(aod_data)
    aod_data = np.where(aod_data > 0, aod_data, np.nan)
    stdev = generic_filter(aod_data, np.nanstd, size=winsize)

    aod_stdv = np.where(stdev > aot_thresh, 1, 0)

    return edge, aod_stdv


def opening_test(data, thresh, winsize=5):
    """Run the opening test"""
    from cv2 import morphologyEx, MORPH_TOPHAT

    kernel = np.ones((winsize, winsize))
    tmp = np.where(data <= 0, 0, data)
    opening = morphologyEx(tmp, MORPH_TOPHAT, kernel)
    test = opening >= (thresh / 255.)
    return test


@jit(nopython=True)
def _sort_window(inidx, winsize, ranger, imsize):
    """Find correct window indices."""
    a1 = inidx - ranger
    a2 = inidx + ranger
    b1 = max(0, a1)
    b2 = min(a2, imsize - 1)
    if b1 == 0:
        a1 = abs(a1)
    else:
        a1 = 0

    if b2 == imsize - 1:
        a2 = winsize - a2 + imsize - 2
    else:
        a2 = winsize - 1

    return int(b1), int(b2), int(a1), int(a2)


def _make_gdal(fname, data):
    from osgeo import gdal
    shaper = data.shape
    output_raster = gdal.GetDriverByName('GTiff').Create(fname, shaper[1], shaper[0], 1, gdal.GDT_Float32)
    output_raster.GetRasterBand(1).WriteArray(data)
    output_raster.FlushCache()
    del output_raster


@jit(nopython=True)
def orac_1km_dist_to_cloud(cld, qual, window=31, nx=False, ny=False):
    """Find the distance of each pixel to the nearest cloud."""
    if not nx:
        nx = cld.shape[0]
    if not ny:
        ny = cld.shape[1]

    dist_arr, ranger = _make_dist_arr(window)

    max_dist = np.nanmax(dist_arr)

    cld_dist = np.zeros(cld.shape)

    # Create QA / CLMK mixed array
    tcld = np.where(qual > 0, 1, cld)
    tcld = np.where(tcld > 0, 1, 0)

    for x in prange(0, nx):
        x1, x2, i1, i2 = _sort_window(x, window, ranger, nx)
        for y in prange(0, ny):
            y1, y2, j1, j2 = _sort_window(y, window, ranger, ny)
            tdist = dist_arr[i1:i2, j1:j2].copy()
            cldi = tcld[x1:x2, y1:y2]
            tdist = np.where(cldi > 0, tdist, np.nan)
            if np.nanmax(cldi > 0):
                cld_dist[x, y] = np.nanmin(tdist)
            else:
                cld_dist[x, y] = max_dist
    return cld_dist


def seviri_additional_cloud_tests(data_dict,
                                  qcselection=284,
                                  dist2cld=3):
    """Apply extra cloud / QC filtering to the ORAC data.
    Inputs:
        -   Stuff
    Outputs:
        -   Stuff
    """

    # Define bit - mask values for each QC test
    bits = {'conv': 0,
            'cost': 1,
            'cld_edge': 2,
            'sd_aod': 3,
            'aod_open': 4,
            'ref_open': 5,
            'ice': 6,
            'ang_open': 7,
            'cld_dist': 8}

    bit_vals = {'conv': 1,
                'cost': 2,
                'cld_edge': 4,
                'sd_aod': 8,
                'aod_open': 16,
                'ref_open': 32,
                'ice': 64,
                'ang_open': 128,
                'cld_dist': 256}

    cld = data_dict['cldmask']

    if dist2cld > 0 and qcselection < 256:
        qcselection += 256

    # This routine masks out pixels which either have cloud in more than 1 / 3
    # of their surrounding pixels, or have(when combined with their neighbouring pixels)
    # have an AOD standard deviation gt 0.1
    if _bit_check(qcselection, bits['cld_edge']) or _bit_check(qcselection, bits['sd_aod']):
        logging.info(f'Applying cloud edge and AOD smoothness tests.')
        edge, aod_stdv = fmi_cloud_qc_1km(cld, data_dict['aot550'])

    # Create and set up new QA datalayer
    qual_arr = np.zeros_like(data_dict['qcflag'])

    if _bit_check(qcselection, bits['conv']):
        tmparr = np.where(data_dict['niter'] > 25, 1, 0)
        qual_arr = qual_arr + (tmparr * bit_vals['conv'])

    if _bit_check(qcselection, bits['cost']):
        tmparr = np.where(data_dict['costjm'] > 3, 1, 0)
        qual_arr = qual_arr + (tmparr * bit_vals['cost'])

    if _bit_check(qcselection, bits['cld_edge']):
        qual_arr = qual_arr + (edge * bit_vals['cld_edge'])

    if _bit_check(qcselection, bits['sd_aod']):
        qual_arr = qual_arr + (aod_stdv * bit_vals['sd_aod'])

    if _bit_check(qcselection, bits['aod_open']):
        logging.info(f'Applying AOD opening test.')
        aod_open = opening_test(data_dict['aot550'], 80)
        qual_arr = qual_arr + (aod_open * bit_vals['aod_open'])

    if _bit_check(qcselection, bits['ref_open']):
        logging.info(f'Applying ref opening test.')
        ref_open = opening_test(data_dict['aer'], 300)
        qual_arr = qual_arr + (ref_open * bit_vals['ref_open'])

    # Compute cloud distance parameter
    if dist2cld > 0:
        clddist = orac_1km_dist_to_cloud(cld, qual_arr, window=2*dist2cld+1)
        clddist = np.where(clddist < dist2cld, clddist, 0)
        clddist = np.where(clddist > 0, 1, 0)
        data_dict['clddist'] = clddist

    if _bit_check(qcselection, bits['cld_dist']):
        qual_arr = qual_arr + (data_dict['clddist'] * bit_vals['cld_dist'])
    data_dict['qcflag'] = data_dict['qcflag'] + 64 * qual_arr
    data_dict['qcflag'] = data_dict['qcflag'].astype(np.uint32)
    return data_dict


def apply_filters(pri_data, opts):
    """Apply filters to data to remove bad pixels"""

    # Filter CTH
    cth_thresh = 0.5  # Min CTH to pass
    pri_data['cth'] = np.where(pri_data['cldmask'] == 1, pri_data['cth'], -999)
    pri_data['cth'] = np.where(pri_data['cth'] >= cth_thresh, pri_data['cth'], -999)

    # Filter dCTH
    cth_thresh2 = 0
    pri_data['cth_uncertainty'] = np.where(pri_data['cth'] > cth_thresh2, pri_data['cth_uncertainty'], -999)
    pri_data['cth_uncertainty'] = np.where(pri_data['illum'] != 2, pri_data['cth_uncertainty'], -999)
    pri_data['cth_uncertainty'] = np.where(pri_data['cldmask'] == 1, pri_data['cth_uncertainty'], -999)

    # Filter COT
    cot_thresh = 0.
    _make_gdal('E:/cot_orig.tif', pri_data['cot'])
    pri_data['cot'] = np.where(pri_data['cldmask'] == 1, pri_data['cot'], -999)
    _make_gdal('E:/cot_1.tif', pri_data['cot'])
    pri_data['cot'] = np.where(pri_data['cot'] > cot_thresh, pri_data['cot'], -999)
    _make_gdal('E:/cot_2.tif', pri_data['cot'])
    pri_data['cot'] = np.where(pri_data['illum'] == 1, pri_data['cot'], -999)
    _make_gdal('E:/cot_3.tif', pri_data['cot'])

    # Filter AOD
    if opts.aerosol_qc is not None:
        aod_thresh = 0
        pri_data['aot550'] = np.where(pri_data['cldmask'] == 0, pri_data['aot550'], -999)
        pri_data['aot550'] = np.where(pri_data['aot550'] > aod_thresh, pri_data['aot550'], -999)
        pri_data['aot550'] = np.where(pri_data['illum'] == 1, pri_data['aot550'], -999)
        pri_data['aot550'] = np.where(np.bitwise_and(pri_data['qcflag'], 3328), -999, pri_data['aot550'])

    return pri_data
