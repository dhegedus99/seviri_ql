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


def fmi_cloud_qc_1km():
    return
    raise NotImplementedError('Currently fmi_cloud_qc_1km is not implemented')


def opening_test():
    return
    raise NotImplementedError('Currently opening_test is not implemented')


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
                                  dist2cld=3,
                                  minclr=5):
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
        fmi_cloud_qc_1km()

    # The 1 km data seems to have a lot of residual cloud "speckle", which
    # looks like so called "salt noise" in the AOD imagery... deal with it
    # using the morphological opening transformation to highlight and
    # remove this noise
    if _bit_check(qcselection, bits['aod_open']):
        logging.info(f'Applying AOD opening test.')
        opening_test()
    # If required, do same for effective radius data
    if _bit_check(qcselection, bits['ref_open']):
        logging.info(f'Applying ref opening test.')
        opening_test()

    # Create and set up new QA datalayer
    qual_arr = np.zeros_like(data_dict['qcflag'])
    if _bit_check(qcselection, bits['conv']):
        qual_arr = qual_arr + (data_dict['niter'] * bit_vals['conv'])
    if _bit_check(qcselection, bits['cost']):
        qual_arr = qual_arr + (data_dict['costjm'] * bit_vals['cost'])
    if _bit_check(qcselection, bits['cld_edge']):
        qual_arr = qual_arr + (data_dict['niter'] * bit_vals['cld_edge'])
    if _bit_check(qcselection, bits['sd_aod']):
        qual_arr = qual_arr + (data_dict['niter'] * bit_vals['sd_aod'])
    if _bit_check(qcselection, bits['aod_open']):
        qual_arr = qual_arr + (data_dict['niter'] * bit_vals['aod_open'])
    if _bit_check(qcselection, bits['ref_open']):
        qual_arr = qual_arr + (data_dict['niter'] * bit_vals['ref_open'])

    # Perform cloud distance parameter
    if dist2cld > 0:
        clddist = orac_1km_dist_to_cloud(cld, qual_arr, window=2*dist2cld+1)
        clddist = np.where(clddist < dist2cld, clddist, 0)
        data_dict['clddist'] = clddist

    if _bit_check(qcselection, bits['cld_dist']):
        qual_arr = qual_arr + (data_dict['clddist'] * bit_vals['cld_dist'])

    data_dict['qcflag'] = data_dict['qcflag'] + 64 * qual_arr
    return data_dict

