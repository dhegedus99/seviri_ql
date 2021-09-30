import numpy as np
import logging


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
    raise NotImplementedError('Currently fmi_cloud_qc_1km is not implemented')


def opening_test():
    raise NotImplementedError('Currently opening_test is not implemented')


def seviri_additional_cloud_tests(data_dict,
                                  qcselection=255,
                                  dist2cld=0,
                                  minclr=5):
    """Apply extra cloud / QC filtering to the ORAC data.
    Inputs:
        -   Stuff
    Outputs:
        -   Stuff
    """

    # Define bit - mask values for each QC test
    bits = {'conv': 0,
            'cost': 1.,
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

    data_dict['qcflag'] = data_dict['qcflag'] + 64 * qual_arr
    return data_dict

