from datetime import datetime

def_pvar = ('lat', 'lon', 'cldmask', 'illum', 'solar_zenith_view_no1',
            'cth', 'cot', 'aot550', 'aer', 'qcflag', 'lsflag', 'niter', 'costjm')

def_svar = ('reflectance_in_channel_no_1', 'reflectance_in_channel_no_2',
            'reflectance_in_channel_no_3', 'brightness_temperature_in_channel_no_9')


class QuickLook_Opts:
    def __init__(self,
                 in_dtstr='202109101100',
                 indir='./',
                 outdir_top='./TEST/',
                 coast_dir='C:/Users/simon/OneDrive/Documents/Shapefiles/',
                 clobber=False,
                 aerosol_qc=None,
                 dist2cloud=10,
                 aerosol_landsea=False,
                 cesium=True,
                 flip_data=True,
                 auto_out=True,
                 out_img_pix=(1700, 597),
                 out_img_ll=(-71.8154, 29.1062, 30.1846, 64.7465),
                 pvar=def_pvar,
                 svar=def_svar):
        # Directory containing the ORAC pri + sec files
        self.indir = indir
        # Top level directory for the output files
        self.outdir_top = outdir_top

        # Directory containing coastline shapefiles
        self.coast_dir = coast_dir

        # Timeslot to search for, YYYYMMDDHHMM
        self.dtstr = in_dtstr
        self.dater = datetime.strptime(self.dtstr, "%Y%m%d%H%M")

        # Whether to overwrite files
        self.clobber = clobber

        # Aerosol quality control flag
        self.aerosol_qc = aerosol_qc

        # Distance to cloud for flagging
        self.dist2cloud = dist2cloud

        # Aerosol land sea mask flag
        self.aerosol_landsea = aerosol_landsea

        # Are we doing cesium or quicklooks
        self.cesium = cesium

        # Should we flip data? ORAC outputs transposed in both directions
        self.flip_data = flip_data

        # Shall we make subdirectories based on date
        self.auto_out = auto_out

        # Set image output size in pixels
        self.out_img_pix = out_img_pix

        # Set lat/lon limits of output image
        self.out_img_ll = out_img_ll

        # List of variables to be read from primary file
        self.pvar = pvar

        # List of variables to be read from secondary file
        self.svar = svar
