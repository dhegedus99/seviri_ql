function seviri_aerosol_additional_cloud_tests, x, qcselection=qcselection, $
   dist2cld=dist2cld

; Define bit-mask values for each QC test
  conv_bit     =   1
  cost_bit     =   2
  cld_edge_bit =   4
  sd_aod_bit   =   8
  aod_open_bit =  16
  ref_open_bit =  32
  ice_bit      =  64
  ang_open_bit = 128
  cld_dist_bit = 256  

; Set a default qcselection mask to apply all tests
  if ~keyword_set(qcselection) then qcselection=255

; Set default behaviour for distance-to-cloud calculation:
;  If keyword isn't set, set the value to 0
;  If the keyword is set to a positive value, use this as the QC
;  minimum distance to cloud allowed, and switch on the test.
  if n_elements(dist2cld) eq 0 then dist2cld = 0
  if dist2cld gt 0 and qcselection lt 256 then qcselection = qcselection + 256

; Default values for various thresholds used to control the output
  if ~keyword_set(minclr) then minclr = 5

; Apply some basic quality control to the input data, and set qual
; flag accordingly
; Extract and combine nadir and forward cloud masks.
  cld = x.cldmask[*,*]
; This routine masks out pixels which either have cloud in more than
; 1/3 of their surrounding pixels, or have (when combined with their
; neighbouring pixesl) have an AOD standard deviation gt 0.1
  if (qcselection and 4) gt 0 or (qcselection and 8) then begin
     print, 'Calculating cloud edge QC parameters'
     fmi_cloud_qc_1km, cld, x.aot550, edge, sd_aot
  endif
; The 1 km data seems to have a lot of residual cloud "speckle", which
; looks like so called "salt noise" in the AOD imagery... deal with it
; using the morphological opening transformation to highlight and
; remove this noise
  if (qcselection and 16) gt 0 then begin
     print, 'Calculating AOD opening'
;    Use a 5x5 pixel kernel
     kern = make_array(5,5,/uint,value=1)
;    Convert AOD to a unsigned integer value, with 0 as the fill value
;    and AOD = 1 set to be 255.
     intaod = x.aot550
     intaod[where(intaod lt 0)] = 0
     intaod = uint(255*intaod)
;    Do the opening transform and check for discrepencies between the
;    two images
     opnaod = dilate(erode(intaod,kern,/gray,/preserve), $
                     kern,/gray,/preserve)
     salt = (intaod - opnaod) ge 80
  endif

; Repeat the openning test, this time on effective radius
; Convert to a unsigned integer value, with 0 as the fill value
; and scale by 255.
  if (qcselection and 32) gt 0 then begin
     print, 'Calculating Reff opening'
;    Use a 5x5 pixel kernel
     kern = make_array(5,5,/uint,value=1)
     intefr = x.aer
     intefr[where(intefr lt 0)] = 0
     intefr = uint(255*intefr)
;    Do the opening transform and check for discrepencies between the
;    two images
     opnefr = dilate(erode(intefr,kern,/gray,/preserve), $
                     kern,/gray,/preserve)
     salter = (intefr - opnefr) ge 300
  endif

; Check for ice/snow contamination 
;  ice = (data.rho_dd_in_channel_no_1 - data.rho_dd_in_channel_no_4) ge 0.07


; Set the quality mask. 0 = pixel is ok.
  qual = make_array(dimension=size(x.lat,/dimen), /byte, value=0)
  if (qcselection and  1) gt 0 then qual = qual + (x.niter ge 25)*conv_bit
  if (qcselection and  2) gt 0 then qual = qual + (x.costjm ge 3)*cost_bit
  if (qcselection and  4) gt 0 then qual = qual + (1-edge) * cld_edge_bit
  if (qcselection and  8) gt 0 then qual = qual + (1-sd_aot) * sd_aod_bit
  if (qcselection and 16) gt 0 then qual = qual + salt * aod_open_bit
  if (qcselection and 32) gt 0 then qual = qual + salter * ref_open_bit
;  if (qcselection and 64) gt 0 then qual = qual + ice * ice_bit

; Now we have our completely determined "cloud mask", we can calculate
; the distance to cloud parameter
  if dist2cld ne 0 then begin
;     cldbits = sd_aod_bit + aod_open_bit + ref_open_bit + ang_open_bit
     cldbits = 0
     clddist = orac_1km_dist_to_cloud(cld, qual, bitmask=cldbits, $
                                      window=2*dist2cld+1)
     if (qcselection and cld_dist_bit) gt 0 then $
        qual = qual + (clddist le dist2cld and clddist gt 0.0) * cld_dist_bit
  endif else clddist = make_array(x.d_across_track, x.d_along_track, $
                                  /float, value=-999.0)

; Modify the qcflag in the product accordingly 
  x.qcflag = x.qcflag + 64 * qual
  
  return, x

end
