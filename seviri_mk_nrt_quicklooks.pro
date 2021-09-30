function get_seviri_map_param, lat, lon, centre=centre, rotation=rotation, $
                               limcoord=limcoord
  xmax = n_elements(lat[*,0])-1
  ymax = n_elements(lon[0,*])-1
  xmid = xmax/2 & ymid = ymax/2
  lat0 = lat[0,ymid]    & lon0 = lon[0,ymid]
  lat1 = lat[xmid,ymax] & lon1 = lon[xmid,ymax]
  lat2 = lat[xmax,ymid] & lon2 = lon[xmax,ymid]
  lat3 = lat[xmid,0]    & lon3 = lon[xmid,0]

  limcoord = [0,0,ymax,xmax]

  if lat0 lt -90. or ~finite(lat0) then begin
     latmp = reform(lat[0,*])
     latok = where(latmp ge -89.9, nlatok)
     if nlatok eq 0 then begin
        lat0 = lat[xmid,ymid]
        wh = where(lon ge -180)
        lon0 = min(lon[wh],wh2,/nan)
        wh3 = array_indices(lon,wh)
        limcoord[1] = wh3[wh2,0]
     endif else begin
        dlat = min(abs(latok - ymid), ii,/nan)
        lat0 = latmp[latok[ii]]
        lon0 = lon[0,latok[ii]]
     endelse
  endif
  if lat1 lt -90. or ~finite(lat1) then begin
     latmp = reform(lat[*,ymax])
     latok = where(latmp ge -89.9, nlatok)
     if nlatok eq 0 then begin
        lat1 = max(lat,wh,/nan)
        lon1 = lon[xmid,ymid]
        wh3 = array_indices(lat,wh)
        limcoord[2] = wh3[1]
     endif else begin
        dlat = min(abs(latok - xmid), ii, /nan)
        lat1 = latmp[latok[ii]]
        lon1 = lon[latok[ii],ymax]
     endelse
  endif
  if lat2 lt -90. or ~finite(lat2) then begin
     latmp = reform(lat[xmax,*])
     latok = where(latmp ge -89.9, nlatok)
     if nlatok eq 0 then begin
        lat2 = lat[xmid,ymid]
        lon2 = max(lon,wh,/nan)
        wh3 = array_indices(lon,wh)
        limcoord[3] = wh3[0]
     endif else begin
        dlat = min(abs(latok - ymid), ii, /nan)
        lat2 = latmp[latok[ii]]
        lon2 = lon[xmax,latok[ii]]
     endelse
  endif
  if lat3 lt -90. or ~finite(lat3) then begin
     latmp = reform(lat[*,0])
     latok = where(latmp ge -89.9, nlatok)
     if nlatok eq 0 then begin
        wh = where(lat ge -90)
        lat3 = min(lat[wh],wh2,/nan)
        lon3 = lon[xmid,ymid]
        wh3 = array_indices(lat,wh)
        limcoord[0] = wh3[wh2,1]
     endif else begin
        dlat = min(abs(latok - xmid), ii, /nan)
        lat3 = latmp[latok[ii]]
        lon3 = lon[latok[ii],0]
     endelse
  endif

  if ~keyword_set(centre) then centre = [0.0, 0.427]
  if ~keyword_set(rotation) then rotation = 0.0

  limit = [lat0,lon0, lat1,lon1, lat2,lon2, lat3,lon3]

  map_param={lat0: centre[0], lon0: centre[1], rotation: rotation, $
             sat_p: [6.6023,0.0,rotation], limit: limit}

  return, map_param
end


pro plot_fc, p, s, geostationary=geostationary, gamma=gamma, $
             title=title, image=img
; Create a false-colour image of the scene
  img = [[[rotate(s.REFLECTANCE_IN_CHANNEL_NO_3,5)]],$
         [[rotate(s.REFLECTANCE_IN_CHANNEL_NO_2,5)]],$
         [[rotate(s.REFLECTANCE_IN_CHANNEL_NO_1,5)]]]
; Scale to 0-1 range and set missing values to 0
  img = img/max(img)
  img[where(img lt 0.0)] = 0.0
; If this scene contains night pixels, we'll fill them with
; imagery taken from the 11 micron BT
  if max(p.solar_zenith_view_no1) ge 70 then begin
     sz = rotate(p.solar_zenith_view_no1,5)
     bt11 = rotate(s.BRIGHTNESS_TEMPERATURE_IN_CHANNEL_NO_9,5)
     ok = where(bt11 gt 0, complement=bad)
;    Invert the image so that cold clouds are bright, and scale to
;    0-1 range
     bt11[ok] = max(bt11[ok]) - bt11[ok]
     bt11[ok] = bt11[ok] / max(bt11)
     bt11[bad] = 0.0
;       Create a scaling image to merge the BT and FC images smoothly
     scale = sz
     scale = (sz - 70.0)/20.0
     scale[where(sz gt 90)] = 1.0
     scale[where(sz lt 70)] = 0.0
     for i=0,2 do img[*,*,i] = img[*,*,i] + scale * bt11
  endif
  if n_elements(gamma) gt 0 then img = img^gamma
  image_plot, img, /sideb, /hide_cbar, xticks=1, yticks=1, $
              xtickname=[' ',' '], ytickname=[' ',' '], $
              title=title, geostationary=geostationary, $
              position=[0.08,0.03,0.97,0.92], plot_pos=pos, $
              barsize=0.015, baroffset=0.5, /hires, foreground=0, $
              background='ffffff'x
;  print, pos
end

pro plot_data, data, filter, range, log=log, title=title, $
               keytitle=keytitle, keyticks=keyticks, $
               geostationary=geostationary
  dat = rotate(data,5)
  min = 0.0
  dat[where(dat) gt 0.99*range[1]] = 0.99*range[1]
  dat[where(dat gt min(dat) and dat lt range[0])] = range[0]
  image_plot, dat, /sidebar, /nocolourt, xticks=1, yticks=1, log=log, $
              xtickname=[' ',' '], ytickname=[' ',' '], range=range, $
              title=title, keytitle=keytitle, keytickname=keyticks, $
              filt=rotate(filter,5), geostationary=geostationary, $
              barsize=0.015, baroffset=0.5, /hires, foreground=0, $
              background='ffffff'x, position=[0.08,0.03,0.97,0.92]
              
end

pro wrap_seviri_mk_nrt_quicklooks, indir=indir, outdir=outdir, $
                                   datetime=datetime, clobber=clobber, $
                                   aerosol_qc=aerosol_qc, $
                                   dist2cld=dist2cld, notrap=notrap, $
                                   aerosol_landsea=aerosol_landsea, $
                                   cesium=cesium

  if ~keyword_set(notrap) then begin
     error=0
     catch, error
     if error ne 0 then begin
        help, /last_message, output=error_message
        catch, /cancel ; prevent looping 
        print, error_message
        return
     endif
  endif

  if ~keyword_set(indir) then indir = getenv('QL_INDIR')
  if ~keyword_set(outdir) then outdir = getenv('QL_OUTDIR')
  if ~keyword_set(datetime) and getenv('QL_DATETIME') ne '' $
  then datetime = getenv('QL_DATETIME')
  if n_elements(clobber) eq 0 then clobber = getenv('QL_CLOBBER')
  if ~keyword_set(aerosol_qc) then aerosol_qc = getenv('QL_AEROSOL_QC')
  if ~keyword_set(dist2cld) then dist2cld = getenv('QL_DIST2CLD')
  if ~keyword_set(aerosol_landsea) then $
     aerosol_landsea = getenv('QL_AEROSOL_LANDSEA')
  if n_elements(cesium) eq 0 then cesium = getenv('QL_CESIUM')

  if ~file_test(indir, /directory, /read) then $
     message,"Specified input dir can't be read: "+indir
  if outdir eq '' then outdir=0
  if clobber eq '' or clobber eq '0' then clobber=0 else clobber=1

  seviri_mk_nrt_quicklooks, indir, outdir=outdir, datetime=datetime, $
                            clobber=clobber, aerosol_qc=aerosol_qc, $
                            dist2cld=dist2cld, cesium=cesium, $
                            aerosol_landsea=aerosol_landsea

  skip_run:
end


pro seviri_mk_nrt_quicklooks, indir, outdir=outdir, datetime=datetime, $
                              p=p, s=s, clobber=clobber, aerosol_qc=aerosol_qc, $
                              dist2cld=dist2cld, cesium=cesium, $
                              aerosol_landsea=aerosol_landsea

; Define variables to read in from the data files
  pvar = ['lat', 'lon', 'cldmask', 'illum', 'solar_zenith_view_no1', $
          'cth', 'cot', 'aot550', 'aer', 'qcflag', 'lsflag']
  svar = ['reflectance_in_channel_no_1', 'reflectance_in_channel_no_2', $
          'reflectance_in_channel_no_3', $
          'brightness_temperature_in_channel_no_9']


; if the outdir hasn't been specified through the keyword, use
; the default. This default directory will be appended with yr/mt/dy
; sub-dirs; a specified outdir will be used as is.
  if n_elements(outdir) eq 0 then begin
     auto_out = 1
;     outdir = '/group_workspaces/cems2/rsgnceo2/scratch/seviri_nrt/quick_look'
     outdir = '/group_workspaces/cems2/rsgnceo2/Data/seviri_msg3/nrt_processing/quick_look'
  endif else begin
     auto_out = 0
  endelse

  if n_elements(aerosol_landsea) eq 0 then aerosol_landsea=0

; Find the primary and secondary output files
  if keyword_set(datetime) then $
     fsearch='*SEVIRI_ORAC_MSG*'+datetime+'*.primary.nc' $
  else $
     fsearch='*SEVIRI_ORAC_MSG*'+datetime+'*.primary.nc'
  fp = file_search(indir,fsearch, count=np)
;  fs = file_search(indir,'*SEVIRI_ORAC_MSG*.secondary.nc', count=ns)
  
; Check to make sure we've actually found a primary file...
  if np eq 0 then message, 'Primary file not found: '+fsearch
; Check the file is okish
  primary_ok = is_ncdf(fp[0])
  if primary_ok ne 1 then message, 'Primary file seems to be bad: ',fp[0]

; Set up colour table
;  set_plot,'ps'
;  device,/color
  print, 'SEVIRI_MK_NRT_QUICKLOOKS: Input and output directories'
  print, indir
  print, outdir
  if keyword_set(cesium) then print, cesium
  print, 'Attempting to read from file: '+fp[0]
  
  stat = read_ncdf(fp[0],p,/no_data)
  set_plot,'z'
  device,set_pixel_depth=24,decomposed=1, $
         set_resolution=[p.d_across_track/0.86, p.d_along_track/0.788]
  !p.multi=0
  erase
;  set_plot,'x'
  !p.font=-1
  !p.charsize=1.5

;  window, 1, xsize=p.d_across_track/0.86, ysize=p.d_along_track/0.86

  cesgrid=0 ; Used for flagging when the Cesium image output grid has
            ; been created

  for f = 0,np-1 do begin
;  for f=4,4 do begin
;    Note that we need to reload the default colour table for each
;    file, as we replace it with something else for the AOD Cesium image
     colour_mycarta, /cube1, /white
;     !p.multi=[0,2,2]
;    Extract date information
     date = strmid(file_basename(fp[f]), 39, 8)
     time = strmid(file_basename(fp[f]), 47, 4)
;    Search for the secondary file...
     fs = file_search(indir,'*SEVIRI_ORAC_MSG*'+date+time+'*.secondary.nc')
     if fs eq '' then message,/info, 'Warning: Secondary data file not found'

     yr = strmid(date,0,4)
     mt = strmid(date,4,2)
     dy = strmid(date,6,2)
;    Define output filename
     fo = file_basename(fp[f], '.primary.nc')
     if auto_out then begin
        file_mkdir, outdir+'/'+yr+'/'+mt+'/'+dy
        fo = outdir+'/'+yr+'/'+mt+'/'+dy+'/'+fo 
     endif else begin
        file_mkdir, outdir
        fo = outdir+'/'+fo
     endelse
;    Do likewise for Cesium imagery, if required
     if n_elements(cesium) ne 0 then begin
        if auto_out then begin
           file_mkdir, cesium+'/'+yr+'/'+mt+'/'+dy
           cfo = cesium+'/'+yr+'/'+mt+'/'+dy+'/'+file_basename(fo)
        endif else begin
           file_mkdir, cesium
           cfo = cesium+'/'+file_basename(fo)
        endelse
     endif

;    If all the plots are present and we're not clobbering,
;    then skip this file right here
     if ~keyword_set(clobber) and file_test(fo+'_FC.png') and $
        file_test(fo+'_CTH.png') and file_test(fo+'_COT.png') and $
        file_test(fo+'_AOD.png') then begin
        if n_elements(cesium) ne 0 then begin
           if file_test(cfo+'_FC.png') and file_test(cfo+'_CTH.png') and $
              file_test(cfo+'_COT.png') and file_test(cfo+'_AOD.png') then $
                 goto, skip_scene
        endif else goto, skip_scene
     endif

     print, 'Plotting data from: ',file_basename(fp[f], '.primary.nc')

;    Read data
     stat = read_ncdf(fp[f], p, variable_l=pvars)
     if fs ne '' then stat = read_ncdf(fs, s, variable_l=svars)

;    If required, apply additional quality control to the products
;    before plotting
     if keyword_set(aerosol_qc) then $
        p = seviri_aerosol_additional_cloud_tests(p, qcsel=aerosol_qc, $
                                                  dist2cld=dist2cld)

;    Extract the qcflag variable and make sure that it is being stored
;    as a long-int (some ORAC output files define a scale_factor and
;    add_offset for this variable, which causes read_ncdf to convert
;    the variable to a floating point type)
     qcflag = long(p.qcflag)

     ptags = tag_names(p)

;    Define image coordinates so image_plot can overplot continents
     map = get_seviri_map_param(rotate(p.lat,5), rotate(p.lon,5))

;     device, filename=fo
;    Produce a false-colour image - always try to regenerate, as the
;    fcimg is needed for the Cesium plotting later
     title_str = 'MSG3-SEVIRI '+yr+'-'+mt+'-'+dy+' '+time
     if fs ne '' then begin
        plot_fc, p, s, geostationary=map, gamma=0.6, title=title_str, $
                 image=fcimg
        if !error_state.code eq 0 then begin
           print, 'Writing false colour QL to :',fo+'_FC.png'
           write_png, fo+'_FC.png', tvrd(/true)
        endif else begin
           print, 'Non-zero error code: ', !error_state.code
           !error_state.code = 0
        endelse
;    If we can't produce the false colour image, create an
;    empty image array, so that the Cesium image production
;    doesn't crash
     endif else fcimg = fltarr(p.d_across_track, p.d_along_track, 3)

;    Now plot some cloud variables
;    CTH
     filt = p.cldmask eq 1 and p.cth gt 0.5
     if (keyword_set(clobber) or ~file_test(fo+'_CTH.png')) $
        and (where(filt))[0] ge 0 then begin
        thistitle = title_str+' CTH'
        plot_data, p.cth, filt, [0.0,15.0], geostationary=map, $
                   title=thistitle, keytitle='Cloud-top height (km)', $
                   keyticks=['0','3','6','9','12','>15']
        if !error_state.code eq 0 then $
           write_png, fo+'_CTH.png', tvrd(/true) $;/true) $
        else !error_state.code = 0
     endif

     filt= p.cldmask eq 1 and p.cth gt 0.0 and p.illum ne 2
     if (keyword_set(clobber) or ~file_test(fo+'_dCTH.png')) $
        and (where(filt))[0] ge 0 then begin
        thistitle = title_str+' del(CTH)'
        plot_data, p.cth_uncertainty, filt, [0.1,10], /log, $
                   geostationary=map, title=thistitle, $
                   keyticks=['0.1', '1.0', '>10.'], $
                   keytitle='Uncertainty in Cloud-top height (km)'
        if !error_state.code eq 0 then $
           write_png, fo+'_dCTH.png', tvrd(/true) $;/true) $
        else !error_state.code = 0
     endif

;    alog(COD)
     filt = p.cldmask eq 1 and p.illum eq 1 and p.cot gt 0.0
     if (keyword_set(clobber) or ~file_test(fo+'_COT.png')) $
        and (where(filt))[0] ge 0 then begin
        thistitle = title_str+' COT'
        plot_data, p.cot, filt, [0.1,100.0], /log, geostationary=map, $
                   title=thistitle, keytitle='Cloud optical depth (550 nm)', $
                   keyticks=['0.1','1','10','>100']
        if !error_state.code eq 0 then $
           write_png, fo+'_COT.png', tvrd(/true) $
        else !error_state.code = 0
     endif

;    And AOD
;    Note we're attempting to use the additional QC flag values I've
;    added for aerosol QC:
;    cld_edge 256
;    aod_open 1024
;    ref_open 2048 => total: 3328
     if (keyword_set(clobber) or ~file_test(fo+'_AOD.png')) and $
        (where(ptags eq 'AOT550'))[0] ge 0 then begin
        filt = p.cldmask eq 0 and p.illum eq 1 and p.aot550 gt 0.0 and $
               (qcflag and 3328) eq 0
        if aerosol_landsea eq 1 then filt = filt and p.lsflag eq 0 $
        else if aerosol_landsea eq 2 then filt = filt and p.lsflag eq 1
        if (where(filt))[0] ge 0 then begin
           thistitle = title_str+' AOD'
           plot_data, p.aot550, filt, [0.0,1.2], geostationary=map, $
                      title=thistitle, $
                      keytitle='Aerosol optical depth (550 nm)', $
                      keyticks=['0.0','0.2','0.4','0.6','0.8','1.0','>1.2']
           if !error_state.code eq 0 then $
              write_png, fo+'_AOD.png', tvrd(/true) $
           else !error_state.code = 0
        endif
     endif
;********************************************************************;
;    Have we been asked to produce images for use with the Cesium    ;
;    viewer?                                                         ;
;********************************************************************;
;    Note that we check for the existance of the Secondary file, as a
;    check for the existence of a real false colour image array. For
;    some reason, the empty fake false colour array is causing the
;    generation of the Cesium imagery to hang indefinitely...
     if n_elements(cesium) gt 0 and fs ne '' then begin
        print,'Producing Cesium imagery...'

        time_setint = 0d0
        time_where_img = 0d0
        time_where_aod = 0d0
        time_where_cth = 0d0
        time_where_cot = 0d0
        time_mean_img = 0d0
        time_mean_aod = 0d0
        time_mean_cth = 0d0
        time_mean_cot = 0d0
        
;       Only generate the output grid once per run
        if ~cesgrid then begin
;          Define the output image-grid. We use the nadir satellite
;          resolution (3 km, ~ 0.03 degree, for SEVIRI)
           res = 0.06
           ok = where(p.lat ne -999.0)
           latrng = [min(p.lat[ok]), max(p.lat[ok])]
           lonrng = [min(p.lon[ok]), max(p.lon[ok])]
           latlon_arrays, lat, lon, 180./res, 360./res, /nglobal, $
                          /b, latrange=latrng, lonrange=lonrng
           lathist = histogram(p.lat, min=latrng[0], binsize=res, $
                               nbins=n_elements(lat)-1, reverse=lati, $
                               locat=latbin)
           lonhist = histogram(p.lon, min=lonrng[0], binsize=res, $
                               nbins=n_elements(lon)-1, reverse=loni, $
                               locat=lonbin)
           nlat = n_elements(lat)-1
           nlon = n_elements(lon)-1
           cesgrid = 1
;           print, min(latbin), max(latbin), res
;           print, min(lonbin), max(lonbin), res
        endif

;       Define output arrays of the rebinning process
        nin  = fltarr(nlon,nlat)
        naer = fltarr(nlon,nlat)
        aod  = fltarr(nlon,nlat)
        ncth = fltarr(nlon,nlat)
        cth  = fltarr(nlon,nlat)
        ncot = fltarr(nlon,nlat)
        cot  = fltarr(nlon,nlat)
        nimg = fltarr(nlon,nlat)
        img  = fltarr(3,nlon,nlat)
        
;       Go through each lat-lon bin, using the reverse-indices from
;       the histogram calls above to populate each in turn
        print,'   Beginning data binning:'
        lastmod = 0.0
        tstart = systime(1)
        for j=0,nlat-1 do if lati[j] ne lati[j+1] then begin
           percent = float(j+1)/float(nlat)
           thismod = percent mod 0.1
           if thismod lt lastmod then begin
              print, string(100*percent, format='(f6.2)')+'% complete after '+$
                     string(systime(1)-tstart, format='(f8.2)')+'s'
              lastmod = 0.0
           endif else lastmod = thismod
           for i=0,nlon-1 do if loni[i] ne loni[i+1] then begin
              t0 = systime(1)
              inbin = setintersection(lati[lati[j]:lati[j+1]-1], $
                                      loni[loni[i]:loni[i+1]-1] )
              time_setint = time_setint + systime(1) - t0

;             inbin holds the indices of the unfiltered image pixels
;             which lie in the current lat-lon bin.
              if inbin[0] ne -1 then begin
                 nin[i,j] = n_elements(inbin)
                 t0 = systime(1)
;                Now we filter for valid reflectances
                 oki = where(s.reflectance_in_channel_no_1[inbin] ge 0, nimg)
                 time_where_img = time_where_img + systime(1) - t0

                 t0 = systime(1)
;                We use the previously generated false colour image,
;                which has the nice transition from colour to thermal
;                for night scenes.
                 for k=0,2 do begin
;                   Note that we have to un-rotate the image for the
;                   indexing to work
                    tmp = rotate(fcimg[*,*,k],5)
                    img[k,i,j] = mean(tmp[inbin[oki]])
                 endfor
;                 img[0,i,j] = mean(s.reflectance_in_channel_no_3[inbin[oki]])
;                 img[1,i,j] = mean(s.reflectance_in_channel_no_2[inbin[oki]])
;                 img[2,i,j] = mean(s.reflectance_in_channel_no_1[inbin[oki]])
                 time_mean_img = time_mean_img + systime(1) - t0
                 
                 t0 = systime(1)
;                Now repeat the procedure for the AOD. Here we also
;                keep track of which pixels lie in the lat-lon bin,
;                but don't have valid data. We need to
;                differentiate between lat-lon bins which have no data
;                and those which have only missing values so the
;                filling will work properly later
                 if aerosol_landsea eq 1 then begin
                    oka = where(p.cldmask[inbin] eq 0 and $
                                p.illum[inbin] eq 1 and $
                                p.aot550[inbin] gt 0.0 and $
                                (qcflag[inbin] and 3328) eq 0 and $
                                p.lsflag[inbin] eq 0, noka, $
                                complement=miss, ncomp=nmiss)
                 endif else begin
                    oka = where(p.cldmask[inbin] eq 0 and $
                                p.illum[inbin] eq 1 and $
                                p.aot550[inbin] gt 0.0 and $
                                (qcflag[inbin] and 3328) eq 0, noka, $
                                complement=miss, ncomp=nmiss)
                 endelse
                 time_where_aod = time_where_aod + systime(1) - t0
                 if noka ne 0 then begin
                    naer[i,j] = noka
                    t0 = systime(1)
                    aod[i,j] = mean(p.aot550[inbin[oka]])
                    time_mean_aod = time_mean_aod + systime(1) - t0
                 endif
;                If the bin has only missing data, then use a fill value
                 if nmiss ne 0 and noka eq 0 then begin
                    aod[i,j] = -999.0
                 endif

                 t0 = systime(1)
;                Cloud-top height is treated exactly the same as AOD
                 okct = where(p.cldmask[inbin] eq 1 and $
                              p.cth[inbin] gt 0.5, nokct, $
                              complement=miss, ncomp=nmiss)
                 time_where_cth = time_where_cth + systime(1) - t0
                 if nokct ne 0 then begin
                    ncth[i,j] = nokct
                    t0 = systime(1)
                    cth[i,j] = mean(p.cth[inbin[okct]])
                    time_mean_cth = time_mean_cth + systime(1) - t0
                 endif
                 if nmiss ne 0 and noka eq 0 then begin
                    cth[i,j] = -999.0
                 endif

                 t0 = systime(1)
;                Finally, cloud-optical depth
                 okcd = where(p.cldmask[inbin] eq 1 and $
                              p.illum[inbin] eq 1 and $
                              p.cot[inbin] gt 0.0, nokcd, $
                              complement=miss, ncomp=nmiss)
                 time_where_cot = time_where_cot + systime(1) - t0
                 if nokcd ne 0 then begin
                    ncot = nokcd
                    t0 = systime(1)
                    cot[i,j] = alog10(mean(p.cot[inbin[okcd]]))
                    time_mean_cot = time_mean_cot + systime(1) - t0
                 endif
                 if nmiss ne 0 and noka eq 0 then begin
                    cot[i,j] = -999.0
                 endif
              endif
           endif
        endif
        
;       In order to speed things up, we generate and save a mask which
;       limits fill_grid to operate only on pixels which lie within
;       the SEVIRI field of view (as defined by the false colour
;       image)
;       NOTE: If the resolution of the output imagery, or the size of
;       the SEVIRI sub-image processed, changes, this file should be
;       deleted (and thus regenerated for the new dimensions)
        have_save=1
        if n_elements(skip) eq 0 then begin
           save_file=getenv('HOME')+'/tmp/seviri_fill_mask.sav'
           have_save = file_test(save_file, /regular, /read)
           if have_save then restore, save_file $
           else skip = make_array([nlon,nlat], /byte, value=1)
        endif
              
        t0 = systime(1)
;       Run fill_grid on the false colour imagery, filling in the
;       missing pixels as we move away from the sub-satellite point
;       (essentially degrading the lat-lon resolution to something
;       resembling what SEVIRI actually sees)
        miss=min(img)
        for i=0,2 do begin
           tmp = reform(img[i,*,*])
           fill_grid, tmp, 0.0, skip, max=1
           img[i,*,*] = tmp
;          Here is where we save the skip mask, if it wasn't
;          read in earlier
           if i eq 0 and ~have_save then begin
              skip[where(tmp eq 0)] = 0
              save, skip, filename=save_file
           endif
        endfor
        miss = min(aod)
;       Now repeat the procedure on the retrieval products
        aods = aod
        fill_grid, aod, 0.0, skip, max=1
        fill_grid, cth, 0.0, skip, max=1
        fill_grid, cot, 0.0, skip, max=1
        time_fill = systime(1) - t0

;       Now create the output images. We set the z-buffer size to
;       match our images, write the data into the buffer, and then
;       write this data into a PNG file.
        device, set_resolution=[nlon,nlat]
        tvscale, img
        write_png, cfo+'_cs_FC.png', tvrd(/true)
        cth[where(cth eq -999.0)] = 0
        cth[where(cth gt 15.0)] = 14.9
        cth[0] = 15.0
        tvscale, cth
        write_png, cfo+'_cs_CTH.png', tvrd(/true)
        cot[where(cot lt -1.0 or cot eq 0.0)] = -1.01
        cot[where(cot gt 2.0)] = 1.95
        cot[0] = 2.0
        tvscale, cot
        write_png, cfo+'_cs_COT.png', tvrd(/true)
;       Use a different colour scheme for AOD, so plotting both cloud
;       and aerosol provides a clear differentiation
        colour_mycarta, red=  2.55*[ 30,  0,100,100], $
                        green=2.55*[  0,  0, 80,  0], $
                        blue= 2.55*[ 30, 80, 30,  0], /custom
        aod[where(aod eq -999.0)] = 0
        aod[where(aod gt 1.2)] = 1.19
        aod[0] = 1.2
        tvscale, aod
        write_png, cfo+'_cs_AOD.png', tvrd(/true)

;       Now use external command, namely Imagemagick's mogrify
;       command, to add transparancy to these images
;       (IDL only supports transparancy for 8-bit images, not 24-bit
;       like ours)
        spawn, 'mogrify -transparent "rgb(0,0,0)" '+cfo+'_cs_FC.png' 
        spawn, 'mogrify -transparent "rgb(255,255,255)" '+cfo+'_cs_FC.png'
        spawn, 'mogrify -transparent "rgb(0,0,0)" -fuzz 1% '+cfo+'_cs_AOD.png'
        spawn, 'mogrify -transparent "rgb(255,255,255)" -fuzz 1% '+cfo+'_cs_CTH.png'
        spawn, 'mogrify -transparent "rgb(255,255,255)" -fuzz 1% '+cfo+'_cs_COT.png'

     endif
;********************************************************************;
;    End of Cesium imagery section                                   ;
;********************************************************************;
     skip_scene:
  endfor
end
