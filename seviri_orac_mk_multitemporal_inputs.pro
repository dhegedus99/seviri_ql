;+
; This procedure takes the ORAC input files for a series of SEVIRI
; scenes (with common geographical bounds and pre-processor config)
; and creates a single set of multi-temporal input files.
;-
function update_substruct, in, data
  tags = tag_names(in)
  out = {}
  for i=0,n_elements(tags)-1 do begin
     if tags[i] ne 'DATA' then out = create_struct(out, tags[i], in.(i)) $
     else out = create_struct(out, 'DATA', data)
  endfor
  return, out
end

pro seviri_orac_mk_multitemporal_inputs, msiroot, outroot

  nscenes = n_elements(msiroot)
  fsuffix = ['.alb.nc', '.clf.nc', '.config.nc', '.geo.nc', '.loc.nc', $
             '.lsf.nc', '.lwrtm.nc', '.msi.nc', '.prtm.nc', '.swrtm.nc']
  nfiles  = n_elements(fsuffix)
; Check that all the required input files are available
  for s=0,nscenes-1 do begin
     check = file_test(msiroot[s]+fsuffix, /read, /regular)
     bad = where(check ne 1)
     if bad[0] ne -1 then begin
        for i=0,n_elements(bad)-1 do $
           message, /info, 'Missing file: '+msiroot[s]+fsuffix[bad[i]]
        message, 'Run aborted'
     endif
  endfor

; Generate time stamp for output files
  caldat, systime(/julian,/utc), MM,DD,YYYY,HH,NN,SS
  tstamp = string(YYYY,MM,DD,HH,NN,format='(i04,i02,i02,i02,i02)')

; Process each filetype in turn...
; Start with the config file, which sets up the channel index numbers,
; data array sizes, etc
  stat = read_ncdf(msiroot[0]+fsuffix[2], config, /attrib, /store_name)
; Define the output dimensions of the config structure first off
  oconfig = { D_NC_CONF: config.d_nc_conf * nscenes, $ ; Total number of channels
              D_NC_ALB : config.d_nc_alb * nscenes,  $ ; Number of sw channels
              D_NC_EMIS: config.d_nc_emis * nscenes, $ ; Number of lw channels
              D_NLAT_CONF: config.d_nlat_conf,       $ ; Preproc grid dimension
              D_NLON_CONF: config.d_nlon_conf,       $ ; Preproc grid dimension
              D_NLEVELS_CONF: config.d_nlevels_conf, $ ; Preproc grid dimension
              D_NX_CONF: config.d_nx_conf,           $ ; Data grid x-dimension
              D_NY_CONF: config.d_ny_conf,           $ ; Data grid y-dimension
              D_NATTRIBUTES: config.d_nattributes }
; Now step through the actual variables in the config structure,
; copying them and, when necessary, redefining variables as we go.
  tags1 = tag_names(config)
  for t=9,n_elements(tags1)-1 do begin
     switch tags1[t] of
        'FILE_NAME': begin
           ostr = update_substruct(config.(t), outroot+fsuffix[2])
           break
        end
        'PRODUCTION_TIME': begin
           ostr = update_substruct(config.(t), tstamp)
           break
        end
        'MSI_INSTR_CH_NUMBERS':
        'MSI_ABS_CH_WL':
        'MSI_CH_SWFLAG':
        'MSI_CH_LWFLAG':
        'ALB_ABS_CH_NUMBERS':
        'EMIS_ABS_CH_NUMBERS': begin
           ntmp = n_elements(config.(t).data)
           tmp = make_array(ntmp*nscenes, type=size(config.(t).data,/type))
           for i=0,nscenes-1 do tmp[i*ntmp:(i+1)*ntmp-1] = config.(t).data
           ostr = update_substruct(config.(t), tmp)
           break
        end
        'MSI_CH_VIEW': begin
           ntmp = n_elements(config.(t).data)
           tmp = make_array(ntmp*nscenes, type=size(config.(t).data,/type))
           for i=0,nscenes-1 do tmp[i*ntmp:(i+1)*ntmp-1] = i+1
           ostr = update_substruct(config.(t), tmp)
           break
        end
        ELSE: ostr = config.(t)
     endswitch
     oconfig = create_struct(oconfig,tags1[t], ostr)
  endfor
  
; Now move on to the ALB file



end
