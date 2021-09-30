;+
; Simple function that extracts primary and secondary ORAC netcdf
; output files from a given directory and builds a data cube from
; them. Note - be careful howmuch data you try to read in!
;-

function seviri_orac_mk_data_cube, dir
  files_p = file_search(dir,'*.primary.nc', count=nf)
;  files_s = file_search(dir,'*.secondary.nc', count=ns)
;  if np ne ns then stop, 'Different numbers of primary ('+ $
;                         strtrim(nf,2)+') and secondary ('+ $
;                         strtrim(ns,2)+') files found.'
  print, 'Loading ',strtrim(nf,2),' scenes'
  
  stat = read_ncdf(files_p[0], p)
;  stat = read_ncdf(files_s[0], s)

  stat = read_ncdf(files_p[1], p1)
;  stat = read_ncdf(files_s[1], s1)

  exceptions = ['D_ACROSS_TRACK','D_ALONG_TRACK','D_VIEWS']

  p = concat_struct_ex(p,p1,/new,except=exceptions)
;  s = concat_struct_ex(s,s1,/new,except=exceptions)
  nf = 10
  for i=2,nf-1 do begin
     stat = read_ncdf(files_p[i], p1)
;     stat = read_ncdf(files_s[i], s1)
     p = concat_struct_ex(p,p1,except=exceptions)
;     s = concat_struct_ex(s,s1,except=exceptions)
  endfor

  return, p
end
