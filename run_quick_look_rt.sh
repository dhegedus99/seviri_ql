#!/bin/bash

idldir=${1}
export QL_INDIR=${2}
export QL_OUTDIR=${3}
thumbdir=${4}
datetime=${5}
if [[ $datetime != 'null' ]]; then
    export QL_DATETIME=${datetime}
fi
export QL_AEROSOL_QC=${6}
export QL_DIST2CLD=${7}
export QL_CESIUM=${8}
export QL_CLOBBER=${9}
export QL_AEROSOL_LANDSEA=${10}
module add idl
idl_exe=`which idl`

echo $QL_INDIR

${idl_exe} -IDL_CPU_TPOOL_NTHREADS 4 -rt=${idldir}'/wrap_seviri_mk_nrt_quicklooks.sav'

chmod -R go+rX ${QL_OUTDIR}

# Note that the following functionality has been moved to the separate
# make_quick_look_rt_thumbnail.sh script now.
# Produce the thumbnail images for linking on the online calendar viewer
# This assumes the input files have a filename format like this:
# NCEO-L2-CLOUD-AEROSOL-SEVIRI_ORAC_MSG3_201709190900_R4769_CTH.png
# and produces filenames like this:
# 201709190900_CTH.png
#for f in ${QL_OUTDIR}/*${datetime}*.png; do
#    fout=${f##*/}
#    fout=${fout:57}
    # False-colour images don't have the colour bar
#    if [[ $fout == *FC* ]]; then
#	convert $f -crop 1318x732+213+51 -resize 360x200 \
#	    $thumbdir/${datetime}${fout}
    # All other images do have the colour bar
#    else
#	convert $f -crop 1491x735+41+50 -resize 406x200 \
#	    $thumbdir/${datetime}${fout}
#    fi
#done