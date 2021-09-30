#!/bin/bash

indir=/group_workspaces/cems2/rsgnceo2/Data/seviri_msg3/nrt_processing/quick_look/2017
thumbdir=/group_workspaces/cems/rsgnceo/public/nrt/nrt_part_seviri_msg3/quick_look/2017

# Produce the thumbnail images for linking on the online calendar viewer
for f in ${indir}/08/*/*.png; do
    echo $f
    fout=${f##*/}
    datetime=${fout:39:12}
    fout=${fout:57}
    mth=${datetime:4:2}
    day=${datetime:6:2}
    echo $thumbdir/${mth}/${day}/${datetime}${fout}
    mkdir -p --mode=0775 $thumbdir/${mth}/${day}
    if [[ $fout == *FC* ]]; then
	convert $f -crop 1318x732+213+51 -resize 360x200 \
	    ${thumbdir}/${mth}/${day}/${datetime}${fout}
#    else
#	convert $f -crop 1478x732+54+51 -resize 404x200 \
#	    ${thumbdir}/${mth}/${day}/${datetime}${fout}
    fi
done