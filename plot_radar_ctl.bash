#!/bin/bash

CWD=`dirname $0`
echo plot_ref
if [ $# -lt 1 ]
then
   echo Usage: $0 ctlfile
   exit 0
fi

rdigsiz=0.01
vdigsiz=0.003
for ctlfile in $@
do
var=`basename $ctlfile|cut -d"." -f4`
lev=`basename $ctlfile|cut -d"." -f5`
tim=`basename $ctlfile|cut -d"." -f3`
 id=`basename $ctlfile|cut -d"." -f2`
cas=`basename $ctlfile|cut -d"." -f1`

rm -f plot_radar.gs
touch plot_radar.gs

line=`grep ^$id radar_info.txt`
lat=`echo $line|cut -d"," -f3`
lon=`echo $line|cut -d"," -f4`

echo $ctlfile
stnmap -i $ctlfile

cat>>plot_radar.gs<<EOF
'open $ctlfile'
lat1=$lat-4
lat2=$lat+4
lon1=$lon-5
lon2=$lon+5
'set lat ' lat1 ' ' lat2
'set lon ' lon1 ' ' lon2
'set mpdset hires cnhimap'
'set gxout shaded'
'set gxout stnmark'
EOF
if [ $var == "ref" ]
then
cat>>plot_radar.gs<<EOF
* color map
'set rgb 71  40 153 196'
'set rgb 72  78 120 177'
'set rgb 73   0 161 247'
'set rgb 74   0 237 237'
'set rgb 75   0 217   0'
'set rgb 76   0 145   0'
'set rgb 77 255 255   0'
'set rgb 78 231 193   0'
'set rgb 79 255 145   0'
'set rgb 80 255   0   0'
'set rgb 81 215   0   0'
'set rgb 82 193   0   0'
'set rgb 83 255   0 241'
'set rgb 84 151   0 181'
'set rgb 85 173 145 241'
'set rgb 86 255 255 255'
'set cint 2'
'set clevs    0 5  10 15 20 25 30 35 40 45 50 55 60 65 70 75'
'set ccols 0 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86'
'set digsiz $rdigsiz'
'set cmark 3'

'd z'
'colmap.gs'
'draw title ID:$id $tim $cas $var LEV:$lev'
'printim $cas.$id.$tim.$var.$lev.gif white x1000 y800'
EOF
elif [ $var == "vel" ]
then
cat >> plot_radar.gs<<EOF
'set rgb 51   0 224 255'
'set rgb 52   0 128 255'
'set rgb 53  50   0 150'
'set rgb 54   0 251 144'
'set rgb 55   0 187   0'
'set rgb 56   0 143   0'
'set rgb 57 205 192 159'
'set rgb 58 118 118 118'
'set rgb 59 248 135   0'
'set rgb 60 255 207   0'
'set rgb 61 255 255   0'
'set rgb 62 174   0   0'
'set rgb 63 201 100   0'
'set rgb 64 255   0   0'
'set cint 2'
'set clevs  -60 -50 -40 -30 -20 -10  0  10  20 30 40 50 60'
'set ccols 51  52  53  54  55 56 57 58 59 60 61 62 63 64 '
'set digsiz $vdigsiz'
'set cmark 3'

'd v'
'colmap.gs'
'draw title ID:$id $tim $cas vel LEV:$lev'
'printim $cas.$id.$tim.vel.$lev.gif white x1000 y800'

'c'
'set cint 2'
'set clevs  0 0.5 1  2  3  4  8 12 16 20 24'
'set ccols 0 57 56 55 54 58 59 60 61 62 63 64 '
'set digsiz $vdigsiz'
'set cmark 3'

'd w'
'colmap.gs'
'draw title ID:$id $tim $cas spw LEV:$lev'
'printim $cas.$id.$tim.spw.$lev.gif white x1000 y800'
EOF
elif [ $var == "zdr" ]
then
cat >> plot_radar.gs<<EOF
'set digsiz $rdigsiz'
'set cmark 3'
* 'set cint 0.1'
'set clevs -5 -4 -3 -2 -1 0 1 2 3 4 5'
'd zdr'
'colmap.gs'
'draw title ID:$id $tim $cas zdr LEV:$lev'
'printim $cas.$id.$tim.zdr.$lev.gif white x1000 y800'

'c'
'set cint 0.1'
'set digsiz $rdigsiz'
'set cmark 3'
'd cc'
'colmap.gs'
'draw title ID:$id $tim $cas cc LEV:$lev'
'printim $cas.$id.$tim.cc.$lev.gif white x1000 y800'

'c'
'set cint 30'
'set digsiz $rdigsiz'
'set cmark 3'
'd fdp'
'colmap.gs'
'draw title ID:$id $tim $cas fdp LEV:$lev'
'printim $cas.$id.$tim.fdp.$lev.gif white x1000 y800'

'c'
'set cint 0.5'
'set digsiz $rdigsiz'
'set cmark 3'
'd kdp'
'colmap.gs'
'draw title ID:$id $tim $cas kdp LEV:$lev'
'printim $cas.$id.$tim.kdp.$lev.gif white x1000 y800'
'c'

'set cint 10'
'set digsiz $rdigsiz'
'set cmark 3'
'd snr'
'colmap.gs'
'draw title ID:$id $tim $cas snr LEV:$lev'
'printim $cas.$id.$tim.snr.$lev.gif white x1000 y800'
EOF

fi
cat >> plot_radar.gs<<EOF
'quit'
EOF
grads -lbc "plot_radar.gs"
done
