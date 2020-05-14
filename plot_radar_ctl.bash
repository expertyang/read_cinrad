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

odigsiz=0.04

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
* vel color map
* 'set rgb 51   0 224 255'
* 'set rgb 52   0 128 255'
* 'set rgb 53  50   0 150'
* 'set rgb 54   0 251 144'
* 'set rgb 55   0 187   0'
* 'set rgb 56   0 143   0'
* 'set rgb 57 205 192 159'
* 'set rgb 58 118 118 118'
* 'set rgb 59 248 135   0'
* 'set rgb 60 255 207   0'
* 'set rgb 61 255 255   0'
* 'set rgb 62 174   0   0'
* 'set rgb 63 201 100   0'
* 'set rgb 64 255   0   0'
'set rgb 41  25  25 112'
'set rgb 42   0   0 139'
'set rgb 43   0   0 205'
'set rgb 44   0   0 238'
'set rgb 45   0   0 255'
'set rgb 46  67 110 238'
'set rgb 47  28 134 238'
'set rgb 48  30 144 255'
'set rgb 49   0 178 238'
'set rgb 50   0 191 255'
'set rgb 51   0 250 154'
'set rgb 52  34 139  34'
'set rgb 53 190 190 190'
'set rgb 54 238 154  73'
'set rgb 55 255 215   0'
'set rgb 56 255 255   0'
'set rgb 57 255 140 105'
'set rgb 58 255  99  71'
'set rgb 59 255  69   0'
'set rgb 60 255   0   0'
'set rgb 61 238   0   0'
'set rgb 62 205  55   0'
'set rgb 63 205   0   0'
'set rgb 64 205  51  51'
'set rgb 65 139  26  26'
'set rgb 66 112  25  25'

* ref color map
* 'set rgb 71  40 153 196'
* 'set rgb 72  78 120 177'
* 'set rgb 73   0 161 247'
* 'set rgb 74   0 237 237'
* 'set rgb 75   0 217   0'
* 'set rgb 76   0 145   0'
* 'set rgb 77 255 255   0'
* 'set rgb 78 231 193   0'
* 'set rgb 79 255 145   0'
* 'set rgb 80 255   0   0'
* 'set rgb 81 215   0   0'
* 'set rgb 82 193   0   0'
* 'set rgb 83 255   0 241'
* 'set rgb 84 151   0 181'
* 'set rgb 85 173 145 241'
* 'set rgb 86 255 255 255'
' set rgb 71   0 100   0'
' set rgb 72  85 107  47'
' set rgb 73  34 139  34'
' set rgb 74   0 205 102'
' set rgb 75  60 179 113'
' set rgb 76 102 205 170'
' set rgb 77 123 104 238'
' set rgb 78   0   0 255 '
' set rgb 79   0   0 139'
' set rgb 80 104  34 139'
' set rgb 81 139  58  98'
' set rgb 82 176  48  96'
' set rgb 83 139  34  82'
' set rgb 84 160  82  45'
' set rgb 85 218 165  32'
' set rgb 86 255 255   0'
' set rgb 87 233 150 122'
' set rgb 88 250 128 114'
' set rgb 89 238  44  44'
' set rgb 90 255  20 147'
' set rgb 91 173 145 241'
EOF
if [ $var == "ref" ]
then
cat>>plot_radar.gs<<EOF
'set cint 2'
* 'set clevs    0 5  10 15 20 25 30 35 40 45 50 55 60 65 70 75'
'set clevs   -15 -10 -6 -3 0   3  6  9 12 15 18 21 24 27 31 35 40 45 50 55 60'
'set ccols 0  71  72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91'
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
'set cint 2'
* 'set clevs  -40 -30 -20 -15 -10 -5  0  5  10 15 20 30 40'
'set clevs   -60 -45 -40 -35 -30 -25 -20 -15 -12 -9 -6 -3 -1  1  3  6  9 12 15 20 25 30 35 40 45 60'
'set ccols 41 41  42  43  44  45  46  47  48  49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 '
'set digsiz $vdigsiz'
'set cmark 3'

'd v'
'colmap.gs'
'draw title ID:$id $tim $cas vel LEV:$lev'
'printim $cas.$id.$tim.vel.$lev.gif white x1000 y800'

'c'
'set cint 2'
'set clevs    0 0.5  1  2  3  4  6  8  10 12 15 20'
*'set ccols 0 53  54 52 55 51 56 50 57 49 58 48'
'set digsiz $vdigsiz'
'set cmark 3'

'd w'
'colmap.gs'
'draw title ID:$id $tim $cas spw LEV:$lev'
'printim $cas.$id.$tim.spw.$lev.gif white x1000 y800'
EOF
elif [ $var == "obs" ]
then
cat>>plot_radar.gs<<EOF
'set cint 2'
* 'set clevs    0 5  10 15 20 25 30 35 40 45 50 55 60 65 70 75'
* 'set ccols 0 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86'
'set clevs   -15 -10 -6 -3 0   3  6  9 12 15 18 21 24 27 31 35 40 45 50 55 60'
'set ccols 0  71  72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91'
'set digsiz $odigsiz'
'set cmark 3'

'd z'
'colmap.gs'
'draw title ID:$id $tim $cas ref LEV:$lev'
'printim $cas.$id.$tim.ref.$lev.gif white x1000 y800'

'c'
'set cint 2'
* 'set clevs  -40 -30 -20 -15 -10 -5  0  5  10 15 20 30 40'
* 'set ccols 51  52  53  54  55 56 57 58 59 60 61 62 63 64 '
'set clevs   -60 -45 -40 -35 -30 -25 -20 -15 -12 -9 -6 -3 -1  1  3  6  9 12 15 20 25 30 35 40 45 60'
'set ccols 41 41  42  43  44  45  46  47  48  49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 '
'set digsiz $odigsiz'
'set cmark 3'

'd v'
'colmap.gs'
'draw title ID:$id $tim $cas vel LEV:$lev'
'printim $cas.$id.$tim.vel.$lev.gif white x1000 y800'

'c'
'set cint 2'
* 'set clevs  0 0.5 1  2  3  4  8 12 16 20 24'
* 'set ccols 0 57 56 55 54 58 59 60 61 62 63 64 '
* 'set clevs    0 0.5  1  2  3  4  6  9 12 15 20'
* 'set ccols 0 53  54 55 56 57 58 59 60 61 62 63'
* 'set ccols 0 53  54 52 55 51 56 50 57 49 58 48'
'set clevs    0 0.5  1  2  3  4  6  8  10 12 15 20'
'set digsiz $odigsiz'
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
