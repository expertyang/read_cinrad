'query gxinfo'
rec2=sublin(result,2)
rec3=sublin(result,3)
rec4=sublin(result,4)
xsiz=subwrd(rec2,4)
ysiz=subwrd(rec2,6)
ylo =subwrd(rec4,4)
xhi =subwrd(rec3,6)
xd=xsiz-xhi
'query shades'
shdinfo =result
if(subwrd(shdinfo,1)='None')
  drawmap=0
else
  drawmap=1
endif
if(drawmap=1)
cnum=subwrd(shdinfo,5)
if(ylo<0.6 xd>1.5)
 xl=xhi+xd/2-0.4
 xr=xl+0.2
 xwid=0.2
 ywid=0.5
 if(ywid*cnum>ysiz*0.8)
  ywid=ysiz*0.8/cnum
 endif
 ymid=ysiz/2
 yb=ymid-ymid*cnum/2
 'set string 1 1 5'
 vert=1
else
 ymid=ylo/2
 yt=ymid+0.2
 yb=ymid
 xmid=xsiz/2
 xwid=0.8
 if(xwid*cnum>xsiz*0.8)
  xwid=xsiz*0.8/cnum
 endif
 xl=xmid-xwid*cnum/2
 'set string 1 tc 5'
 vert = 0
endif  
'set strsiz 0.12 0.13'
num=0
while(num<cnum)
 rec=sublin(shdinfo,num+2)
 col=subwrd(rec,1)
 hi=subwrd(rec,3)
 low=subwrd(rec,2)
 'set line 'col
 if(vert)
  yt=yb+ywid
 else
  xr=xl+xwid
 endif
 'draw recf 'xl' 'yb' 'xr' 'yt
 if(num<cnum-1)
   if(vert)
     'draw string '%(xr+0.05)%' 'yt' 'hi
   else
     'draw string 'xr' '%(yb-0.05)%' 'hi
   endif
 endif
 num=num+1
 if(vert);yb=yt;
 else;xl=xr;endif;
endwhile 
endif
'set string 1 bl 1'
