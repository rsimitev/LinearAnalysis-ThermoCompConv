; visualization of stream function in equatorial plane
; input: datafiles from galerkin program, converted by 'lara'
; the program is started under idl by entering:   .run lara
; input: the files idl.z, idl.x, idl.y, created by lara
; ----- pro make_levels ---------------------------------

pro make_levels, ma,mi,v,vl,vls,levelnum

  delta=(ma-mi)/(levelnum-1)
  v=findgen(levelnum)*delta + ceil(mi/delta)*delta
  for i=0,levelnum-1 do begin
  
;    if(v(i) lt 0.0) then vls(i)=2
;    if(v(i) eq 0.0) then vls(i)=1
;    if(v(i) gt 0.0) then vls(i)=0

  ;---Radostin
  ; black&white
  ;
  ;  if(v(i) lt 0.0) then vls(i)=2
  ;  if(v(i) eq 0.0) then vls(i)=1
  ;  if(v(i) gt 0.0) then vls(i)=0
 
  ;  color
  ; 
   if(v(i) lt 0.0) then vls(i)=4
   if(v(i) eq 0.0) then vls(i)=3
   if(v(i) gt 0.0) then vls(i)=2

  endfor
end

; ----- MAIN PROGRAM ---------------------------------------

runname=strarr(30)
runname=''

;------ dimensions are fixed for idl.z,... from lara.f:
nmx=65      ; ncos   +1    always check the number of points in 
nmy=128     ; nangle +1     the output file and adjust.

;------ device:
ps=1
;------ arrange scales so that circles actually are circular on paper 
if (ps eq 1) then posv=[0.3,0.1,0.85,0.9]
if (ps ne 1) then posv=[0.2,0.1,0.9,0.9]
if(ps eq 1) then begin
set_plot, 'PS'
device,filename='plot.ps' ,/color  ;,XSIZE=8.0, YSIZE=12.0
end

;------ dont't plot the axes:
  !x.style=4
  !y.style=4

ncos=nmx
nangle=nmy
z=dblarr(ncos,nangle)
y=dblarr(ncos,nangle)
x=dblarr(ncos,nangle)
ncos=ncos-1
nangle=nangle-1

;------ Equatorial plot (cut_type=2) or Spherical plot (cut_type=3)
cut_type=2
 

; =====-----------------=====
; ====-  equatorial cut  -====
; =====-----------------=====

if(cut_type eq 2) then begin     

;------ read data
print, 'reading data from idl.z, idl,x, idl.y ...'
openr, 1, 'idl.z'
for i=0,ncos do for j=0,nangle do begin 
readf,1,tz
z(i,j)=tz
end
close,1

openr,2,'idl.y'
for i=0,ncos do for j=0,nangle do begin 
readf,2,ty
y(i,j)=ty
end
close,2

openr,3,'idl.x'
for i=0,ncos do for j=0,nangle do begin
readf,3,tx
x(i,j)=tx
end
close,3

print, '...done.'

;------ contour,z,x,y,nlevels=8,/follow,pos=posv
; arrays for the contour defs:
levelnum=19
v=dblarr(levelnum)
vl=intarr(levelnum)
vls=intarr(levelnum)
ma=max(z)
mi=min(z)
print, ma, mi
;stop
make_levels, ma,mi,v,vl,vls,levelnum

  red =   [0,1,1,0,0,1] 
  green = [0,1,0,1,0,1] 
  blue =  [0,1,0,0,1,0] 
  tvlct, 255*red, 255*green, 255*blue 
 
 ;------Radostin 
 ; black and white plots 
 ; 
 ; contour, z,x,y, /noerase, color=plot_color, clip=clipvec, levels=v, c_linestyle=vls, c_colors=0, c_labels=vl, /follow, pos=posv 
 
 ; color plots  
 ; 
   contour, z,x,y, /noerase, color=plot_color, clip=clipvec, levels=v, c_linestyle=0, c_colors=vls, c_labels=vl, /follow, pos=posv 

 ; contour, z,x,y, /noerase, color=plot_color, clip=clipvec, levels=v, c_linestyle=vls, c_colors=0, c_labels=vl, /follow, pos=posv


;------ draw lines at inner and outer radius
for i=0,ncos do for j=0,nangle do z(i,j)=0
i=0
for j=0,nangle do z(i,j)=1
i=ncos
for j=0,nangle do z(i,j)=1
;j=nangle
;for i=0,ncos do z(i,j)=1
;j=0
;for i=0,ncos do z(i,j)=1
; contour, z,x,y, /noerase, color=plot_color, clip=clipvec, levels=v, c_linestyle=vls, c_colors=0, c_labels=vl, /follow, pos=posv,/overplot
contour,z,x,y ,levels=1,c_thick= [1,1],/noerase,pos=posv,c_charsize=0,/overplot
if(ps eq 1) then begin
   device,/close
   set_plot,'X'
end
ncos=ncos+1

endif


; =====-----------------=====
; ====-  spherical cut  -====
; =====-----------------=====

if(cut_type eq 3) then begin


nangle=nmx-1
lmax=nmy-1
z=dblarr(nangle+1,lmax+1)
lats=dblarr(lmax+1)
lons=dblarr(nangle+1)

openr,3,'idl.x'     
  for i=0,nangle do begin
    readf,3,tx
    lons(i)=tx
  end
  close,3

openr,2,'idl.y'
  for i=0,lmax do begin
    readf,2,ty
    lats(i)=ty
  end
  close,2

openr,1,'idl.z'
  for i=0,nangle do for j=0,lmax do begin
    readf,1,tz
    z(i,j)=tz
  end
  close,1

; sphere
  map_set, 20, -90, 0, color=plot_color, /noerase, latdel=45, /orthographic, glinestyle=2, glinethick=1, /grid, /noborder, pos=posv

; ellipsoid
;   map_set, 0, -90, 0, color=plot_color, /isotropic, /horizon, /noerase, latdel=45, /aitoff, glinestyle=2, glinethick=1, /grid, /noborder, pos=posv

;define levels for contourplot (including one zero contour!):
;  make_levels, abs_max,abs_min,v,vl,vls,levelnum

;------ contour,z,x,y,nlevels=8,/follow,pos=posv
; arrays for the contour defs:
levelnum=19
v=dblarr(levelnum)
vl=intarr(levelnum)
vls=intarr(levelnum)
ma=max(z)
mi=min(z)
print, ma, mi
;v=findgen(levelnum)*(ma-mi)/(levelnum-1)+mi
;vl=intarr(levelnum)
;lev=0.0
make_levels, ma,mi,v,vl,vls,levelnum

  red =   [0,1,1,0,0,1]
  green = [0,1,0,1,0,1]
  blue =  [0,1,0,0,1,0]
  tvlct, 255*red, 255*green, 255*blue                         


  ;----Radostin
  ;  black and white plots
  ;
  ;  contour, z,lons,lats, /noerase, color=plot_color, levels=v, c_linestyle=vls, c_colors=0 ,c_labels=vl, /overplot, pos=posv

  ;  color plots
  ;
   contour, z,lons,lats, /noerase, color=plot_color, levels=v, c_linestyle=0, c_colors=vls ,c_labels=vl, /overplot, pos=posv

endif

;---------------------------end of program------------

end
