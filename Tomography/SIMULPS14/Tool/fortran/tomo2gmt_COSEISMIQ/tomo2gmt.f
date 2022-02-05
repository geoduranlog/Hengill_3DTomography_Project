c
	program tomo2gmt
c
c 
c
c  This program is a merged version of 'sim2gmt' (Graeber, Haslinger, Husen)
c  and 'fatomo2gmt' (Husen). It has been set up to read the main output file  
c  (named 'output') written by SIMULPS13Q or by FATOMO and to
c  extract
c       - horizontal depth slices 
c       - vertical depth section along constant lat/lon centered
c         around the the origin and spaced as desired
c       - depth sections along user defined profiles
c  for further use in GMT.
c  Velocities (Vp & Vp/Vs) and model changes are extrapolated from
c  the grid they are defined on in the tomography (called seismic grid)
c  on a finer grid with a spacing defined by the user. Interpolation
c  can be done either by trilinear interpolation or by nearest neighbor.
c  Resolution estimates such as KHIT,DWS,RDE are not interpolated.
c
c  FATOMO uses a seismic and an inversion grid. Velocities are defined
c  on the former and it is made of grid nodes with linear interpolation. 
c  Resolution estimates are computed for the inversion grid, which consists
c  of cells. Grid nodes of the seismic grid are stored in the variables 
c  sei_xp,yp,zp in increasing order and block boundaries of the inversion 
c  in the variables inv_xp,yp,zp also increasing. Since SIMULPS is using 
c  only one grid (here called seismic grid) only the variables sei_xp,yp,zp
c  are used. 
c   
c  Velocities and resolution estimates are stored in 3D arrays, which are
c  defined in the 'tomo2gmt_common.inc' file. Thereby, the first element of
c  these arrays corresponds to the grid point with the lowest x,y,z-value (z is
c  counted positively to greater depth). So if you want to use tomo2gmt
c  for your output, just write your own routines to fill in the different
c  arrays correctly. 
c
c  tomo2gmt also support plotting of ray density tensors computed by 
c  TOMEK (E. Kissling, kiss@tomo.ig.erdw.ethz.ch). This option requires
c  two outputfiles of TOMEK: raydens.out & resol_main.out. As a matter of fact
c  the first two blocks of each layer in TOMEK must be combined when computing
c  ray density tensors. tomo2gmt analyzes the ratio of the eigenvalues  of the 
c  ray density tensor to determine four classes of ray illumination within each 
c  block. Thereby, class 1 corresponds to an equal ray distribution in all three
c  dimensions, class 2 to two main ray direction, and class 3&4 to one ray direction.   
c              
c
C       outputfiles of tomo2gmt (in the directory VP):
c	
c       final hypocenter location       -> 'hypo.xyz' (only for SIMULPS)
c       station location                -> 'sta.xyz'
c       coordinates (lat/lon) of grid nodes in 
c        horizontal plane               -> 'plane_nodes.xyz'
c        vertical sections              -> 'vert_nodes.xyz'
c       coordinates (lat/lon) of inversion block boundaries in  (only FATOMO) 
c        horizontal plane               -> 'plane_invblo.xyz'
c        vertical sections              -> 'vert_invblo.xyz'
c	values of each plane (regrided)
c         -> planeXXXXXX.xyz  XXXXXX=depth of plane
c       values for each longitudinal section 
c         -> lonsecNO.xyz  NO=number of section 
c       values for each latitude section
c	  -> latsecNO.xyz  NO=number of section
c       values for profiles
c	  -> profileNO.xyz NO=number of profile
c
c                                                        Aug. 99, S. Husen
c                                                (stephan@tomo.ig.erdw.ethz.ch)
c
c       insert option to output an ascii-file for MATLAB use
c       file is called 'matlab.xyz' and is ordered
c         x  y  z  Vp/VpVs
c       with x as fastest index.
c                                                        Jul. 12 00, S.Husen
c
c       modify subroutine "read_fatomo" to handle modified output of FATOMO,
c       Ver. Oct., 2000. This version of FATOMO allows a more flexible 
c       gridding of the inversion grid.
c						        Oct., 2000, S. Husen
c
c       change code to deal with more than 20 velocity layers 
c       (nz>20)
c                                                       May, 2001, S. Husen
c
c       bugfix: replaced nz with sei_nz in the section where spread
c               values are computed
c                                                       11/26/01 S. Husen
c
c       if all Vp grid nodes are fixed - f.e. inverting only for Vp/Vs -
c       no final Vp, KHIT, DWS, RDE are printed out -> changed tomo2gmt
c       to handle this; final Vp velocities are set to initial Vp velocities
c
c                                                       02/12/02 S. Husen
c
c       made changes to cope with simul2000 output
c        -> no more than 13 values of RDE are printed per row in simul2000
c       allow spacing between lat and lon depth sections to be either in 
c         degrees or kilometers; this was introduced by Steve Sherburn
c
c                                                       02/10/04 S. Husen
c
c      change in format to read DWS for simulps14 output
c                                                       08/09/04 S. Husen
       implicit none
c
	include 'tomo2gmt_common.inc'
c
	real hitct,cutoff
	real vp_tmp(mx),vs_tmp
	real g(1,mx),gsum(mx),gsum2(mx)
	real gfact(mx),dist(1,mx),tmp,tmp2
	real sum,radd
	real xg(mx),yg(mx),zg(mx)
	real distx,disty,distz
        integer ialt
	real x,y,z,alt,lat,lon,dlat,dlon,in
	real xsh(2,20),ysh(2,20),xs,ys,zstep,xystep
	real wzmax,wzmin
	real rota
	real dx,dy,dxy
	real v_no,v_dist
	real ddist,proflnght,intpoints
	real vtot,vstot,area,sarea,vavr(nzmax),vsavr(nzmax)
	real rv,rvps
	real xll(2,20),yll(2,20)
	real drad
	real pz(10)
	real xmin,xmax,ymin,ymax,zmin,zmax
	parameter (drad=1.7453292d-2)
	integer i,ii,j,jj,k,kk,l,nx,ny,nz,ix,iy,iz,kz,switch,k1,k2
	real ilat,ilon
	integer slat,slon
	integer profno,nop
	integer nyfrst,nxfrst
	integer irow(mx),ig(mx),jg(mx),kg(mx)
	integer ndexfx(mx),mdexfx(mx)
	integer nxy,nx2,nxy2,nz2
	integer nparp, nparcur
	integer npar,nparv,npari,nparvi,nparpi,nparsi,nrowp,nrows,nrowst
	integer nrow,nn,jy,i_ii,j_ii,k_ii
	integer iuses,avp
	integer prog_ver
	character*7 verno
	character line*132,depth*6,cn(25)*2,sta*4
	character*1 ap,av,al,as,am,ano,asp
        character*1 interp,aref,ao,ap2,raydens,ab,zoom
	character*1 amat
csh 10/02/04  new parameter dk; specifies if distance between lat or lon sections are in deg or km
	character*1 dk
	character*8 alat,alon
	character ns,ew
	character*20 cvref
	character*6  fname
c
 	data cn/'01','02','03','04','05','06','07','08','09','10','11',
     &  '12','13','14','15','16','17','18','19','20','21','22','23',
     &  '24','25'/

	verno='072704'
c
	write(6,*)' '
	write(6,*)'           tomo2gmt Ver. ',verno
	write(6,*)'-------------------------------------'
	write(6,*)' '
	write(6,*)'This program reads the output file "output"'
        write(6,*)'of SIMULPS13Q'
	write(6,*)'or FATOMO and performs regridding and interpolating '
	write(6,*)'of p-velocity, velocity change, vp/vs in'
	write(6,*)'  - plane sections, '
	write(6,*)'  - depth sections along longitude/latitude, '
	write(6,*)'  - depth sections along specified profiles'
	write(6,*)' for further use with GMT to produce b/w or colour coded plots.'
	write(6,*)' Outputfiles will be written in the directory VP/.'
	write(6,*)' '
	write(6,*)'which output is read in SIMULPS or FATOMO (s/f) or none (n)?'
	write(6,*)' the latter option can be used to plot only ray density tensors'
	read(5,'(a1)')ao
	if (ao.eq.'S') ao='s'
	if (ao.eq.'F') ao='f'
	if (ao.eq.'N') ao='n'
c  next questions only if to analyze SIMULPS or FATOMO output
	if (ao.ne.'n') then
	   write(6,*)'what kind of interpolation for velocities'
	   write(6,*)'NEAREST NEIGHBOUR(n) or LINEAR INTERPOLATION(i) ?'
	   read(5,'(a1)')interp
	   write(6,*)'ENTER DESIRED GRID SPACING FOR REGRIDDING IN XY-DIRECTION'
	   write(6,*)'  (for horizontal depth sections)'
	   write(6,*)'  (for SIMULPS grid spacing must be in bld-units!)'
	   read(5,*)xystep
	   write(6,*)'ENTER DESIRED GRID SPACING FOR REGRIDDING IN Z-DIRECTION'
	   write(6,*)'  (for vertical depth sections)'
	   write(6,*)'  (for SIMULPS grid spacing must be in bld-units!)'
	   read(5,*)zstep
	endif
	write(6,*)'Do you want to create HORIZONTAL DEPTH SECTIONS (y/n)?'
	read(5,*)ap
 	write(6,*)'  VERTICAL DEPTH SECTIONS along CONSTANT LATITUDE (y/n)?'
	read(5,*)av
	write(6,*)'  VERTICAL DEPTH SECTIONS along CONSTANT LONGITUDE (y/n)?'
	read(5,*)al
	write(6,*)'  VERTICAL DEPTH SECTIONS along USER DEFINED PROFILES (y/n)?'
	write(6,*)'    ray density tensors for such profiles are not supported yet!'
	read(5,*)as
	write(6,*)' '
	write(6,*)'  an MATLAB inputfile (ascii, xyzv) (y/n)?'
	write(6,*)' '
	read(5,*)amat
c  next questions only if to analyze SIMULPS or FATOMO output
	if (ao.ne.'n') then
	   write(6,*)'Changes (abs. and %) of model parameters can be calculated'
	   write(6,*)'relative to:'
	   write(6,*)'  (1)  initial reference 1D model'
	   write(6,*)'  (2)  initial reference 2D/3D model'
	   write(6,*)'  (3)  1D model of final average layer velocity'
	   write(6,*)'  (4)  user defined, external 1D model'
	   write(6,*)'Default is (1).'
	   read(5,'(a1)')aref
	   if (aref.eq.' ') aref='1'
	   write(6,*)' '
	   write(6,*)'Do you want to calculate SPREAD VALUES [n]?'
	   write(6,*)' This has not been implemented for FATOMO'
	   write(6,*)' and flexible gridding, yet!!'
	   write(6,*)' SPREAD VALUES requires the prior computation of the '
	   write(6,*)' FULL resolution matrix stored in the file '
	   write(6,*)' resol.out'
	   read(5,'(a1)')asp
	   if (asp.eq.' ') asp='n'
	   write(6,*)' '
	   write(6,*)'Do you want to have 1D models of absolute'
	   write(6,*)'velocities at user defined lon/lat positions'
	   write(6,*)'to be calculated [n]'
	   read(5,'(a1)')am
	   if (am.eq.' ') am='n'
	   write(6,*)'Do you want a normalization KHIT/DWS relative'
	   write(6,*)'to the maximum [n]?'
	   read(5,'(a1)')ano
	   if (ano.eq.' ') ano='n'
	   write(6,*)' '
	endif
	write(6,*)'Do you want to plot ray density tensors [n]?'
	write(6,*)' Up to now Ray Density Tensors are only computed'
	write(6,*)' by TOMEK and plotting of Ray Density Tensors'
	write(6,*)' requires TOMEK outputfiles "raydens.out" and'
	write(6,*)' "resol_main.out". This means you need to have'
	write(6,*)' an extra run of TOMEK along with SIMULPS or FATOMO'
	write(6,*)' since inversion results by TOMEK are not supported.'
	read(5,'(a1)')raydens
	if (raydens.eq.' ') raydens='n'
c open output file (if SIMULPS or FATOMO)
	if (ao.ne.'n') open(2,file='output',status='old',err=9988)
c
	switch=0
	hitct=0
	csr=1
	snr=0
csh 29/03/01 
c  what about people using bld=0.1?????	
c	bld=1.0
c
c initialization
c
	pobsmax=1
	phitmax=1
	pdwsmax=1.0
	presmax=0.1
	sobsmax=1
	sdwsmax=1.0
	sresmax=0.1
c
	do iz=1,nzmax
	 do iy=1,nymax
	  do ix=1,nxmax
	   pobs(ix,iy,iz)=0
	   phit(ix,iy,iz)=0
	   pdws(ix,iy,iz)=0.01
	   pres(ix,iy,iz)=0.0
	   sobs(ix,iy,iz)=0
	   sdws(ix,iy,iz)=0.01
	   sres(ix,iy,iz)=0.0
	  enddo
	  enddo
	enddo
c
c  compute ray density tensor first if required
c
	if (raydens.eq.'y') then
	   if (ao.eq.'n') then 
	     call ray_dens(ap,av,al)
	   else
	     write(6,*)' '
	     write(6,*)' ERROR: '
	     write(6,*)'  To plot ray density tensors you need to have'
	     write(6,*)'  output files from TOMOEK !'
	     write(6,*)'  Specify "none" when asked for output file format'
	     goto 7777
	   endif
	endif
c
c  read cut off value hitct for DWS
c   (only for SIMULPS)
c
	do j=1,50000
c
	  if (ao.eq.'s') then
 1	  read(2,'(a132)',end=99)line
c csh change to cope with simulps2000
c
c  determining which version of simulps
c
	    if(line(2:8).eq.'Program') then
	     if(line(10:18).eq.'Simulps13') then
	      prog_ver=13
	      write(6,*)' '
	      write(6,*)' You were running Simulps13....'
	      write(6,*)' '
	     endif
	     if(line(10:18).eq.'Simulps14') then
	      prog_ver=14
	      write(6,*)' '
	      write(6,*)' You were running Simulps14....'
	      write(6,*)' '
	     endif
	     if(line(10:18).eq.'simul2000') then
	      prog_ver=2000
	      write(6,*)' '
	      write(6,*)' You were running simul2000....'
	      write(6,*)' '
	     endif
	     if (prog_ver.ne.13.and.prog_ver.ne.14.and.prog_ver.ne.2000) then
	      write(6,*)' '
	      write(6,*)' Error: You used a different version of Simulps'
	      write(6,*)'        that is not supported by tomo2gmt'
	      write(6,*)' tomo2gmt only supports Simulps14 or simul2000'
	      write(6,*)' terminating....'
	      stop
	     endif
	    else
	     goto 1
	   endif
	  endif
c
	  if (ao.eq.'s') then
 2      read(2,'(a132)',end=99)line
 	   if (line(2:5).eq.'neqs') then
 	    read(2,'(t58,f5.0)')hitct
 	    write (6,1111)hitct
 1111	    format('  nodes with a DWS less than ',f5.0,' are fixed\n')
 	    write (6,*)' '
 	   else
 	    goto 2
 	   endif
c
	  endif
c
c  reading iuses
c   (only for SIMULPS)
c
	 if (ao.eq.'s') then
 3	 read(2,'(a132)',end=99)line
	  if(line(2:5).eq.'Only') then
	   iuses=0
	  else if(line(2:5).eq.'Both') then
           iuses=1
	  else
	   goto 3
        endif
	  if (iuses.eq.0) then
	     write(6,*)' Inversion was not for Vp/Vs'
	     write(6,*)' '   
         else
	     write(6,*)' Inversion was for Vp/Vs'
	     write(6,*)' '   
        endif
	 endif

c
c  reading origin and short distance conversion
c
	  write(6,*)'READING ORIGIN AND SHORT DISTANCE PARAMETER'
 4      read(2,'(a132)',end=99)line
	  if (ao.eq.'s') then   ! for SIMULPS
	   if(line(3:10).eq.'origin :')then
	    read(2,*)ilat,dlat,ilon,dlon,rota
	    if (ilat.lt.0) then
	       olat=ilat-dlat/60
	    else
	       olat=ilat+dlat/60
	    endif
	    if (ilon.lt.0) then
	     olon=-1.*(ilon-dlon/60.)
            else
	     olon=-1.*(ilon+dlon/60.)
	    endif
	    write(6,*)'ORIGIN  :',olat,olon
	   else
	    goto 4
	   endif
	  else   ! for FATOMO
	   if(line(6:26).eq.'short distance origin')then
	    read(2,*)ilat,dlat,ilon,dlon

	    if (ilat.lt.0) then
	       olat=ilat-dlat/60
	    else
	       olat=ilat+dlat/60
	    endif
	    if (ilon.lt.0) then
	     olon=-1.*(ilon-dlon/60.)
            else
	     olon=-1.*(ilon+dlon/60.)
	    endif
	    write(6,*)'ORIGIN  :',olat,olon
	   else
	    goto 4
	   endif
	  endif ! short distance reading
c
 5      read(2,'(a132)',end=99)line
	  if (ao.eq.'s') then   ! for SIMULPS
	   if(line(11:24).eq.'one min lat   ')then
	    read(line,'(t25,f6.4)')latkm
	    read(2,'(t25,f6.4)')lonkm
	    latkm=latkm*60.
	    lonkm=lonkm*60.
	    write(6,*)'SHORT DISTANCE CONVERSION :'
	    write(6,*)'latkm= ',latkm
	    write(6,*)'lonkm= ',lonkm
	    write(6,*)'ROTATION : ',rota,' deg'
	    rota=rota*drad
	    csr=cos(rota)
	    snr=sin(rota)
	   else
	    goto 5
	   endif
	  else  ! for FATOMO
	   if(line(11:24).eq.'one min lat   ')then
	    read(line,'(t25,f6.4)')latkm
	    read(2,'(t25,f6.4)')lonkm
	    latkm=latkm*60.
	    lonkm=lonkm*60.
	    write(6,*)'SHORT DISTANCE CONVERSION :'
	    write(6,*)'latkm= ',latkm
	    write(6,*)'lonkm= ',lonkm
	   else
	    goto 5
	   endif
	  endif
c
c  reading stations
cfh read N,S E,W from station format to get correct lat,lon
cfh  for GMT convention (S,W  < 0 ; N,E > 0) 
cfh  assuming that all degrees in simul are unsigned
cfh  z in depth > 0 !!! (so < 0 above sea level)
c
 6	  read(2,'(a132)',end=99)line
 1005     format(t9,a4,1x,i3,a1,f5.2,3x,i3,a1,f5.2,3x,i4,1x,3f7.2)
	  if(line(1:10).eq.'   station'.or.line(1:8).eq.' station') then
	    write(6,*)' '
	    write(6,*)'READING STATION COORDINATES'
	    open(1,file='VP/sta.xyz')
	    do i=1,1000
	      read(2,1005)sta,slat,ns,dlat,slon,ew,dlon,ialt,x,y,z
              alt = float(ialt)
	      if(sta.eq.' ')goto 10
            if(ns.eq.'S'.or.ns.eq.'s') then
              lat=-1.*(slat+dlat/60.)
            else
	          lat=(slat+dlat/60.)
            endif
            if(ew.eq.'W'.or.ew.eq.'w') then
	          lon=-1.*(slon+dlon/60.)
            else
	           lon=(slon+dlon/60.)
	       endif
	     write(1,'(2f10.2,3f8.1,2x,a4)')lon,lat,x,y,z,sta
	    enddo
	  else
	    goto 6
	  endif
 10       close(1)
c
c  now reading velocities and resolution estimates
c  depending on tomography code used
c
	  if (ao.eq.'s') then
	   call read_simulps(iuses,aref,ab,hitct,prog_ver)
	   goto 99
	  else
	   call read_fatomo(aref)
	   goto 99
	  endif
c
        enddo		! end of big reading do-loop
c
 99     continue
   	close(77)
c
	write(6,*)' '
	write(6,*)'  MAXIMUM OF '
	write(6,*) '  POBSMAX =',pobsmax
	write(6,*) '  PHITMAX =',phitmax
	write(6,*) '  PDWSMAX =',pdwsmax
	write(6,*) '  PRESMAX =',presmax
	write(6,*) '  TRACE of resolution matrix =',trace
       if (iuses.eq.1) then
	write(6,*) '  SOBSMAX =',sobsmax
	write(6,*) '  SDWSMAX =',sdwsmax
	write(6,*) '  SRESMAX =',sresmax
       endif
c
c
c normalisation of POBS(KHIT) and PDWS 
C
      if (ano.eq.'y') then
	write(6,*)'NORMALISATION OF POBS(KHIT), PHIT(HIT) AND PDWS '
	write(6,*)' and SOBS and SDWS'
c defining nx,ny,nzout
	if (ao.eq.'s') then 
	   nx=sei_nx    ! SIMULPS
	   ny=sei_ny
	else
	   nzout=inv_nz
	endif
	do iz=1,nzout
	  if (ao.eq.'f') then
	   nx=inv_nx(iz)-1  ! FATOMO: KHIT & DWS are defined on inv. grid
	   ny=inv_ny(iz)-1
	  endif
	  do iy=1,ny
	    do ix=1,nx
	      pobs(ix,iy,iz)=1.-pobs(ix,iy,iz)/pobsmax
	      phit(ix,iy,iz)=1.-phit(ix,iy,iz)/phitmax
	      pdws(ix,iy,iz)=1.-pdws(ix,iy,iz)/pdwsmax  
	      sobs(ix,iy,iz)=1.-sobs(ix,iy,iz)/pobsmax
	      sdws(ix,iy,iz)=1.-sdws(ix,iy,iz)/pdwsmax  
	    enddo
	  enddo
	enddo
      endif
c
c calculation of velocity change in absolut values and percentage
c also for percentage change for vpvs-ratio
c
c if change rel. to average layer velocity compute first
c
c for SIMULPS edge-nodes are ot included in the inversion,
c velocities are defined on seismic grid for FATOMO
	if (ao.eq.'s') then
	   nxfrst=2      ! SIMULPS
	   nx=sei_nx-2
	   nyfrst=2
	   ny=sei_ny-2
	   nz=nzout-1
	else             ! FATOMO
	   nxfrst=1
	   nx=sei_nx
	   nyfrst=1
	   ny=sei_ny
	   nzfrst=1
	   nz=sei_nz
	endif
	if (aref.eq.'3') then
	 write(6,*)'Average layer velocities after inversion:'
	 do k=nzfrst,nzout
 	  vtot=0
	  area=0
	  sarea=0
	  vstot=0
	  do j=nyfrst,ny
	   do i=nxfrst,nx
	    if (pdws(i,j,k).ge.hitct) then    ! take only nodes which are included in inversion
	     dx=abs(sei_xp(i+1)-sei_xp(i))
	     dy=abs(sei_yp(j+1)-sei_yp(j))
	     dxy=dx*dy
	     vtot=vtot+vpm(i,j,k)*dxy
	     area=area+dxy
	    endif
	    if (iuses.eq.1) then    ! Vp/Vs part
	     vstot=vstot+vpvs(i,j,k)*dxy
	     sarea=sarea+dxy
	    endif
           enddo
          enddo
	  if (area.gt.0.00001) then	! only then it makes sense
	   vavr(k)=vtot/area
	   if (iuses.eq.1) then 
	     vsavr(k)=vstot/sarea
	   write(6,*)vavr(k),vsavr(k)
           else
	    write(6,*)vavr(k)
           endif
          else 
           vavr(k)=0.00         ! to have a value in there
	   vsavr(k)=0.00
          endif
         enddo
	 do iz=nzfrst,nzout
	  do iy=nyfrst,ny
	   do ix=nxfrst,nx
	    refvp(ix,iy,iz)=vavr(iz)
	    if (iuses.eq.1) then
	      refvpvs(ix,iy,iz)=vsavr(iz)
            endif
	   enddo
	  enddo
	 enddo
	endif
c now compute change         	 
	do iz=nzfrst,nz
	  do iy=nyfrst,ny
	    do ix=nxfrst,nx
	      dvp(ix,iy,iz)=vpm(ix,iy,iz)-refvp(ix,iy,iz)
	      pvp(ix,iy,iz)=(100.*dvp(ix,iy,iz))/refvp(ix,iy,iz)
	      if (iuses.eq.1) then
                pvpvs(ix,iy,iz)=(100.*(vpvs(ix,iy,iz)-refvpvs(ix,iy,iz))
     .              )/refvpvs(ix,iy,iz)
              endif
	    enddo
	  enddo
	enddo
c
9998    continue
c
c
	wzmin=0.
	wzmax=0.
c
c  if desired calculation of spread value S(ijk)
C   it is defined as 
c    S(ijk)=Log {1/norm(row(ijk)) sum {dis(ijk:ijk)*res(ijk)**2}}
c
C  the computation requires the calculation of the FULL resolution
c  matrix in simul (ires=2) and printed in file=resol.out
c
c  resolution matrix consists of nrowp rows of Vp parameters
c   and nrows rows of Vp/Vs parameters
c  each row has nparvi elements (nparvi = nparpi + nparsi)
c
csh has been not implemented for FATOMO yet
c
	if (asp.eq.'y') then
c  nodes file output from simulps (file 18, f18)
	 if (ao.eq.'s') then
          open(unit=10,file='nodes.out',status='old',form='formatted')
          rewind (10)
	 endif
c  resolution file (file 17, resol)
         open(unit=12,file='resol.out',status='old',form='formatted')
         rewind (12)
	 write(6,*)' '
	 write(6,*)' compute SPREAD VALUES '
	 write(6,*)'  this may take some time depending on your grid size'
	 write(6,*)' '
c  set nx,ny
	 if (ao.eq.'s') then
	    nx=sei_nx   ! seismic grid of SIMULPS
	    ny=sei_ny
	 else
	    nzout=inv_nz
	 endif
c
c  set all spread values to 9999
c
	 do iz=1,nzout
	  if (ao.eq.'f') then
	   nx=inv_nx(iz)-1  ! FATOMO: reso estimates are defined on inv. grid
	   ny=inv_ny(iz)-1
	  endif
	  do iy=1,ny
	   do ix=1,nx
	     spred_p(ix,iy,iz)=9999
	     spred_s(ix,iy,iz)=9999
	     spred_vpvs(ix,iy,iz)=9999
	   enddo
	  enddo
	 enddo
c
c  compute total number of gridpoints (nodes)
c
         nxy=nx*ny
	 if (ao.eq.'s') then
            nx2=nx-2  ! number non-edge nodes in row
	    nxy2=nx2*(sei_ny-2)	! number non-edge nodes in layer
	    nz2=sei_nz-2
	 endif
c         if (ao.eq.'f') then            
c            nx2=inv_nx-1         ! number of inv. blocks in row
c	    nxy2=nx2*(inv_ny-1)  ! number of inv. blocks in layer
c	    nz2=inv_nz-1
c	 endif
c         nodes2=nz2*nxy2
c
c define the number of columns and rows for the resolution matrix
c
cfh For now we assume that station parameters are not included in the resolution
cfh   matrix! check that some day...
c
	 if (ao.eq.'s') then 
	   read(10,*)npar,nparv,npari,nparvi,nparpi,nparsi,
     &               nrowp,nrows,nrowst,radd
	 endif
	 if (ao.eq.'f') then
	   npar=nxy2*nz2  ! number of parameters
	   nparv=npar     ! in FATOMO no parameters can be apriori fixed
	   npari=npar
	   nparvi=npar
	   nparpi=npar
	   nrowp=npar     
	   nrows=0        ! no Vp/Vs in FATOMO
	   radd=0.0001
	 endif
c  npar = number of parameters
c  nparv = number of velocity parameters (npar - station corrections parameter)
c  npari = number of parameters inverted for
c  nparvi = number of velocity parameters inverted for (Vp + Vp/Vs)
c  nparpi =  number of Vp parameters inverted for (without taking into account the deliberatly fixed ones!)
c  nparsi =  number of Vp/Vs parameters inverted for ( " )
c  nrowp = number of Vp rows  (now taking into account the deliberatly fixed ones)
c  nrows = number of Vp/Vs rows ( " )
c  nrowst = number of station corrections rows
c 
	 nrow=nrowp+nrows ! total number of rows (Vp + Vp/Vs)
	 nparp = nparv/2  ! number of p-velo parameters in grid (=nodes2 ?)
	 if(npar.eq.npari) goto 940
         do 920 l=1,npar
          read(10,3020) i,ndexfx(i)
          if(ndexfx(l).eq.0) goto 920
          mdexfx(ndexfx(i))=i
 3020     format(2i7)
  920    continue
  940    continue
cfh         nnodes=nparv ! nnodes is never used for anything
cfh 	 nodes_p=nnodes ! nodes_p is never used for anything...
c
c compute x,yz, positions in the absolute grid
c of each node within a row of the resolution matrix
c
	do i=1,nparvi
	   if (ao.eq.'s') then  ! SIMULPS
	      if (npar.ne.npari) then
                 nn=mdexfx(i)   ! grid nodes apriori fixed
	      else
		 nn=i           !no grid nodes apriori fixed
	      endif
	      kz=(nn-1)/nxy2+2
	      jy=2+(nn-1+(2-kz)*nxy2)/nx2
	      ix=1+nn+nx2*(2-jy)+nxy2*(2-kz)
	      ig(i)=ix
	      xg(i)=sei_xp(ig(i))
	      jg(i)=jy
	      yg(i)=sei_yp(jg(i))
	      if(kz.ge.sei_nz) kz=kz-nz2 ! if s-velocity node
	      kg(i)=kz
	      zg(i)=sei_zp(kg(i))
c	   else   ! FATOMO
c	      nn=i  ! no apriori fixed grid nodes
c	      kz=(nn-1)/nxy2+1
c	      jy=1+((nn-1+(1-kz)*nxy2)/nx2)
c	      ix=nx2-(nn-1)+nx2*(jy-1)+nxy2*(kz-1)
c	      ig(i)=ix
c	      xg(i)=(abs(inv_xp(ix))+abs(inv_xp(ix+1)))/2
c	      jg(i)=jy
c	      yg(i)=(abs(inv_yp(jy))+abs(inv_yp(jy+1)))/2
c	      kg(i)=kz
c	      zg(i)=(abs(inv_zp(kz))+abs(inv_zp(kz+1)))/2
	   endif
	enddo
c
c read input data file (fort.17 format)
c    process one row after the other (i-index)
c
      i=1
  950 read(12,3720,end=777) irow(i)
 3720 format(/,4x,i5,/)
csh 02/10/04 change to cope with simul2000
 3721   format(18f7.4)
 3722   format(20f7.4)
	if (prog_ver.eq.14.or.prog_ver.eq.13) then
        read(12,3721)(g(i,j),j=1,npari)
 	else ! simul2000
 	 read(12,3722)(g(i,j),j=1,npari)
       endif
       do 945 j=1,npari
  945   g(i,j)=g(i,j)+radd
c  gridpoint number in array of all nodes
      if (npar.ne.npari) nn=mdexfx(irow(i))
      if (npar.eq.npari) nn=i  ! no apriori fixed nodes in FATOMO
c
c calculate the averaging factor for a given row
c if i<=nrowp -> do Vp spread calculation
c if i>nrowp && i<=nrows do Vp/Vs calculation
c
c for Vp spread function take only narpi nodes
cfh 
cfh here come the changes: do pure Vp-spread if iuses==0 (no S-data used
cfh   in inversion), in all other cases do Vp and Vp/Vs spread as DMEP:
cfh   use full row for averaging vector and scaling, and compute Vp/Vs-spread
cfh   as spread (vp_elements + vp/vs_elements)
c
c     use nparcur as the valid number of parameters (=number of elements per row) 
       if (iuses.eq.0) nparcur = nparpi
       if (iuses.eq.1) nparcur = nparpi + nparsi
       tmp=0.
       tmp2=0.
c       do j=1,nparpi
       do j=1,nparcur
          tmp=g(i,j)+tmp
          tmp2=g(i,j)**2 +tmp2
       enddo
       gsum(i)=tmp
       gsum2(i)=tmp2
c normalize such that the sum of the square of g(i,j) for each row
c is equal to 1
c
cfh I don't understand the next part: if tmp is ~ 0 , set g(i,j) to 0 for 
cfh  every element and then compute tmp again??????????
cfh  why not just set tmp to 1 e-14?
        if(tmp2.LT.1.0e-5) then
c	 do j=1,nparpi
	 do j=1,nparcur
	   g(i,j)=1.0e-7
	 enddo
	 tmp=0.
c         do j=1,nparpi
         do j=1,nparcur
	    tmp2=g(i,j)**2 + tmp
         enddo
        endif
	sum=sqrt(tmp2)
c
c	do j=1,nparpi
	do j=1,nparcur
	  g(i,j)=g(i,j)/sum
	enddo
cfh - for rows with sums close to zero g(i,j) now is of order 10 e4 ?!?!?!?!
c
	gfact(i)=(1./sum)
c
c determine the position in the 3D grid and the distances
c
c	do jj=1,nparpi
	do jj=1,nparcur
c calculate distance
          distx=xg(irow(i))-xg(jj)
          disty=yg(irow(i))-yg(jj)
          distz=zg(irow(i))-zg(jj)
c
          dist(i,jj)=sqrt(distx*distx + disty*disty + distz*distz)
	enddo
c
c calculate the spread function for point irow which is at the
c position i_ii,j_ii,k_jj in the 3D grid
c
        i_ii=ig(irow(i))
        j_ii=jg(irow(i))
        k_ii=kg(irow(i))
c	
cfh here split into Ponly, Pfrom P+PS, PSfrom P+PS
c
cfh P only
	if (iuses.eq.0) then
          tmp=0.
	  do jj=1,nparpi
             tmp=dist(i,jj)*g(i,jj)**2+tmp
	  enddo
c
csh
c	  write(6,7878)irow(i),i_ii,j_ii,k_ii
c 7878	  format('node ',i5,i3,i3,i3)
csh
	  spred_p(i_ii,j_ii,k_ii)= log10(gfact(i)*tmp)  ! Vp node
        else
c
c  now for Vp + Vp/Vs spread 
c
c  if the current row is a P-row use P-elements only:
          if (i.le.nrowp) then
            tmp=0.
            do jj=1,nparpi
              tmp= dist(i,jj)*g(i,jj)**2+tmp
            enddo
            if (gsum2(i).lt.1.0e-5) then
	       vp_tmp(nn)=10000		!need vp_tmp as array for all nodes now
	       spred_p(i_ii,j_ii,k_ii)= 5.		! should not be neccessary ?
            else
	       vp_tmp(nn)=gfact(i)*tmp
	       spred_p(i_ii,j_ii,k_ii)=log10(gfact(i)*tmp)
            endif
          else
c  for the S-rows, use lower half of elements (correspond to S-S)
            tmp=0.
cfh dist(i,jj) should be the same for jj and nparpi+jj, as the P and PS grid are equal, 
cfh  and nparpi = nparsi. spred_s is the SxS spread,and spred_vpvs is computed from 
cfh  vp_tmp(nn-nparp)*vs_tmp, where nn-nparp is the index of the p-value for the
cfh  current s-value (= same grid node)
           do jj=1,nparsi
             tmp= dist(i,nparpi+jj)*g(i,nparpi+jj)**2+tmp
           enddo
           if (gsum2(i).lt.1.0e-5) then
	      vs_tmp=10000
	      spred_s(i_ii,j_ii,k_ii)= 5.		! should not be neccessary ?
           else
	      vs_tmp=gfact(i)*tmp
	      spred_s(i_ii,j_ii,k_ii)=log10(gfact(i)*tmp)
           endif
cfh  the vp/vs spread - now coupled
cfh  (nn-nparp) should be the index of the P-row for this S-row 
	 spred_vpvs(i_ii,j_ii,k_ii)=log10(vp_tmp(nn-nparp)+vs_tmp)
c
          endif	!if p or s element
        endif     !iuses?
c
	i=i+1
	goto 950
	endif
  777	close(10)
	close(12)
c
csh basically, output structure for planes, sections and profiles is identical
csh for SIMULPS & FATOMO. except: 
csh   - first grid points to be output for SIMULPS are the second ones
csh   - HIT,KHIT,DWS,RDE are defined on the inversion grid for FATOMO.
csh     hence we need a different interpolation, resp. no interpolation
c    
cfh common output format for planes and sections:
cfh 2 headerlines: first one describing the contents of the file (plane or section, which one, 
cfh   which kind of reference velocity (and the value for planes)
cfh  second one describing the format, which is:
cfh  x y lat lon z vp pvp dvp pkhit pdws presol spread_p vpvs pvpvs skhit sdws sresol spread_s
cfh x,y,z are the cartesian coordinates (simul convention, x increasing to W, y increasing to N, 
cfh   z increasing to depth)
cfh lat, lon are the geographic coordinates (GMT-convention, N,E positive)
cfh vp, pvp, dvp are final P-velocity, percentage change against reference, absolute change against
cfh   reference (km/s)
cfh pkhit, pdws, presol, spread_p are hitcount, derivative weight sum, resolution diagonal and
cfh   spread value for Vp
cfh vpvs, pvpvs  are final Vp/Vs, percentage change against reference
cfh skhit, sdws, sresol, spread_s are the same values for Vp/Vs as the ones above for Vp
cfh !! watch out - for now spread_s is maybe not correct !!
cfh for profiles the inline-coordinate is included as first value
c
        if (aref.eq.'1') then 
	   cvref = 'min 1D initial model'
        else if (aref.eq.'2') then
	   cvref = '2D/3D initial model'
        else if (aref.eq.'3') then
	   cvref = '1D layer average'
        else if (aref.eq.'4') then
	   cvref = 'ext. 1D model'
        endif
c  regriding of plane-depth-sections
c
        if(ap.eq.'y')then
c
	if (ao.eq.'s') then
	 if (ab.eq.'n') then  ! do not output boundary nodes
	   nxfrst=2    ! first grid node is the second one (SIMULPS)
	   nx=sei_nx-1
	   nyfrst=2
	   ny=sei_ny-1
	 else
	   nxfrst=1    
	   nx=sei_nx
	   nyfrst=1
	   ny=sei_ny
	 endif	   
	else  ! FATOMO
	   nxfrst=1
	   nx=sei_nx
	   nyfrst=1
	   ny=sei_ny
	   nzfrst=1
	   nzout=sei_nz
	endif
c
c calculation of coordinates for original nodes in
c horizontal plane
c
	open(33,file='VP/plane_nodes.xyz',err=9999)
	write(6,*)' '
	write(6,*)'--writing plane nodes in file--'
	write(6,*)'   VP/plane_nodes.xyz'
	do i=nyfrst,ny
	  do j=nxfrst,nx
	    y=sei_yp(i)
	    x=sei_xp(j)
            call km2deg(x,y,lat,lon)
	    write(33,'(2f9.2,2f10.4)')x,y,lat,lon
	  enddo
	enddo
	close(33)
c now for inversion block edges (only FATOMO)
	if (ao.eq.'f') then
	   write(6,*)'--writing inversion block boundaries in file--'
	 do iz=1,inv_nz
	  write(fname,'(f6.1)')inv_zp(iz)
	  do l=1,6
	    if(fname(l:l).eq.' ')fname(l:l)='0'
	  enddo
	  open(34,file='VP/plane_invblo'//fname//'.xyz')  ! file for inversion block edges
	  write(6,*)'   VP/plane_invblo'//fname//'.xyz'
	  do j=1,inv_ny(iz)
	   x=inv_xp(iz,1)
	   y=inv_yp(iz,j)
	   call km2deg(x,y,lat,lon)
	   write(34,'(6f10.4,1x)')lat,lon,y,x,inv_zp(iz),inv_zp(iz+1)
	   x=inv_xp(iz,inv_nx(iz))
	   call km2deg(x,y,lat,lon)
	   write(34,'(6f10.4,1x)')lat,lon,y,x,inv_zp(iz),inv_zp(iz+1)
	   write(34,*)'>'
	  enddo
	  do i=1,inv_nx(iz)
	    y=inv_yp(iz,1)
	    x=inv_xp(iz,i)
            call km2deg(x,y,lat,lon)
	    write(34,'(6f10.4,1x)')lat,lon,y,x,inv_zp(iz),inv_zp(iz+1)
	    y=inv_yp(iz,inv_ny(iz))
	    x=inv_xp(iz,i)
            call km2deg(x,y,lat,lon)
	    write(34,'(6f10.4,1x)')lat,lon,y,x,inv_zp(iz),inv_zp(iz+1)
	    write(34,*)'>'
	  enddo
	 enddo
	 close(34)
	endif
c
c now output the planes
c
c  choose between layers at seismic gridpoints or input put user
c
	write(6,*)' '
	write(6,*)' Horizontal depth sections can be printed out'
	write(6,*)' either at the predefined seismic grid (as in the'
	write(6,*)' MOD file) or at user defined grid points.'
	write(6,*)' '
	write(6,*)' Do you want horizontal depth sections at the pre-'
	write(6,*)' defined seismic grid [y] ?'
	read(5,'(a1)')ap2
	if (ap2.eq.' ') ap2='y'
	if (ap2.eq.'y') then
c now define boundary of horizontal section
	  write(6,*)' '
	  write(6,*)'Do you want to plot just a part of the model between'
	  write(6,*)'xmin,xmax and ymin,ymax [n] ? '
	  write(6,*)'Otherwise the entire model will be plotted.'
	  read(5,'(a1)')zoom
	  if (zoom.eq.' ') zoom='n'
	  if (zoom.eq.'y') then
	     write(6,*)'  please enter xmin,xmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)xmin,xmax
	     write(6,*)'  please enter ymin,ymax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)ymin,ymax
	  else
	     xmin=sei_xp(nxfrst)
	     xmax=sei_xp(nx)
	     ymin=sei_yp(nyfrst)
	     ymax=sei_yp(ny)
	  endif
	 write(6,*)' '
	 write(6,*)'--creating files with PLANE SECTIONS--'
	 do iz=nzfrst,nzout-1
	  depth=azp(iz)
	  z=sei_zp(iz)
	  do l=1,6
	    if(depth(l:l).eq.' ')depth(l:l)='0'
	  enddo
	  write(6,*)' VP/plane'//depth//'.xyz'
	  open(3,file='VP/plane'//depth//'.xyz')
c add one descriptive headerline above
          if(aref.eq.'2') then	! dummy values for reference velo 
	     rv = 99.99
	     rvps = 99.99
          else 
             rv = refvp(1,1,iz)
             rvps = refvpvs(1,1,iz)
          endif
            write(3,1999)iz,z,rv,rvps,cvref
 1999     format('#> plane ',i3,' depth ',f5.2,' ref. Vp ',f5.2,
     &     ' ref. VpVs ',f5.2,3x, a20)
c headerline in output to identify the values easylier
	  write(3,'(a,a)')
     &'# x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,',
     &'presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s'
	  do y=ymin,ymax,xystep
	    do x=xmin,xmax,xystep
	      call vel3mod(vpm,x,y,z,v,interp)
	      call vel3mod(dvp,x,y,z,dv,interp)
	      call vel3mod(pvp,x,y,z,pv,interp)
	      call vel3mod(vpvs,x,y,z,vps,interp)
	      call vel3mod(pvpvs,x,y,z,pvps,interp)
c resolution tools
	      if (ao.eq.'s') then  ! SIMULPS
	       call vel3mod(pobs,x,y,z,po,'n')
	       call vel3mod(phit,x,y,z,ph,'n')
	       call vel3mod(pdws,x,y,z,pd,'n')
	       call vel3mod(pres,x,y,z,pr,'n')
	       call vel3mod(spred_p,x,y,z,sprp,'n')
	       call vel3mod(sres,x,y,z,sr,'n')
               call vel3mod(sobs,x,y,z,so,'n')
               call vel3mod(sdws,x,y,z,sd,'n')
	       call vel3mod(spred_vpvs,x,y,z,sprs,'n')
	      else  ! FATOMO on inversion grid
	       call getinvcell(x,y,z,ii,jj,kk) ! get inversion cell
	       po=pobs(ii,jj,kk)
	       ph=phit(ii,jj,kk)
	       pd=pdws(ii,jj,kk)
	       pr=pres(ii,jj,kk)
	       sprp=spred_p(ii,jj,kk)
	      endif
	      call km2deg(x,y,lat,lon)
c cutoff of results that are unreliable
c            if(pr.ge.cutoff)write(3,2000)x,y,lon,lat,z,v,po,pd,pr,
c     &                        sprp,pv,dv,vps,pvps,sr,so,sd,sprs
c or take all values
	     write(3,2000)x,y,lon,lat,z,v,pv,dv,ph,po,pd,pr,sprp,vps,pvps,
     .             so,sd,sr,sprs
	    enddo
	  enddo
	  close(3)
   	 enddo
	else  ! planes are input by user
	 write(6,*)' '
	 write(6,*)' Horizontal depth sections at user defined'
	 write(6,*)' seismic grid'
	 write(6,*)' Please specify no. of grid points in depth (max=10):'
	 read(5,*)nop
	 write(6,*)' For each depth section specify depth in km'
	 write(6,*)'  (hit return after each value)'
	 do i=1,nop
	    read(5,*)pz(i)
	 enddo
c now define boundary of horizontal section
	  write(6,*)' '
	  write(6,*)'Do you want to plot just a part of the model between'
	  write(6,*)'xmin,xmax and ymin,ymax [n] ? '
	  write(6,*)'Otherwise the entire model will be plotted.'
	  read(5,'(a1)')zoom
	  if (zoom.eq.' ') zoom='n'
	  if (zoom.eq.'y') then
	     write(6,*)'  please enter xmin,xmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)xmin,xmax
	     write(6,*)'  please enter ymin,ymax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)ymin,ymax
	  else
	     xmin=sei_xp(nxfrst)
	     xmax=sei_xp(nx)
	     ymin=sei_yp(nyfrst)
	     ymax=sei_yp(ny)
	  endif
	 do iz=1,nop
	  z=pz(iz)
	  write(depth,'(f6.1)')z
	  do l=1,6
	    if(depth(l:l).eq.' ')depth(l:l)='0'
	  enddo
	  write(6,*)' VP/plane'//depth//'.xyz'
	  open(3,file='VP/plane'//depth//'.xyz')
c add one descriptive headerline above
          if(aref.eq.'2') then	! dummy values for reference velo 
	     rv = 99.99
	     rvps = 99.99
          else 
             rv = refvp(1,1,iz)
             rvps = refvpvs(1,1,iz)
          endif
            write(3,1999)iz,z,rv,rvps,cvref
c headerline in output to identify the values easylier
	  write(3,'(a,a)')
     & '# x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,',
     &'presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s'
	  do y=ymin,ymax,xystep
	    do x=xmin,xmax,xystep
	      call vel3mod(vpm,x,y,z,v,interp)
	      call vel3mod(dvp,x,y,z,dv,interp)
	      call vel3mod(pvp,x,y,z,pv,interp)
	      call vel3mod(vpvs,x,y,z,vps,interp)
	      call vel3mod(pvpvs,x,y,z,pvps,interp)
c resolution tools
	      if (ao.eq.'s') then  ! SIMULPS
	       call vel3mod(pobs,x,y,z,po,'n')
	       call vel3mod(phit,x,y,z,ph,'n')
	       call vel3mod(pdws,x,y,z,pd,'n')
	       call vel3mod(pres,x,y,z,pr,'n')
	       call vel3mod(spred_p,x,y,z,sprp,'n')
	       call vel3mod(sres,x,y,z,sr,'n')
               call vel3mod(sobs,x,y,z,so,'n')
               call vel3mod(sdws,x,y,z,sd,'n')
	       call vel3mod(spred_vpvs,x,y,z,sprs,'n')
	      else  ! FATOMO on inversion grid
	       call getinvcell(x,y,z,ii,jj,kk) ! get inversion cell
	       po=pobs(ii,jj,kk)
	       ph=phit(ii,jj,kk)
	       pd=pdws(ii,jj,kk)
	       pr=pres(ii,jj,kk)
	      endif
	      call km2deg(x,y,lat,lon)
c cutoff of results that are unreliable
c            if(pr.ge.cutoff)write(3,2000)x,y,lon,lat,z,v,po,pd,pr,
c     &                        sprp,pv,dv,vps,pvps,sr,so,sd,sprs
c or take all values
	     write(3,2000)x,y,lon,lat,z,v,pv,dv,ph,po,pd,pr,sprp,vps,pvps,
     .             so,sd,sr,sprs
	    enddo
	  enddo
	  close(3)
   	 enddo
        endif
c
	goto 999
	endif
cfh @@@
 2000	format(2(f7.2,1x),2(f9.4,1x),f6.2,2x,f5.2,1x,2(f7.2,1x),3(f9.1,1x)
     .     ,f7.4,1x,f8.3,2x,f5.2,1x,f7.2,1x,2(f8.1,1x),f7.4,1x,f8.3)
c
c
 999	continue

c
c regriding of vertical latitude sections centered around the origin line 
c
        if(av.eq.'y')then
	write(6,*)' '
	write(6,*)'-- calculation of latitude depth sections--'
	write(6,*)'latitude depth sections are centered around the origin'
	write(6,*)'and spaced at user defined distance (in degrees)'
	write(6,*)' ' 
        write(6,*)' DISTANCE BETWEEN DEPTH SECTIONS (in degrees)'
csh 02/10/04 changes to allow distance to be in deg or km
csh  modifications were originally done by Steve Sherburn
	  write(6,*)' measured in degrees (d) or km (k):'
	  read(5,'(a1)')dk
csh end change
	  write(6,*)' DISTANCE BETWEEN DEPTH SECTIONS '
        write(6,*)' along constant latitude (1 float):'
        read(5,*)v_dist
        write(6,*)' NUMBER OF SECTIONS (odd number)'
        read(5,*)v_no
c  calculation of vertical nodes and output to plot the nodes later on
	open(88,file='VP/latvert_nodes.xyz',err=9999)
	open(89,file='VP/latvert_invblo.xyz')
	write(6,*)' '
	write(6,*)'writing vertical nodes to file'
	write(6,*)' VP/latvert_nodes.xyz'
	write(6,*)' VP/latvert_invblo.xyz'
c
 	y=0
	if (ao.eq.'s') then
	 if (ab.eq.'n') then  ! do not output boundary nodes
	   nxfrst=2    ! first grid node is the second one (SIMULPS)
	   nx=sei_nx-1
	   nyfrst=2
	   ny=sei_ny-1
	 else  ! inlcude boundary nodes
	   nxfrst=1    
	   nx=sei_nx
	   nyfrst=1
	   ny=sei_ny
	 endif	   
	else  ! FATOMO
	   nxfrst=1
	   nx=sei_nx
	   nzfrst=1
	   nz=sei_nz
	endif
c
	do j=nzfrst,nz
	  do k=nxfrst,nx
	    x=sei_xp(k)
	    z=sei_zp(j)
	    call km2deg(x,y,lat,lon)
	    write(88,'(2f9.4,1x)')lon,z
	  enddo
	enddo
c for inversion grid write out block boundaries (FATOMO)
	if (ao.eq.'f') then
	 do k=1,inv_nz
	    x=inv_xp(k,1)
	    z=inv_zp(k+1)
	    call km2deg(x,y,lat,lon)
	    write(89,*)'>'
	    write(89,'(3f9.3,1x)')lon,x,z
	    x=inv_xp(k,inv_nx(k))
	    call km2deg(x,y,lat,lon)
	    write(89,'(3f9.3,1x)')lon,x,z
	 enddo
	 do k=1,inv_nz
	  do i=1,inv_nx(k)
	    x=inv_xp(k,i)
	    z=inv_zp(k)
	    call km2deg(x,y,lat,lon)
	    write(89,*)'>'
	    write(89,'(3f9.3,1x)')lon,x,z
	    z=inv_zp(k+1)
	    call km2deg(x,y,lat,lon)
	    write(89,'(3f9.3,1x)')lon,x,z
	  enddo
	 enddo
	endif
	close(88)
	close(89)
c now define boundary of horizontal section
	  write(6,*)' '
	  write(6,*)'Do you want to plot just a part of the model between'
	  write(6,*)'xmin,xmax and zmin,zmax [n] ? '
	  write(6,*)'Otherwise the entire model will be plotted.'
	  read(5,'(a1)')zoom
	  if (zoom.eq.' ') zoom='n'
	  if (zoom.eq.'y') then
	     write(6,*)'  please enter xmin,xmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)xmin,xmax
	     write(6,*)'  please enter zmin,zmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)zmin,zmax
	  else
	     xmin=sei_xp(nxfrst)
	     xmax=sei_xp(nx)
	     zmin=sei_zp(nzfrst)
	     zmax=sei_zp(nz)
	  endif
c  output model dimension on screen
	write(6,*)' '
	write(6,*)' Dimension of Depth Section: '
	write(6,*)'   xmin = ',xmin,' xmax = ',xmax
	write(6,*)'   zmin = ',zmin,' zmin = ',zmax
c  now calculate latitude depth section
	do i=1,v_no
	  open(8,file='VP/latsec'//cn(i)//'.xyz')
	  write(6,*)'VP/latsec'//cn(i)//'.xyz'
csh 02/10/04 change for distance between latitude depth sections to be in deg or km
csh   see comments above
	  if (dk.eq.'d') then
	    y=(i-(aint(v_no/2)+1))*latkm*v_dist
	  else
	    y=(i-(aint(v_no/2)+1))*v_dist
	  endif
csh end change
cfh additional descriptive headerline
          write(8,1888)y,cvref
 1888     format('#> cross section along fixed y-coordinate, y = ',f7.2,
     .          4x,a20)
csh headerline in output to identify the values easier
	  write(8,'(a,a)')
     & '# x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,',
     & 'presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s'
	  do z=zmin,zmax,zstep  !  for depth sections use entire model
	    do x=xmin,xmax,xystep
	      call vel3mod(vpm,x,y,z,v,interp)
	      call vel3mod(dvp,x,y,z,dv,interp)
	      call vel3mod(pvp,x,y,z,pv,interp)
	      call vel3mod(vpvs,x,y,z,vps,interp)
	      call vel3mod(pvpvs,x,y,z,pvps,interp)
c resolution tools
	      if (ao.eq.'s') then  ! SIMULPS
	       call vel3mod(pobs,x,y,z,po,'n')
	       call vel3mod(phit,x,y,z,ph,'n')
	       call vel3mod(pdws,x,y,z,pd,'n')
	       call vel3mod(pres,x,y,z,pr,'n')
	       call vel3mod(spred_p,x,y,z,sprp,'n')
	       call vel3mod(sres,x,y,z,sr,'n')
               call vel3mod(sobs,x,y,z,so,'n')
               call vel3mod(sdws,x,y,z,sd,'n')
	       call vel3mod(spred_vpvs,x,y,z,sprs,'n')
	      else  ! FATOMO on inversion grid
	       call getinvcell(x,y,z,ii,jj,kk) ! get inversion cell
	       po=pobs(ii,jj,kk)
	       ph=phit(ii,jj,kk)
	       pd=pdws(ii,jj,kk)
	       pr=pres(ii,jj,kk)
	      endif
	      call km2deg(x,y,lat,lon)
c cut off values which are unreliable	      
c	      if(pr.ge.cutoff)write(8,2000)x,y,lon,lat,z,v,po,pd,pr,
c     &          	      sprp,dv,pv,vps,pvps,sr,so,sd,sprs
c or take all values
	     write(8,2000)x,y,lon,lat,z,v,pv,dv,ph,po,pd,pr,sprp,vps,pvps,
     .             so,sd,sr,sprs
	    enddo
	  enddo
	  close(8)
	enddo
	endif
c
c regriding of vertical longitude sections centered around the origin line
c
        if(al.eq.'y')then
	write(6,*)' '
	write(6,*)'-- calculation of longitude depth sections--'
	write(6,*)'longitude depth sections are centered around the origin'
	write(6,*)'and spaced at user defined distance (in degrees)' 
	write(6,*)' '
        write(6,*)'DISTANCE BETWEEN DEPTH SECTIONS (in degrees)'
csh 02/10/04 changes to allow distance to be in deg or km
csh  modifications were originally done by Steve Sherburn
	  write(6,*)' measured in degrees (d) or km (k):'
	  read(5,'(a1)')dk
csh end change
	  write(6,*)' DISTANCE BETWEEN DEPTH SECTIONS '
        write(6,*)' along constant longitude (1 float):'
        read(5,*)v_dist
        write(6,*)'NUMBER OF SECTIONS:'
        read(5,*)v_no
c  calculation of vertical nodes and output to plot the nodes later on
	open(88,file='VP/lonvert_nodes.xyz',err=9999)
	open(89,file='VP/lonvert_invnodes.xyz')
	write(6,*)' '
	write(6,*)'writing vertical nodes to file'
	write(6,*)' VP/lonvert_nodes.xyz'
	write(6,*)' VP/lonvert_invnodes.xyz'
c
 	x=0
	if (ao.eq.'s') then
	 if (ab.eq.'n') then  ! do not output boundary nodes
	   nxfrst=2    ! first grid node is the second one (SIMULPS)
	   nx=sei_nx-1
	   nyfrst=2
	   ny=sei_ny-1
	 else   ! include boundary nodes
	   nxfrst=1    
	   nx=sei_nx
	   nyfrst=1
	   ny=sei_ny
	 endif	   
	else  ! FATOMO
	   nyfrst=1
	   ny=sei_ny
	   nzfrst=1
	   nz=sei_nz
	endif
c
	do j=nzfrst,nz
	  do k=nyfrst,ny
	    y=sei_yp(k)
	    z=sei_zp(j)
	    call km2deg(x,y,lat,lon)
	    write(88,'(2f9.4,1x)')lat,z
	  enddo
	enddo
c inversion edges (only FATOMO)
	if (ao.eq.'f') then
	 do k=1,inv_nz
	    y=inv_yp(k,1)
	    z=inv_zp(k+1)
	    call km2deg(x,y,lat,lon)
	    write(89,*)'>'
	    write(89,'(3f9.4,1x)')lat,y,z
	    y=inv_yp(k,inv_ny(k))
	    call km2deg(x,y,lat,lon)
	    write(89,'(3f9.4,1x)')lat,y,z
	 enddo
	 do k=1,inv_nz
	  do i=1,inv_ny(k)
	    y=inv_yp(k,i)
	    z=inv_zp(k)
	    call km2deg(x,y,lat,lon)
	    write(89,*)'>'
	    write(89,'(3f9.4,1x)')lat,y,z
	    z=inv_zp(k+1)
	    call km2deg(x,y,lat,lon)
	    write(89,'(3f9.4,1x)')lat,y,z
	  enddo
	 enddo
	endif
	close(88)
	close(89)
c now define boundary of horizontal section
	  write(6,*)' '
	  write(6,*)'Do you want to plot just a part of the model between'
	  write(6,*)'ymin,ymax and zmin,zmax [n] ? '
	  write(6,*)'Otherwise the entire model will be plotted.'
	  read(5,'(a1)')zoom
	  if (zoom.eq.' ') zoom='n'
	  if (zoom.eq.'y') then
	     write(6,*)'  please enter ymin,ymax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)ymin,ymax
	     write(6,*)'  please enter zmin,zmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)zmin,zmax
	  else
	     ymin=sei_yp(nyfrst)
	     ymax=sei_yp(ny)
	     zmin=sei_zp(nzfrst)
	     zmax=sei_zp(nz)
	  endif
c now the sections
	do i=1,v_no
	  open(10,file='VP/lonsec'//cn(i)//'.xyz')
	  write(6,*)'VP/lonsec'//cn(i)//'.xyz'
csh 02/10/04 change for distance between latitude depth sections to be in deg or km
csh   see comments above
	  if (dk.eq.'d') then
	    x=(i-(aint(v_no/2)+1))*lonkm*v_dist
	  else
	    x=(i-(aint(v_no/2)+1))*v_dist
	  endif
csh end change
cfh additional descriptive header
          write(10,1110)x,cvref
 1110     format('#> cross section along fixed x-coordinate, x = ',f7.2,
     .           4x,a20)
csh headerline in output to identify the values easylier
	  write(10,'(a,a)')
     & '# x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,',
     & 'presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s'
	  do z=zmin,zmax,zstep  !  for depth sections use entire model
	    do y=ymin,ymax,xystep
	      call vel3mod(vpm,x,y,z,v,interp)
	      call vel3mod(dvp,x,y,z,dv,interp)
	      call vel3mod(pvp,x,y,z,pv,interp)
	      call vel3mod(vpvs,x,y,z,vps,interp)
	      call vel3mod(pvpvs,x,y,z,pvps,interp)
c resolution tools
	      if (ao.eq.'s') then  ! SIMULPS
	       call vel3mod(pobs,x,y,z,po,'n')
	       call vel3mod(phit,x,y,z,ph,'n')
	       call vel3mod(pdws,x,y,z,pd,'n')
	       call vel3mod(pres,x,y,z,pr,'n')
	       call vel3mod(spred_p,x,y,z,sprp,'n')
	       call vel3mod(sres,x,y,z,sr,'n')
               call vel3mod(sobs,x,y,z,so,'n')
               call vel3mod(sdws,x,y,z,sd,'n')
	       call vel3mod(spred_vpvs,x,y,z,sprs,'n')
	      else  ! FATOMO on inversion grid
	       call getinvcell(x,y,z,ii,jj,kk) ! get inversion cell
	       po=pobs(ii,jj,kk)
	       ph=phit(ii,jj,kk)
	       pd=pdws(ii,jj,kk)
	       pr=pres(ii,jj,kk)
	      endif
	      call km2deg(x,y,lat,lon)
c cut off values which are unreliable
c	      if(pr.ge.cutoff)write(10,2000)x,y,lon,lat,z,v,po,pd,pr,dv,
c    &           pv,vps,pvpvs,sr,so,sd,sprs
c or take all values
             write(10,2000)x,y,lon,lat,z,v,pv,dv,ph,po,pd,pr,sprp,vps,
     &             pvps,so,sd,sr,sprs
	    enddo
	  enddo
	  close(10)
	enddo
	endif
c
c  calculate profiles
c
        if(as.eq.'y')then
	 write(6,*)' '
	 write(6,*)'----calculating user defined profiles ----'
	 write(6,*)' ! endpoint of profile will be adjusted if'
	 write(6,*)'   profile lenght is not a multiple of spacing'
	 write(6,*)'   of interpolation points along profile !' 
	 write(6,*)' '
         write(6,*)'HOW MANY PROFILES ?'
	 read(5,*)profno
c
	 if (ao.eq.'s') then
	  if (ab.eq.'n') then  ! do not output boundary nodes
	   nxfrst=2    ! first grid node is the second one (SIMULPS)
	   nx=sei_nx-1
	   nyfrst=2
	   ny=sei_ny-1
	  else   ! include boundary nodes
	   nxfrst=1    
	   nx=sei_nx
	   nyfrst=1
	   ny=sei_ny
	  endif	   
	 else
	    nxfrst=1
	    nx=sei_nx
	    nyfrst=1
	    ny=sei_ny
	    nzfrst=1
	    nz=sei_nz
	 endif
c
         do i=1,profno
	  write(6,*)' '
	  write(6,*)'--working on profile ',i,'--'
	  write(6,*)' '
	  open(9,file='VP/profile'//cn(i)//'.xyz',err=9999)
	  write(6,*)'spacing of interpolation points along profile:'
	  write(6,*)'  (for SIMULPS must be in bld-units)'
	  read(5,*)ddist
          write(6,*)'start point of profile (lon lat) in GMT-coord'
          read(5,*)xll(1,i),yll(1,i)
          write(6,*)'end point of profile (lon lat) in GMT-coord'
          read(5,*)xll(2,i),yll(2,i)
c now define boundary 
	  write(6,*)' '
	  write(6,*)'Do you want to plot just a part of the model between'
	  write(6,*)'xmin,xmax and zmin,zmax [n] ? '
	  write(6,*)'Otherwise the entire model will be plotted.'
	  read(5,'(a1)')zoom
	  if (zoom.eq.' ') zoom='n'
	  if (zoom.eq.'y') then
	     write(6,*)'  please enter xmin,xmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)xmin,xmax
	     write(6,*)'  please enter zmin,zmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)zmin,zmax
	  else
	     xmin=sei_xp(nxfrst)
	     xmax=sei_xp(nx)
	     zmin=sei_zp(nzfrst)
	     zmax=sei_zp(nz)
	  endif
c
	  do j=1,2
	    call deg2km(yll(j,i),xll(j,i),xsh(j,i),ysh(j,i))
	  enddo
          write(9,1988)xll(1,i),yll(1,i),xll(2,i),yll(2,i),cvref
 1988     format('#> prof. start: ',2(f9.4,2x),'  prof. end: ',
     &    2(f9.4,2x),' (probably adjusted) ',4x,a20)
	  write(9,'(a,a,a)')
     & '# in(=km along profile),x,y,lon,lat,z,vp,pvp(%),dvp(absol.),',
     & 'phit,pkhit,pdws,presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,',
     & 'sresol,spread_s'
          proflnght=sqrt((xsh(2,i)-xsh(1,i))**2+(ysh(2,i)-ysh(1,i))**2)
	  intpoints=proflnght/ddist
	  xs=(xsh(2,i)-xsh(1,i))/intpoints
	  ys=(ysh(2,i)-ysh(1,i))/intpoints
	  intpoints=aint(intpoints)+1
	  write(6,*)'profile lenght = ',proflnght,' no. points = ',
     &            intpoints
	  write(6,*)'delta_x= ',xs,'  delta_y= ',ys
          write(6,*)' '
c get starting and ending point if plotting just part of the model
	  if (zoom.eq.'y') then
	     k1=aint(xmin/ddist)+1
	     k2=aint(xmax/ddist)+1
	  else
	     k1=1
	     k2=intpoints
	  endif
c
	  write(6,*)'--writing output to file--'
	  write(6,*)' VP/profile'//cn(i)//'.xyz'
c
	  do z=zmin,zmax,zstep
            do k=k1,k2
	      x=(k-1)*xs+xsh(1,i)
	      y=(k-1)*ys+ysh(1,i)
	      if (x.ge.sei_xp(nxfrst).and.x.le.sei_xp(nx).and.
     &            y.ge.sei_yp(nyfrst).and.y.le.sei_yp(ny)) then
		 call vel3mod(vpm,x,y,z,v,interp)
		 call vel3mod(dvp,x,y,z,dv,interp)
		 call vel3mod(pvp,x,y,z,pv,interp)
		 call vel3mod(vpvs,x,y,z,vps,interp)
		 call vel3mod(pvpvs,x,y,z,pvps,interp)
c resolution tools
		 if (ao.eq.'s') then ! SIMULPS
		    call vel3mod(pobs,x,y,z,po,'n')
		    call vel3mod(phit,x,y,z,ph,'n')
		    call vel3mod(pdws,x,y,z,pd,'n')
		    call vel3mod(pres,x,y,z,pr,'n')
		    call vel3mod(spred_p,x,y,z,sprp,'n')
		    call vel3mod(sres,x,y,z,sr,'n')
		    call vel3mod(sobs,x,y,z,so,'n')
		    call vel3mod(sdws,x,y,z,sd,'n')
		    call vel3mod(spred_vpvs,x,y,z,sprs,'n')
		 else		! FATOMO on inversion grid
		    call getinvcell(x,y,z,ii,jj,kk) ! get inversion cell
		    po=pobs(ii,jj,kk)
		    ph=phit(ii,jj,kk)
		    pd=pdws(ii,jj,kk)
		    pr=pres(ii,jj,kk)
		 endif
		 call km2deg(x,y,lat,lon)
		 in=sqrt((x-xsh(1,i))**2+(y-ysh(1,i))**2)
		 write(9,4010)in,x,y,lon,lat,z,v,pv,dv,ph,po,pd,pr,sprp,vps,
     .             pvps,so,sd,sr,sprs
	      else ! if x,y are outside model skip
		 write(6,*)' '
		 write(6,4011)x,y,z
	      endif
	    enddo
	  enddo
	  write(6,*)' '
	  write(6,*)'ADJUSTED ENDPOINT of profile (lon lat) :',lon,lat
	  close(9)
 4010   format(3(f7.2,1x),2(f9.4,1x),f6.2,2x,f5.2,1x,2(f7.2,1x),
     .     3(f9.1,1x),f7.4,1x,f8.3,1x,f5.2,1x,f7.2,1x,2(f8.1,1x),f7.4,
     .     1x,f8.3)
 4011	format('point x=',f5.1,' y=',f5.1,' z=',f5.1,' is outside model -> 
     &          skipped!')
	enddo
	endif
c
c
c 1-D section
c
        if(am.eq.'y')then
 200	  write(6,*)'ENTER LAT/LON OF DESIRED 1D-MODEL(f8.4,1x,f8.4),
     &    >RETURN< MEANS ABORT'
	  read(5,'(a8,1x,a8)')alat,alon
	  if(alat.eq.' '.or.alon.eq.' ')goto 300
          write(6,*)'DATA  WILL BE WRITTEN TO
     &    1D_'//alat//'_'//alon//'.DAT'
	  read(alat,'(f8.4)')lat
	  read(alon,'(f8.4)')lon
	  open(20,file='1D_'//alat//'_'//alon//'.dat')
	  call deg2km(lat,lon,x,y)
c	  do z=-5.,150.,5.    !Original by SH
          do z=-5.,50.,1.     !Changed by tdiehl for CH
	     call vel3mod(vpm,x,y,z,v,interp)
	     write(20,2000)x,y,lon,lat,z,v
	  enddo
	  close(20)
	  goto 200
 300	  continue
	  write(6,*)' '
	endif
csh 07/12/00
c
c  insert part to export x,y,z,value file for MATLAB use
c
	if(amat.eq.'y')then
	   write(6,*)' '
	   write(6,*)'--- output MATLAB file  ---'
	   write(6,*)' '
	   write(6,*)' which parameter:?'
	   write(6,*)'  (1) abs. Vp'
	   write(6,*)'  (2) % Vp change'
	   write(6,*)'  (3) KHIT (Vp)'
	   write(6,*)'  (4) DWS (Vp)'
	   write(6,*)'  (5) RDE (Vp)'
	   write(6,*)'  (6) spread (Vp)'
	   write(6,*)'  (7) abs. Vp/Vs'
	   write(6,*)'  (8) % Vp?Vs change'
	   write(6,*)'  (9) KHIT (Vp/Vs)'
	   write(6,*)'  (10) DWS (Vp/Vs)'
	   write(6,*)'  (11) RDE (Vp/Vs)'
	   write(6,*)'  (12) spread (Vp/Vs)'
	   read(5,*)avp
c open corresponding output file
	   if (avp.eq.1) open(3,file='VP/vp3D.xyz')
	   if (avp.eq.2) open(3,file='VP/vpertu3D.xyz')
	   if (avp.eq.3) open(3,file='VP/pkhit3D.xyz')
	   if (avp.eq.4) open(3,file='VP/pdws3D.xyz')
	   if (avp.eq.5) open(3,file='VP/prde3D.xyz')
	   if (avp.eq.6) open(3,file='VP/psprd3D.xyz')
	   if (avp.eq.7) open(3,file='VP/vpvs3D.xyz')
	   if (avp.eq.8) open(3,file='VP/vpvspertu3D.xyz')
	   if (avp.eq.9) open(3,file='VP/skhit3D.xyz')
	   if (avp.eq.10) open(3,file='VP/sdws3D.xyz')
	   if (avp.eq.11) open(3,file='VP/srde3D.xyz')
	   if (avp.eq.12) open(3,file='VP/ssprd3D.xyz')
c enter cutoff for resolution
	   write(6,*)' '
	   write(6,*)' enter cutoff for resolution (RDE): '
	   read(5,*)cutoff
c
	   if (ao.eq.'s') then
	      nxfrst=2  ! first grid node is the 2nd one (SIMULPS)
	      nx=sei_nx-1
	      nyfrst=2
	      ny=sei_ny-1
	      nzfrst=2
	      nz=sei_nz-1
	   else ! FATOMO
	      nxfrst=1
	      nx=sei_nx
	      nyfrst=1
	      ny=sei_ny
	      nzfrst=1
	      nz=sei_nz
	   endif
c
 	   write(6,*)' '
	   write(6,*)'Do you want to plot just a part of the model between'
	   write(6,*)'xmin,xmax and zmin,zmax [n] ? '
	   write(6,*)'Otherwise the entire model will be plotted.'
	   read(5,'(a1)')zoom
	   if (zoom.eq.' ') zoom='n'
	   if (zoom.eq.'y') then
	     write(6,*)'  please enter xmin,xmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)xmin,xmax
	     write(6,*)'  please enter ymin,ymax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)ymin,ymax	     
	     write(6,*)'  please enter zmin,zmax [in km]:'
	     write(6,*)'   (must be multiple of grid spacing)'
	     read(5,*)zmin,zmax
	   else
	     xmin=sei_xp(nxfrst)
	     xmax=sei_xp(nx)
	     ymin=sei_yp(nyfrst)
	     ymax=sei_yp(ny)
	     zmin=sei_zp(nzfrst)
	     zmax=sei_zp(nz)
	   endif
c
	   do z=zmin,zmax,zstep
	      do y=ymin,ymax,xystep
		 do x=xmin,xmax,xystep
		    call vel3mod(vpm,x,y,z,v,interp)
		    call vel3mod(pvp,x,y,z,pv,interp)
		    call vel3mod(vpvs,x,y,z,vps,interp)
		    call vel3mod(pvpvs,x,y,z,pvps,interp)
c resolution tools
		 if (ao.eq.'s') then ! SIMULPS
		    call vel3mod(phit,x,y,z,ph,'n')
		    call vel3mod(pdws,x,y,z,pd,'n')
		    call vel3mod(pres,x,y,z,pr,'n')
		    call vel3mod(spred_p,x,y,z,sprp,'n')
		    call vel3mod(sres,x,y,z,sr,'n')
		    call vel3mod(sobs,x,y,z,so,'n')
		    call vel3mod(sdws,x,y,z,sd,'n')
		    call vel3mod(spred_vpvs,x,y,z,sprs,'n')
		 else		! FATOMO on inversion grid
		    call getinvcell(x,y,z,ii,jj,kk) ! get inversion cell
		    po=pobs(ii,jj,kk)
		    ph=phit(ii,jj,kk)
		    pd=pdws(ii,jj,kk)
		    pr=pres(ii,jj,kk)
		 endif
		 call km2deg(x,y,lat,lon)
c output depending on which parameter
		 if (avp.eq.1) then
		    if (pr.ge.cutoff) write(3,2222)lon,lat,x,y,z,v
		    if (pr.lt.cutoff) write(3,2223)lon,lat,x,y,z,'NaN'
		 endif
		 if (avp.eq.2) then
		    if (pr.ge.cutoff) write(3,2222)lon,lat,x,y,z,pv
		    if (pr.lt.cutoff) write(3,2223)lon,lat,x,y,z,'0'
		 endif
		 if (avp.eq.3) write(3,2222)lon,lat,x,y,z,ph
		 if (avp.eq.4) write(3,2222)lon,lat,x,y,z,pd
		 if (avp.eq.5) write(3,2222)lon,lat,x,y,z,pr
		 if (avp.eq.6) write(3,2222)lon,lat,x,y,z,sprp
		 if (avp.eq.7) write(3,2222)lon,lat,x,y,z,vps
		 if (avp.eq.8) write(3,2222)lon,lat,x,y,z,pvps
		 if (avp.eq.9) write(3,2222)lon,lat,x,y,z,so
		 if (avp.eq.10) write(3,2222)lon,lat,x,y,z,sd
		 if (avp.eq.11) write(3,2222)lon,lat,x,y,z,sr
		 if (avp.eq.12) write(3,2222)lon,lat,x,y,z,sprs		 
c		    call km2deg(x,y,lat,lon)
c		    if (avp.eq.'p') then ! Vp
c		       call vel3mod(vpm,x,y,z,v,interp)
c		       call vel3mod(pres,x,y,z,pr,'n')
c		       call vel3mod(pdws,x,y,z,pd,'n')
c		       if (pd.lt.hitct) v=0    ! mask unresolved values
c		       write(3,2222)lon,lat,z,v
c		    else ! VpVs
c		       call vel3mod(vpvs,x,y,z,v,interp)
c		       call vel3mod(sdws,x,y,z,sd,'n')
c		       if (sd.lt.hitct) v=0    ! mask unresolved values
c		       write(3,2222)x,y,lon,lat,z,v
c		    endif
		 enddo
	      enddo
	    enddo
	    close(3)
	endif
 2222	format(2f10.4,3f8.1,f7.2)
 2223	format(2f10.4,3f8.1,a3)
csh end
	goto 7777
c
c error messages
c
 9988	continue  ! error on opening file 'output'
	write(6,*)' '
	write(6,*)'***  error on opening file main output file ***'
	write(6,*)'  Did you run tomo2gmt in the directory where'
	write(6,*)'  you copied the file "output" ?'
	stop
 9999	continue ! error on opening output files of tomo2gmt
	write(6,*)' '
	write(6,*)'***  error on opening output file of tomo2gmt ***'
	write(6,*)'  Did you create the subdirectory VP ?'
	stop
c
 7777	continue
	write(6,*)'FINISHED TOMO2GMT'
c	
	end
c
c######### end main ######################################################
c
c
      subroutine getinvcell(x,y,z,ii,jj,kk)
c
c  for a given x,y,z find ii,jj,kk of invcell
c
	include 'tomo2gmt_common.inc'
c
	integer ii,jj,kk
	real x,y,z
c
	do 300 k=1,inv_nz
	   if (z.ge.inv_zp(k).and.z.lt.inv_zp(k+1)) then
	      kk=k
	      goto 301
	   endif
 300	continue
 301	continue
	do 100 i=1,inv_nx(kk)-1
	   if (x.ge.inv_xp(kk,i).and.x.lt.inv_xp(kk,i+1)) then
	      ii=i
	      goto 101
	   endif
 100	continue
 101	continue
	do 200 j=1,inv_ny(kk)-1
	   if (y.ge.inv_yp(kk,j).and.y.lt.inv_yp(kk,j+1)) then
	      jj=j
	      goto 201
	   endif
 200	continue
 201	continue
c
	return
	end
c
c ------------------------------------------------------------------
c
      subroutine vel3mod(vel,x,y,z,val,interp)
c 
c  Interpolation routine taken from SIMULPS. It has been completed
c  to perform a nearest neighbor interpolation apart from the standard
c  trilinear velocity interpolation.
c
c
        implicit none
	include 'tomo2gmt_common.inc'
c
      integer i,ix,iy,iz,ip,jp,kp,ip1,jp1,kp1,wvno
      real vel(nxmax,nymax,nzmax)
      real val,x,y,z,xf,yf,zf,xf1,yf1,zf1,wv(8),wvmax
      character*1 interp

c
      do ix=1,sei_nx-1
        if(x.ge.sei_xp(ix).and.x.lt.sei_xp(ix+1))then
          ip=ix
	  goto 100
        endif
      enddo
 100  continue
      do iy=1,sei_ny-1
        if(y.ge.sei_yp(iy).and.y.lt.sei_yp(iy+1))then
          jp=iy
	  goto 200
        endif
      enddo
 200  continue
      do iz=1,sei_nz-1
        if(z.ge.sei_zp(iz).and.z.lt.sei_zp(iz+1))then
          kp=iz
	  goto 300
        endif
      enddo
 300  continue
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c	write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
c100	format(3(2f7.3,i3))
c
      xf=(x-sei_xp(ip))/(sei_xp(ip1)-sei_xp(ip))
      yf=(y-sei_yp(jp))/(sei_yp(jp1)-sei_yp(jp))
      zf=(z-sei_zp(kp))/(sei_zp(kp1)-sei_zp(kp))
c
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c 
c  calculating weigthing for each neighbor
c
c
      wv(1)=xf1*yf1*zf1
      wv(2)=xf*yf1*zf1
      wv(3)=xf1*yf*zf1
      wv(4)=xf*yf*zf1
      wv(5)=xf1*yf1*zf
      wv(6)=xf*yf1*zf
      wv(7)=xf1*yf*zf
      wv(8)=xf*yf*zf
c
c finding maximum if nearest neighbor interpolation is wanted
c
      if (interp.eq.'n') then
       wvmax=wv(1)
       wvno=1
       do i=2,8
	if (wv(i).gt.wvmax) then
	  wvmax=wv(i)
	  wv(i-1)=0
	  wv(wvno)=0
	  wvno=i
        else
          wv(i)=0
	endif
       enddo
       wv(wvno)=1
      endif
c
c  calculate interpolated velocity
c
      val=wv(1)*vel(ip,jp,kp)+wv(2)*vel(ip1,jp,kp)
     & +wv(3)*vel(ip,jp1,kp)+wv(4)*vel(ip1,jp1,kp)
     & +wv(5)*vel(ip,jp,kp1)+wv(6)*vel(ip1,jp,kp1)
     & +wv(7)*vel(ip,jp1,kp1)+wv(8)*vel(ip1,jp1,kp1)
      return
c
c***** end of subroutine vel3mod *****
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	subroutine km2deg(x,y,lat,lon)
c
cfh 980804 added rotation of coords, sine and cosine of rotation
cfh        angle are in snr and csr. km2deg assumes input x,y coords
cfh        are probably rotated and converts back to unrotated fx,fy before
cfh        converting to degrees. as simul (and sim2gmt, for that) has
cfh        x (=lon) positive towards left (=W), the rotation matrix is inverted.
cfh        with respect to general mathematics.
cfh        as x,y are in units of bld, and the conversion expects km, one has
cfh        to account for that also
csh 28/03/01
csh        x,y are not in units of bld, but in km!!
c
        implicit none
	include 'tomo2gmt_common.inc'
c
c	real latkm,lonkm,olat,olon,snr,csr
	real  x,y,lat,lon,fx,fy
c	common/conv/ latkm,lonkm,olat,olon,snr,csr,bld
c
c convert x,y to unrotated coordinates and get km out of it (multiply with bld)
c	fx=(snr*y+csr*x) * bld
c	fy=(csr*y-snr*x) * bld
        fx=(snr*y+csr*x)
	fy=(csr*y-snr*x)
csh	character*1 hemi
c
csh        if (hemi.eq.'s') then
csh	 lat=-1*(y/latkm-olat)
csh        else
csh         lat=-1*(-y/latkm-olat)
csh	endif
	if (olat.lt.0) then
	   lat=fy/latkm+olat
	else
	   lat=-1*(-fy/latkm-olat)
	endif
	lon=-1.*(fx/lonkm-olon)
c	lon=-1.*(x/latkm/cos(lat*3.141592654/180.)-olon)
c
	return
c
	end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	subroutine deg2km(lat,lon,x,y)
c
cfh     see comments on subroutine km2deg
c       here input degrees are first converted to unrotated
c       coordinates and then rotated with snr, csr
c       bld is applied to get x,y in units of bld instead of km
csh 28/03/01
c      x,y are in km and not in bld-units!!!
        implicit none
	include 'tomo2gmt_common.inc'
c
	real lat,lon,x,y,fx,fy
c	real latkm,lonkm,olat,olon,snr,csr,bld
c	common/conv/ latkm,lonkm,olat,olon,snr,csr,bld
csh	character*1 hemi
c
csh	if (hemi.eq.'s') then
csh	 y=-1*(latkm*(lat-olat))
csh	else
csh         y=latkm*(lat-olat)
csh	endif
	if (olat.lt.0) then
	   fy=latkm*(lat+olat)
	else
	   fy=latkm*(lat-olat)
	endif
	fx=-1*lonkm*(lon-olon)
c
c now rotate coordinates
c        y = (csr*fy+snr*fx) / bld
c	x = (csr*fx-snr*fy) / bld
	y = (csr*fy+snr*fx)
	x = (csr*fx-snr*fy)
c
	return
	end
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	subroutine read_simulps(iuses,aref,ab,hitct,prog_ver)
c
c  this subroutines reads grid, final velocities, and 
c  resolution estimates out of the SIMULPS13Q output file
c
	implicit none
	include 'tomo2gmt_common.inc'
	integer     iuses
	real        hitct
	character*1 aref,ab
	integer     prog_ver
c
	integer   i,j,k,nlay,itd
	integer   ix,iy,iz
	integer   ilat,ilon
	real      rfvp,rfvpvs
	real      lat,lon,dlat,dlon,dep
	real      x,y,z
	character line*132,extmod*32,ns*1,ew*1
c
 10	  read(2,'(a132)',end=9999)line
 	  write(6,*)line
c
c  reading grid planes 
c
	  if(line(1:20).eq.' velocity grid size:')then
	    write(6,*)' '
	    write(6,*)'READING GRID PLANES'
c
	    write(6,*)' '
	    write(6,*)'do you want to include (fixed) boundary nodes [n]'
	    read(5,'(a1)')ab
	    if (ab.eq.' ') ab='n'
c
	    read(2,1010)bld,sei_nx,sei_ny,sei_nz
c 	    write(6,'(f4.1i3i3i3)')bld,sei_nx,sei_ny,sei_ny
 1010	    format(t7,f3.2,t21,i2,t33,i2,t45,i2)
 1020	    format(3x,20f7.1)
            nzout = sei_nz
cfh 1212   cope with lost layers ...
	    write(6,'(a26,i2)')'nz from grid definition = ',nzout
	    write(6,*)'Simulps does not output layers with no hits,'
	    write(6,*)'so if the highest layer number printed in your'
	    write(6,*)'output has less than nz-1 layers please'
	    write(6,*)'specify your highest layer number (else return):'
	    read(5,'(i2)')nlay
	    if (nlay.ne.0) nzout = nlay+1
csh 1402
	    write(6,*)'if you our output does not start with first layer'
	    write(6,*)' please enter no. of first plane (else return)'
	    read(5,'(i2)')nzfrst
	    if (nzfrst.eq.0) nzfrst=2
c  
	    write(6,*)'number of output layers:',nzout-1
	    write(6,*)'starting with layer: ',nzfrst
c
	    read(2,'(/a132)')

ctd         Test if reading routine is compatable with grid dimensions:
            if(sei_nx.gt.60) then
               write(*,'(a,a)')
     &         '--> number of nodes in X not supported by read-routine',
     &         ' - adjust read routine to read more than 60 grid nodes'
               stop
            endif
            if(sei_ny.gt.60) then
               write(*,'(a,a)')
     &         '--> number of nodes in Y not supported by read-routine',
     &         ' - adjust read routine to read more than 60 grid nodes'
               stop
            endif
            if(sei_nz.gt.40) then
               write(*,'(a,a)')
     &         '--> number of nodes in Z not supported by read-routine',
     &         ' - adjust read routine to read more than 60 grid nodes'
               stop
            endif


ctd      increased for nx>40, before only ny>40 considered...
cfh 0807 if-block for nx,ny > 40
            if (sei_nx.gt.40) then
               read(2,1020)(sei_xp(k),k=1,20)
               read(2,'(/a132)')
               read(2,1020)(sei_xp(k),k=21,40)
               read(2,'(/a132)')
               read(2,1020)(sei_xp(k),k=41,sei_nx)
            else
	     if (sei_nx.gt.20) then
	        read(2,1020)(sei_xp(k),k=1,20)
	        read(2,'(/a132)')
	        read(2,1020)(sei_xp(k),k=21,sei_nx)
             else 
	        read(2,1020)(sei_xp(k),k=1,sei_nx)	      
             endif
            endif
            do itd=1,sei_nx
               write(*,*)'x:',sei_xp(itd)
            enddo
	    read(2,'(/a132)')
csh 2303 if-block for ny >40
	    if (sei_ny.gt.40) then
	       read(2,1020)(sei_yp(k),k=1,20)
	       read(2,'(/a132)')
	       read(2,1020)(sei_yp(k),k=21,40)
	       read(2,'(/a132)')
	       read(2,1020)(sei_yp(k),k=41,sei_ny)
	    else
	     if (sei_ny.gt.20) then
	       read(2,1020)(sei_yp(k),k=1,20)	       
	       read(2,'(/a132)') 
	       read(2,1020)(sei_yp(k),k=21,sei_ny)
	     else
	       read(2,1020)(sei_yp(k),k=1,sei_ny)	       
             endif
	    endif
            do itd=1,sei_ny
               write(*,*)'y:',sei_yp(itd)
            enddo
chf 0807 end
csh 0514 if block for nz>20
	    if (sei_nz.gt.20) then
	       read(2,'(//a132)')line
	       read(line,1020)(sei_zp(k),k=1,20)
	       read(line,1021)(azp(k),k=1,20)
	       read(2,'(///a132)')line
	       read(line,1020)(sei_zp(k),k=21,sei_nz)
	       read(line,1021)(azp(k),k=21,sei_nz)
	    else
	       read(2,'(//a132)')line
	       read(line,1020)(sei_zp(k),k=1,sei_nz)
	       read(line,1021)(azp(k),k=1,sei_nz)
	    endif
csh
	    write(6,1021)(azp(k),k=1,sei_nz)
c
c azp(k) contains in character the original depth nodes as
c  zp(k) contains this in real format 
 1021	    format(3x,20a7)
	    goto 14
	  else
	    goto 10
	  endif
c
c reading initial velocities
c
 14	 if (aref.eq.'1') then   ! minimum 1D input model
	  read(2,'(a132)',end=9999)line
	  if (line(2:19).eq.'velocity values on') then  
		write(6,*)' '
		write(6,*)'READING INITIAL VELOCITY REFVP(X,Y,Z)'
		write(6,*)' '
		do iz=1,sei_nz
	       	 read(2,'(a132)',end=9999)line ! read 2 lines until input starts
		 read(2,'(a132)',end=9999)line
		  do iy=1,sei_ny
c                  write(*,*)sei_nx
c                  This has been extended to read sei_nx<=60 
c                  Before it was limited to 40
		   if (sei_nx.gt.40) then
		        read(2,'(f6.2)')rfvp
c                       write(*,*)rfvp
			read(2,'(a132)')line
                        read(2,'(a132)')line
                   else
                      if (sei_nx.gt.20) then
                          read(2,'(f6.2)')rfvp
                          read(2,'(a132)')line
		      else
		          read(2,'(f6.2)')rfvp
		      endif
                   endif
		   do ix=1,sei_nx
		     refvp(ix,iy,iz)=rfvp
		   enddo
		  enddo
		 enddo
csh  now for vpvs-ratio 
		read(2,'(a132)',end=9999)line
		read(2,'(a132)',end=9999)line
		if (iuses.ne.0) then
	 	 do iz=1,sei_nz
		  do iy=1,sei_ny
		    if (sei_nx.gt.20) then
			read(2,'(f6.2)')rfvpvs
			read(2,'(a132)')line
		    else
			read(2,'(f6.2)')rfvpvs
		    endif
		    do ix=1,sei_nx
		      refvpvs(ix,iy,iz)=rfvpvs
	            enddo
		  enddo
	       	  read(2,'(a132)',end=9999)line ! read 2 lines until input starts
		  read(2,'(a132)',end=9999)line
		 enddo
	        endif
	       goto 15
	  else
	   goto 14
	  endif
	 endif
	 if (aref.eq.'2') then  ! 2D/3D input model
	  read(2,'(a132)',end=99)line
	  if (line(2:19).eq.'velocity values on'.or.
     &        line(1:17).eq.' initial velocity') then
		write(6,*)' '
		write(6,*)'READING INITIAL VELOCITY REFVP(X,Y,Z)'
		write(6,*)' '
		do iz=1,sei_nz
	       	 read(2,'(a132)',end=99)line ! read 2 lines until input starts		 
		 read(2,'(a132)',end=99)line
		 do iy=1,sei_ny
		    if (sei_nx.gt.20) then
			read(2,1030)(refvp(ix,iy,iz),ix=1,20)
			read(2,1030)(refvp(ix,iy,iz),ix=21,sei_nx)
		    else
			read(2,1030)(refvp(ix,iy,iz),ix=1,sei_nx)
		    endif
		 enddo
		enddo
c  now for vpvs-ratio 
		read(2,'(a132)',end=99)line
		read(2,'(a132)',end=99)line
		if (iuses.ne.0) then
		do iz=1,sei_nz
		  do iy=1,sei_ny
		    if (sei_nx.gt.20) then
			read(2,1030)(refvpvs(ix,iy,iz),ix=1,20)
			read(2,1030)(refvpvs(ix,iy,iz),ix=21,sei_nx)
		    else
			read(2,1030)(refvpvs(ix,iy,iz),ix=1,sei_nx)
		    endif
		  enddo
	       	 read(2,'(a132)',end=99)line ! read 2 lines until input starts
		 read(2,'(a132)',end=99)line
		enddo
		endif
	       goto 15
	   else
		goto 14
	   endif
	 endif
c read in external 1D velocity model
	if (aref.eq.'4') then
	 write(6,*)' '
	 write(6,*)'enter filename (32 char. max) of external 1D model'
	 write(6,*)' file must consist of nz lines, one value per line (Vp,Vp/Vs (opt.))'
	 read(5,'(a32)')extmod
	 open(11,file=extmod)
	 do iz=1,sei_nz
	  read(11,'(f5.2)',end=15)rfvp
	  do iy=1,sei_ny
	   do ix=1,sei_nx
	    refvp(ix,iy,iz)=rfvp
	   enddo
	  enddo
	 enddo
c now VpVs 
	 do iz=1,sei_nz
	  read(11,'(f5.2)',end=15)rfvpvs
	  do iy=1,sei_ny
	   do ix=1,sei_nx
	    refvpvs(ix,iy,iz)=rfvpvs
	   enddo
	  enddo
	 enddo
	endif
c
c reading final hypocenters (only for SIMULPS)
cfh same change as for stations - read N,S ; E,W for lat,lon
c
 15	  read(2,'(a132)',end=99)line
csh 02/10/04 change to cope with simul2000
 1025     format(23x,2(i3,a1,f6.2,1x),f6.2)
 1026     format(22x,2(i3,a1,f6.2,1x),f6.2)
csh 092500 change in SIMULPS14 of eq-no to i4, so change here 1025
c 1025     format(23x,2(i3,a1,f6.2,1x),f6.2)
	   if(line(16:30).eq.'FINAL LOCATIONS')then
	    write(6,*)' '
	    write(6,*)'READING HYPOCENTERS'
	    open(4,file='VP/hypo.xyz')
	    read(2,'(a132)')line
	    read(2,'(a132)')line
	    do i=1,20000
	      read(2,'(a132)')line
	      if(line.ne.' '.and.line(1:6).ne.'RMSALL')then
	        if (prog_ver.eq.14) then
	          read(line,1025)ilat,ns,dlat,ilon,ew,dlon,dep
	        else ! simul2000
	          read(line,1026)ilat,ns,dlat,ilon,ew,dlon,dep
	        endif
              if (ns.eq.'S'.or.ns.eq.'s') then
                  lat=-1.*(ilat+dlat/60.)
              else
	          lat=(ilat+dlat/60.)
              endif
              if (ew.eq.'W'.or.ew.eq.'w') then
	          lon=-1.*(ilon+dlon/60.)
              else
	          lon=(ilon+dlon/60.)
	        endif
csh	      write(4,'(2f10.4,f10.2)')lon,lat,dep
csh 082300 now output x,y,lon,lat,dep to hypo.xyz
	      call deg2km(lat,lon,x,y)
	      write(4,'(2f8.2,2f10.4,f8.2)')x,y,lon,lat,dep
	      else
	        goto 20
	      endif
	    enddo
	    close(4)
	    goto 20
	  else
	    goto 15
	  endif
c
c now read in Vp,KHIT,DWS,RDE 
c
c reading final P-velocity
c
 20	  read(2,'(a132)',end=9999)line
 1030	  format(20f6.2)
	  if(line(1:20).eq.' FINAL P-VELOCITY MO')then
	    write(6,*)' '
	    write(6,*)'READING FINAL VELOCITY VP(X,Y,Z)'
	    do iz=nzfrst,nzout-1
csh 02/12/02 if all Vp nodes are fixed no output is written
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      if (line(1:20).eq.' OBSERVATION MATRIX ') then
		 do k=nzfrst,nzout-1
		    do j=1,sei_ny
		       do i=1,sei_nx
			  vpm(i,j,k)=refvp(i,j,k) ! set final velo to initial velos
		       enddo
		    enddo
		  enddo
		  goto 35
	       endif
	      do iy=1,sei_ny
		if (sei_nx .gt. 20) then
		   read(2,1030)(vpm(ix,iy,iz),ix=1,20)
		   read(2,1030)(vpm(ix,iy,iz),ix=21,sei_nx)
                else
	          read(2,1030)(vpm(ix,iy,iz),ix=1,sei_nx)
                endif
	      enddo
	    enddo
csh 1502
c  set velocities of planes not printed in the output to initial values
c   this is necessary for depth sections
	    if (nzfrst.ne.2) then
	       do iz=2,nzfrst-1
		  do iy=1,sei_ny
		     do ix=1,sei_nx
			vpm(ix,iy,iz)=refvp(ix,iy,iz)
		     enddo
		  enddo
	       enddo
	    endif
	    if (nlay.ne.0) then
	       do iz=nzout,sei_nz-1
		  do iy=1,sei_ny
		     do ix=1,sei_nx
			vpm(ix,iy,iz)=refvp(ix,iy,iz)
		     enddo
		  enddo
		enddo
	    endif
	    goto 30
	  else
	    goto 20
	  endif
c
c reading numbers of rays passed by a mode (KHIT)
c
 30	  read(2,'(a132)',end=9999)line
1040	  format(21f6.0)
1041      format(T7,20f6.0)
 35	  if(line(1:20).eq.' OBSERVATION MATRIX ')then
	    write(6,*)' '
	    write(6,*)'READING KHIT MATRIX POBS(X,Y,Z)'
	    pobsmax=0
	    do iz=nzfrst,nzout-1
csh 02/12/02 if all Vp nodes are fixed no output is written
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      if (line(1:20).eq.' DERIVATIVE WEIGHT S') goto 45
	      do iy=2,sei_ny-1
                if (sei_nx.gt.21) then
		 read(2,1040)(pobs(ix,iy,iz),ix=1,21)
		 read(2,1041)(pobs(ix,iy,iz),ix=22,sei_nx)
                else
	         read(2,1040)(pobs(ix,iy,iz),ix=1,sei_nx)
                endif
		do ix=1,sei_nx
		  if(pobs(ix,iy,iz).gt.pobsmax)pobsmax=pobs(ix,iy,iz)
		enddo
	      enddo
c set values corresponding to boundary nodes to zero
	      do ix=1,sei_nx
	         pobs(ix,1,iz)=777
		 pobs(ix,sei_ny,iz)=777
	      enddo			      
	    enddo	    
	    goto 40
	  else
	    goto 30
	  endif
c
c reading DWS
c
 40	  read(2,'(a132)',end=9999)line
csh 09082004
c    change in Simulps14 to output format of DWS
 1050	  format(f3.0,18f8.0)
 1051	  format(3x,18f8.0)
c 1050	  format(f3.0,18f7.0)
c 1051	  format(3x,18f7.0)
csh 02/19/04
c   in simul2000 output format of DWS has been changed
 1052	  format(f3.0,20f7.0)
 1053    format(3x,20f7.0)
 45	  if(line(1:20).eq.' DERIVATIVE WEIGHT S')then
	    write(6,*)' '
	    write(6,*)'READING DWS MATRIX PDWS(X,Y,Z)'
	    pdwsmax=0.0
	    open(77,file='VP/fixnodes.dat')
	    open(78,file='VP/fixnodes.xyz')
	    do iz=nzfrst,nzout-1
csh 02/12/02 if all Vp nodes are fixed no output is written
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      if (line(1:20).eq.' RESOLUTION : GRIDOI') goto 55
	      do iy=2,sei_ny-1
	       if (prog_ver.eq.14.or.prog_ver.eq.13) then
	        if (sei_nx.gt.19) then
	           read(2,1050)(pdws(ix,iy,iz),ix=1,19)
	           read(2,1051)(pdws(ix,iy,iz),ix=20,sei_nx)
	        else 
	          read(2,1050)(pdws(ix,iy,iz),ix=1,sei_nx)	 
                       endif
                      else ! simul2000
                       if (sei_nx.gt.21) then
	           read(2,1052)(pdws(ix,iy,iz),ix=1,21)
	           read(2,1053)(pdws(ix,iy,iz),ix=22,sei_nx)
	        else 
	          read(2,1052)(pdws(ix,iy,iz),ix=1,sei_nx)	 
                       endif
                      endif
	       do ix=1,sei_nx
	         if(pdws(ix,iy,iz).gt.pdwsmax) pdwsmax=pdws(ix,iy,iz)
c  do not output x-edge nodes
	         if (ix.ne.1.and.ix.ne.sei_nx) then
	            if(pdws(ix,iy,iz).le.hitct) then
		       write(77,'(3i3)')ix,iy,iz
csh output x,y,z,lon,lat of fixed node
		       x=sei_xp(ix)
		       y=sei_yp(iy)
		       z=sei_zp(iz)
		       call km2deg(x,y,lat,lon)
		       write(78,'(3f7.1,1x,2(f9.4,1x))')x,y,z,lon,lat
	            endif
	           endif
	        enddo
	      enddo
c set values corresponding to boundary nodes to zero
	      do ix=1,sei_nx
	         pdws(ix,1,iz)=777
		 pdws(ix,sei_ny,iz)=777
	      enddo			       
	    enddo
	    close(77)
	    close(78)
	    goto 50
	  else
	    goto 40
	  endif
c
c  reading resolution diagonal elements
c
 50	  read(2,'(a132)',end=9999)line
 1060	  format(1x,20(6x,f7.4))
 1061	  format(20(6x,f7.4))
csh chnange for simul2000
 1065   format(1x,13(6x,f7.4))
 1066   format(20(6x,f7.4))
 55	  if(line(1:20).eq.' RESOLUTION : GRIDOI')then
	    write(6,*)' '
	    write(6,*)'READING RESOLUTION DIAGONAL ELEMENTS PRES(X,Y,Z)'
	    presmax=0.0
	    do iz=nzfrst,nzout-1
csh 02/12/02 if all Vp nodes are fixed no output is written
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      if (line.eq.' ') goto 60
	      do iy=2,sei_ny-1
csh change for simul2000
	       if (prog_ver.eq.14.or.prog_ver.eq.13) then
	        if (sei_nx-2.gt.20) then
	           read(2,1060)(pres(ix,iy,iz),ix=2,21)
	           read(2,1061)(pres(ix,iy,iz),ix=22,sei_nx-1)
	        else
	           read(2,1060)(pres(ix,iy,iz),ix=2,sei_nx-1)
	        endif
	       else ! simul2000
	        if (sei_nx-2.gt.13) then
	         read(2,1065)(pres(ix,iy,iz),ix=2,13)
	         read(2,1066)(pres(ix,iy,iz),ix=14,sei_nx-1)
                       else
	         read(2,1065)(pres(ix,iy,iz),ix=2,sei_nx-1)
                       endif
                      endif
c set values corresponding to boundary nodes to -1.0
	        pres(1,iy,iz)=-1.0
	        pres(sei_nx,iy,iz)=-1.0
c				
		do ix=1,sei_nx
		  if(pres(ix,iy,iz).gt.presmax)presmax=pres(ix,iy,iz)
		enddo
	      enddo
c set values corresponding to boundary nodes to -1.0
	      do ix=1,sei_nx
	         pres(ix,1,iz)=-1.0
		 pres(ix,sei_ny,iz)=-1.0
	      enddo			            
	    enddo
csh 2302
c  compute trace of the resolution matrix
c
	  trace=0
	  do iz=nzfrst,nzout-1
	     do iy=2,sei_ny-1
		do ix=2,sei_nx-1
		 if(pres(ix,iy,iz).gt.0.0) trace=trace+pres(ix,iy,iz)
		enddo
	     enddo
	  enddo
	    goto 60
	  else
	    goto 50
	  endif
c
cfh next section only if iuses = 1
 60	 if (iuses.eq.1) then
c
c  reading final Vp/Vs model
c
 	  read(2,'(a132)',end=99)line
	  if(line(1:20).eq.' FINAL S-VELOCITY MO')then
	    write(6,*)' '
	    write(6,*)'READING VP/VS RATIO VPVS(X,Y,Z)'
	    write(6,*)' '
	    do iz=nzfrst,nzout-1
	      read(2,'(/a132)')line
	      if(line(1:5).eq.' OBSE')then
	        write(6,*)'NO S MODEL FOUND'
	        do iy=1,sei_ny
	          do ix=1,sei_nx
	            vpvs(ix,iy,iz)=1.73
	          enddo
	        enddo
	        goto 99
	      else
	        do iy=1,sei_ny
		if (sei_nx.gt.20) then
		   read(2,1030)(vpvs(ix,iy,iz),ix=1,20)
		   read(2,1030)(vpvs(ix,iy,iz),ix=21,sei_nx)
                else
	          read(2,1030)(vpvs(ix,iy,iz),ix=1,sei_nx)	
                endif
	        enddo
		read(2,'(a132)')line
		read(2,'(a132)')line
	        do iy=1,sei_ny
	          read(2,'(a132)')line
		  if (sei_nx.gt.20) then 
		   read(2,'(a132)')line
		  endif
	        enddo
	      endif
	    enddo
	    goto 70
	  else
	    goto 60
	  endif
c
c   reading KHIT for Vp/Vs
c
 70       read(2,'(a132)',end=9999)line
          if(line(1:20).eq.' OBSERVATION MATRIX ')then
            write(6,*)'READING SOBS(X,Y,Z)'
            write(6,*)' '
            sobsmax=0
            do iz=nzfrst,nzout-1
              read(2,'(/a132)')line
              do iy=2,sei_ny-1
	if (sei_nx.gt.21) then
	  read(2,1040)(sobs(ix,iy,iz),ix=1,21)
	  read(2,1041)(sobs(ix,iy,iz),ix=22,sei_nx)
                else
	  read(2,1040)(sobs(ix,iy,iz),ix=1,sei_nx)	
                endif
                do ix=1,sei_nx
                  if(sobs(ix,iy,iz).gt.sobsmax)sobsmax=sobs(ix,iy,iz)
		  if(sobs(ix,iy,iz).lt.hitct)write(77,'(3i3)')ix,iy,iz+sei_nz
                enddo
              enddo
            enddo
            goto 80
          else
            goto 70
          endif
c
c  reading DWS for Vp/Vs
c
 80       read(2,'(a132)',end=9999)line
          if(line(1:20).eq.' DERIVATIVE WEIGHT S')then
            write(6,*)'READING SDWS(X,Y,Z)'
            write(6,*)' '
            do iz=nzfrst,nzout-1
              read(2,'(/a132)')line
              do iy=2,sei_ny-1
csh 02/19/2004 change for simul2000
	  if (prog_ver.eq.14.or.prog_ver.eq.13) then
	      if (sei_nx.gt.19) then
	           read(2,1050)(sdws(ix,iy,iz),ix=1,19)
	           read(2,1051)(sdws(ix,iy,iz),ix=20,sei_nx)
	      else 
	          read(2,1050)(sdws(ix,iy,iz),ix=1,sei_nx)	 
                     endif
                  else ! simul2000
                      if (sei_nx.gt.21) then
	           read(2,1052)(sdws(ix,iy,iz),ix=1,21)
	           read(2,1053)(sdws(ix,iy,iz),ix=22,sei_nx)
	       else 
	          read(2,1052)(sdws(ix,iy,iz),ix=1,sei_nx)	 
                      endif
                  endif
                do ix=1,sei_nx
                  if(sdws(ix,iy,iz).gt.sdwsmax)sdwsmax=sdws(ix,iy,iz)
                enddo
              enddo
            enddo
            goto 90
          else
            goto 80
          endif
c
c  reading reso. diag. element for Vp/Vs
c
 90       read(2,'(a132)',end=9999)line
          if(line(1:20).eq.' RESOLUTION : GRIDOI')then
            write(6,*)'READING SRES(X,Y,Z)'
            write(6,*)' '
            sresmax=0.0
            do iz=nzfrst,nzout-1
              read(2,'(/a132)')line
              do iy=2,sei_ny-1
csh change for simul2000
	 if (prog_ver.eq.14.or.prog_ver.eq.13) then
	    if (sei_nx-2.gt.20) then
	         read(2,1060)(sres(ix,iy,iz),ix=2,21)
	         read(2,1061)(sres(ix,iy,iz),ix=22,sei_nx-1)
	    else
	         read(2,1060)(sres(ix,iy,iz),ix=2,sei_nx-1)
	    endif
	  else ! simul2000
	    if (sei_nx-2.gt.13) then
	         read(2,1065)(sres(ix,iy,iz),ix=2,13)
	         read(2,1066)(sres(ix,iy,iz),ix=14,sei_nx-1)
                   else
	         read(2,1065)(sres(ix,iy,iz),ix=2,sei_nx-1)
                   endif
                endif
c
                do ix=1,sei_nx
                  if(sres(ix,iy,iz).gt.sresmax)sresmax=sres(ix,iy,iz)
                enddo
              enddo
            enddo
            goto 99
          else
            goto 90
          endif
c
         endif		! end of if block for vp/vs
c
c
 99	  continue
	  return  
 9999	  continue
	  write(6,*)' WARNING:  CHECk OUTPUT OF SIMULPS'
c
	  end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
	subroutine read_fatomo(aref)
c
c  this subroutines reads grid, final velocities, and 
c  resolution estimates out of the FATOMO output file
c
	implicit none
	include 'tomo2gmt_common.inc'
	character aref*1
c
	integer   i,k,ix,iy,iz
	integer   ilat,ilon
	real     dlat,dlon,lat,lon,dep,x,y
	real      rfvp
	character*1 ew,ns
	character line*132,extmod*32
c
c
c  reading grid planes 
c
 10	read(2,'(a132)',end=9999)line
	if(line(2:15).eq.'INVERSION GRID')then
	    read(2,'(a132)')line
	    write(6,*)' '
	    write(6,*)'READING GRID PLANES'
	    read(2,2010)inv_nz
	    write(6,2010)inv_nz
 2010	    format(3x,'number of layers: ',i3)
	    read(2,2015)(inv_zp(k),k=1,inv_nz+1)
	    write(6,'(20f5.1)')(inv_zp(k),k=1,inv_nz+1)
 2015	    format(25x,20f5.1)
	    goto 11
	 else
	     goto 10
	 endif
 11	 continue
	 read(2,'(a132)',end=9999)line
	 if (line(4:12).eq.'inversion') then
	    do k=1,inv_nz
	       read(2,2016)inv_nx(k),inv_ny(k)
c	       write(6,'(2i3)')inv_nx(k),inv_ny(k)
	       read(2,2017)(inv_xp(k,i),i=inv_nx(k),1,-1)
	       read(2,2017)(inv_yp(k,i),i=1,inv_ny(k))
c	       write(6,2017)(inv_xp(k,i),i=inv_nx(k),1,-1)
c	       write(6,2017)(inv_yp(k,i),i=inv_ny(k),1,-1)
 2016	       format(21x,i2,9x,i2)
 2017	       format(6x,20f6.1)
	    enddo
	    goto 12
	  else
	    goto 11
	  endif
 12	  continue
	  read(2,'(a132)')line
	  if (line(4:20).eq.'no. of gridnodes:') then
c now read seismic grid
 2011	    format(3x,'sei_nx = ',i3,
     &        ' sei_ny = ',i3,' sei_nz = ',i3)
 2020	    format(20f6.1)
	    read(2,2011)sei_nx,sei_ny,sei_nz
	    write(6,2011)sei_nx,sei_ny,sei_nz
	    read(2,'(///a132)')line
cfh 0807 if-block for nx,ny > 20
	    if (sei_nx.gt.20) then
	       read(2,2020)(sei_xp(k),k=sei_nx,sei_nx-19,-1)
	       read(2,'(/a132)')
	       read(2,2020)(sei_xp(k),k=sei_nx-20,1,-1)
c	       write(6,2020)(sei_xp(k),k=sei_nx,1,-1)
            else 
	      read(2,2020)(sei_xp(k),k=sei_nx,1,-1)	      
            endif
	    read(2,'(/a132)')line
	    if (sei_ny .gt. 20) then
	       read(2,2020)(sei_yp(k),k=1,20)
	       read(2,'(/a132)')
	       read(2,2020)(sei_yp(k),k=21,sei_ny)
c	       write(6,2020)(sei_yp(k),k=1,sei_ny)
	    else
	       read(2,2020)(sei_yp(k),k=1,sei_ny)	       
            endif
chf 0807 end
	    if (sei_nz.le.20) then
	     read(2,'(/a132)')
	     read(2,'(a132)')line
	     read(line,2020)(sei_zp(k),k=1,sei_nz)
 	     write(6,2020)(sei_zp(k),k=1,sei_nz)
	     read(line,2021)(azp(k),k=1,sei_nz)
	     goto 14
	    else   ! nz >20
	     read(2,'(/a132)')
	     read(2,'(a132)')line
	     read(line,2020)(sei_zp(k),k=1,20)
 	     write(6,2020)(sei_zp(k),k=1,20)
	     read(line,2021)(azp(k),k=1,20)
	     read(2,'(/a132)')
	     read(2,'(a132)')line
	     read(line,2020)(sei_zp(k),k=21,sei_nz)
 	     write(6,2020)(sei_zp(k),k=21,sei_nz)
	     read(line,2021)(azp(k),k=21,sei_nz)
	    endif
	  else
	    goto 12
	  endif
 2021	  format(20a6)
c
c reading initial velocities
c
 14	 if (aref.eq.'1') then   ! minimum 1D input model
	  read(2,'(a132)',end=9999)line
	  if (line(1:17).eq.' initial velocity') then  
		write(6,*)' '
		write(6,*)'READING INITIAL VELOCITY REFVP(X,Y,Z)'
		write(6,*)' '
		do iz=1,sei_nz
	       	 read(2,'(a132)',end=99)line ! read 2 lines until input starts
		 read(2,'(a132)',end=99)line
		 do iy=1,sei_ny
		   if (sei_nx.gt.20) then
		      read(2,'(f5.2)')rfvp
		      read(2,'(a132)')line
		   else
		      read(2,'(f5.2)')rfvp
		   endif
		   do ix=1,sei_nx
		     refvp(ix,iy,iz)=rfvp
		   enddo
		 enddo
		enddo
	        goto 220
	  else
	   goto 14
	  endif
	 endif
	 if (aref.eq.'2') then  ! 2D/3D input model
	  read(2,'(a132)',end=9999)line
	  if (line(2:19).eq.'velocity values on'.or.
     &        line(1:17).eq.' initial velocity') then
		write(6,*)' '
		write(6,*)'READING INITIAL VELOCITY REFVP(X,Y,Z)'
		write(6,*)' '
		do iz=1,sei_nz
	       	 read(2,'(a132)',end=99)line ! read 2 lines until input starts		 
		 read(2,'(a132)',end=99)line
		  do iy=sei_ny,1,-1
		    if (sei_nx.gt.20) then
			read(2,2030)(refvp(ix,iy,iz),ix=sei_nx,sei_nx-19,-1)
			read(2,2030)(refvp(ix,iy,iz),ix=sei_nx-20,1,-1)
c			write(6,2030)(refvp(ix,iy,iz),ix=sei_nx,1,-1)
		    else
			read(2,2030)(refvp(ix,iy,iz),ix=sei_nx,1,-1)
c			write(6,2030)(refvp(ix,iy,iz),ix=1,sei_nx)
		    endif		    
		  enddo
	       enddo
	       goto 220
	   else
		goto 14
	   endif
	 endif
c read in external 1D velocity model
	if (aref.eq.'4') then
	 write(6,*)' '
	 write(6,*)'enter filename (32 char. max) of external 1D model'
	 write(6,*)' file must consist of nz lines, one value per line (Vp,Vp/Vs (opt.))'
	 read(5,'(a32)')extmod
	 open(11,file=extmod)
	 do iz=1,sei_nz
	  read(11,'(f5.2)',end=220)rfvp
	  do iy=1,sei_ny
	   do ix=1,sei_nx
	    refvp(ix,iy,iz)=rfvp
	   enddo
	  enddo
	 enddo
	endif
c
c
c reading numbers of rays passed by a mode (KHIT)
c
 220	  read(2,'(a132)',end=9999)line
2040	  format(20f9.0)
	  open(77,file='VP/fixnodes.dat')
	  if(line(3:12).eq.'HIT MATRIX')then
	    write(6,*)' '
	    write(6,*)'READING HIT MATRIX POBS(X,Y,Z)'
	    phitmax=0
	    read(2,'(a132)')line
	    do iz=1,inv_nz
	      read(2,'(a132)')
	      read(2,'(a132)')
	      read(2,'(a132)')		! 3 lines to next layer
	      do iy=inv_ny(iz)-1,1,-1
                if (inv_nx(iz).gt.20) then
		 read(2,2040)(phit(ix,iy,iz),ix=inv_nx(iz)-1,inv_nx(iz)-20,-1)
		 read(2,2040)(phit(ix,iy,iz),ix=inv_nx(iz)-21,1,-1)
                else
	         read(2,2040)(phit(ix,iy,iz),ix=inv_nx(iz)-1,1,-1)
                endif
		do ix=1,inv_nx(iz)-1
		  if (phit(ix,iy,iz).gt.phitmax) phitmax=phit(ix,iy,iz)
c  put out nodes which are not hit
		  if(phit(ix,iy,iz).eq.0)write(77,'(3i3)')ix,iy,iz
		enddo
	      enddo
	    enddo
	    goto 225
	  else
	    goto 220
	  endif
c
c reading numbers of rays passed by a mode (KHIT)
c
 225	  read(2,'(a132)',end=9999)line
	  if(line(3:13).eq.'KHIT MATRIX')then
	    write(6,*)' '
	    write(6,*)'READING KHIT MATRIX KHIT(X,Y,Z)'
	    pobsmax=0
	    read(2,'(a132)')line
	    do iz=1,inv_nz
	      read(2,'(a132)')
	      read(2,'(a132)')
	      read(2,'(a132)')		! 3 lines to next layer
	      do iy=inv_ny(iz)-1,1,-1
                if (inv_nx(iz).gt.20) then
		 read(2,2040)(pobs(ix,iy,iz),ix=inv_nx(iz)-1,inv_nx(iz)-20,-1)
		 read(2,2040)(pobs(ix,iy,iz),ix=inv_nx(iz)-21,1,-1)
                else
	         read(2,2040)(pobs(ix,iy,iz),ix=inv_nx(iz)-1,1,-1)
                endif
		do ix=1,inv_nx(iz)-1
		  if(pobs(ix,iy,iz).gt.pobsmax)pobsmax=pobs(ix,iy,iz)
		enddo
	      enddo
	    enddo
	    goto 230
	  else
	    goto 225
	  endif
c
c reading DWS
c
 230	  read(2,'(a132)',end=9999)line
 2050	  format(20f7.0)
	  if(line(3:12).eq.'DERIVATIVE')then
	    write(6,*)' '
	    write(6,*)'READING DWS MATRIX PDWS(X,Y,Z)'
	    pdwsmax=0.0
	    read(2,'(a132)')line
	    do iz=1,inv_nz
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      do iy=inv_ny(iz)-1,1,-1
		if (inv_nx(iz).gt.20) then
	           read(2,2050)(pdws(ix,iy,iz),ix=inv_nx(iz)-1,inv_nx(iz)-20,-1)
	           read(2,2050)(pdws(ix,iy,iz),ix=inv_nx(iz)-21,1,-1)
                else 
	           read(2,2050)(pdws(ix,iy,iz),ix=inv_nx(iz)-1,1,-1)
                endif
		do ix=1,inv_nx(iz)-1
		  if(pdws(ix,iy,iz).gt.pdwsmax)pdwsmax=pdws(ix,iy,iz)
		enddo
	      enddo
	    enddo
	    goto 240
	  else
	    goto 230
	  endif
c
c reading final P-velocity
c
 240	  read(2,'(a132)',end=9999)line
 2030	  format(20f5.2)
	  if(line(3:18).eq.'FINAL VELOCITIES')then
	    write(6,*)' '
	    write(6,*)'READING FINAL VELOCITY VP(X,Y,Z)'
	    read(2,'(a132)')line
	    do iz=1,sei_nz   ! velocities are defined on seismic grid!
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      read(2,'(a132)')line
	      do iy=sei_ny,1,-1
		if (sei_nx .gt. 20) then
		   read(2,2030)(vpm(ix,iy,iz),ix=sei_nx,sei_nx-19,-1)
		   read(2,2030)(vpm(ix,iy,iz),ix=sei_nx-20,1,-1)
                else
	          read(2,2030)(vpm(ix,iy,iz),ix=sei_nx,1,-1)	!original code
                endif
	      enddo
	    enddo
	    goto 250 
	  else
	    goto 240
	  endif
c
c  reading resolution diagonal elements
c
 250	  read(2,'(a132)',end=9999)line
 2060	  format(1x,10(6x,f7.4))
 2061	  format(10(6x,f7.4))
	  if(line(1:20).eq.' RESOLUTION : GRIDOI')then
	    write(6,*)' '
	    write(6,*)'READING RESOLUTION DIAGONAL ELEMENTS PRES(X,Y,Z)'
	    write(6,*)' '
	    presmax=0.0
	    do iz=1,inv_nz
	      read(2,'(//a132)')line
	      do iy=inv_ny(iz)-1,1,-1
		if(inv_nx(iz)-1.gt.10) then
		  read(2,2060)(pres(ix,iy,iz),ix=inv_nx(iz)-1,inv_nx(iz)-10,-1)
		  read(2,2061)(pres(ix,iy,iz),ix=inv_nx(iz)-11,1,-1)
                else
	          read(2,2060)(pres(ix,iy,iz),ix=inv_nx(iz)-1,1,-1)
                endif
		do ix=1,inv_nx(iz)-1
		  if(pres(ix,iy,iz).gt.presmax)presmax=pres(ix,iy,iz)
		enddo
	      enddo
	    enddo
c
c  compute trace of the resolution matrix
c
	    trace=0
	    do iz=1,inv_nz
	     do iy=1,inv_ny(iz)-1
		do ix=1,inv_nx(iz)-1
		 if(pres(ix,iy,iz).gt.0.0) trace=trace+pres(ix,iy,iz)
		enddo
	     enddo
	    enddo
	    goto 300
	  else
	    goto 250
	  endif

c  reading final hypocenter locations and output to file 'hypo.xyz'
  300	  read(2,'(a132)',end=9999)line
 2070     format(24x,i2,a1,f5.2,1x,i3,a1,f5.2,1x,f6.2)
	  if(line(16:30).eq.'FINAL LOCATIONS') then
	   open(4,file='VP/hypo.xyz')
	   write(6,*)' '
	   write(6,*)' READING FINAL HYPOCENTER '
	   read(2,'(a132)',end=9999)line
	   read(2,'(a132)',end=9999)line
	   do i=1,5000
	    read(2,'(a132)',end=9999)line
	    if (line.ne.' ') then
	     read(line,2070)ilat,ns,dlat,ilon,ew,dlon,dep
	     if (ns.eq.'s'.or.ns.eq.'S') then
	      lat=-1*(ilat+dlat/60.)
	     else
	      lat=(ilat+dlat/60.)
	     endif
	     if (ew.eq.'W'.or.ew.eq.'w') then
	      lon=-1*(ilon+dlon/60.)
	     else
	      lon=(ilon+dlon/60.)
	     endif
	     call deg2km(lat,lon,x,y)
	     write(4,'(2f9.4,3f8.2)')lon,lat,x,y,dep
	    else
	     close(4)
	     goto 99
	    endif
	   enddo   
	  else
	    goto 300
	  endif
c
 99	  continue
	  return  
 9999	  continue
	  write(6,*)' WARNING:  CHECK OUTPUT OF FATOMO'
c
	  end
c
c----------------------------------------------------------------
c
	subroutine ray_dens(ap,av,al)
c
c  this subroutine prepares the ray density tensor to plot with GMT
c  
c  ray density tensor is computed by the TOMEK code of Edi Kissling
c  this means to compute ray density tensors you must run the TOMOEK program
c
c  two outputfiles are needed:
c    - raydens.out   ! outputfile of RESOL program
c    - resol.out     ! main print out of RESOL
c
c  this subroutine has been set up to a grid design of TOMEK in which
c  the first two seismic blocks (i.e. upper left corner) of each layer
c  are combined to one inversion block. 
c
c  the file raydens.out contains for each inversion block the three
c  eigenvalues and eigenvectors of the ray density tensor. according to the
c  ratio of the eigenvalues four different quality clases are constructed 
c  resembling the uniformity of the ray illumination within the inversion 
c  block. the projection of the eigenvectors on the xy-plane (for plane sections)
c  or on the x(y)z-plane (for depth sections) are computed to represent the main 
c  ray direction within the inversion block. if two main directions exist (quality
c  class = 2) both projections are computed.
c
c  short distance and grid information are read from the file resol.out.
c
c
c                                                   S. Husen (09/07/99)
c
c
c===============================================================================
c                   DEFINITION OF VARIABLES 
c===============================================================================
c======================
	implicit none
c======================
	include 'tomo2gmt_common.inc'
c
	character*1 ap,av,al
c
	real RADDEG,PI
	parameter(PI=3.14159,RADDEG=180./PI)
	integer blocmax
	parameter(blocmax=10000) ! max. no. of inv. blocks
c
	real rtio1,rtio2
	real x,y
	real dx,dy,dz,v_no,v_dist
	real ev(blocmax,3),vec(3,blocmax,3)
	real projx,projy,projz,projlon,projlat,lat,lon
	real blwray(blocmax)
	real sei_xpb(nxmax),sei_ypb(nymax)
	integer no_invlay,invlay,nxy,nx1
	integer ix,jx,iy,jy,iz,ze,j,k,i,l,m
        integer top_inv(nzmax)
	integer inv_bloc(nzmax)
	integer qual,bloc
	character*6 depth
	character line*132   
	character*1 cew,cns
	character*2 cn(10)
c
 	data cn/'01','02','03','04','05','06','07','08','09','10'/
c
c open input files (as output from tomo programs)
c
C===============================================================================
C           open inputfiles
C===============================================================================
	open(2,file='resol_main.out',status='old',err=800)
	open(4,file='raydens.out',status='old',err=801)

C===============================================================================
C                FIRST READING, FROM FILE 2 (resol.out)
C===============================================================================
c                    ==================================
c               reading origin and short distance conversion
c                    ==================================
c
	write(6,*)'    '
	write(6,*)'READING ORIGIN AND SHORT DISTANCE PARAMETER'
	do j=1,90000
 2	  read(2,'(a132)',end=99,err=802)line
	  if(line(6:21).eq.'NUMBER OF LAYERS')then
	   read(line,'(t28,i2)')sei_nz ! no. of seis. layers
          else
	   goto 2
	  endif
 3        read(2,'(a132)',end=99,err=802)line	  
	  if(line(2:7).eq.'Origin')then
	    read(2,*)olat,cns,olon,cew
	    if (cew.eq.'W') then
	     olon=-1.*olon
	    endif
	    write(6,*)'    '
	    write(6,*)'                       ORIGIN  :',olat,olon
	    write(6,*)'    '
	  else
	    goto 3
	  endif
 4        read(2,'(a132)',end=99,err=802)line	  
	  if(line(11:21).eq.'one min lat')then
	    read(line,'(t25,f6.4)',err=802)latkm
	    read(2,'(t25,f6.4)')lonkm
	    latkm=latkm*60.
	    lonkm=lonkm*60.
	    write(6,*)'    '
	    write(6,*)'                   SHORT DISTANCE CONVERSION :'
	    write(6,*)'                      latkm= ',latkm
	    write(6,*)'                      lonkm= ',lonkm
	    write(6,*)'    '
	  else
	    goto 4
	  endif
c
c                    ==================================
c                   reading grid planes & inversion blocks
c                    ==================================
c
 5	  read(2,'(a132)',end=99,err=803)line
	  if(line(3:18).eq.'NUMBER OF BLOCKS')then
	    write(6,*)' '
	    write(6,*)'           READING GRID PLANES'
	    read(line,'(t56,i3,t77,i3)',err=803)sei_nx,sei_ny
	    write(6,*)'number of blocks in x,y: ',sei_nx,sei_ny
	  else
	    goto 5
	  endif
c
	  nxy=(sei_nx*sei_ny)-1  ! no. of blocks per plane; first two are joined
	  nx1=sei_nx-1           ! no. of blocks in first x-row
c read inv. block boundaries
 61	  read(2,'(a132)',end=99)line
1020      format(10f8.1)
	  if(line(3:16).eq.'X-DIR. BLOCK-B')then
	   ze=NINT((REAL(sei_nx+1)+5)/10)
	   do i=0,ze-1
	    read(2,1020)(sei_xpb(k+i*10),k=1,10)
	   enddo
	  else 
	    goto 61
	  endif
 71	  read(2,'(a132)',end=99)line
	  if(line(3:16).eq.'Y-DIR. BLOCK-B')then
	   ze=NINT((REAL(sei_ny+1)+5)/10)
	   do i=0,ze
	    read(2,1020)(sei_ypb(k+i*10),k=1,10)
	   enddo
	  else 
	    goto 71
	  endif	 
c now read block center 
6	  read(2,'(a132)',end=99)line
	  if(line(3:16).eq.'X-DIR. CENTERS')then
	   read(2,'(a132)')
	   ze=NINT((REAL(sei_nx)+5)/10)
	   do i=0,ze-1
	    read(2,1020)(sei_xp(k+i*10),k=1,10)
	   enddo
	  else 
	    goto 6
	  endif
7	  read(2,'(a132)',end=99)line
	  if(line(3:16).eq.'Y-DIR. CENTERS')then
	   read(2,'(a132)')
	   ze=NINT((REAL(sei_nx)+5)/10)
	   do i=0,ze
	    read(2,1020)(sei_yp(k+i*10),k=1,10)
	   enddo
	  else 
	    goto 7
	  endif
c reading first block per inv. layer
8	 read(2,'(a132)',end=99)line
	 if(line(3:17).eq.'NUMBER OF BLOCK') then
          write(6,*)' '
	  write(6,*)'reading INVERSION BLOCKS'
	  write(6,*)' '
	 else
	  goto 8
	 endif
	 read(2,'(a132)')line
	 read(2,'(a132)')line
	 invlay=1
81	 read(2,'(4x,i6,6x,i4)')inv_bloc(invlay),invlay
	 if (inv_bloc(invlay).ne.0) then
	   invlay=invlay+1
	   goto 81
	 endif
c reading inversion layers
9	  read(2,'(a132)',end=99)line
	  if (line(6:28).eq.'TOT. NR. OF INV. LAYER;')then
	   read(line,'(t39,i2)')no_invlay ! no. of inver. layers
	   write(6,*)' '
	   write(6,*)'number of inv. layers',no_invlay
	   write(6,*)' '
          else
	   goto 9
	  endif
c
 20	  read(2,'(a132)',end=99)line
	  if (line(3:12).eq.'SEIS-LAYER')then
	   read(2,'(a132)')
           k=1
	   do i=1,sei_nz
	    read(2,'(a132)')line
	    read(line,'(t45,f6.2)')sei_zp(i)
	    read(line,'(t45,a6)')azp(i) ! reading depth as character for filename
	    read(2,'(a132)')line
            read(line,'(t24,i3)')invlay
            if (i.eq.1) then  ! first layer
              top_inv(k)=i
              k=k+1
             else
              if (invlay.eq.k) then
                top_inv(k)=i
                k=k+1
              endif
            endif
	   enddo
	   goto 99
	   else
	    goto 20
	   endif
c
	 enddo  ! finish with inputfile 2
 99	 continue
	 close(2)
c
c==============================================================================
c           read RAY DENSITY EIGENVECTORS and VALUES (file 4)
c==============================================================================
c
	 write (6,*)' '
	 write (6,*)'reading RAY DENSITY TENSOR (EW & EV)'
	 write (6,*)' ' 
80	 read(4,*,end=888)bloc,blwray(bloc)
        read(4,'(3(f8.3,1x))') (ev(bloc,k),k=1,3) ! eigenvalues
      	 read(4,'(3(f6.3,1x))') (vec(1,bloc,k),k=1,3) ! 1. eigenvector
      	 read(4,'(3(f6.3,1x))') (vec(2,bloc,k),k=1,3) ! 2. eigenvector
      	 read(4,'(3(f6.3,1x))') (vec(3,bloc,k),k=1,3) ! 3. eigenvector
         bloc=bloc+1
         goto 80
888	 continue
	 close(4)
	 write(6,*)'Number of inversions block: ',bloc
c
c=============================================================================
c           determine quality classes and projection of RAY DENSITY
c=============================================================================
c
	 write(6,*)' '
	 write(6,*)'Value of rtio1 (=high ratio of eigenvalues) for quality class [0.8]: '
	 read(5,'(f3.1)')rtio1
	 if (rtio1.eq.0) then
	    rtio1=0.8
	    write(6,'(f3.1)')rtio1
	 endif
	 write(6,*)'Value of rtio2 (=low ratio of eigenvalues) for quality class [0.4]: '
	 read(5,'(f3.1)')rtio2
	 if (rtio2.eq.0) then
	    rtio2=0.4
	    write(6,'(f3.1)')rtio2
	 endif
	 write(6,*)' '
	 if (ap.eq.'y') then
c
c  for plane sections
c
	 write(6,*)'writing inversion grid to file'
	 write(6,*)'  VP/invgrd.xy'
	 open(3,file='VP/invgrd.xy')
	 do ix=1,sei_nx    ! first two blocks are combined
	  call km2deg(sei_xp(ix)-5,sei_yp(1)+5,lat,lon)
	  write(3,'(2(f8.4,1x))')lon,lat
	  call km2deg(sei_xp(ix)-5,sei_yp(sei_ny)-5,lat,lon)
	  write(3,'(2(f8.4,1x))')lon,lat
	  write(3,'(a1)')'>'
	 enddo
	 do iy=1,sei_ny
	  call km2deg(sei_xp(1)+5,sei_yp(iy)-5,lat,lon)
	  write(3,'(2(f8.4,1x))')lon,lat
	  call km2deg(sei_xp(sei_nx)-5,sei_yp(iy)-5,lat,lon)
	  write(3,'(2(f8.4,1x))')lon,lat
	  write(3,'(a1)')'>'
	 enddo
	 call km2deg(sei_xp(1)-5,sei_yp(2)+5,lat,lon)
	 write(3,'(2(f8.4,1x))')lon,lat
	 call km2deg(sei_xp(1)-5,sei_yp(sei_ny)-5,lat,lon)
	 write(3,'(2(f8.4,1x))')lon,lat
	 close (3)
c
         write(6,*)' '
	 write(6,*)'--creating files for ray density--'
	 k=1    ! first block
         do iz=1,no_invlay 
	  depth=azp(top_inv(iz))  ! sorted after inversion layers
	  do l=1,6
	   if(depth(l:l).eq.' ')depth(l:l)='0'
	  enddo
	  write(6,*)'  VP/rayden'//depth//'.xyz'
	  open(3,file='VP/rayden'//depth//'.xyz') ! open outputfile
	  dx=sei_xp(1)-sei_xp(2)
	  dy=sei_yp(1)-sei_yp(2)
	  do iy=1,sei_ny
	   do ix=1,sei_nx
c now determine 4 clases for ray distribution
	   qual=0
	   if ((ev(k,2)/ev(k,3)).gt.rtio1.and.(ev(k,1)/ev(k,3)).gt.rtio1) then
            qual=1
	    goto 90
	   endif
	   if ((ev(k,2)/ev(k,3)).gt.rtio1.and.(ev(k,1)/ev(k,3)).lt.rtio1) then
	    qual=2
	    goto 90
	   endif
	   if ((ev(k,2)/ev(k,3)).gt.rtio2.and.(ev(k,1)/ev(k,3)).lt.rtio1) then
	    qual=3
	    goto 90
	   endif
	   if ((ev(k,2)/ev(k,3)).lt.rtio2.and.(ev(k,1)/ev(k,3)).lt.rtio2) then
	    qual=4
	   endif
90         continue
c calculate projection of vector onto xy-plane
	   if (qual.eq.2) then
	    do l=2,3     ! two main directions, print both on output
	     projx=sei_xp(ix)+(vec(l,k,1)*dx/2) 
	     projy=sei_yp(iy)+(vec(l,k,2)*dy/2)
	     call km2deg(sei_xp(ix),sei_yp(iy),lat,lon)
	     call km2deg(projx,projy,projlat,projlon)
	     write(3,1050)qual,k,blwray(k),lon,lat,projlon,projlat,
     &       sei_xp(ix),sei_yp(iy),projx,projy
	    enddo
	   else
	    projx=sei_xp(ix)+(vec(3,k,1)*dx/2)
	    projy=sei_yp(iy)+(vec(3,k,2)*dy/2)
	    call km2deg(sei_xp(ix),sei_yp(iy),lat,lon)
	    call km2deg(projx,projy,projlat,projlon)
	    write(3,1050)qual,k,blwray(k),lon,lat,projlon,projlat,
     &      sei_xp(ix),sei_yp(iy),projx,projy
	   endif  !qual.eq.2
	   do m=1,invlay  ! if first block of layer skip next
	      if (k.eq.inv_bloc(m).and.iy.eq.1) goto 667
	   enddo
	   k=k+1  ! next block
1050       format(i2,1x,i4,1x,f8.2,1x,2(f8.4,1x),2(f8.4,1x),4(f7.2,1x))
 667       enddo !ix
	  enddo !iy
	  close(3)
	 enddo !iz
	 endif  ! calculate plane sections
c
c  for depth sections
c     latitude depth sections
c
	 if (av.eq.'y') then   ! compute depth sections
	  write(6,*)' '
	  write(6,*)'-- calculation of latitude depth sections--'
	  write(6,*)'latitude depth sections are centered around the origin'
	  write(6,*)'and spaced in the desired distance (in degrees)'
	  write(6,*)' ' 
          write(6,*)' DISTANCE BETWEEN VERTICAL SECTIONS (in degrees)'
          read(5,*)v_dist
          write(6,*)' NUMBER OF VERTICAL SECTIONS (odd number)'
          read(5,*)v_no
c  create file with borders of inversions blocks
	  open(8,file='VP/latsec_invgrd.xz')
	  iy=1  ! dummy y value
	  do iz=2,no_invlay-1   ! borders in z-direction
	   write(8,'(a1)')'>'  ! seperator for GMT
	   call km2deg(sei_xp(1)+5,sei_yp(iy),lat,lon)
	   write(8,'(2(f8.4,1x))')lon,sei_zp(top_inv(iz))   
	   call km2deg(sei_xp(sei_nx)-5,sei_yp(iy),lat,lon)
	   write(8,'(2(f8.4,1x))')lon,sei_zp(top_inv(iz))
	  enddo
	  do ix=1,sei_nx   ! borders in x-direction
	   write(8,'(a1)')'>'  ! seperator for GMT
	   call km2deg(sei_xp(ix)+5,sei_yp(iy),lat,lon)
	   write(8,'(2(f8.4,1x))')lon,sei_zp(top_inv(1))   
	   call km2deg(sei_xp(ix)+5,sei_yp(iy),lat,lon)
	   write(8,'(2(f8.4,1x))')lon,sei_zp(top_inv(no_invlay))
	  enddo
	  close(8)
c now create latitude depth sections
	  do i=1,v_no
	   open(8,file='VP/latsec.rayden.'//cn(i)//'.xyz')
	   write(6,*)'writing in VP/latsec.rayden.'//cn(i)//'.xyz'
c
	   y=(i-aint(v_no/2+1))*latkm*v_dist   ! y value of section
	   call km2deg(0.0,y,lat,lon)
cfh additional descriptive headerline
          write(8,1888)lat
 1888     format('#> cross section along fixed latitude lat = ',f7.2)
	   do iy=1,sei_ny
	    if (y.lt.sei_ypb(iy).and.y.ge.sei_ypb(iy+1))then
	     jy=iy       ! determine iy value
	     goto 91
	    endif
	   enddo !iy
91	   continue
           do iz=1,no_invlay-1
	    dz=sei_zp(top_inv(iz+1))-sei_zp(top_inv(iz))
	    do ix=1,sei_nx
	     dx=sei_xp(ix+1)-sei_xp(ix)
	     k=nx1+ix+(jy-2)*sei_nx+(iz-1)*nxy ! determine inv. block number
c	     write(6,*) top_inv(iz),k,jy
c now determine 4 clases for ray distribution
	     qual=0
	     if ((ev(k,2)/ev(k,3)).gt.rtio1.and.(ev(k,1)/ev(k,3)).gt.rtio1) then
              qual=1
	      goto 92
	     endif
	     if ((ev(k,2)/ev(k,3)).gt.rtio1.and.(ev(k,1)/ev(k,3)).lt.rtio1) then
	      qual=2
	      goto 92
	     endif
	     if ((ev(k,2)/ev(k,3)).gt.rtio2.and.(ev(k,1)/ev(k,3)).lt.rtio1) then
	      qual=3
	      goto 92
	     endif
	     if ((ev(k,2)/ev(k,3)).lt.rtio2.and.(ev(k,1)/ev(k,3)).lt.rtio2) then
	      qual=4
	      endif
92           continue
c calculate projection of vector onto xz-plane
	     if (qual.eq.2) then
	      do l=2,3     ! two main directions, print both on output
	       projx=sei_xp(ix)+(vec(l,k,1)*dx/2) 
	       projz=sei_zp(top_inv(iz))+dz/2+(vec(l,k,3)*dz/2)
	       call km2deg(sei_xp(ix),sei_yp(jy),lat,lon)
	       call km2deg(projx,sei_yp(jy),lat,projlon)
	       write(8,1050)qual,k,blwray(k),lon,sei_zp(top_inv(iz))+dz/2,
     &         projlon,projz,sei_xp(ix),sei_zp(top_inv(iz))+dz/2,
     &         projx,projz
	      enddo
	     else
	      projx=sei_xp(ix)+(vec(3,k,1)*dx/2)
	      projz=sei_zp(top_inv(iz))+dz/2+(vec(3,k,3)*dz/2)
	      call km2deg(sei_xp(ix),sei_yp(jy),lat,lon)
	      call km2deg(projx,sei_yp(jy),lat,projlon)
	      write(8,1050)qual,k,blwray(k),lon,sei_zp(top_inv(iz))+dz/2,
     &         projlon,projz,sei_xp(ix),sei_zp(top_inv(iz))+dz/2,
     &         projx,projz
	     endif ! qual.eq.2
            enddo  ! ix
	   enddo  !iz
	   close(8)
	  enddo  ! i=1,v_no
	 endif ! create latitude depth sections
c 
c     longitude depth sections
c
	 if (al.eq.'y') then   ! compute depth sections
	  write(6,*)' '
	  write(6,*)'-- calculation of longitude depth sections--'
	  write(6,*)'longitude depth sections are centered around the origin'
	  write(6,*)'and spaced in the desired distance (in degrees)'
	  write(6,*)' ' 
          write(6,*)' DISTANCE BETWEEN VERTICAL SECTIONS (in degrees)'
          read(5,*)v_dist
          write(6,*)' NUMBER OF VERTICAL SECTIONS (odd number)'
          read(5,*)v_no
c  create file with borders of inversions blocks
	  open(8,file='VP/lonsec_invgrd.xz')
	  ix=1  ! dummy x value
	  do iz=2,no_invlay-1   ! borders in z-direction
	   write(8,'(a1)')'>'  ! seperator for GMT
	   call km2deg(sei_xpb(1),sei_ypb(1),lat,lon)
	   write(8,'(2(f8.4,1x))')lat,sei_zp(top_inv(iz))   
	   call km2deg(sei_xpb(1),sei_ypb(sei_ny+1),lat,lon)
	   write(8,'(2(f8.4,1x))')lat,sei_zp(top_inv(iz))
	  enddo
	  do iy=2,sei_ny   ! borders in y-direction
	   write(8,'(a1)')'>'  ! seperator for GMT
	   call km2deg(sei_xpb(1),sei_ypb(iy),lat,lon)
	   write(8,'(2(f8.4,1x))')lon,sei_zp(top_inv(1))   
	   call km2deg(sei_xpb(sei_nx+1),sei_ypb(iy),lat,lon)
	   write(8,'(2(f8.4,1x))')lon,sei_zp(top_inv(no_invlay))
	  enddo
	  close(8)
c now create latitude depth sections
	  do i=1,v_no
	   open(8,file='VP/lonsec.rayden.'//cn(i)//'.xyz')
	   write(6,*)'writing in VP/lonsec.rayden.'//cn(i)//'.xyz'
c
	   x=(i-aint(v_no/2+1))*lonkm*v_dist   ! x value of section
	   call km2deg(x,0.0,lat,lon)
cfh additional descriptive headerline
          write(8,1889)lon
 1889     format('#> cross section along fixed longitude lon = ',f7.2)
	   do ix=1,sei_nx
	    if (x.lt.sei_xpb(ix).and.x.ge.sei_xpb(ix+1))then
	     jx=ix-1       ! determine ix value
	     if (jx.eq.2) jx=1  ! first two blocks in x are combined
	     goto 991
	    endif
	   enddo !ix
 991	   continue
           do iz=1,no_invlay-1
	    dz=sei_zp(top_inv(iz+1))-sei_zp(top_inv(iz))
	    do iy=1,sei_ny
	     dy=sei_yp(iy+1)-sei_yp(iy)
	     k=jx+(iy-1)*sei_nx+(iz-1)*nxy ! determine inv. block number
c	     write(6,*) top_inv(iz),k,jy
c now determine 4 clases for ray distribution
	     qual=0
	     if ((ev(k,2)/ev(k,3)).gt.rtio1.and.(ev(k,1)/ev(k,3)).gt.rtio1) then
              qual=1
	      goto 992
	     endif
	     if ((ev(k,2)/ev(k,3)).gt.rtio1.and.(ev(k,1)/ev(k,3)).lt.rtio1) then
	      qual=2
	      goto 992
	     endif
	     if ((ev(k,2)/ev(k,3)).gt.rtio2.and.(ev(k,1)/ev(k,3)).lt.rtio1) then
	      qual=3
	      goto 992
	     endif
	     if ((ev(k,2)/ev(k,3)).lt.rtio2.and.(ev(k,1)/ev(k,3)).lt.rtio2) then
	      qual=4
	      endif
992           continue
c calculate projection of vector onto yz-plane
	     if (qual.eq.2) then
	      do l=2,3     ! two main directions, print both on output
	       projy=sei_yp(iy)+(vec(l,k,2)*dy/2) 
	       projz=sei_zp(top_inv(iz))+dz/2+(vec(l,k,3)*dz/2)
	       call km2deg(sei_xp(jx),sei_yp(iy),lat,lon)
	       call km2deg(sei_xp(jx),projy,projlat,lon)
	       write(8,1050)qual,k,blwray(k),lat,sei_zp(top_inv(iz))+dz/2,
     &         projlat,projz,sei_yp(iy),sei_zp(top_inv(iz))+dz/2,
     &         projy,projz
	      enddo
	     else
	      projy=sei_yp(iy)+(vec(3,k,2)*dy/2)
	      projz=sei_zp(top_inv(iz))+dz/2+(vec(3,k,3)*dz/2)
	      call km2deg(sei_xp(jx),sei_yp(iy),lat,lon)
	      call km2deg(sei_xp(jx),projy,projlat,lon)
	      write(8,1050)qual,k,blwray(k),lat,sei_zp(top_inv(iz))+dz/2,
     &         projlat,projz,sei_yp(iy),sei_zp(top_inv(iz))+dz/2,
     &         projy,projz
	     endif ! qual.eq.2
            enddo  ! iy
	   enddo  !iz
	   close(8)
	  enddo  ! i=1,v_no
	 endif ! create depth sections	     
c
	 return
c
c==============================================================================
c           ERROR MESSAGES FOR 999 OUTPUT OF SUBROUTINE
c==============================================================================
800     write(6,*)'Error when opening RESOL file'
	stop
801     write(6,*)'Error when opening RAY DENSITY FILE'
	stop
802     write(6,*)'Error when reading SHORT DISTANCE'
	stop
803     write(6,*)'Error when reading BLOCK GEOMETRY'
	stop
c
	end



