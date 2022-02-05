c
c  global variables for tomo2gmt
c
c  12. aug. 1999
c
c
c --------------------------------------------------------------------
c    MAXIMUM NO. OF NODES/CELLS
c
	integer nxmax,nymax,nzmax
	parameter(nxmax=29,nymax=29,nzmax=9) 
	integer mx
	parameter(mx=8000) ! maximal number of nodes
c
c --------------------------------------------------------------------
c
c     SHORT DISTANCE CONVERSION
c
c
	real latkm,lonkm,olat,olon,bld,snr,csr
c
	common/conv/ latkm,lonkm,olat,olon,snr,csr,bld
c
c --------------------------------------------------------------------
c    GRID
c       inv_nx,ny,nz   =  no. of inversion block boundaries
c       sei_nx,ny,nz   =  no. of seismic gridpoints
c       inv_xp(),inv_yp(),inv_zp() = location of inversion block boundaries [km]
c	sei_xp(),sei_yp(),sei_zp() = location of seismic gridpoints [km]
c
	integer inv_nx(nzmax),inv_ny(nzmax),inv_nz
	integer sei_nx,sei_ny,sei_nz
	integer nzfrst,nzout
	real inv_xp(nzmax,nxmax),inv_yp(nzmax,nymax),inv_zp(nzmax)
	real sei_xp(nxmax),sei_yp(nymax),sei_zp(nzmax)
	character azp(nzmax)*6
c
	common/invgrid/ inv_nx,inv_ny,inv_nz,inv_xp,inv_yp,inv_zp
	common/seigrid/ sei_nx,sei_ny,sei_nz,sei_xp,sei_yp,sei_zp,
     &                  nzfrst,nzout,azp
c
c --------------------------------------------------------------------
c    3D ARRAYS 
c
c	refvp()  =  3D matrix of initial p-velocity
c	vpm()    =  3D matrix of final p-velocity  [km/s]
c       dvp()    =  3D matrix of absolut velocity change 
c       pvp()	 =  3D matrix of percentage velocity change
c       phit()   =  3D matrix of HIT count (only FATOMO)
c       pobs()   =  3D matrix of KHIT values for each node of p-wave model
c       pdws()   =  3D matrix of DWS values for each node of p-wave model
c       pres()   =  3D matrix of resol. diag. elements of p-wave model
c	spred_p,s()=  3D matrix of spread function values
c       refvpvs()  =  3D matrix of initial vp/vs ratio
c       vpvs()   =  3D matrix of final vp/vs ratio
c       pvpvs()  =  3D matrix of precentage vp/vs change
c       sobs()   =  3D matrix of KHIT values for each node of vp/vs model
c       sdws()   =  3D matrix of DWS values for each node of vp/vs model
c       sres()   =  pres but for vp/vs model
c
	real vpm(nxmax,nymax,nzmax),vpvs(nxmax,nymax,nzmax)
       	real pres(nxmax,nymax,nzmax),sres(nxmax,nymax,nzmax)
       	real dvp(nxmax,nymax,nzmax),pvp(nxmax,nymax,nzmax)
	real pvpvs(nxmax,nymax,nzmax)
        real pdws(nxmax,nymax,nzmax),sdws(nxmax,nymax,nzmax)
       	real refvp(nxmax,nymax,nzmax),refvpvs(nxmax,nymax,nzmax)
	real spred_p(nxmax,nymax,nzmax),spred_s(nxmax,nymax,nzmax)
	real spred_vpvs(nxmax,nymax,nzmax)
	real pobs(nxmax,nymax,nzmax),sobs(nxmax,nymax,nzmax)
	real phit(nxmax,nymax,nzmax)
c
	common/vpgrid/ vpm,dvp,pvp,refvp,pobs,phit,pdws,pres,spred_p
	common/vpvsgrid/ vpvs,pvpvs,sobs,sdws,sres,spred_s,spred_vpvs,refvpvs 
c
c --------------------------------------------------------------------
c   MAXIMUM VALUES
c
	real pdwsmax,presmax,sdwsmax,sresmax,trace
        integer pobsmax,phitmax,sobsmax
c
	common/maxval/ pdwsmax,presmax,sdwsmax,sresmax,trace,
     &                 pobsmax,phitmax,sobsmax
c
c --------------------------------------------------------------------
c   INTERPOLATED VALUES OF 3D ARRAYS
c
	real v,dv,pv,pd,pr,sprp
	real vps,pvps,sd,sr,sprs
	real po,ph,so	
c
	common/interpo/ v,dv,pv,pd,pr,sprp,vps,pvps,sd,sr,sprs,po,ph,so
