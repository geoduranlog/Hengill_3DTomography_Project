c     PROGRAM make3DsimulMOD.f
c     ========================
c     Generates a SIMULPS14 3D input model from a 1D model with option to add
c     checkerboard velocity pertubations as described in Husen et al. 2004
c     "Evidence for gas and magmatic sources beneath the Yellowstone volcanic field..." 
c
c     T. Diehl 2007/07

      implicit none
      character*500 cline
      character*500 simulhdr1
      character*500 simulhdr2
      character*500 simulhdr3
      character*500 simulhdr4
      character*200 modfi
      character*18  cmnfi
      character*2   modtp
      character*1   chkflg

      integer trimlen
      integer i,j
      integer iz,ix,iy
      integer zgap,ygap,xgap
      integer zlen,ylen,xlen
      integer zflg,yflg,xflg
      integer zpol,ypol,xpol
      integer nx,ny,nz,vlr
      integer pertfir
      integer pertpol
      integer pertlon,pertlat
      integer pertgapH,pertgapV

      real grid,pertamp,chkbrdamp

c     Array:
      real velmod(50,4)
      real tmpvel(300)

c     Initial values:
      cmnfi='make3DsimulMOD.cmn'
      vlr=0
      
      do i = 1,50
         velmod(i,1)=0.0
         velmod(i,2)=0.0
         velmod(i,3)=0.0
         velmod(i,4)=0.0
      enddo

c     Load control parameters:
      open(10,file=cmnfi,status='old',err=900)
      read(10,*,end=901)                                      !Ignore comment
      read(10,'(a)',end=901) simulhdr1                        !Read 1. line of simul-header
      read(10,'(a)',end=901) simulhdr2                        !Read 2. line of simul-header
      read(10,'(a)',end=901) simulhdr3                        !Read 3. line of simul-header
      read(10,'(a)',end=901) simulhdr4                        !Read 4. line of simul-header
      read(10,*,end=901)                                      !Ignore comment
      read(10,'(a)',end=901) modfi                            !Name of model file
      read(10,*,end=901)                                      !Ignore comment
      read(10,'(a)',end=901) modtp                            !Read model type (P or PS)
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)                                      !Ignore comment
      read(10,'(a)',end=901) chkflg                           !Read checkerboard flag (y/n)
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)     pertfir                          !First layer with pertubation (1/2)
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)     pertpol                          !Polarity of first checker board pertubation
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)     pertlon                          !Pertubation length in longitude (number of nodes)
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)     pertlat                          !Pertubation length in latitude  (number of nodes)
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)     pertgapH                         !Horizontal gap
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)     pertgapV                         !Vertical   gap
      read(10,*,end=901)                                      !Ignore comment
      read(10,*,end=901)                                      !Ignore comment

c     Now read numerical grid spaceing, nx, ny, nz from simulhdr1:
      read(simulhdr1,*) grid,nx,ny,nz

c     Control:
      write(*,'(a,i4,a,i4,a,i4)')
     &   'Setup model nx x ny x nz: ',
     &   nx,' x ',ny,' x ',nz

c     Now read 1D velocity model:
11    read(10,'(a)',end=20) cline
      if(cline(1:5).ne.'PLANE') goto 11
      vlr=vlr+1

c     PS with pertubation
      if(chkflg.eq.'y'.AND.
     &   modtp.eq.'PS') then
         read(cline,'(9x,f5.2,11x,f5.1,
     &                3x,f5.2,11x,f5.1)')
     &    velmod(vlr,1),
     &    velmod(vlr,2),
     &    velmod(vlr,3),
     &    velmod(vlr,4)
      endif

c     Only P with pertubation
      if(chkflg.eq.'y'.AND.
     &   modtp.eq.'P') then
         read(cline,'(9x,f5.2,11x,f5.1)')
     &    velmod(vlr,1),
     &    velmod(vlr,2)
      endif

c     Only P without pertubation:
      if(chkflg.ne.'y'.AND.
     &   modtp.eq.'P') then
         read(cline,'(9x,f5.2)')
     &    velmod(vlr,1)
      endif

      goto 11                                                 !Next layer   

20    close(10)

c     Check if model is complete:
      if(nz.ne.vlr) goto 902

c     Control:
      do i = 1,vlr
         write(*,'(a,i2,a,f5.2,a,f5.1,
     &                  a,f5.2,a,f5.1,
     &                  a)')
     &   ' LAYER # ',i,
     &   ' vp: ',velmod(i,1),
     &   ' Pertubation-P: ',velmod(i,2),
     &   '% | vs: ',velmod(i,3),
     &   ' Pertubation-S: ',velmod(i,4),
     &   '% |'
      enddo
      write(*,*)

c     Write SIMULPS model file:
      open(10,file=modfi,status='unknown')
c     Write header:
      write(10,'(a)') simulhdr1(1:trimlen(simulhdr1))
      write(10,'(a)') simulhdr2(1:trimlen(simulhdr2))
      write(10,'(a)') simulhdr3(1:trimlen(simulhdr3))
      write(10,'(a)') simulhdr4(1:trimlen(simulhdr4))
      write(10,'(a)') '  0  0  0'

c-----Loop over horizontal planes for Vp:
      if(pertfir.eq.1) then
         zgap=1
      endif
      if(pertfir.eq.2) then
         zgap=0
      endif
      zpol=pertpol*(-1)
      do iz = 1,vlr

c        Calculate P-wave pertubation amplitude:
         pertamp=(velmod(iz,1)/100)*velmod(iz,2)    

c        Control:
         write(*,'(a,i2,a,f5.2,a,f5.2)')
     &         ' Create vp Layer # ',iz,
     &         ' vp = ',velmod(iz,1),
     &         ' +/- ',pertamp
 
c        Initial pertubation for this layer:
         chkbrdamp=0.0

c        Check checkerboard-Z:
c        (verical gaps)
         if(zgap.gt.pertgapV) then
            zflg=1                                            !If gap-counter for Z is greater    than number of gaps ->    pertubation layer
c           Swap first polarity of first pertubation:
            zpol=zpol*(-1)
         else
            zflg=0                                            !If gap-counter for Z is less-equal than number of gaps -> no pertubation layer
         endif
         if(zflg.eq.1) zgap=0                                 !Reset  gap-counter for Z in case of pertubation layer
         zgap=zgap+1                                          !Update gap-counter for next layer

c--------Loop over y-nodes:
c        ygap=pertgapH
         ygap=0
         ylen=0
         ypol=zpol
         do iy = 1,ny

c           Check for checkerboard-Y
c           (horizontal gaps, longitude-width):
            if(ygap.gt.pertgapH.AND.
     &         ylen.le.pertlat) then
               yflg=1                                         !If gap-counter for Y is greater    than number of gaps ->    pertubation Y-column
            else 
               yflg=0                                         !If gap-counter for Y is less-equal than number of gaps -> no pertubation Y-column
            endif

c           In case of    pertubation y-column:
            if(yflg.eq.1.AND.
     &         ylen+1.le.pertlat) then
               ylen=ylen+1                                    !Update pertubation length for Y-column in case of    pertubation Y-column
            endif                                             !Next y-column is less-equal than the pertubation length in latitude -> pertubation continious in next y-column
            if(yflg.eq.1.AND.
     &         ylen+1.gt.pertlat) then
               ygap=1                                         !Reset  gap-counter        for Y-column in case of    pertubation Y-column
               ylen=ylen+1                                    !Update pertubation length for Y-column in case of    pertubation Y-column
            endif                                             !Next y-column is greater    than the pertubation length in latitude -> gap starts in next y-column 

c           In case of no pertubation y-column:
            if(yflg.eq.0) then
               if(ylen.ne.0) ypol=ypol*(-1)                   !Swap polarity for first pertubation of next          pertubation Y-column
               ylen=0                                         !Reset  pertubation length for Y-column in case of no pertubation Y-column
               ygap=ygap+1                                    !Update gap-counter        for Y-column in case of no pertubation Y-column
            endif 
            
c-----------Loop over x-nodes:
c           xgap=pertgapH
            xgap=0
            xlen=0
            xpol=ypol
            do ix = 1,nx

c              Check for checkerboard-X
c              (horizontal gaps, longitude-width)
c              xflg=1

               if(xgap.gt.pertgapH.AND.
     &            xlen.le.pertlon) then
                  xflg=1                                      !If gap-counter for X is greater    than number of gaps ->    pertubation X-node
               else
                  xflg=0                                      !If gap-counter for X is less-equal than number of gaps -> no pertubation X-node
               endif

c              In case of    pertubation x-node:
               if(xflg.eq.1.AND.
     &            xlen+1.le.pertlon) then
                  xlen=xlen+1                                    !Update pertubation length for Y-column in case of    pertubation X-node
               endif                                             !Next X-column is less-equal than the pertubation length in longitude -> pertubation continious in next X-node
               if(xflg.eq.1.AND.
     &            xlen+1.gt.pertlon) then
                  xgap=1                                         !Reset  gap-counter        for Y-column in case of    pertubation X-node
                  xlen=xlen+1                                    !Update pertubation length for Y-column in case of    pertubation X-node
               endif                                             !Next X-node is greater    than the pertubation length in longitude -> gap starts in next y-node

c              In case of no pertubation x-node:
               if(xflg.eq.0) then
                  if(xlen.ne.0) xpol=xpol*(-1)                   !Swap polarity for first pertubation of next        pertubation X-node
                  xlen=0                                         !Reset  pertubation length for X-node in case of no pertubation X-node
                  xgap=xgap+1                                    !Update gap-counter        for X-node in case of no pertubation X-node
               endif

c              Apply pertubation if all x-y-zflgs are 1,
c              checkerboard-flag is "y",
c              and node is not located near edge of model (at least 2 nodes undisturbed) 
               if(chkflg.eq.'y'.AND.
     &            zflg.eq.1.AND.
     &            yflg.eq.1.AND.
     &            xflg.eq.1.AND.
     &            iy.lt.ny-1.AND.
     &            ix.lt.nx-1) then
                  chkbrdamp=pertamp*xpol
               else
                  chkbrdamp=0.0
               endif

c              Copy velocity (+ possible pertubation) to X-vector     
               tmpvel(ix)=velmod(iz,1)+chkbrdamp
            enddo
            
c           Write x-column for current y node:
            write(10,'(300(1x,f4.2))') (tmpvel(ix),ix=1,nx)
         enddo            

      enddo      
 
      close(10)

c     Summary:
      write(*,*)
      write(*,'(a,a)') 'Model file written to: ',
     &      modfi(1:trimlen(modfi))

c     Regular end of main code:
      stop
c
c     Start of error messages:
900   write(*,'(a,a)') 'ERROR 900: Cannot open file : ',
     &      cmnfi(1:trimlen(cmnfi))
      close(10)
      stop
901   write(*,'(a,a)') 'ERROR 901: CMN file incomplete: ',
     &      cmnfi(1:trimlen(cmnfi))
      close(10)
      stop
902   write(*,'(a,a)') 'ERROR 902: Velocity model incomplete: ',
     &      'nz != number of specified layers...'
      stop

c     End of main code:
      end

c     Start of subroutines:
c
      integer function TRIMLEN(t)   ! Urs Kradolfer, June 1986
c     ===========================
c
c     Call:    nc=TRIMLEN(char)
c
c          --> nc says, how many characters the input-string has
c              (ignoring trailing blanks!).
c
      implicit none
      character t*(*)
      do 1 trimlen=LEN(t),1,-1
    1    if(t(trimlen:trimlen).ne.' ')RETURN
      trimlen=1
      end ! of integer function trimlen
c
