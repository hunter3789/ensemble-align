!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE RSHIFT2DGRD                 ######
!######                                                      ######
!######                Copyright (c) 1995                    ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rshift2dgrd(nx,ny,nvar,ipass,mxpass,mxzone,                 &
           max_elem_send,applyshft,posdef,time,xs,ys,                  &
           fcsta,fcstb,fcsta_sm,fcstb_sm,fcstc,                        &
           istart,jstart,ifinish,jfinish,                              &
           ibkshift,jbkshift,xshift,yshift,xsum,ysum,wgtsum,           &
           dxfld,dyfld,rdxfld,rdyfld,                                  &
           slopey,alphay,betay,                                        &
           tem2d,recv_buf,lvldbg)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Identify a set of horizontal translation vectors (xshift and yshift)
!  that when applied to fcsta best matches fcstb.
!  Apply translation to the original gridded data, output as fcstc.
!
!  NOTE: The calling program is responsible for
!        making sure any necessary boundary conditions are satisfied
!        after adjustments are made by this subroutine.
!
!  AUTHOR:
!  Keith Brewster, CAPS, March, 2017
!  After rshift3d
!
!  MODIFICATION HISTORY:
!  Feb 2018
!  Message Passing implementation 
!  Keith Brewster, CAPS
!
!  Modification:
!    ChangJae Lee, KMA, October 2022
!    - Apply Nested Parallelizaion
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mpif.h'
!
  INTEGER, INTENT(IN) :: nx,ny,nvar
  INTEGER, INTENT(IN) :: mxpass,ipass
  INTEGER, INTENT(IN) :: mxzone
  INTEGER, INTENT(IN) :: max_elem_send
  INTEGER, INTENT(IN) :: applyshft
  REAL, INTENT(IN) :: time
!
!-----------------------------------------------------------------------
!
!  Model variables
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(IN)    :: posdef(nvar)
  REAL, INTENT(IN)    :: xs(nx)
  REAL, INTENT(IN)    :: ys(ny)
!
  REAL, INTENT(IN)    :: fcsta(nx,ny,nvar)
  REAL, INTENT(IN)    :: fcstb(nx,ny,nvar)
  REAL, INTENT(OUT)   :: fcsta_sm(nx,ny,nvar)
  REAL, INTENT(OUT)   :: fcstb_sm(nx,ny,nvar)
  REAL, INTENT(OUT)   :: fcstc(nx,ny,nvar)
  INTEGER, INTENT(IN) :: ibkshift(nx,ny), jbkshift(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Shift arrays
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(OUT)    :: xshift(nx,ny)
  REAL, INTENT(OUT)    :: yshift(nx,ny)
  REAL, INTENT(OUT)    :: xsum(nx,ny)
  REAL, INTENT(OUT)    :: ysum(nx,ny)
  REAL, INTENT(OUT)    :: wgtsum(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Interpolation work arrays
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(OUT) :: dxfld(nx)
  REAL, INTENT(OUT) :: dyfld(ny)
  REAL, INTENT(OUT) :: rdxfld(nx)
  REAL, INTENT(OUT) :: rdyfld(ny)
  REAL, INTENT(OUT) :: slopey(nx,ny)
  REAL, INTENT(OUT) :: alphay(nx,ny)
  REAL, INTENT(OUT) :: betay(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(OUT) :: tem2d(nx,ny)
  REAL, INTENT(OUT) :: recv_buf(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Debug Flag
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: lvldbg
!
!-----------------------------------------------------------------------
!
!  Shift zone parameters
!
!-----------------------------------------------------------------------
!
!  REAL, PARAMETER :: hrzlap = 0.50   ! fraction overlapping in horizontal
  INTEGER, PARAMETER :: iplot = 0
!
!-----------------------------------------------------------------------
!
!  Other parameters
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: s = 0.5
!
!-----------------------------------------------------------------------
!
!  Zone storage vectors
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: istart(mxzone)
  INTEGER, INTENT(OUT) :: jstart(mxzone)
  INTEGER, INTENT(OUT) :: ifinish(mxzone)
  INTEGER, INTENT(OUT) :: jfinish(mxzone)
! 
!-----------------------------------------------------------------------
!
! Message Passiing Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: mpstatus(MPI_STATUS_SIZE)
  INTEGER, PARAMETER :: root = 0
  INTEGER, PARAMETER :: send_xsum_tag = 10000
  INTEGER, PARAMETER :: send_ysum_tag = 20000
  INTEGER, PARAMETER :: send_wgt_tag  = 30000
  INTEGER :: buf_size,iproc,ierr
  INTEGER :: color, procspg, ngroup, row_rank, row_size, row_comm
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: xmin,xmax,ymin,ymax,vmin,vmax
  REAL :: dxinv,dyinv,dsdr,riloc,rjloc
  REAL :: deltar,flkdat,rmin,xoptshft,yoptshft,slnratio
  INTEGER :: ilap,jlap,nzone
  INTEGER :: ioffset,joffset
  INTEGER :: istep,jstep
  INTEGER :: i,j,iz,jz,mz,ivar
  INTEGER :: inear,jnear
  INTEGER :: send_xsum,send_ysum,send_wgt
  REAL :: cput1, cput2
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'align.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  buf_size=nx*ny
  dxinv =1./(xs(2)-xs(1))
  dyinv =1./(ys(2)-ys(1))
  xshift(:,:) = 0.
  yshift(:,:) = 0.
!
!-----------------------------------------------------------------------
!
!  Smooth the gridded fields to avoid being misled by small scale
!  features
!
!-----------------------------------------------------------------------

! WRITE(6,'(a,i5)') 'Inside rshift2dgrd, myproc=',myproc!
  IF( myproc == root .and. lvldbg > 10 ) THEN 
      print *, ' nx,ibgn,iend = ',nx,ibgn,iend
      print *, ' ny,jbgn,jend = ',ny,jbgn,jend
      print *, ' nvar = ',nvar
      print *, ' nbaksmth = ',nbaksmth,' s: ',s
  END IF
  fcsta_sm(:,:,:) = fcsta(:,:,:)
  fcstb_sm(:,:,:) = fcstb(:,:,:)
  DO i=1,nbaksmth
    DO ivar=1,nvar
      CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,s,                       &
                     fcsta_sm(1,1,ivar),tem2d,fcsta_sm(1,1,ivar))
      CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,s,                       &
                     fcstb_sm(1,1,ivar),tem2d,fcstb_sm(1,1,ivar))
    END DO
  END DO
  IF( myproc == root .and. lvldbg > 10 ) THEN
    CALL a2dmax0(fcsta_sm(1,1,1),1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
    WRITE(6,'(1x,2(a,f13.4))')                                         &
          'Smoothed field a min = ', vmin,',  max =',vmax
    CALL a2dmax0(fcstb_sm(1,1,1),1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
    WRITE(6,'(1x,2(a,f13.4))')                                         &
          'Smoothed field b min = ', vmin,',  max =',vmax
  END IF
!
!!  DO ipass=1,nshfpass

    nzone=nizone(ipass)*njzone(ipass)
    IF(nzone > mxzone) THEN
      IF( myproc == root ) THEN
        WRITE(6,'(//a,i4)') ' Pass number:',ipass
        WRITE(6,'(a,i4,a,i4,a,i4)')                                    &
            ' nizone= ',nizone(ipass),' njzone= ',njzone(ipass)
        WRITE(6,'(a,i4,a,i4)') ' nzone=',nzone,' mxzone=',mxzone
        WRITE(6,'(a)') ' Increase mxzone sent to rshift3d'
        STOP
        !CALL ARPSSTOP('rshiftgrd')
      END IF
    END IF
!
!-----------------------------------------------------------------------
!
!  Find the appropriate zone size based on the number
!  of zones and overlap required.
!
!-----------------------------------------------------------------------
!
    !izsize=                                                            &
    !    IFIX((FLOAT(nx-3)/((nizone(ipass)-1)*(1.-hrzlap) + 1.))+1.0)
    !jzsize=                                                            &
    !    IFIX((FLOAT(ny-3)/((njzone(ipass)-1)*(1.-hrzlap) + 1.))+1.0)

    !IF (MOD(izsize,2).ne.0) izsize = izsize+1
    !IF (MOD(jzsize,2).ne.0) jzsize = jzsize+1
!
    ilap=MAX(IFIX((izsize(ipass)*hrzlap)+0.5),1)
    jlap=MAX(IFIX((jzsize(ipass)*hrzlap)+0.5),1)
!
!-----------------------------------------------------------------------
!
!  slen is the distance parameter (in grid lengths)
!  used to prevent aliasing -- shifting more than one wavelength
!
!-----------------------------------------------------------------------
!
    !slnratio48h = 0.5
    !slnratio0h = 0.1
    IF( time >= 48. ) THEN
      slnratio=slnratio48h
    ELSE
      slnratio=slnratio0h+((time/48.)*(slnratio48h-slnratio0h))
    END IF
    IF( myproc == root .and. lvldbg > 10 ) WRITE(6,'(/f10.5,f10.5,f10.5)') slnratio, slnratio0h, slnratio48h

    !slnratio = 0.5
    slen=slnratio*sqrt(float(izsize(ipass)*izsize(ipass)+jzsize(ipass)*jzsize(ipass)))
!
    IF( myproc == root ) THEN
      WRITE(6,'(//a,i4)') ' Pass number:',ipass
      WRITE(6,'(a,i4,a,i4)') ' Shift zone dimensions: i=',             &
                         izsize(ipass),'  j=',jzsize(ipass)
      WRITE(6,'(a,i4,a,i4//)') ' Shift zone overlaps  : i=',           &
                         ilap,'  j=',jlap
      WRITE(6,'(a,f10.1)') ' Scaling length ratio: ',slnratio
      WRITE(6,'(a,f10.1)') ' Scaling length: ',slen
    END IF
!
!-----------------------------------------------------------------------
!
    istep=izsize(ipass)-ilap
    ioffset=ibgn+1
    jstep=jzsize(ipass)-jlap
    joffset=jbgn+1

    DO jz=1,njzone(ipass)
      DO iz=1,nizone(ipass)
        mz=(jz-1)*nizone(ipass) + iz
        istart(mz)=((iz-1)*istep) + ioffset
        IF(istart(mz) >= (iend-1)) istart(mz)=iend-1
        ifinish(mz)=istart(mz)+izsize(ipass)
        ifinish(mz)=MIN(ifinish(mz),iend)
        jstart(mz)=((jz-1)*jstep) + joffset
        IF(jstart(mz) >= (jend-1)) jstart(mz)=jend-1
        jfinish(mz)=jstart(mz)+jzsize(ipass)
        jfinish(mz)=MIN(jfinish(mz),jend)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Initialize shift and weight grids
!
!-----------------------------------------------------------------------
!
    xsum = 0.
    ysum = 0.
    wgtsum = 0.
!
!-----------------------------------------------------------------------
!
!  Loop through shift zones, finding appropriate shift for each
!
!-----------------------------------------------------------------------
!
    IF( myproc == root) WRITE(6,'(a)') ' Beginning zone loop'
    IF( myproc == root) CALL CPU_TIME(cput1)


    IF (ipass.eq.1) THEN
      IF (nprocs.ge.8) THEN
        procspg = 8
      ELSE IF (nprocs.ge.4) THEN
        procspg = 4
      ELSE
        procspg = 1
      END IF
    ELSE
      IF (nprocs.ge.2) THEN
        procspg = 2
      ELSE
        procspg = 1
      END IF
    END IF

    ngroup = nprocs/procspg
    color = myproc/procspg
    CALL MPI_Comm_split(MPI_COMM_WORLD, color, myproc, row_comm, ierr)

    CALL MPI_Comm_rank(row_comm, row_rank, ierr)
    CALL MPI_Comm_size(row_comm, row_size, ierr)

    !print *, myproc, color, row_rank, row_size

    DO mz=1,nzone
      IF( mod(mz,ngroup) == color ) THEN
        IF(istart(mz) > 0 .AND. jstart(mz) > 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Find shift for this zone
!
!-----------------------------------------------------------------------
!
!      PRINT *, ' Finding shift for zone: ',mz
!      PRINT *, ' imin,imax,jmin,jmax: ',                             &
!                 istart(mz),ifinish(mz),jstart(mz),jfinish(mz)
!
          CALL rsh2dgrd(nx,ny,nvar,nvarshf,ipass,mz,                   &
                    ibgn,iend,jbgn,jend,fcsta_sm,fcstb_sm,             &
                    istart(mz),ifinish(mz),jstart(mz),jfinish(mz),     &
                    iplot,minkdat,minkdratio,slen,wgtvar,              &
                    ibkshift,jbkshift,izsize(ipass),jzsize(ipass),     &
                    xoptshft,yoptshft,rmin,deltar,lvldbg,              &
                    row_comm,row_rank,row_size)

!       WRITE(6,'(a,f6.1,a,f6.1)')                                     &
!             ' Shift increment  X: ',xoptshft,'  Y: ',yoptshft
!       WRITE(6,'(a,f14.6,a,f14.6)')                                   &
!             '   best r: ',rmin, '  delta r: ',deltar
!
!-----------------------------------------------------------------------
!
!  Fill zone with shift info
!
!-----------------------------------------------------------------------
!
          IF (row_rank == root) THEN
            DO j=jstart(mz),jfinish(mz)
              DO i=istart(mz),ifinish(mz)
                xsum(i,j)=xsum(i,j)+xoptshft
                ysum(i,j)=ysum(i,j)+yoptshft
                wgtsum(i,j)=wgtsum(i,j)+1.0
              END DO
            END DO
          END IF

        END IF ! valid zone
      END IF ! my turn
    END DO ! zone loop

    IF( nprocs > 1 ) THEN
      recv_buf=0.
!     print *, 'myproc: ',myproc,'  calling xsum reduce'
      CALL LGARRAY_REDUCER(xsum,recv_buf,buf_size,max_elem_send, &
                         myproc,root,MPI_SUM,MPI_COMM_WORLD,ierr)
!     print *, 'myproc: ',myproc,'  mpi_reduce status: ',ierr
      IF( myproc == root ) THEN
         xsum(:,:) = recv_buf(:,:)
      END IF

      recv_buf=0.
!     print *, 'myproc: ',myproc,'  calling ysum reduce'
      CALL LGARRAY_REDUCER(ysum,recv_buf,buf_size,max_elem_send, &
                         myproc,root,MPI_SUM,MPI_COMM_WORLD,ierr)
      IF( myproc == root ) THEN
         ysum(:,:) = recv_buf(:,:)
      END IF

      recv_buf=0.
!     print *, 'myproc: ',myproc,'  calling wgtsum reduce'
      CALL LGARRAY_REDUCER(wgtsum,recv_buf,buf_size,max_elem_send, &
                         myproc,root,MPI_SUM,MPI_COMM_WORLD,ierr)
      IF( myproc == root ) THEN
         wgtsum(:,:) = recv_buf(:,:)
      END IF
    END IF

    CALL MPI_Comm_free(row_comm, ierr)

    IF( myproc == root ) WRITE(6,'(a)') ' End of zone loop...'
    IF( myproc == root) CALL CPU_TIME(cput2)
    IF( myproc == root) WRITE(6,'(a,f10.2,a)') ' CPU time:',(cput2-cput1),' seconds'
!
!-----------------------------------------------------------------------
!
!  Normalize shift vectors
!
!-----------------------------------------------------------------------
!
    IF( myproc == root ) THEN
      DO j=jbgn,jend
        DO i=ibgn,iend
          IF(wgtsum(i,j) > 0.) THEN
            xshift(i,j)=xshift(i,j)+(xsum(i,j)/wgtsum(i,j))
            yshift(i,j)=yshift(i,j)+(ysum(i,j)/wgtsum(i,j))
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Set shift to zero-gradient at boundaries
!
!-----------------------------------------------------------------------
!
      DO i=ibgn,iend
        xshift(i,jbgn) = xshift(i,(jbgn+1))
        xshift(i,jend) = xshift(i,(jend-1))
        yshift(i,jbgn) = yshift(i,(jbgn+1))
        yshift(i,jend) = yshift(i,(jend-1))
      END DO

      DO j=jbgn,jend
        xshift(ibgn,j) = xshift((ibgn+1),j)
        xshift(iend,j) = xshift((iend-1),j)
        yshift(ibgn,j) = yshift((ibgn+1),j)
        yshift(iend,j) = yshift((iend-1),j)
      END DO
!
!-----------------------------------------------------------------------
!
!  Smooth the shift vectors
!
!-----------------------------------------------------------------------
!
      DO i=1,nshfsmth
        CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,s,                     &
                          xshift,tem2d,xshift)
        CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,s,                     &
                          yshift,tem2d,yshift)
      END DO


!
!-----------------------------------------------------------------------
!
!  Report min and max of shift vactors
!
!-----------------------------------------------------------------------
!
      IF( myproc == root .and. lvldbg > 10 ) THEN
        CALL a2dmax0(xshift,1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
        WRITE(6,'(1x,2(a,f13.4))')                                     &
              'xshift min = ', vmin,',  xshift max =',vmax
        CALL a2dmax0(yshift,1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
        WRITE(6,'(1x,2(a,f13.4))')                                     &
              'yshift min = ', vmin,',  yshift max =',vmax
      END IF
!
!-----------------------------------------------------------------------
!
!  Apply the shifts, one level at a time
!  Report min and max of shift vactors and variable values.
!
!-----------------------------------------------------------------------
!
      IF( applyshft > 0 ) THEN
        DO ivar=1,nvar
          IF( myproc == root .and. lvldbg > 10 ) THEN
            WRITE(6,'(a,i6)') 'shifting variable:',ivar
            CALL a2dmax0(fcsta(1,1,ivar),1,nx,ibgn,iend,               &
                         1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))')                                 &
              'Pre-shift  min = ', vmin,',  max =',vmax
          END IF
          CALL movegr(nx,ny,fcsta(1,1,ivar),tem2d,                     &
                    xshift,yshift,fcstc(1,1,ivar),                     &
                    ibgn,iend,jbgn,jend,                               &
                    dxfld,dyfld,rdxfld,rdyfld,                         &
                    slopey,alphay,betay)
          IF( myproc == root .and. lvldbg > 10 ) THEN
            CALL a2dmax0(fcstc(1,1,ivar),1,nx,ibgn,iend,               &
                         1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))')                                 &
              'Post-shift min = ', vmin,',  max =',vmax
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!  Make sure the final water fields are positive.
!
!-----------------------------------------------------------------------
!
      DO ivar=1,nvar
        IF( posdef(ivar) == 1) THEN
          DO j=jbgn,jend
            DO i=ibgn,iend
              fcstc(i,j,ivar)=MAX(0.,fcstc(i,j,ivar))
            END DO
          END DO
        END IF
      END DO
    END IF 
!!  END DO
  RETURN
END SUBROUTINE rshift2dgrd
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE RSH2dGRD                    ######
!######                                                      ######
!######                Copyright (c) 1995                    ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rsh2dgrd(nx,ny,nvar,nvarshf,ipass,mz,                       &
                    ibgn,iend,jbgn,jend,fcsta,fcstb,                   &
                    iminz,imaxz,jminz,jmaxz,                           &
                    iplot,minkdat,minkdratio,slen,wgtvar,              &
                    ibkshift,jbkshift,izsize,jzsize,                   &
                    xoptshft,yoptshft,rmin,deltar,lvldbg,              &
                    row_comm,row_rank,row_size)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Match two 2d grids by shifting first grid
!
!  AUTHOR:
!  Keith Brewster, CAPS, March, 2017
!  Based on rshft2dint
!
!  Modification:
!    ChangJae Lee, KMA, October 2022
!    - Apply Nested Parallelizaion
!    - Penalize rshift Using Threshold Value
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mpif.h'
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: nx,ny,nvar,nvarshf,ipass,mz
  INTEGER, INTENT(IN) :: ibgn,iend,jbgn,jend
!
  REAL, INTENT(IN)    :: fcsta(nx,ny,nvar)
  REAL, INTENT(IN)    :: fcstb(nx,ny,nvar)
!
  INTEGER, INTENT(IN) :: iminz,imaxz,izsize
  INTEGER, INTENT(IN) :: jminz,jmaxz,jzsize
  INTEGER, INTENT(IN) :: iplot
  INTEGER, INTENT(IN) :: minkdat
  INTEGER, INTENT(IN) :: ibkshift(nx,ny), jbkshift(nx,ny)
  REAL, INTENT(IN)    :: minkdratio
  REAL, INTENT(IN)    :: slen
  REAL, INTENT(IN)    :: wgtvar(nvarshf)
!
  REAL, INTENT(OUT)   :: xoptshft,yoptshft
  REAL, INTENT(OUT)   :: rmin,deltar
  INTEGER, INTENT(IN) :: lvldbg
  INTEGER, INTENT(IN) :: row_comm,row_rank,row_size
!
!-----------------------------------------------------------------------
!
!  Sector search arrays
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: root = 0
  INTEGER :: nishift, njshift
  
  CHARACTER (LEN=3), ALLOCATABLE :: chpr(:)
  CHARACTER (LEN=5) :: lengthr
  REAL, ALLOCATABLE :: r(:,:)
  REAL, ALLOCATABLE :: kdratio(:,:)
  REAL :: npts, kpts, sumsqx, recv_buf
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,ii,jj,is,js
  INTEGER :: izero,jzero,ioff,joff
  INTEGER :: iis,jjs,iioff,jjoff
  INTEGER :: imin,jmin,ierr
  REAL :: rthresh,fnpts0inv,nvalid,npts0,kpts0
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Intialize r array to 9.E30.
!
!-----------------------------------------------------------------------
!
  !nishift = izsize
  !njshift = jzsize
  nishift = int(slen)
  njshift = int(slen)

  IF (MOD(nishift,2).eq.0) then
    nishift = nishift + 1
  END IF
  IF (MOD(njshift,2).eq.0) then
    njshift = njshift + 1
  END IF
  IF (nishift.lt.11) then
    nishift = 11
  END IF
  IF (njshift.lt.11) then
    njshift = 11
  END IF
  IF (nishift.gt.51) then
    nishift = 51
  END IF
  IF (njshift.gt.51) then
    njshift = 51
  END IF
  !print *, myproc, nishift, njshift
  allocate(r(nishift,njshift), kdratio(nishift,njshift))
  r(:,:) = 9.0E30
  kdratio(:,:)=0.0
  allocate(chpr(nishift))

!
!-----------------------------------------------------------------------
!
!  Compute RMS difference in the center of the
!  shift grid (i.e., for zero displacement)
!  Search by ones in the vicinity of zero offset.
!  Find local minimum.
!
!
!-----------------------------------------------------------------------
!
  izero=(nishift/2)+1
  jzero=(njshift/2)+1
  rmin=9.0E30

  ioff=0
  joff=0
  is=ioff+izero
  js=joff+jzero
  CALL grdoffstrms(nx,ny,nvar,nvarshf,ipass,mz,                        &
                   ibgn,iend,jbgn,jend,                                &
                   ibkshift,jbkshift,fcsta,fcstb,                      &
                   iminz,imaxz,jminz,jmaxz,                            &
                   minkdat,slen,wgtvar,ioff,joff,npts0,r(is,js),       &
                   row_comm,row_rank,row_size)

  IF (r(is,js) .GE. 0.) THEN
    rmin=r(is,js)
  END IF
  IF (row_rank == root) THEN
    fnpts0inv=1.0/(npts0)
    kdratio(is,js)=1.0
  END IF

  DO ioff=-1,1
    DO joff=-1,1
      IF( ioff == 0 .AND. joff == 0) CYCLE
      is=ioff+izero
      js=joff+jzero
      CALL grdoffstrms(nx,ny,nvar,nvarshf,ipass,mz,                    &
                       ibgn,iend,jbgn,jend,                            &
                       ibkshift,jbkshift,fcsta,fcstb,                  &
                       iminz,imaxz,jminz,jmaxz,                        &
                       minkdat,slen,wgtvar,ioff,joff,nvalid,r(is,js),  &
                       row_comm,row_rank,row_size)
      IF (r(is,js) .GE. 0.) THEN
        rmin=AMIN1(rmin,r(is,js))
      END IF
      IF (row_rank == root) THEN
        kdratio(is,js)=(nvalid)*fnpts0inv
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Search entire shift-possibility matrix by threes.
!  Search by ones around each of those shift points
!  1) If the RMS found is less than 2.0 times min RMS
!     from near-zero shift
!  2) If the distance is within a block of 5 points from zero.
!
!  Note rmin is from the near-zero search and is not replaced
!  by any lesser rms's found elsewhere.
!
!-----------------------------------------------------------------------
!
  rthresh=2.0*rmin
  DO js=1,njshift,3
    DO is=1,nishift,3
      ioff=is-izero
      joff=js-jzero
      CALL grdoffstrms(nx,ny,nvar,nvarshf,ipass,mz,                    &
                       ibgn,iend,jbgn,jend,                            &
                       ibkshift,jbkshift,fcsta,fcstb,                  &
                       iminz,imaxz,jminz,jmaxz,                        &
                       minkdat,slen,wgtvar,ioff,joff,nvalid,r(is,js),  &
                       row_comm,row_rank,row_size)
      IF (row_rank == root) THEN
        kdratio(is,js)=(nvalid)*fnpts0inv
      END IF
      IF( (ABS(ioff) < 5 .AND. ABS(joff) < 5) .OR. &
          (r(is,js) < rthresh .AND. r(is,js) .GE. 0.) ) THEN
        DO ii=-1,1,1
          DO jj=-1,1,1
            IF( ii /= 0 .OR. jj /= 0 ) THEN
              iis=is+ii
              jjs=js+jj
              IF( iis > 0 .AND. iis <= nishift .AND.                   &
                    jjs > 0 .AND. jjs <= njshift ) THEN
                iioff=ioff+ii
                jjoff=joff+jj
                CALL grdoffstrms(nx,ny,nvar,nvarshf,ipass,mz,          &
                       ibgn,iend,jbgn,jend,                            &
                       ibkshift,jbkshift,fcsta,fcstb,                  &
                       iminz,imaxz,jminz,jmaxz,                           &
                       minkdat,slen,wgtvar,iioff,jjoff,nvalid,r(iis,jjs), &
                       row_comm,row_rank,row_size)
                IF (row_rank == root) THEN
                  kdratio(iis,jjs)=(nvalid)*fnpts0inv
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Search for minimum
!
!-----------------------------------------------------------------------
!
  IF (row_rank == root) THEN
    rmin=r(izero,jzero)
    imin=izero
    jmin=jzero
    IF (rmin .GE. 0) THEN
      DO js=1,njshift
        DO is=1,nishift
          IF(kdratio(is,js) >= minkdratio .AND. r(is,js) .GE. 0. .AND. &
             r(is,js) < rmin) THEN
            imin=is
            jmin=js
            rmin=r(is,js)
          END IF
        END DO
      END DO
    END IF
    xoptshft=FLOAT(imin-izero)
    yoptshft=FLOAT(jmin-jzero)
    deltar=r(izero,jzero)-rmin
  END IF

  RETURN
END SUBROUTINE rsh2dgrd
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GRDOFFSTRMS                 ######
!######                                                      ######
!######                   Developed by                       ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grdoffstrms(nx,ny,nvar,nvarshf,ipass,mz,ibgn,iend,jbgn,jend,&
                       ibkshift,jbkshift,fcsta,fcstb,                  &
                       iminz,imaxz,jminz,jmaxz,                        &
                       minkdat,slen,wgtvar,ioff,joff,npts,rshift,      &
                       row_comm,row_rank,row_size)
!
!-----------------------------------------------------------------------
!
!
!  PURPOSE:
!
!  Compute normalized root-mean-square error statistic for
!  a single offset grid shift (ioff,joff).
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, Sep, 1997
!
!  MODIFICATION HISTORY:
!
!  Keith Brewster, CAPS, May, 2021
!  Updated to argument list to allow for variable valid grid limits
!  (ibgn,iend,jbgn,jend)
!
!  Modification:
!    ChangJae Lee, KMA, October 2022
!    - Apply Nested Parallelizaion
!    - Penalize rshift Using Threshold Value
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mpif.h'
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN)  :: nx,ny,nvar,nvarshf,ipass,mz
  INTEGER, INTENT(IN)  :: ibgn,iend,jbgn,jend
  INTEGER, INTENT(IN)  :: ibkshift(nx,ny), jbkshift(nx,ny)
!
  REAL, INTENT(IN)     :: fcsta(nx,ny,nvar)
  REAL, INTENT(IN)     :: fcstb(nx,ny,nvar)
!
  INTEGER, INTENT(IN)  :: iminz,imaxz
  INTEGER, INTENT(IN)  :: jminz,jmaxz

  INTEGER, INTENT(IN)  :: minkdat
  REAL, INTENT(IN)     :: slen
  REAL, INTENT(IN)     :: wgtvar(nvarshf)
  INTEGER, INTENT(IN)  :: ioff
  INTEGER, INTENT(IN)  :: joff
  REAL, INTENT(OUT)    :: npts,rshift
  INTEGER, INTENT(IN)  :: row_comm,row_rank,row_size
  REAL :: kpts,sumsqx,recv_buf

  INTEGER :: ierr
  INTEGER, PARAMETER :: root = 0
!
!-----------------------------------------------------------------------
!
!  External functions
!
!-----------------------------------------------------------------------
!
  REAL :: soar
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: thresh,kratio,nratio
  REAL :: vardiff,mserr,dist
  INTEGER :: kdat,halfk
  INTEGER :: i,j,ii,jj,ivar,iii,jjj
  REAL  :: cput1, cput2
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  !CALL CPU_TIME(cput1)
  thresh = 0.25
  !thresh = 0.01
  kdat=((imaxz-iminz)+1)*((jmaxz-jminz)+1)
  sumsqx=0.
  !dist=0.
  npts=0
  kpts=0
  halfk=(kdat/2)+1
  DO j=jminz,jmaxz
    IF( mod(j,row_size) == row_rank) THEN
      DO i=iminz,imaxz
!      print *, ' OFFSTRMS, ioff,joff,kdat: ',kdat,ioff,joff
        ii=i+ioff
        jj=j+joff
!
!-----------------------------------------------------------------------
!
!    Do not use points that go outside the boundaries.
!
!-----------------------------------------------------------------------
!
        IF( ii <= iend .AND.  ii >= ibgn  .AND. jj <= jend .AND.  jj >= jbgn) THEN
          iii = ii+ibkshift(i,j)
          jjj = jj+jbkshift(i,j)

          IF( iii <= iend .AND.  iii >= ibgn  .AND. jjj <= jend .AND.  jjj >= jbgn ) THEN
            npts=npts+1
            DO ivar=1,nvar
              IF (fcsta(iii,jjj,ivar) < thresh .AND. fcstb(i,j,ivar) < thresh) THEN
                kpts=kpts+1
              ENDIF
              vardiff = (fcsta(iii,jjj,ivar) - fcstb(i,j,ivar))
              sumsqx = sumsqx + wgtvar(ivar)*(vardiff*vardiff)
            END DO     
          ELSE
            npts=npts+1
            DO ivar=1,nvar
              IF (fcstb(i,j,ivar) < thresh) THEN
                kpts=kpts+1
              ENDIF
              vardiff = fcstb(i,j,ivar)
              sumsqx = sumsqx + wgtvar(ivar)*(vardiff*vardiff)
            END DO
          END IF

        ELSE
          npts=npts+1
          DO ivar=1,nvar
            IF (fcstb(i,j,ivar) < thresh) THEN
              kpts=kpts+1
            ENDIF
            vardiff = fcstb(i,j,ivar)
            sumsqx = sumsqx + wgtvar(ivar)*(vardiff*vardiff)
          END DO
        END IF

      END DO
    END IF
  END DO

  IF( row_size > 1 ) THEN
    recv_buf=0.
    CALL MPI_REDUCE(sumsqx,recv_buf,1, &
                    MPI_REAL,MPI_SUM,root,row_comm,ierr)
    IF (row_rank == root) THEN
      sumsqx = recv_buf
    END IF

    recv_buf=0.
    CALL MPI_REDUCE(npts,recv_buf,1, &
                    MPI_REAL,MPI_SUM,root,row_comm,ierr)
    IF (row_rank == root) THEN
      npts = recv_buf
    END IF

    recv_buf=0.
    CALL MPI_REDUCE(kpts,recv_buf,1, &
                    MPI_REAL,MPI_SUM,root,row_comm,ierr)
    IF (row_rank == root) THEN
      kpts = recv_buf
    END IF
  END IF

  IF (row_rank == root) THEN
    kratio = (kpts)/(npts*nvar+1)

    IF( npts >= minkdat  .AND. npts >= halfk .AND. npts*nvar .NE. kpts) THEN
      mserr=sumsqx/(npts)
      dist=SQRT(FLOAT(ioff*ioff + joff*joff))
      rshift=mserr/soar(dist,slen*(1-kratio))
    ELSE
      !rshift=99999.
      rshift=-1.
      !rshift=9.0E30
    END IF
  END IF

  IF( row_size > 1 ) THEN
    CALL MPI_BCAST(rshift,1,MPI_REAL,root,row_comm,ierr)
  END IF

  RETURN
END SUBROUTINE grdoffstrms
!

FUNCTION soar(dist,s)
  IMPLICIT NONE
!
  REAL :: soar,dist,s
  REAL :: ratio
!
  ratio=dist/s
  soar=(1. + ratio)*EXP(-ratio)
  RETURN
END FUNCTION soar
