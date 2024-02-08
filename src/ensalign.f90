 PROGRAM ENSALIGN
!
! Program to create ensemble mean after spatially aligning results
!
! Author: Keith Brewster
! CAPS/University of Oklahoma
! March, 2017
!
! Modification:
! MPI Version
!    Keith Brewster, CAPS, October 2017
!
!
! Modification:
!    ChangJae Lee, KMA, October 2022
!    - Use eccodes to read Grib file
!    - Write shifted means (also include PM and LPM) and shift vectors
!      into NETCDF file
!

  USE eccodes
  USE netcdf

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'netcdf.inc'
  INCLUDE 'mpif.h'
  INCLUDE 'align.inc'
!
!-----------------------------------------------------------------------
!
! Forecast Fields
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  CHARACTER(LEN=512), ALLOCATABLE :: ensfname(:)
  REAL, ALLOCATABLE :: ensfcst(:,:,:)
  REAL, ALLOCATABLE :: ensfcst_smt(:,:,:,:)
  REAL, ALLOCATABLE :: ensfcst_shf(:,:,:,:)
  REAL, ALLOCATABLE :: ensmean(:,:)
  REAL, ALLOCATABLE :: ensmax(:,:)
  REAL, ALLOCATABLE :: ensshfmn(:,:,:)
  REAL, ALLOCATABLE :: ensshfmax(:,:,:)
  REAL, ALLOCATABLE :: enslpm(:,:), enspm(:,:), enspm_smt(:,:,:)
  REAL, ALLOCATABLE :: ensshflpm(:,:,:), ensshfpm(:,:,:)
  REAL, ALLOCATABLE :: lats(:,:), lons(:,:)
!
!-----------------------------------------------------------------------
!
! Spatial Alignment Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: posdef(:)
  INTEGER, ALLOCATABLE :: ibkshift(:,:,:,:), jbkshift(:,:,:,:)
  REAL, ALLOCATABLE :: xshift(:,:)
  REAL, ALLOCATABLE :: yshift(:,:)
  REAL, ALLOCATABLE :: xshiftmn(:,:,:,:)
  REAL, ALLOCATABLE :: yshiftmn(:,:,:,:)
  REAL, ALLOCATABLE :: xsum(:,:)
  REAL, ALLOCATABLE :: ysum(:,:)
  REAL, ALLOCATABLE :: wgtsum(:,:)
  REAL, ALLOCATABLE :: tem2d1(:,:)
  REAL, ALLOCATABLE :: tem2d2(:,:)
  REAL, ALLOCATABLE :: tem2d3(:,:) 
  REAL, ALLOCATABLE :: recv_buf(:,:)
  INTEGER, ALLOCATABLE :: zero2d(:,:) 
  CHARACTER(LEN=1), ALLOCATABLE :: mpi_work_buf(:)
!
!-----------------------------------------------------------------------
!
!  Interpolation work arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: istart(:)
  REAL, ALLOCATABLE :: jstart(:)
  REAL, ALLOCATABLE :: ifinish(:)
  REAL, ALLOCATABLE :: jfinish(:)

  REAL, ALLOCATABLE :: dxfld(:)
  REAL, ALLOCATABLE :: dyfld(:)
  REAL, ALLOCATABLE :: rdxfld(:)
  REAL, ALLOCATABLE :: rdyfld(:)
  REAL, ALLOCATABLE :: slopey(:,:)
  REAL, ALLOCATABLE :: alphay(:,:)
  REAL, ALLOCATABLE :: betay(:,:)
  REAL, ALLOCATABLE :: gs_weight(:)
!
!-----------------------------------------------------------------------
!
! Dimensions 
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nxny,nelements,nx_grib,ny_grib,cx,cy
  INTEGER :: unum,ilcc,acctime,min_nmem
!
!-----------------------------------------------------------------------
!
!  Misc Internal Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: root = 0
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: max_elem_send = 10000
  REAL, PARAMETER    :: smt_param = 0.5

  INTEGER :: i,j,k,k1,k2
  INTEGER :: ictr,jctr,realdata,missing
  INTEGER :: lvldbg,ipass,mxzone,applyshft
  INTEGER :: itime,ifhr,ntstep,itstep,rtime
  INTEGER :: narg,membknt,nelem2d,nelem3d
  INTEGER :: lstring,rsize,bufelems,mpibufsize
  INTEGER :: istat,istatus,inqstat,ierr,ireturn
  INTEGER :: istat_nx,istat_ny
  INTEGER :: year,month,day,hour,minute
  INTEGER :: ilap,jlap,istep,jstep
  INTEGER :: odate(4),pdate(4)
  REAL :: dx,dy
  REAL :: ctrx0,ctry0,ctrx1,ctry1,ctrx,ctry
  REAL :: delx,dely
  REAL :: refmax,radx,rady
  REAL :: piov2,radx2inv,rady2inv,radius
  REAL :: rninv,fnorm,fcstmin,fcstmax,favgmin,favgmax
  REAL :: time,vmin,vmax
  REAL :: cput0,cput1,cput2,cput3,cput4,cput5,cput6,cput7
  REAL, PARAMETER :: rmisg_data = 9.0E36

  INTEGER             :: ifile,idx,igrib,fsout_opt
  INTEGER             :: iret,ncid,ndim,count,nvarout
  INTEGER,ALLOCATABLE :: dimids(:),vardims(:),chunks(:),varids(:)
  INTEGER,ALLOCATABLE :: vardims2(:),chunks2(:),vardims3(:)
  REAL,ALLOCATABLE    :: values(:),projx(:),projy(:)

  CHARACTER(LEN=4)    :: fhrstr
  CHARACTER(LEN=512)  :: fnamelist
  CHARACTER(LEN=256)  :: finname
  CHARACTER(LEN=512)  :: foutname, shfsoutdir
  CHARACTER(LEN=512)  :: command_string
  CHARACTER(LEN=256)  :: fout_head, fsout_tail

  LOGICAL :: file_exists,dir_exists,lposdef 
!
!-----------------------------------------------------------------------
!
!  Shift namelists
!
!-----------------------------------------------------------------------
!
  NAMELIST /ens_data/ realdata,ilcc,basedir,varname,nmembers,membname, &
                      acctime,min_nmem
  NAMELIST /shift_const/ alignopt,calcmean,iposdef,nshfpass,           &
                 nbaksmth,nshfsmth,noutsmth,minkdat,minkdratio,        &
                 slnratio0h,slnratio48h,applyshft
  NAMELIST /shift_zone/ hrzlap,iborder,jborder,izsize,jzsize,          &
                        loopstep,procspg
  NAMELIST /shift_wgt/ wgtvar,thresh_flag,threshvar
  NAMELIST /lpm_const/ patch_nx,patch_ny,ovx,ovy,               &
                       gauss_sigma,filt_min
  NAMELIST /shift_out/ shfoutdir,fout_head,fsout_opt,shfsoutdir,fsout_tail
  NAMELIST /debug/ lvldbg
!
!-----------------------------------------------------------------------
!
! Intializations
! For now these are assigned.
! Some should come from NAMELIST
!
!-----------------------------------------------------------------------
! 
  CALL MPI_INIT( ierr )
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, myproc, ierr)

  IF(myproc == root) THEN
    WRITE(6,'(a,i5)') ' Number of processors: ',nprocs
    CALL CPU_TIME(cput0)
  END IF
 
  membknt = 0

  narg=COMMAND_ARGUMENT_COUNT()
  IF( narg .eq. 3 ) THEN
    call getarg(1,datehrstr)
    call getarg(2,fhrstr)
    call getarg(3,fnamelist)
    read(fhrstr,*) ifhr
    READ(datehrstr,'(i4,i2,i2,i2)') odate(1), odate(2), odate(3), odate(4)
    CALL GETH_NEWDATE(pdate, odate, ifhr*3600)
    call calc_hour(pdate, itime)
    call calc_hour(odate, rtime)
  ELSE
    IF(myproc == root) THEN
      WRITE(stdout,'(1x,a)') 'needs three argument: ${YYYYMMDDHH} ${fhr} ${fnamelist}'
    END IF
    STOP
  END IF
!
!-------------------------------------------------------------------------
!
! Get NAMELIST variables
!
!-----------------------------------------------------------------------
!
  CALL GET_UNIT(unum)
  INQUIRE ( file = TRIM(fnamelist) , exist = file_exists )

  if ( file_exists ) then
    OPEN(unum,FILE=TRIM(fnamelist),STATUS='OLD',FORM='FORMATTED')

    READ(unum, ens_data)
    READ(unum, shift_const)
    READ(unum, shift_zone)
    READ(unum, shift_wgt)
    READ(unum, lpm_const)
    READ(unum, shift_out)
    READ(unum, debug)

    CLOSE(unum)
  else
    print '(A,A)',' Namelist could not be founded : ', TRIM(fnamelist)
    stop
  endif

!
!-------------------------------------------------------------------------
!
! Get Grid information (nx, ny, lats, lons)
! realdata - 0 : analytic case
! realdata - 1 : CAPS Ensemble / 2 : HREF / 3 : GEFS
!
!-----------------------------------------------------------------------
!
  IF ( realdata == 0 ) THEN
    nx=75
    ny=75
    nmembers = 5
  ELSE  IF ( realdata == 1 .OR. realdata == 2 ) THEN
    nx_grib = 0
    ny_grib = 0

    IF (myproc == root) THEN
      DO k=1,nmembers
        IF ( realdata == 1 ) THEN
          WRITE(finname,'(a,i4.4,3(i2.2),a,a,a,i4.4,3(i2.2),a,i3.3,a)') trim(basedir), &
              odate(1), odate(2), odate(3), odate(4), '/', trim(membname(k)), '_', odate(1), odate(2), odate(3), odate(4), &
              'f', ifhr, '.grib2'
        ELSE IF ( realdata == 2 ) THEN
          WRITE(finname,'(a,i4.4,2(i2.2),a,a,a,i4.4,3(i2.2),a,i3.3,a)') trim(basedir), &
              odate(1), odate(2), odate(3), '/', trim(membname(k)), '_', odate(1), odate(2), odate(3), odate(4), &
              'f', ifhr, '.grib2'
        END IF
        WRITE(6,'(2a)') ' read file: ', trim(finname)

        call codes_open_file(ifile,finname,'r',istat)
        if (istat .eq. 0) then
          call codes_index_create(idx,finname,'shortName,level,lengthOfTimeRange',istat)

          call codes_index_select(idx,'shortName',varname,istat)
          call codes_index_select(idx,'level',0,istat)
          IF ( realdata == 1 ) THEN
            call codes_index_select(idx,'lengthOfTimeRange',acctime,istat)
          ELSE IF ( realdata == 2 ) THEN
            call codes_index_select(idx,'lengthOfTimeRange',1,istat)
          END IF

          call codes_new_from_index(idx,igrib,istat)

          call codes_get(igrib,'Nx',nx,istat_nx)
          call codes_get(igrib,'Ny',ny,istat_ny)

          if (istat_nx .eq. 0 .and. istat_ny .eq. 0) then
            allocate(values(nx*ny), lats(nx,ny), lons(nx,ny))
            call codes_get(igrib,'latitudes',values,istat)

            if (istat .eq. 0) then
              do i=1,nx
                do j=1,ny
                  lats(i,j) = values((j-1)*nx+i)
                enddo
              enddo

              nx_grib = 1
            end if

            call codes_get(igrib,'longitudes',values,istat)

            if (istat .eq. 0) then
              do i=1,nx
                do j=1,ny
                  lons(i,j) = values((j-1)*nx+i)
                enddo
              enddo

              ny_grib = 1
            end if
          end if

          if (allocated(values)) then
            deallocate(values)
          end if
          call codes_release(igrib,istat)
        END IF

        if (nx_grib .eq. 1 .and. ny_grib .eq. 1) then
          exit
        end if
      END DO
    END IF

    CALL MPI_BCAST(nx,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ny,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  END IF

  dx = 3000.
  dy = 3000.

  ibgn = 1
  jbgn = 1
  iend = nx
  jend = ny

  lposdef=(iposdef > 0)
!
!#######################################################################
!
! Shift Zones
!
! nizone: number of zones in x direction for each shift iteration
! njzone: number of zones in y direction for each shift iteration
!
!#######################################################################

  mxzone=0
  DO ipass=1,nshfpass
    ilap=MAX(IFIX((izsize(ipass)*hrzlap)+0.5),1)
    jlap=MAX(IFIX((jzsize(ipass)*hrzlap)+0.5),1)
    istep=izsize(ipass)-ilap
    jstep=jzsize(ipass)-jlap
    !nizone(ipass)=(nx-(2*iborder(ipass)))/istep
    !njzone(ipass)=(ny-(2*jborder(ipass)))/jstep
    nizone(ipass)=MAX((nx-(2*iborder(ipass)))/istep,1)
    njzone(ipass)=MAX((ny-(2*jborder(ipass)))/jstep,1)
    mxzone=MAX(mxzone,(nizone(ipass)*njzone(ipass)))
    IF (myproc == root) THEN
      IF (lvldbg > 0) THEN
        WRITE(6,'(4(a,i0))') 'nizone(', ipass, ') : ', nizone(ipass), &
            ', njzone(', ipass, '): ', njzone(ipass)
      END IF
    END IF
  END DO

!
!#######################################################################
!
!  weighting for each variable (in addition to normalization to account
!  for standard error of observation).
!  index   variable
!    1     Precipitation
!
!#######################################################################

!
!-----------------------------------------------------------------------
!
!     Check on status of output directories
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Allocate Arrays
!
!-----------------------------------------------------------------------
!
  ALLOCATE(posdef(nvarshf))
  ALLOCATE(xs(nx))
  ALLOCATE(ys(ny))
  ALLOCATE(ensfname(nmembers))
  ALLOCATE(ensmean(nx,ny))
  ALLOCATE(ensmax(nx,ny))
  ALLOCATE(xshift(nx,ny))
  ALLOCATE(yshift(nx,ny))
  ALLOCATE(ibkshift(nx,ny,nmembers,nshfpass))
  ALLOCATE(jbkshift(nx,ny,nmembers,nshfpass))
  ALLOCATE(xshiftmn(nx,ny,nmembers,nshfpass))
  ALLOCATE(yshiftmn(nx,ny,nmembers,nshfpass))
  ALLOCATE(xsum(nx,ny))
  ALLOCATE(ysum(nx,ny))
  ALLOCATE(wgtsum(nx,ny))
  ALLOCATE(ensshfmn(nx,ny,nshfpass),ensshfmax(nx,ny,nshfpass))
  ALLOCATE(enslpm(nx,ny),enspm(nx,ny),enspm_smt(nx,ny,nshfpass))
  ALLOCATE(ensshflpm(nx,ny,nshfpass),ensshfpm(nx,ny,nshfpass))
  ALLOCATE(gs_weight(nx*ny))

  ALLOCATE(ensfcst(nx,ny,nmembers))
  ALLOCATE(ensfcst_smt(nx,ny,nmembers,nshfpass))
  ALLOCATE(ensfcst_shf(nx,ny,nmembers,nshfpass))

  ALLOCATE(istart(mxzone))
  ALLOCATE(jstart(mxzone))
  ALLOCATE(ifinish(mxzone))
  ALLOCATE(jfinish(mxzone))

  ALLOCATE(tem2d1(nx,ny))
  ALLOCATE(tem2d2(nx,ny))
  ALLOCATE(tem2d3(nx,ny))
  ALLOCATE(recv_buf(nx,ny))
  ALLOCATE(dxfld(nx))
  ALLOCATE(dyfld(ny))
  ALLOCATE(rdxfld(nx))
  ALLOCATE(rdyfld(ny))
  ALLOCATE(slopey(nx,ny))
  ALLOCATE(alphay(nx,ny))
  ALLOCATE(betay(nx,ny))
  ALLOCATE(zero2d(nx,ny))

  if (.not. allocated(lats)) then
    ALLOCATE(lats(nx,ny))
  end if
  if (.not. allocated(lons)) then
    ALLOCATE(lons(nx,ny))
  end if
!
!-----------------------------------------------------------------------
!
!  Initialize Arrays
!
!-----------------------------------------------------------------------
!
  IF(myproc == root) THEN
    IF (lvldbg > 0) THEN
      WRITE(6,'(a)') 'Initialize Arrays' 
    END IF
  END IF
  CALL CPU_TIME(cput0)
  minute=0
  ensfname(:)='NULL'
  ensmean(:,:) = 0.
  ensmax(:,:) = 0.
  ensshfmn(:,:,:) = 0.
  ensshfmax(:,:,:) = 0.
  enslpm(:,:) = 0.
  ensshflpm(:,:,:) = 0.
  enspm(:,:) = 0.
  ensshfpm(:,:,:) = 0.
  posdef(:)=iposdef
  dxfld(:) = 0.
  dyfld(:) = 0.
  rdxfld(:) = 0.
  rdyfld(:) = 0.

! bufelems=nx*ny
! CALL MPI_PACK_SIZE(bufelems,MPI_REAL,MPI_COMM_WORLD,rsize,ierr)
! print *, ' MPI PACK REAL size: ',rsize
! mpibufsize=nprocs*(rsize+MPI_BSEND_OVERHEAD)
! print *, ' mpibufsize: ',mpibufsize
! ALLOCATE(mpi_work_buf(mpibufsize),STAT=ierr)
! print *, ' mpi_work_buf status: ',ierr
  IF(myproc == root) THEN
    CALL gaussian_weight(nx,ny,gauss_sigma,gs_weight)
  END IF

  IF (gauss_sigma > 0) THEN
    IF ( nprocs > 1 ) THEN
      nelem2d=nx*ny
      IF( nelem2d < max_elem_send) THEN
        CALL MPI_BCAST(gs_weight,nelem2d,MPI_REAL,root, &
                   MPI_COMM_WORLD,ierr)
      ELSE
        CALL LGARRAY_BCASTR(gs_weight,nelem2d,max_elem_send, &
                        myproc,root,MPI_COMM_WORLD,ierr)
      END IF
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Get Ensemble Forecast Data
!
!-----------------------------------------------------------------------
!
  IF(myproc == root) WRITE(6,'(a)') ' Get Ensemble Forecast Data' 
  DO i=1,nx
    xs(i)=(float(i-2)+0.5)*dx
  END DO
  DO j=1,ny
    ys(j)=(float(j-2)+0.5)*dy
  END DO

  CALL setdxdy(nx,ny,ibgn,iend,jbgn,jend, &
           xs,ys,dxfld,dyfld,rdxfld,rdyfld)

  IF(myproc == root) THEN
    IF (lvldbg > 0) THEN
      WRITE(6,'(a)') 'Set dx dy' 
    END IF
  ENd IF

  IF ( realdata .eq. 0 ) THEN
    membknt=min(nmembers,5)
  END IF

  IF( myproc == root ) THEN
    IF (lvldbg > 0) THEN
      WRITE(6,'(a,i5,a,i5)') ' MPI nprocs:',nprocs,'  myproc:',myproc
    END IF
    CALL CPU_TIME(cput1)
    WRITE(6,'(a,f10.2,a)') ' Start-up CPU time: ',(cput1-cput0),' seconds'
  END IF

  IF( myproc == root) CALL CPU_TIME(cput1)
  ensfcst(:,:,:) = 0.

  IF( myproc == root) THEN
!
!---------------------------
!  realdata - 0 : analytic case
!  realdata - 1 : CAPS Ensemble
!  realdata - 2 : HREF
!  realdata - 3 : GEFS
!---------------------------
!
    IF ( realdata == 0 ) THEN
      piov2 = 2.0*atan(1.0)
      refmax=50.
      ctrx0=60.0E03
      ctry0=60.0E03
      ctrx1=90.0E03
      ctry1=40.0E03
      radx=30.0E03
      rady=27.0E03
  
      radx2inv=1.0/(radx*radx)
      rady2inv=1.0/(rady*rady)
  
      DO k=1,nmembers
        IF(k==1) THEN
          ctrx=ctrx0-12.0E03
          ctry=ctry0
        ELSE IF(k==2) THEN
          ctrx=ctrx0+12.0E03
          ctry=ctry0 
        ELSE IF(k==3) THEN
          ctrx=ctrx0
          ctry=ctry0+12.0E03
        ELSE IF(k==4) THEN
          ctrx=ctrx1
          ctry=ctry1-12.0E03
        ELSE IF(k==5) THEN
          ctrx=ctrx1+12.0E03
          ctry=ctry1+12.0E03
        END IF
        DO j=1,ny-1
          DO i=1,nx-1
            delx=xs(i)-ctrx
            dely=ys(j)-ctry
            radius=SQRT(radx2inv*(delx*delx) + rady2inv*(dely*dely))
            IF(radius < 1.0) THEN
              ensfcst(i,j,k)=refmax*(COS(piov2*radius)**2)
            END IF
          END DO
        END DO
      END DO

    ELSE IF (realdata == 1) THEN
      DO k=1,nmembers
        WRITE(finname,'(a,i4.4,3(i2.2),a,a,a,i4.4,3(i2.2),a,i3.3,a)') trim(basedir), &
          odate(1), odate(2), odate(3), odate(4), '/', trim(membname(k)), '_', odate(1), odate(2), odate(3), odate(4), &
          'f', ifhr, '.grib2'
        WRITE(6,'(2a)') ' read file: ', trim(finname)

        call codes_open_file(ifile,finname,'r',istat)
        if (istat .eq. 0) then
          call codes_index_create(idx,finname,'shortName,level,lengthOfTimeRange',istat)

          call codes_index_select(idx,'shortName',varname,istat)
          call codes_index_select(idx,'level',0,istat)
          call codes_index_select(idx,'lengthOfTimeRange',acctime,istat)

          call codes_new_from_index(idx,igrib,istat)

          allocate(values(nx*ny))
          call codes_get(igrib,'values',values,istat)

          if (istat .eq. 0) then
            membknt = membknt + 1
            call codes_get(igrib,'missingValue',missing,istat)

            do i=1,nx
              do j=1,ny
                if (values((j-1)*nx+i).eq.missing) then
                  ensfcst(i,j,membknt) = 0.
                else
                  ensfcst(i,j,membknt) = values((j-1)*nx+i)
                end if
              enddo
            enddo

            ensfname(membknt) = membname(k)
          end if

          if (allocated(values)) then
            deallocate(values)
          end if
          call codes_release(igrib,istat)
          call codes_close_file(ifile,istat)
        end if
      END DO
    ELSE IF (realdata == 2) THEN
      DO k=1,nmembers
        WRITE(finname,'(a,i4.4,2(i2.2),a,a,a,i4.4,3(i2.2),a,i3.3,a)') trim(basedir), &
          odate(1), odate(2), odate(3), '/', trim(membname(k)), '_', odate(1), odate(2), odate(3), odate(4), &
          'f', ifhr, '.grib2'
        WRITE(6,'(2a)') ' read file: ', trim(finname)

        call codes_open_file(ifile,finname,'r',istat)
        if (istat .eq. 0) then
          IF (membname(k)(1:3) .NE. 'nam') THEN
            call codes_index_create(idx,finname,'shortName,level,lengthOfTimeRange',istat)
          ELSE
            call codes_index_create(idx,finname,'shortName,level',istat)
          END IF

          call codes_index_select(idx,'shortName',varname,istat)
          call codes_index_select(idx,'level',0,istat)
          IF (membname(k)(1:3) .NE. 'nam') THEN
            call codes_index_select(idx,'lengthOfTimeRange',ifhr,istat)
          END IF

          call codes_new_from_index(idx,igrib,istat)

          allocate(values(nx*ny))
          call codes_get(igrib,'values',values,istat)

          if (istat .eq. 0) then
            membknt = membknt + 1
            call codes_get(igrib,'missingValue',missing,istat)

            do i=1,nx
              do j=1,ny
                if (values((j-1)*nx+i).eq.missing) then
                  ensfcst(i,j,membknt) = 0.
                else
                  ensfcst(i,j,membknt) = values((j-1)*nx+i)
                end if
              enddo
            enddo

            ensfname(membknt) = membname(k)
          end if

          if (allocated(values)) then
            deallocate(values)
          end if
          call codes_release(igrib,istat)
          call codes_close_file(ifile,istat)

          IF (membname(k)(1:3) .NE. 'nam') THEN
            IF (ifhr-acctime .gt. 0) THEN
              WRITE(finname,'(a,i4.4,2(i2.2),a,a,a,i4.4,3(i2.2),a,i3.3,a)') trim(basedir), &
                odate(1), odate(2), odate(3), '/', trim(membname(k)), '_', odate(1), odate(2), odate(3), odate(4), &
                'f', ifhr-acctime, '.grib2'
              WRITE(6,'(2a)') ' read file: ', trim(finname)

              call codes_open_file(ifile,finname,'r',istat)
              if (istat .eq. 0) then
                call codes_index_create(idx,finname,'shortName,level,lengthOfTimeRange',istat)
    
                call codes_index_select(idx,'shortName',varname,istat)
                call codes_index_select(idx,'level',0,istat)
                call codes_index_select(idx,'lengthOfTimeRange',ifhr-acctime,istat)
   
                call codes_new_from_index(idx,igrib,istat)
 
                allocate(values(nx*ny))
                call codes_get(igrib,'values',values,istat)

                if (istat .eq. 0) then
                  call codes_get(igrib,'missingValue',missing,istat)

                  do i=1,nx
                    do j=1,ny
                      if (values((j-1)*nx+i).ne.missing) then
                        ensfcst(i,j,membknt) = ensfcst(i,j,membknt) - values((j-1)*nx+i)
                      end if
                    enddo
                  enddo
                end if

                if (allocated(values)) then
                  deallocate(values)
                end if
                call codes_release(igrib,istat)
                call codes_close_file(ifile,istat)
              end if
            END IF

          !! NAM Nest
          ELSE IF (acctime .EQ. 6) THEN

            WRITE(finname,'(a,i4.4,2(i2.2),a,a,a,i4.4,3(i2.2),a,i3.3,a)') trim(basedir), &
              odate(1), odate(2), odate(3), '/', trim(membname(k)), '_', odate(1), odate(2), odate(3), odate(4), &
              'f', ifhr-3, '.grib2'
            WRITE(6,'(2a)') ' read file: ', trim(finname)

            call codes_open_file(ifile,finname,'r',istat)
            if (istat .eq. 0) then
              call codes_index_create(idx,finname,'shortName,level',istat)
    
              call codes_index_select(idx,'shortName',varname,istat)
              call codes_index_select(idx,'level',0,istat)
    
              call codes_new_from_index(idx,igrib,istat)
 
              allocate(values(nx*ny))
              call codes_get(igrib,'values',values,istat)

              if (istat .eq. 0) then
                call codes_get(igrib,'missingValue',missing,istat)

                do i=1,nx
                  do j=1,ny
                    if (values((j-1)*nx+i).ne.missing) then
                      ensfcst(i,j,membknt) = ensfcst(i,j,membknt) + values((j-1)*nx+i)
                    end if
                  enddo
                enddo
              end if

              if (allocated(values)) then
                deallocate(values)
              end if
              call codes_release(igrib,istat)
              call codes_close_file(ifile,istat)
            end if

          END IF !! NAM Nest

        end if !! istat

        IF (k.eq.nmembers/2) THEN
          ifhr = ifhr + 12
          CALL GETH_NEWDATE(odate, odate, -12*3600)
        END IF

      END DO !! k=1,nmembers

      ifhr = ifhr - 12
    END IF !! realdata

    DO ipass=1,nshfpass
      ensfcst_smt(:,:,:,ipass) = ensfcst(:,:,:)
      DO i=1,nbaksmth(ipass)
        DO k=1,membknt
          CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,smt_param,          &
                        ensfcst(:,:,k),tem2d1,ensfcst_smt(:,:,k,ipass))
        END DO
      END DO
    END DO

  END IF !! myproc == root

  IF ( nprocs > 1 ) THEN
    CALL MPI_BCAST(membknt,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    nelem3d=nx*ny*nmembers
    nelem2d=nx*ny
    DO ipass=1,nshfpass
      IF( nelem3d < max_elem_send) THEN
        CALL MPI_BCAST(ensfcst_smt(:,:,:,ipass),nelem3d,MPI_REAL,root,MPI_COMM_WORLD,ierr)
      ELSE IF( nelem2d < max_elem_send) THEN
        DO k=1,nmembers
          CALL MPI_BCAST(ensfcst_smt(:,:,k,ipass),nelem2d,MPI_REAL,root, &
                     MPI_COMM_WORLD,ierr)
        END DO
      ELSE
        DO k=1,nmembers
          CALL LGARRAY_BCASTR(ensfcst_smt(:,:,k,ipass),nelem2d,max_elem_send, &
                          myproc,root,MPI_COMM_WORLD,ierr)
        END DO
      END IF
    END DO

    IF( myproc == root) THEN
      WRITE(6,'(a)') ' ensfcst broadcast complete'
    END IF
  END IF

  IF ( membknt .lt. min_nmem) THEN
    IF( myproc == root) THEN
      WRITE(6,'(a,i3)') ' membknt:', membknt
      WRITE(6,'(a)') ' not enough ensemble member'
    END IF

    STOP
  ELSE
    IF( myproc == root) THEN
      WRITE(6,'(a,i3)') ' membknt:', membknt
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate simple mean and max
!
!-----------------------------------------------------------------------
!
  IF( myproc == root ) THEN
    CALL CPU_TIME(cput2)
    WRITE(6,'(a,f10.2,a)') ' Set-up CPU time: ',(cput2-cput1),' seconds'
    IF (lvldbg > 0) THEN
      WRITE(6,'(a,i5)') ' Align Option: ',alignopt
    END IF

    ensmean(:,:) = rmisg_data
    ensmean(ibgn:iend,jbgn:jend) = 0.
    DO k=1,membknt
      fcstmax=-999.
      fcstmin=999.
      DO j=jbgn,jend
        DO i=ibgn,iend
          fcstmin=min(fcstmin,ensfcst(i,j,k))
          fcstmax=max(fcstmax,ensfcst(i,j,k))
          ensmean(i,j)=ensmean(i,j)+ensfcst(i,j,k)
          ensmax(i,j)=max(ensmax(i,j),ensfcst(i,j,k))
        END DO
      END DO
    END DO
    rninv=1.0/float(membknt)
    fcstmax=-999.
    fcstmin=999.
    DO j=jbgn,jend
      DO i=ibgn,iend
        ensmean(i,j)=rninv*ensmean(i,j)
        IF(lposdef) ensmean(i,j)=max(0.,ensmean(i,j))
        favgmin=min(favgmin,ensmean(i,j))
        favgmax=max(favgmax,ensmean(i,j))
      END DO
    END DO

    IF( alignopt == 2 ) THEN
      CALL pm_mean(nx,ny,membknt,RESHAPE(ensfcst(:,:,1:membknt),(/nx*ny*membknt/)), &
              ensmean,enspm(:,:))

      enspm_smt(:,:,1) = enspm(:,:)
      DO i=1,nbaksmth(1)
        CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,smt_param, &
                  enspm(:,:),tem2d1,enspm_smt(:,:,1))
      END DO
    END IF

    CALL CPU_TIME(cput3)
    IF (lvldbg > 0) THEN
      WRITE(6,'(a,f10.2,a)') ' Simple mean CPU time: ',              &
           (cput3-cput2),' seconds'
    END IF
  END IF

!! Broadcast PM mean for alignopt = 2 (align to ensemble mean location)
  IF( alignopt == 2 ) THEN
    IF( nprocs > 1 ) THEN
      nelem2d=nx*ny

      IF( nelem2d < max_elem_send) THEN
        CALL MPI_BCAST(enspm_smt(:,:,1),nelem2d,MPI_REAL,root, &
                 MPI_COMM_WORLD,ierr)
      ELSE
        CALL LGARRAY_BCASTR(enspm_smt(:,:,1),nelem2d,max_elem_send, &
                    myproc,root,MPI_COMM_WORLD,ierr)
      END IF
    END IF
  END IF

!
!-----------------------------------------------------------------------
!
!   Calculate shift vectors
!   alignopt = 1   Method A : Align to all individual members
!   alignopt = 2   Method B : Align to ensemble PM mean location
!
!-----------------------------------------------------------------------
!
  IF( alignopt == 1) THEN
    ensfcst_shf(:,:,:,:) = 0.
    ibkshift(:,:,:,:) = 0
    jbkshift(:,:,:,:) = 0
    fnorm=1.0/float(membknt)

    DO ipass=1,nshfpass
      DO k1=1,membknt
        !IF( myproc == root ) WRITE(6,'(a,i4)') ' Processing ensemble member: ',k1
        xshiftmn(:,:,k1,ipass) = 0.
        yshiftmn(:,:,k1,ipass) = 0.

        time = real(ifhr)
        !! HREF (realdata == 2) uses time lagged members
        IF (realdata == 2) THEN
          IF (k1.gt.membknt/2) THEN
            time = real(ifhr) + 12.
          END IF
        END IF

        IF( myproc == root ) THEN
           WRITE(6,'(a,i4,a,i4)') ' Processing ensemble member: ',k1, ', pass number: ', ipass
        END IF

        DO k2=1,membknt
          IF( k2 /= k1 ) THEN
            IF ( myproc == root ) THEN
              IF (lvldbg > 0) THEN
                print *, ''
                print *, k1, k2          
              END IF
            END IF
            xshift(:,:) = 0.
            yshift(:,:) = 0.
            CALL rshift2dgrd(nx,ny,nvarshf,ipass,mx_shfpass,mxzone,    &
               max_elem_send,applyshft,posdef,time,xs,ys,              &
               ensfcst_smt(:,:,k1,ipass),ensfcst_smt(:,:,k2,ipass),    &
               ensfcst_shf(:,:,k1,ipass),                              &
               istart,jstart,ifinish,jfinish,                          &
               ibkshift(:,:,k1,ipass),jbkshift(:,:,k1,ipass),          &
               ibkshift(:,:,k2,ipass),jbkshift(:,:,k2,ipass),          &
               xshift,yshift,xsum,ysum,wgtsum,                         &
               dxfld,dyfld,rdxfld,rdyfld,                              &
               slopey,alphay,betay,                                    &
               tem2d3,recv_buf,lvldbg)
             
             IF( myproc == root ) THEN
               xshiftmn(:,:,k1,ipass)=xshiftmn(:,:,k1,ipass)+xshift(:,:)*fnorm
               yshiftmn(:,:,k1,ipass)=yshiftmn(:,:,k1,ipass)+yshift(:,:)*fnorm
             END IF

          END IF
        END DO
!
!-----------------------------------------------------------------------
!
! Calculate mean shift for ensemble member k1. 
! Note that zero shift is implied for k1,k1, so denominator is membknt.
! 
! + Save backshift for next pass
!
!-----------------------------------------------------------------------
!
        IF( myproc == root ) THEN
          IF (ipass .ne. nshfpass) THEN
            ibkshift(:,:,k1,ipass+1)=ibkshift(:,:,k1,ipass)+NINT(xshiftmn(:,:,k1,ipass))
            jbkshift(:,:,k1,ipass+1)=jbkshift(:,:,k1,ipass)+NINT(yshiftmn(:,:,k1,ipass))
          END IF

          IF (ipass .gt. 1) THEN
            xshiftmn(:,:,k1,ipass)=xshiftmn(:,:,k1,ipass)+xshiftmn(:,:,k1,ipass-1)
            yshiftmn(:,:,k1,ipass)=yshiftmn(:,:,k1,ipass)+yshiftmn(:,:,k1,ipass-1)
          END IF

          DO i=1,nshfsmth
            CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                        xshiftmn(:,:,k1,ipass),tem2d3,xshiftmn(:,:,k1,ipass))
            CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                        yshiftmn(:,:,k1,ipass),tem2d3,yshiftmn(:,:,k1,ipass))
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
! Apply shift to ensemble using mean shift and sum 
! to calculate ensemble mean (and broadcast ensemble mean)
!
!-----------------------------------------------------------------------
!
        IF( myproc == root ) THEN
          IF (lvldbg > 0) THEN
            print *, ''
            print *, ''   
          END IF
          IF (lvldbg > 10) THEN
            CALL a2dmax0(ensfcst(:,:,k1),1,nx,ibgn,iend,           &
                         1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))')                             &
                  'Pre-shift min = ', vmin,',  max =',vmax
          END IF
          CALL movegr(nx,ny,ensfcst(:,:,k1),tem2d3,                    &
                    xshiftmn(:,:,k1,ipass),yshiftmn(:,:,k1,ipass),     &
                    ensfcst_shf(:,:,k1,ipass),                         &
                    ibgn,iend,jbgn,jend,                               &
                    dxfld,dyfld,rdxfld,rdyfld,                         &
                    slopey,alphay,betay)

          !! do smooth to shifted fields
          DO i=1,noutsmth
            CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                        ensfcst_shf(:,:,k1,ipass),tem2d3,ensfcst_shf(:,:,k1,ipass))
          END DO

          !! ** IMPORTANT **
          !! rescale shifted field with its original PDF
          CALL pm_mean(nx,ny,1,RESHAPE(ensfcst(:,:,k1),(/nx*ny/)), &
            ensfcst_shf(:,:,k1,ipass),ensfcst_shf(:,:,k1,ipass))

          IF (lvldbg > 10) THEN
            CALL a2dmax0(ensfcst_shf(:,:,k1,ipass),1,nx,ibgn,iend,       &
                       1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))')                                   &
                  'Post-shift min = ', vmin,',  max =',vmax
          END IF

          ensshfmn(:,:,ipass)=ensshfmn(:,:,ipass)+ensfcst_shf(:,:,k1,ipass)
        END IF

        IF( nprocs > 1 ) THEN
          nelem2d=nx*ny

          IF (ipass .ne. nshfpass) THEN
            IF( nelem2d < max_elem_send) THEN
              CALL MPI_BCAST(ibkshift(:,:,k1,ipass+1),nelem2d,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
            ELSE
              CALL LGARRAY_INT_BCASTR(ibkshift(:,:,k1,ipass+1),nelem2d,max_elem_send, &
                          myproc,root,MPI_COMM_WORLD,ierr)
            END IF

            IF( nelem2d < max_elem_send) THEN
              CALL MPI_BCAST(jbkshift(:,:,k1,ipass+1),nelem2d,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
            ELSE
              CALL LGARRAY_INT_BCASTR(jbkshift(:,:,k1,ipass+1),nelem2d,max_elem_send, &
                          myproc,root,MPI_COMM_WORLD,ierr)
            END IF

            IF( myproc == root) THEN
              IF (lvldbg > 0) THEN
                WRITE(6,'(a)') ' back-shift broadcast complete'
              END IF
            END IF
          END IF

        END IF
      END DO  ! k1

!
!-----------------------------------------------------------------------
!
! Calculate shifted ensemble mean and max
!
!-----------------------------------------------------------------------
!
      IF( myproc == root ) THEN
        fnorm=1.0/float(membknt)
        ensshfmn(:,:,ipass)=fnorm*ensshfmn(:,:,ipass)
        IF(lposdef) THEN
          DO j=jbgn,jend
            DO i=ibgn,iend
              ensshfmn(i,j,ipass)=max(0.,ensshfmn(i,j,ipass))
              DO k=1,membknt
                ensshfmax(i,j,ipass)=max(ensshfmax(i,j,ipass),ensfcst_shf(i,j,k,ipass))
              END DO
            END DO
          END DO
        END IF
      END IF

    END DO  ! ipass

    IF( myproc == root ) THEN
      CALL CPU_TIME(cput4)
      WRITE(6,'(a,f10.2,a)') ' Ensemble aligned mean CPU time: ',(cput4-cput3),' seconds'
    END IF
   
  ELSE IF( alignopt == 2) THEN
    ensfcst_shf(:,:,:,:) = 0.
    ibkshift(:,:,:,:) = 0
    jbkshift(:,:,:,:) = 0
    fnorm=1.0/float(membknt)
    zero2d(:,:) = 0

    DO ipass=1,nshfpass
      DO k1=1,membknt
        !IF( myproc == root ) WRITE(6,'(a,i4)') ' Processing ensemble member: ',k1
        xshiftmn(:,:,k1,ipass) = 0.
        yshiftmn(:,:,k1,ipass) = 0.

        time = real(ifhr)
        !! HREF (realdata == 2) uses time lagged members
        IF (realdata == 2) THEN
          IF (k1.gt.membknt/2) THEN
            time = real(ifhr) + 12.
          END IF
        END IF

        IF( myproc == root ) THEN
          !WRITE(6,'(a,i4,a,i4)') ' Processing ensemble member: ',k1, ', leadtime: ', nint(time)
           WRITE(6,'(a,i4,a,i4)') ' Processing ensemble member: ',k1, ', pass number: ', ipass
        END IF

        xshift(:,:) = 0.
        yshift(:,:) = 0.
        CALL rshift2dgrd(nx,ny,nvarshf,ipass,mx_shfpass,mxzone,    &
               max_elem_send,applyshft,posdef,time,xs,ys,              &
               ensfcst_smt(:,:,k1,ipass),enspm_smt(:,:,ipass),         &
               ensfcst_shf(:,:,k1,ipass),                              &
               istart,jstart,ifinish,jfinish,                          &
               ibkshift(:,:,k1,ipass),jbkshift(:,:,k1,ipass),          &
               zero2d(:,:),zero2d(:,:),                                &
               xshift,yshift,xsum,ysum,wgtsum,                         &
               dxfld,dyfld,rdxfld,rdyfld,                              &
               slopey,alphay,betay,                                    &
               tem2d3,recv_buf,lvldbg)
             
!
!-----------------------------------------------------------------------
!
! Calculate mean shift for ensemble member k1. 
! Note that zero shift is implied for k1,k1, so denominator is membknt.
! 
! + Save backshift for next pass
!
!-----------------------------------------------------------------------
!
        IF( myproc == root ) THEN
           xshiftmn(:,:,k1,ipass)=xshift(:,:)
           yshiftmn(:,:,k1,ipass)=yshift(:,:)

          IF (ipass .ne. nshfpass) THEN
            ibkshift(:,:,k1,ipass+1)=ibkshift(:,:,k1,ipass)+NINT(xshiftmn(:,:,k1,ipass))
            jbkshift(:,:,k1,ipass+1)=jbkshift(:,:,k1,ipass)+NINT(yshiftmn(:,:,k1,ipass))
          END IF

          IF (ipass .gt. 1) THEN
            xshiftmn(:,:,k1,ipass)=xshiftmn(:,:,k1,ipass)+xshiftmn(:,:,k1,ipass-1)
            yshiftmn(:,:,k1,ipass)=yshiftmn(:,:,k1,ipass)+yshiftmn(:,:,k1,ipass-1)
          END IF

          DO i=1,nshfsmth
            CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                        xshiftmn(:,:,k1,ipass),tem2d3,xshiftmn(:,:,k1,ipass))
            CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                        yshiftmn(:,:,k1,ipass),tem2d3,yshiftmn(:,:,k1,ipass))
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
! Apply shift to ensemble using mean shift and sum 
! to calculate ensemble mean (and broadcast ensemble mean)
!
!-----------------------------------------------------------------------
!
        IF( myproc == root ) THEN
          IF (lvldbg > 0) THEN
            print *, ''
            print *, ''   
          END IF
          IF (lvldbg > 10) THEN
            CALL a2dmax0(ensfcst(:,:,k1),1,nx,ibgn,iend,           &
                         1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))')                             &
                  'Pre-shift min = ', vmin,',  max =',vmax
          END IF
          CALL movegr(nx,ny,ensfcst(:,:,k1),tem2d3,                    &
                    xshiftmn(:,:,k1,ipass),yshiftmn(:,:,k1,ipass),     &
                    ensfcst_shf(:,:,k1,ipass),                         &
                    ibgn,iend,jbgn,jend,                               &
                    dxfld,dyfld,rdxfld,rdyfld,                         &
                    slopey,alphay,betay)

          !! do smooth to shifted fields
          DO i=1,noutsmth
            CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                        ensfcst_shf(:,:,k1,ipass),tem2d3,ensfcst_shf(:,:,k1,ipass))
          END DO

          !! ** IMPORTANT **
          !! rescale shifted field with its original PDF
          CALL pm_mean(nx,ny,1,RESHAPE(ensfcst(:,:,k1),(/nx*ny/)), &
            ensfcst_shf(:,:,k1,ipass),ensfcst_shf(:,:,k1,ipass))

          IF (lvldbg > 10) THEN
            CALL a2dmax0(ensfcst_shf(:,:,k1,ipass),1,nx,ibgn,iend,       &
                       1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))')                                   &
                  'Post-shift min = ', vmin,',  max =',vmax
          END IF

          ensshfmn(:,:,ipass)=ensshfmn(:,:,ipass)+ensfcst_shf(:,:,k1,ipass)
        END IF

        IF( nprocs > 1 ) THEN
          nelem2d=nx*ny

          IF (ipass .ne. nshfpass) THEN
            IF( nelem2d < max_elem_send) THEN
              CALL MPI_BCAST(ibkshift(:,:,k1,ipass+1),nelem2d,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
            ELSE
              CALL LGARRAY_INT_BCASTR(ibkshift(:,:,k1,ipass+1),nelem2d,max_elem_send, &
                          myproc,root,MPI_COMM_WORLD,ierr)
            END IF

            IF( nelem2d < max_elem_send) THEN
              CALL MPI_BCAST(jbkshift(:,:,k1,ipass+1),nelem2d,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
            ELSE
              CALL LGARRAY_INT_BCASTR(jbkshift(:,:,k1,ipass+1),nelem2d,max_elem_send, &
                          myproc,root,MPI_COMM_WORLD,ierr)
            END IF

            IF( myproc == root) THEN
              IF (lvldbg > 0) THEN
                WRITE(6,'(a)') ' back-shift broadcast complete'
              END IF
            END IF
          END IF

        END IF
      END DO  ! k1

!
!-----------------------------------------------------------------------
!
! Calculate shifted ensemble mean and PM
!
!-----------------------------------------------------------------------
!
      IF( myproc == root ) THEN
        fnorm=1.0/float(membknt)
        ensshfmn(:,:,ipass)=fnorm*ensshfmn(:,:,ipass)
        IF(lposdef) THEN
          DO j=jbgn,jend
            DO i=ibgn,iend
              ensshfmn(i,j,ipass)=max(0.,ensshfmn(i,j,ipass))
              DO k=1,membknt
                ensshfmax(i,j,ipass)=max(ensshfmax(i,j,ipass),ensfcst_shf(i,j,k,ipass))
              END DO
            END DO
          END DO
        END IF

        CALL pm_mean(nx,ny,membknt,RESHAPE(ensfcst_shf(:,:,1:membknt,ipass),(/nx*ny*membknt/)), &
                ensshfmn(:,:,ipass),ensshfpm(:,:,ipass))
        IF (lvldbg > 10) THEN
          CALL a2dmax0(ensshfpm(:,:,ipass),1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
          WRITE(6,'(1x,2(a,f13.4))') 'Shift PM output min = ', vmin,',  max=',vmax
        END IF

        IF (ipass .ne. nshfpass) THEN
          enspm_smt(:,:,ipass+1) = ensshfpm(:,:,ipass)

          DO i=1,nbaksmth(ipass+1)
            CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,smt_param,          &
                      ensshfpm(:,:,ipass),tem2d1,enspm_smt(:,:,ipass+1))
          END DO
        END IF
      END IF

      IF (ipass .ne. nshfpass) THEN
        IF( nprocs > 1 ) THEN
          nelem2d=nx*ny

          IF( nelem2d < max_elem_send) THEN
            CALL MPI_BCAST(enspm_smt(:,:,ipass+1),nelem2d,MPI_REAL,root, &
                     MPI_COMM_WORLD,ierr)
          ELSE
            CALL LGARRAY_BCASTR(enspm_smt(:,:,ipass+1),nelem2d,max_elem_send, &
                        myproc,root,MPI_COMM_WORLD,ierr)
          END IF
        END IF
      END IF

    END DO  ! ipass

    IF( myproc == root ) THEN
      CALL CPU_TIME(cput4)
      WRITE(6,'(a,f10.2,a)') ' Ensemble aligned mean CPU time: ',(cput4-cput3),' seconds'
    END IF

  END IF

!
!-----------------------------------------------------------------------
!
! Calculate LPM and PM
!
!-----------------------------------------------------------------------
!
  IF( myproc == root ) THEN
    CALL CPU_TIME(cput5)
  END IF

  IF( nprocs > 1 ) THEN
    nelem3d=nx*ny*nmembers
    nelem2d=nx*ny

    IF( nelem3d < max_elem_send) THEN
      CALL MPI_BCAST(ensfcst(:,:,:),nelem3d,MPI_REAL,root,MPI_COMM_WORLD,ierr)
    ELSE IF( nelem2d < max_elem_send) THEN
      DO k=1,nmembers
        CALL MPI_BCAST(ensfcst(:,:,k),nelem2d,MPI_REAL,root, &
                   MPI_COMM_WORLD,ierr)
      END DO
    ELSE
      DO k=1,nmembers
        CALL LGARRAY_BCASTR(ensfcst(:,:,k),nelem2d,max_elem_send, &
                        myproc,root,MPI_COMM_WORLD,ierr)
      END DO
    END IF
 
    IF( nelem2d < max_elem_send) THEN
      CALL MPI_BCAST(ensmean(:,:),nelem2d,MPI_REAL,root, &
               MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ensmax(:,:),nelem2d,MPI_REAL,root, &
               MPI_COMM_WORLD,ierr)
    ELSE
      CALL LGARRAY_BCASTR(ensmean(:,:),nelem2d,max_elem_send, &
                  myproc,root,MPI_COMM_WORLD,ierr)
      CALL LGARRAY_BCASTR(ensmax(:,:),nelem2d,max_elem_send, &
                  myproc,root,MPI_COMM_WORLD,ierr)
    END IF
  END IF

  CALL lpm_mean_la(nx,ny,ibgn,iend,jbgn,jend,      &
            membknt,patch_nx,patch_ny,ovx,ovy,     &
            filt_min,float(gauss_sigma),gs_weight, &
            ensfcst(:,:,1:membknt),ensmean,ensmax, &
            enslpm,myproc,nprocs)

  IF( myproc == root ) THEN
    CALL CPU_TIME(cput6)
    WRITE(6,'(a,f10.2,a)') ' Ensemble LPM CPU time: ',(cput6-cput5),' seconds'
    CALL CPU_TIME(cput5)
  END IF

  IF( nprocs > 1 ) THEN
    nelem3d=nx*ny*nmembers
    nelem2d=nx*ny
    DO ipass=1,nshfpass
      IF( nelem3d < max_elem_send) THEN
        CALL MPI_BCAST(ensfcst_shf(:,:,:,ipass),nelem3d,MPI_REAL,root,MPI_COMM_WORLD,ierr)
      ELSE IF( nelem2d < max_elem_send) THEN
        DO k=1,nmembers
          CALL MPI_BCAST(ensfcst_shf(:,:,k,ipass),nelem2d,MPI_REAL,root, &
                     MPI_COMM_WORLD,ierr)
        END DO
      ELSE
        DO k=1,nmembers
          CALL LGARRAY_BCASTR(ensfcst_shf(:,:,k,ipass),nelem2d,max_elem_send, &
                          myproc,root,MPI_COMM_WORLD,ierr)
        END DO
      END IF

      IF( nelem2d < max_elem_send) THEN
        CALL MPI_BCAST(ensshfmn(:,:,ipass),nelem2d,MPI_REAL,root, &
                 MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(ensshfmax(:,:,ipass),nelem2d,MPI_REAL,root, &
                 MPI_COMM_WORLD,ierr)
      ELSE
        CALL LGARRAY_BCASTR(ensshfmn(:,:,ipass),nelem2d,max_elem_send, &
                    myproc,root,MPI_COMM_WORLD,ierr)
        CALL LGARRAY_BCASTR(ensshfmax(:,:,ipass),nelem2d,max_elem_send, &
                    myproc,root,MPI_COMM_WORLD,ierr)
      END IF
    END DO
  END IF

  DO ipass=1,nshfpass
    CALL lpm_mean_la(nx,ny,ibgn,iend,jbgn,jend, &
            membknt,patch_nx,patch_ny,ovx,ovy,  &
            filt_min,float(gauss_sigma),gs_weight,        &
            ensfcst_shf(:,:,1:membknt,ipass),ensshfmn(:,:,ipass),ensshfmax(:,:,ipass), &
            ensshflpm(:,:,ipass),myproc,nprocs)

  END DO

  IF( myproc == root ) THEN
    CALL CPU_TIME(cput6)
    WRITE(6,'(a,f10.2,a)') ' Ensemble aligned LPM CPU time:',(cput6-cput5),' seconds'
  END IF

  IF( myproc == root ) THEN
    IF( alignopt == 1) THEN
      CALL pm_mean(nx,ny,membknt,RESHAPE(ensfcst(:,:,1:membknt),(/nx*ny*membknt/)), &
                ensmean,enspm(:,:))
      IF (lvldbg > 10) THEN
        CALL a2dmax0(enspm(:,:),1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
        WRITE(6,'(1x,2(a,f13.4))') 'PM output min = ', vmin,',  max =',vmax
      END IF

      DO ipass=1,nshfpass
        CALL pm_mean(nx,ny,membknt,RESHAPE(ensfcst(:,:,1:membknt),(/nx*ny*membknt/)), &
                ensshfmn(:,:,ipass),ensshfpm(:,:,ipass))
        IF (lvldbg > 10) THEN
          CALL a2dmax0(ensshfpm(:,:,ipass),1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
          WRITE(6,'(1x,2(a,f13.4))') 'Shift PM output min = ', vmin,',  max=',vmax
        END IF
      END DO
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
! Write output file (NetCDF)
!
!-----------------------------------------------------------------------
!
  IF( myproc == root ) THEN

! Write file(netcdf)
    IF ( realdata == 0 ) THEN
      WRITE(foutname,'(a)') './analytic_case.nc'
    ELSE
      WRITE(foutname,'(5a,i3.3,a)') trim(shfoutdir), datehrstr(1:10), trim(fout_head), &
          datehrstr, '_f', ifhr, '.nc'
    END IF
    iret = nf_create(foutname, NF_NETCDF4, ncid)
    WRITE(6,'(2a)') ' write file: ', trim(foutname)

! Define dimensions
    ndim = 6
    nvarout = 18
    allocate(dimids(ndim), vardims(ndim), chunks(ndim), varids(nvarout))
    allocate(vardims2(ndim), chunks2(ndim), vardims3(2))

    cx = int(real(nx)/2)
    cy = int(real(ny)/2)
    allocate(projx(nx), projy(ny))
    do i=1,nx
      projx(i) = 3.*(i-1-cx)
    enddo
    do j=1,ny
      projy(j) = 3.*(j-1-cy)
    enddo

    iret = nf_def_dim(ncid, 'x', nx, dimids(1))
    iret = nf_def_dim(ncid, 'y', ny, dimids(2))
    iret = nf_def_dim(ncid, 'nmembers', membknt, dimids(3))
    iret = nf_def_dim(ncid, 'nshfpass', nshfpass, dimids(4))
    iret = nf_def_dim(ncid, 'time', 1, dimids(5))  
    iret = nf_def_dim(ncid, 'strlen', len(ensfname(1)), dimids(6))

    vardims(1) = dimids(1)
    vardims(2) = dimids(2)
    vardims(3) = dimids(3)
    vardims(4) = dimids(4)

    vardims2(1) = dimids(1)
    vardims2(2) = dimids(2)
    vardims2(3) = dimids(4)
    vardims2(4) = dimids(3)

    vardims3(1) = dimids(6)
    vardims3(2) = dimids(3)

    chunks(1) = nx
    chunks(2) = ny
    chunks(3) = membknt
    chunks(4) = nshfpass

    chunks2(1) = nx
    chunks2(2) = ny
    chunks2(3) = nshfpass
    chunks2(4) = membknt

! Define variables
    iret = nf_def_var(ncid, 'lat', nf_float, 2, vardims, varids(1))
    iret = nf_def_var(ncid, 'lon', nf_float, 2, vardims, varids(2))
    iret = nf_def_var(ncid, 'ensfcst', nf_float, 3, vardims, varids(3))
    iret = nf_def_var(ncid, 'shiftvec_u', nf_float, 4, vardims, varids(4))
    iret = nf_def_var(ncid, 'shiftvec_v', nf_float, 4, vardims, varids(5))
    iret = nf_def_var(ncid, 'ensshffcst', nf_float, 4, vardims, varids(6))
    iret = nf_def_var(ncid, 'ensmean', nf_float, 2, vardims, varids(7))
    iret = nf_def_var(ncid, 'ensshfmean', nf_float, 3, vardims2, varids(8))
    iret = nf_def_var(ncid, 'enslpm', nf_float, 2, vardims, varids(9))
    iret = nf_def_var(ncid, 'ensshflpm', nf_float, 3, vardims2, varids(10))
    iret = nf_def_var(ncid, 'enspm', nf_float, 2, vardims, varids(11))
    iret = nf_def_var(ncid, 'ensshfpm', nf_float, 3, vardims2, varids(12))
    iret = nf_def_var(ncid, 'time', nf_int, 1, dimids(5), varids(13))
    iret = nf_def_var(ncid, 'forecast_reference_time', nf_int, 1, dimids(5), varids(14))
    iret = nf_def_var(ncid, 'grid_mapping', nf_int, 0, dimids(5), varids(15))
    iret = nf_def_var(ncid, 'x', nf_float, 1, dimids(1), varids(16))
    iret = nf_def_var(ncid, 'y', nf_float, 1, dimids(2), varids(17))
    iret = nf_def_var(ncid, 'members', nf_char, 2, vardims3, varids(18))

! Turn on chunking.
    iret = nf_def_var_chunking(ncid, varids(1), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(2), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(3), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(4), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(5), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(6), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(7), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(8), nf_chunked, chunks2)
    iret = nf_def_var_chunking(ncid, varids(9), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(10), nf_chunked, chunks2)
    iret = nf_def_var_chunking(ncid, varids(11), nf_chunked, chunks)
    iret = nf_def_var_chunking(ncid, varids(12), nf_chunked, chunks2)

! Turn on deflate compression, fletcher32 checksum.
    do k = 1, nvarout
      iret = nf_def_var_deflate(ncid, varids(k), 1, 1, 4)
      iret = nf_def_var_fletcher32(ncid, varids(k), nf_fletcher32)
    enddo

! Put attributes
    iret = nf_put_att_text(ncid, varids(1), 'long_name', len_trim('latitude'), &
                           'latitude')
    iret = nf_put_att_text(ncid, varids(1), 'standard_name', len_trim('latitude'), &
                           'latitude')
    iret = nf_put_att_text(ncid, varids(1), 'units', len_trim('degrees_north'), 'degrees_north')
    iret = nf_put_att_real(ncid, varids(1), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(1), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(1), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    iret = nf_put_att_text(ncid, varids(2), 'long_name', len_trim('longitude'), &
                           'longitude')
    iret = nf_put_att_text(ncid, varids(2), 'standard_name', len_trim('longitude'), &
                           'longitude')
    iret = nf_put_att_text(ncid, varids(2), 'units', len_trim('degrees_east'), 'degrees_east')
    iret = nf_put_att_real(ncid, varids(2), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(2), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(2), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')
  
    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(3), 'long_name', &
                             len_trim('ensemble forecasts of precipitation'), &
                             'ensemble forecasts of precipitation')
      iret = nf_put_att_text(ncid, varids(3), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(3), 'long_name', &
                             len_trim('ensemble forecasts of snowfall'), &
                             'ensemble forecasts of snowfall')
      iret = nf_put_att_text(ncid, varids(3), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(3), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(3), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(3), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(3), 'grid_mapping',len_trim('grid_mapping'), &
                           'grid_mapping')

    iret = nf_put_att_text(ncid, varids(4), 'long_name', len_trim('shift u_vector'), &
                           'shift u_vector')
    iret = nf_put_att_text(ncid, varids(4), 'units', len_trim('grid'), 'grid')
    iret = nf_put_att_real(ncid, varids(4), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(4), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(4), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    iret = nf_put_att_text(ncid, varids(5), 'long_name', len_trim('shift v_vector'), &
                           'shift v_vector')
    iret = nf_put_att_text(ncid, varids(5), 'units', len_trim('grid'), 'grid')
    iret = nf_put_att_real(ncid, varids(5), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(5), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(5), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(6), 'long_name', &
                             len_trim('ensemble forecasts of shifted precipitation field'), &
                             'ensemble forecasts of shifted precipitation field')
      iret = nf_put_att_text(ncid, varids(6), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(6), 'long_name', &
                             len_trim('ensemble forecasts of shifted snowfall field'), &
                             'ensemble forecasts of shifted snowfall field')
      iret = nf_put_att_text(ncid, varids(6), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(6), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(6), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(6), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(6), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(7), 'long_name', &
                             len_trim('ensemble mean of precipitation'), &
                             'ensemble mean of precipitation')
      iret = nf_put_att_text(ncid, varids(7), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(7), 'long_name', &
                             len_trim('ensemble mean of snowfall'), &
                             'ensemble mean of snowfall')
      iret = nf_put_att_text(ncid, varids(7), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(7), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(7), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(7), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(7), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(8), 'long_name', &
                             len_trim('ensemble shifted mean of precipitation'), &
                             'ensemble shifted mean of precipitation')
      iret = nf_put_att_text(ncid, varids(8), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(8), 'long_name', &
                             len_trim('ensemble shifted mean of snowfall'), &
                             'ensemble shifted mean of snowfall')
      iret = nf_put_att_text(ncid, varids(8), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(8), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(8), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(8), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(8), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(9), 'long_name', &
                             len_trim('ensemble LPM mean of precipitation'), &
                             'ensemble LPM mean of precipitation')
      iret = nf_put_att_text(ncid, varids(9), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(9), 'long_name', &
                             len_trim('ensemble LPM mean of snowfall'), &
                             'ensemble LPM mean of snowfall')
      iret = nf_put_att_text(ncid, varids(9), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(9), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(9), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(9), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(9), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(10), 'long_name', &
                             len_trim('ensemble shifted LPM mean of precipitation'), &
                             'ensemble shifted LPM mean of precipitation')
      iret = nf_put_att_text(ncid, varids(10), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(10), 'long_name', &
                             len_trim('ensemble shifted LPM mean of snowfall'), &
                             'ensemble shifted LPM mean of snowfall')
      iret = nf_put_att_text(ncid, varids(10), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(10), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(10), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(10), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(10), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(11), 'long_name', &
                             len_trim('ensemble PM mean of precipitation'), &
                             'ensemble PM mean of precipitation')
      iret = nf_put_att_text(ncid, varids(11), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(11), 'long_name', &
                             len_trim('ensemble PM mean of snowfall'), &
                             'ensemble PM mean of snowfall')
      iret = nf_put_att_text(ncid, varids(11), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(11), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(11), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(11), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(11), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')
 
    if (varname == 'tp') then
      iret = nf_put_att_text(ncid, varids(12), 'long_name', &
                             len_trim('ensemble shifted PM mean of precipitation'), &
                             'ensemble shifted PM mean of precipitation')
      iret = nf_put_att_text(ncid, varids(12), 'units', len_trim('kg m-2'), 'kg m-2')
    else if (varname == 'asnow') then
      iret = nf_put_att_text(ncid, varids(12), 'long_name', &
                             len_trim('ensemble shifted PM mean of snowfall'), &
                             'ensemble shifted PM mean of snowfall')
      iret = nf_put_att_text(ncid, varids(12), 'units', len_trim('m'), 'm')
    end if
    iret = nf_put_att_real(ncid, varids(12), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_int(ncid, varids(12), 'lengthOfTimeRange', nf_int, 1, acctime)
    iret = nf_put_att_text(ncid, varids(12), 'coordinates', len_trim('lat lon'), 'lat lon')
    iret = nf_put_att_text(ncid, varids(12), 'grid_mapping', len_trim('grid_mapping'), &
                           'grid_mapping')

    iret = nf_put_att_text(ncid, varids(13), 'long_name', len_trim('time'), &
                           'time')
    iret = nf_put_att_text(ncid, varids(13), 'units', &
           len_trim('hours since 1970-01-01 00:00:00 0:00'), 'hours since 1970-01-01 00:00:00 0:00')
    iret = nf_put_att_text(ncid, varids(13), 'standard_name', len_trim('time'), &
                           'time')
    iret = nf_put_att_text(ncid, varids(13), 'axis', len_trim('T'), &
                           'T')

    iret = nf_put_att_text(ncid, varids(14), 'long_name', len_trim('forecast_reference_time'), &
                           'forecast_reference_time')
    iret = nf_put_att_text(ncid, varids(14), 'units', &
           len_trim('hours since 1970-01-01 00:00:00 0:00'), 'hours since 1970-01-01 00:00:00 0:00')
    iret = nf_put_att_text(ncid, varids(14), 'standard_name', len_trim('forecast_reference_time'), &
                           'forecast_reference_time')

    iret = nf_put_att_text(ncid, varids(15), 'grid_mapping_name', len_trim('lambert_conformal_conic'), &
                           'lambert_conformal_conic')
    iret = nf_put_att_real(ncid, varids(15), 'standard_parallel', nf_double, 1, 38.5)
    iret = nf_put_att_real(ncid, varids(15), 'longitude_of_central_meridian', nf_double, 1, 262.5)
    iret = nf_put_att_real(ncid, varids(15), 'latitude_of_projection_origin', nf_double, 1, 38.5)
    iret = nf_put_att_int(ncid, varids(15), 'false_easting', nf_int, 1, 0)
    iret = nf_put_att_int(ncid, varids(15), 'false_northing', nf_int, 1, 0)
    iret = nf_put_att_text(ncid, varids(15), 'GRIB_earth_shape', len_trim('spherical'), &
                           'spherical')
    iret = nf_put_att_int(ncid, varids(15), 'GRIB_earth_shape_code', nf_int, 1, 6)

    iret = nf_put_att_text(ncid, varids(16), 'units', len_trim('km'), &
                           'km')
    iret = nf_put_att_text(ncid, varids(16), 'long_name', len_trim('x coordinate of projection'), &
                           'x coordinate of projection')
    iret = nf_put_att_text(ncid, varids(16), 'standard_name', len_trim('projection_x_coordinate'), &
                           'projection_x_coordinate')
    iret = nf_put_att_text(ncid, varids(16), 'grid_spacing', len_trim('3 km'), &
                           '3 km')
    iret = nf_put_att_text(ncid, varids(16), 'axis', len_trim('X'), &
                           'X')

    iret = nf_put_att_text(ncid, varids(17), 'units', len_trim('km'), &
                           'km')
    iret = nf_put_att_text(ncid, varids(17), 'long_name', len_trim('y coordinate of projection'), &
                           'y coordinate of projection')
    iret = nf_put_att_text(ncid, varids(17), 'standard_name', len_trim('projection_y_coordinate'), &
                           'projection_y_coordinate')
    iret = nf_put_att_text(ncid, varids(17), 'grid_spacing', len_trim('3 km'), &
                           '3 km')
    iret = nf_put_att_text(ncid, varids(17), 'axis', len_trim('Y'), &
                           'Y')

    iret = nf_put_att_text(ncid, varids(18), 'long_name', len_trim('ensemble name'), &
                           'ensemble name')

    iret = nf_put_att_text(ncid, NCGLOBAL, 'Conventions', len_trim('CF-1.6'), 'CF-1.6')

    if (alignopt == 1) then
      iret = nf_put_att_text(ncid, NCGLOBAL, 'Alignment_Method', len_trim('Method A'), 'Method A')
    elseif (alignopt == 2) then
      iret = nf_put_att_text(ncid, NCGLOBAL, 'Alignment_Method', len_trim('Method B'), 'Method B')
    endif

    iret = nf_enddef(ncid)

! Put variables
    iret = nf_put_var_real(ncid, varids(1), lats(:,:))
    iret = nf_put_var_real(ncid, varids(2), lons(:,:))
    iret = nf_put_var_real(ncid, varids(3), ensfcst(:,:,:))
    iret = nf_put_var_real(ncid, varids(4), xshiftmn(:,:,:,:))
    iret = nf_put_var_real(ncid, varids(5), yshiftmn(:,:,:,:))
    iret = nf_put_var_real(ncid, varids(6), ensfcst_shf(:,:,:,:))
    iret = nf_put_var_real(ncid, varids(7), ensmean(:,:))
    iret = nf_put_var_real(ncid, varids(8), ensshfmn(:,:,:))
    iret = nf_put_var_real(ncid, varids(9), enslpm(:,:))
    iret = nf_put_var_real(ncid, varids(10), ensshflpm(:,:,:))
    iret = nf_put_var_real(ncid, varids(11), enspm(:,:))
    iret = nf_put_var_real(ncid, varids(12), ensshfpm(:,:,:))
    iret = nf_put_var_int(ncid, varids(13), itime)
    iret = nf_put_var_int(ncid, varids(14), rtime)
    iret = nf_put_var_int(ncid, varids(15), 0)
    iret = nf_put_var_real(ncid, varids(16), projx)
    iret = nf_put_var_real(ncid, varids(17), projy)
    iret = nf_put_var_text(ncid, varids(18), ensfname(:))  

! Close Netcdf File
    iret = nf_close(ncid)

    deallocate(dimids, vardims, chunks, varids, vardims2, chunks2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (fsout_opt == 1 .and. realdata .ne. 0) THEN
! Write file(netcdf) for transfer (extracted)
      WRITE(foutname,'(5a,i3.3,2a)') trim(shfsoutdir), datehrstr(1:10), trim(fout_head), &
            datehrstr, '_f', ifhr, trim(fsout_tail), '.nc'
      iret = nf_create(foutname, NF_NETCDF4, ncid)
      WRITE(6,'(2a)') ' write file: ', trim(foutname)

! Define dimensions
      ndim = 5
      nvarout = 10
      allocate(dimids(ndim), vardims(ndim), chunks(ndim), varids(nvarout), vardims2(ndim))

      iret = nf_def_dim(ncid, 'x', nx, dimids(1))
      iret = nf_def_dim(ncid, 'y', ny, dimids(2))
      iret = nf_def_dim(ncid, 'nmembers', membknt, dimids(3))
      iret = nf_def_dim(ncid, 'time', 1, dimids(4))  
      iret = nf_def_dim(ncid, 'strlen', len(ensfname(1)), dimids(5))

      vardims(1) = dimids(1)
      vardims(2) = dimids(2)

      vardims2(1) = dimids(5)
      vardims2(2) = dimids(3)

      chunks(1) = nx
      chunks(2) = ny

! Define variables
      iret = nf_def_var(ncid, 'lat', nf_float, 2, vardims, varids(1))
      iret = nf_def_var(ncid, 'lon', nf_float, 2, vardims, varids(2))
      iret = nf_def_var(ncid, 'ensshfmean', nf_float, 2, vardims, varids(3))
      iret = nf_def_var(ncid, 'ensshflpm', nf_float, 2, vardims, varids(4))
      iret = nf_def_var(ncid, 'time', nf_int, 1, dimids(4), varids(5))
      iret = nf_def_var(ncid, 'forecast_reference_time', nf_int, 1, dimids(4), varids(6))
      iret = nf_def_var(ncid, 'grid_mapping', nf_int, 0, dimids(4), varids(7))
      iret = nf_def_var(ncid, 'x', nf_float, 1, dimids(1), varids(8))
      iret = nf_def_var(ncid, 'y', nf_float, 1, dimids(2), varids(9))
      iret = nf_def_var(ncid, 'members', nf_char, 2, vardims2, varids(10))

! Turn on chunking.
      iret = nf_def_var_chunking(ncid, varids(1), nf_chunked, chunks)
      iret = nf_def_var_chunking(ncid, varids(2), nf_chunked, chunks)
      iret = nf_def_var_chunking(ncid, varids(3), nf_chunked, chunks)
      iret = nf_def_var_chunking(ncid, varids(4), nf_chunked, chunks)

! Turn on deflate compression, fletcher32 checksum.
      do k = 1, nvarout
        iret = nf_def_var_deflate(ncid, varids(k), 1, 1, 4)
        iret = nf_def_var_fletcher32(ncid, varids(k), nf_fletcher32)
      enddo

! Put attributes
      iret = nf_put_att_text(ncid, varids(1), 'long_name', len_trim('latitude'), &
                             'latitude')
      iret = nf_put_att_text(ncid, varids(1), 'standard_name', len_trim('latitude'), &
                             'latitude')
      iret = nf_put_att_text(ncid, varids(1), 'units', len_trim('degrees_north'), 'degrees_north')
      iret = nf_put_att_real(ncid, varids(1), '_FillValue', nf_float, 1, -999.)
      iret = nf_put_att_text(ncid, varids(1), 'coordinates', len_trim('lat lon'), 'lat lon')
      iret = nf_put_att_text(ncid, varids(1), 'grid_mapping', len_trim('grid_mapping'), &
                             'grid_mapping')

      iret = nf_put_att_text(ncid, varids(2), 'long_name', len_trim('longitude'), &
                             'longitude')
      iret = nf_put_att_text(ncid, varids(2), 'standard_name', len_trim('longitude'), &
                             'longitude')
      iret = nf_put_att_text(ncid, varids(2), 'units', len_trim('degrees_east'), 'degrees_east')
      iret = nf_put_att_real(ncid, varids(2), '_FillValue', nf_float, 1, -999.)
      iret = nf_put_att_text(ncid, varids(2), 'coordinates', len_trim('lat lon'), 'lat lon')
      iret = nf_put_att_text(ncid, varids(2), 'grid_mapping', len_trim('grid_mapping'), &
                             'grid_mapping')
  
      if (varname == 'tp') then
        iret = nf_put_att_text(ncid, varids(3), 'long_name', &
                               len_trim('ensemble shifted mean of precipitation'), &
                               'ensemble shifted mean of precipitation')
        iret = nf_put_att_text(ncid, varids(3), 'units', len_trim('kg m-2'), 'kg m-2')
      else if (varname == 'asnow') then
        iret = nf_put_att_text(ncid, varids(3), 'long_name', &
                               len_trim('ensemble shifted mean of snowfall'), &
                               'ensemble shifted mean of snowfall')
        iret = nf_put_att_text(ncid, varids(3), 'units', len_trim('m'), 'm')
      end if
      iret = nf_put_att_real(ncid, varids(3), '_FillValue', nf_float, 1, -999.)
      iret = nf_put_att_int(ncid, varids(3), 'lengthOfTimeRange', nf_int, 1, acctime)
      iret = nf_put_att_text(ncid, varids(3), 'coordinates', len_trim('lat lon'), 'lat lon')
      iret = nf_put_att_text(ncid, varids(3), 'grid_mapping', len_trim('grid_mapping'), &
                             'grid_mapping')

      if (varname == 'tp') then
        iret = nf_put_att_text(ncid, varids(4), 'long_name', &
                               len_trim('ensemble shifted LPM mean of precipitation'), &
                               'ensemble shifted LPM mean of precipitation')
        iret = nf_put_att_text(ncid, varids(4), 'units', len_trim('kg m-2'), 'kg m-2')
      else if (varname == 'asnow') then
        iret = nf_put_att_text(ncid, varids(4), 'long_name', &
                               len_trim('ensemble shifted LPM mean of snowfall'), &
                               'ensemble shifted LPM mean of snowfall')
        iret = nf_put_att_text(ncid, varids(4), 'units', len_trim('m'), 'm')
      end if
      iret = nf_put_att_real(ncid, varids(4), '_FillValue', nf_float, 1, -999.)
      iret = nf_put_att_int(ncid, varids(4), 'lengthOfTimeRange', nf_int, 1, acctime)
      iret = nf_put_att_text(ncid, varids(4), 'coordinates', len_trim('lat lon'), 'lat lon')
      iret = nf_put_att_text(ncid, varids(4), 'grid_mapping', len_trim('grid_mapping'), &
                             'grid_mapping')

      iret = nf_put_att_text(ncid, varids(5), 'long_name', len_trim('time'), &
                             'time')
      iret = nf_put_att_text(ncid, varids(5), 'units', &
             len_trim('hours since 1970-01-01 00:00:00 0:00'), 'hours since 1970-01-01 00:00:00 0:00')
      iret = nf_put_att_text(ncid, varids(5), 'standard_name', len_trim('time'), &
                             'time')
      iret = nf_put_att_text(ncid, varids(5), 'axis', len_trim('T'), &
                             'T')

      iret = nf_put_att_text(ncid, varids(6), 'long_name', len_trim('forecast_reference_time'), &
                             'forecast_reference_time')
      iret = nf_put_att_text(ncid, varids(6), 'units', &
             len_trim('hours since 1970-01-01 00:00:00 0:00'), 'hours since 1970-01-01 00:00:00 0:00')
      iret = nf_put_att_text(ncid, varids(6), 'standard_name', len_trim('forecast_reference_time'), &
                             'forecast_reference_time')

      iret = nf_put_att_text(ncid, varids(7), 'grid_mapping_name', len_trim('lambert_conformal_conic'), &
                             'lambert_conformal_conic')
      iret = nf_put_att_real(ncid, varids(7), 'standard_parallel', nf_double, 1, 38.5)
      iret = nf_put_att_real(ncid, varids(7), 'longitude_of_central_meridian', nf_double, 1, 262.5)
      iret = nf_put_att_real(ncid, varids(7), 'latitude_of_projection_origin', nf_double, 1, 38.5)
      iret = nf_put_att_int(ncid, varids(7), 'false_easting', nf_int, 1, 0)
      iret = nf_put_att_int(ncid, varids(7), 'false_northing', nf_int, 1, 0)
      iret = nf_put_att_text(ncid, varids(7), 'GRIB_earth_shape', len_trim('spherical'), &
                             'spherical')
      iret = nf_put_att_int(ncid, varids(7), 'GRIB_earth_shape_code', nf_int, 1, 6)

      iret = nf_put_att_text(ncid, varids(8), 'units', len_trim('km'), &
                             'km')
      iret = nf_put_att_text(ncid, varids(8), 'long_name', len_trim('x coordinate of projection'), &
                             'x coordinate of projection')
      iret = nf_put_att_text(ncid, varids(8), 'standard_name', len_trim('projection_x_coordinate'), &
                             'projection_x_coordinate')
      iret = nf_put_att_text(ncid, varids(8), 'grid_spacing', len_trim('3 km'), &
                             '3 km')
      iret = nf_put_att_text(ncid, varids(8), 'axis', len_trim('X'), &
                             'X')

      iret = nf_put_att_text(ncid, varids(9), 'units', len_trim('km'), &
                             'km')
      iret = nf_put_att_text(ncid, varids(9), 'long_name', len_trim('y coordinate of projection'), &
                             'y coordinate of projection')
      iret = nf_put_att_text(ncid, varids(9), 'standard_name', len_trim('projection_y_coordinate'), &
                             'projection_y_coordinate')
      iret = nf_put_att_text(ncid, varids(9), 'grid_spacing', len_trim('3 km'), &
                             '3 km')
      iret = nf_put_att_text(ncid, varids(9), 'axis', len_trim('Y'), &
                             'Y')

      iret = nf_put_att_text(ncid, varids(10), 'long_name', len_trim('ensemble name'), &
                             'ensemble name')

      iret = nf_put_att_text(ncid, NCGLOBAL, 'Conventions', len_trim('CF-1.6'), 'CF-1.6')

      if (alignopt == 1) then
        iret = nf_put_att_text(ncid, NCGLOBAL, 'Alignment_Method', len_trim('Method A'), 'Method A')
      elseif (alignopt == 2) then
        iret = nf_put_att_text(ncid, NCGLOBAL, 'Alignment_Method', len_trim('Method B'), 'Method B')
      endif

      iret = nf_put_att_int(ncid, NCGLOBAL, 'Total_iteration_pass_of_alignment', nf_int, 1, nshfpass)

      iret = nf_enddef(ncid)

! Put variables
      iret = nf_put_var_real(ncid, varids(1), lats(:,:))
      iret = nf_put_var_real(ncid, varids(2), lons(:,:))
      iret = nf_put_var_real(ncid, varids(3), ensshfmn(:,:,nshfpass))
      iret = nf_put_var_real(ncid, varids(4), ensshflpm(:,:,nshfpass))
      iret = nf_put_var_int(ncid, varids(5), itime)
      iret = nf_put_var_int(ncid, varids(6), rtime)
      iret = nf_put_var_int(ncid, varids(7), 0)
      iret = nf_put_var_real(ncid, varids(8), projx)
      iret = nf_put_var_real(ncid, varids(9), projy)
      iret = nf_put_var_text(ncid, varids(10), ensfname(:))  

! Close Netcdf File
      iret = nf_close(ncid)
      deallocate(dimids, varids, vardims, chunks, vardims2)
    END IF

    deallocate(projx, projy)
  END IF

  CALL MPI_FINALIZE(ierr)

  IF( myproc == root ) THEN
    CALL CPU_TIME(cput7)
    WRITE(6,'(a,f10.2,a)') ' Total Execution CPU time:',(cput7-cput0),' seconds'
  END IF

  STOP 
END PROGRAM ENSALIGN

SUBROUTINE a2dmax0(array,idim0,idim1,ibgn,iend,jdim0,jdim1,jbgn,jend,vmin,vmax)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: idim0,idim1,ibgn,iend
  INTEGER, INTENT(IN) :: jdim0,jdim1,jbgn,jend
  REAL, INTENT(IN) :: array(idim0:idim1,jdim0:jdim1)
  REAL, INTENT(OUT) :: vmin
  REAL, INTENT(OUT) :: vmax

  INTEGER :: i,j

  vmin=1.0E30
  vmax=-1.0E30
  DO j=jbgn,jend
    DO i=ibgn,iend
      vmin=MIN(vmin,array(i,j))
      vmax=MAX(vmax,array(i,j))
    END DO
  END DO
  RETURN
END SUBROUTINE a2dmax0

SUBROUTINE LGARRAY_BCASTR( array, nelem, max_elem_send, &
                             myproc, srcproc, comm_channel, ierr)
!
! Broadcasts data to other processors from a large array using MPI_BCAST.
! Solves size limit problem of MPI_BCAST by sending the data in
! smaller chunks.
!
! array: REAL array to be broadcasted
! nelem: Number of total elements in array (can be nx*ny in calling program)
! max_elem_send: Limit to MPI_BCAST (Intel is about 500,000)
! myproc: Processor number of calling processor
! srcproc: Source processor number for broadcast
! comm_channel: MPI Communications channel, e.g., MPI_COMM_WORLD
! ierr: Output status
!
! Keith Brewster, CAPS/Univ of Oklahoma
! March, 2018
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nelem
  REAL, INTENT(INOUT) :: array(nelem)
  INTEGER, INTENT(IN) :: max_elem_send
  INTEGER, INTENT(IN) :: myproc
  INTEGER, INTENT(IN) :: srcproc
  INTEGER, INTENT(IN) :: comm_channel
  INTEGER, INTENT(OUT) :: ierr
!
! Misc local variables
!
  INTEGER :: k,ncall,nelemcall,idxbgn,idxend,nelembcast
  INCLUDE 'mpif.h'
  
  ierr = 0
  ncall=nelem/max_elem_send
  IF( (ncall*max_elem_send) < nelem ) ncall=ncall+1
! IF( myproc == srcproc ) WRITE(6,'(a,i8,a,i8,a,i5)') &
!   ' LGARRAY_BCASTR: nelm=',nelem,' max elem:',max_elem_send,' ncall:',ncall 
  nelemcall=nelem/ncall
  IF( (nelemcall*ncall) < nelem ) nelemcall=nelemcall+1
! IF( myproc == srcproc ) THEN
!   WRITE(6,'(a,i8)') ' LGARRAY_BCASTR: nelem=',nelem
!   WRITE(6,'(a,i5,a,i8)') &
!   ' LGARRAY_BCASTR: ncall=',ncall,' nelem/call:',nelemcall 
! END IF
! WRITE(6,'(a,i5,a)') 'Processor ',myproc,' here 1.'
! FLUSH(6)
  DO k=1,ncall
    IF(k == 1) THEN
      idxbgn=1
    ELSE
      idxbgn=idxend+1
    END IF
    idxend=min((idxbgn+nelemcall-1),nelem)
    nelembcast=(idxend-idxbgn)+1
!   IF( myproc == srcproc ) THEN
!     WRITE(6,'(a,i5,a,i8,a,i8,a,i8)') &
!     '   k=',k,' idxbgn=',idxbgn,' idxend=',idxend,' nelembcast:',nelembcast
!     FLUSH(6)
!   END IF
    CALL MPI_BCAST(array(idxbgn),nelembcast,MPI_REAL,srcproc, &
                       comm_channel,ierr)
    !CALL MPI_BARRIER(comm_channel,ierr)
  END DO
  RETURN
END SUBROUTINE LGARRAY_BCASTR

SUBROUTINE LGARRAY_INT_BCASTR( iarray, nelem, max_elem_send, &
                             myproc, srcproc, comm_channel, ierr)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nelem
  INTEGER, INTENT(INOUT) :: iarray(nelem)
  INTEGER, INTENT(IN) :: max_elem_send
  INTEGER, INTENT(IN) :: myproc
  INTEGER, INTENT(IN) :: srcproc
  INTEGER, INTENT(IN) :: comm_channel
  INTEGER, INTENT(OUT) :: ierr
!
! Misc local variables
!
  INTEGER :: k,ncall,nelemcall,idxbgn,idxend,nelembcast
  INCLUDE 'mpif.h'

  ierr = 0
  ncall=nelem/max_elem_send
  IF( (ncall*max_elem_send) < nelem ) ncall=ncall+1
! IF( myproc == srcproc ) WRITE(6,'(a,i8,a,i8,a,i5)') &
!   ' LGARRAY_BCASTR: nelm=',nelem,' max elem:',max_elem_send,' ncall:',ncall
  nelemcall=nelem/ncall
  IF( (nelemcall*ncall) < nelem ) nelemcall=nelemcall+1
! IF( myproc == srcproc ) THEN
!   WRITE(6,'(a,i8)') ' LGARRAY_BCASTR: nelem=',nelem
!   WRITE(6,'(a,i5,a,i8)') &
!   ' LGARRAY_BCASTR: ncall=',ncall,' nelem/call:',nelemcall
! END IF
! WRITE(6,'(a,i5,a)') 'Processor ',myproc,' here 1.'
! FLUSH(6)
  DO k=1,ncall
    IF(k == 1) THEN
      idxbgn=1
    ELSE
      idxbgn=idxend+1
    END IF
    idxend=min((idxbgn+nelemcall-1),nelem)
    nelembcast=(idxend-idxbgn)+1
!   IF( myproc == srcproc ) THEN
!     WRITE(6,'(a,i5,a,i8,a,i8,a,i8)') &
!     '   k=',k,' idxbgn=',idxbgn,' idxend=',idxend,' nelembcast:',nelembcast
!     FLUSH(6)
!   END IF
    CALL MPI_BCAST(iarray(idxbgn),nelembcast,MPI_INTEGER,srcproc, &
                       comm_channel,ierr)
    !CALL MPI_BARRIER(comm_channel,ierr)
  END DO
  RETURN
END SUBROUTINE LGARRAY_INT_BCASTR

SUBROUTINE LGARRAY_REDUCER( array,recv_buf,nelem,max_elem_send, &
                            myproc,root,reduce_op,comm_channel,ierr)
!
! Reduces array of data to recv_buf for a large array using MPI_REDUCE.
! Solves size limit problem of MPI_REDUCE by sending the data in
! smaller chunks.
!
! array: REAL array to be reduced
! recv_buf: REAL array to hold reduced result in root processor
! nelem: Number of total elements in array (can be nx*ny in calling program)
! max_elem_send: Limit to MPI_BCAST (Intel is about 500,000)
! myproc: Processor number of calling processor
! root: Root processor number, where reduce result is collected
! reduce_op: Reduce operation (such as MPI_SUM)
! comm_channel: MPI Communications channel, e.g., MPI_COMM_WORLD
! ierr: Output status
!
! Keith Brewster, CAPS/Univ of Oklahoma
! March, 2018
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nelem
  REAL, INTENT(INOUT) :: array(nelem)
  REAL, INTENT(OUT) :: recv_buf(nelem)
  INTEGER, INTENT(IN) :: max_elem_send
  INTEGER, INTENT(IN) :: myproc
  INTEGER, INTENT(IN) :: root
  INTEGER, INTENT(IN) :: reduce_op
  INTEGER, INTENT(IN) :: comm_channel
  INTEGER, INTENT(OUT) :: ierr
!
! Misc local variables
!
  INTEGER :: k,ncall,nelemcall,idxbgn,idxend,nelemreduce
  INCLUDE 'mpif.h'
  
  ierr = 0
  ncall=nelem/max_elem_send
  IF( (ncall*max_elem_send) < nelem ) ncall=ncall+1
! IF( myproc == srcproc ) WRITE(6,'(a,i8,a,i8,a,i5)') &
!   ' LGARRAY_BCASTR: nelm=',nelem,' max elem:',max_elem_send,' ncall:',ncall 
  nelemcall=nelem/ncall
  IF( (nelemcall*ncall) < nelem ) nelemcall=nelemcall+1
! IF( myproc == srcproc ) THEN
!   WRITE(6,'(a,i8)') ' LGARRAY_BCASTR: nelem=',nelem
!   WRITE(6,'(a,i5,a,i8)') &
!   ' LGARRAY_BCASTR: ncall=',ncall,' nelem/call:',nelemcall 
! END IF
! WRITE(6,'(a,i5,a)') 'Processor ',myproc,' here 1.'
! FLUSH(6)
  DO k=1,ncall
    IF(k == 1) THEN
      idxbgn=1
    ELSE
      idxbgn=idxend+1
    END IF
    idxend=min((idxbgn+nelemcall-1),nelem)
    nelemreduce=(idxend-idxbgn)+1
!   IF( myproc == srcproc ) THEN
!     WRITE(6,'(a,i5,a,i8,a,i8,a,i8)') &
!     '   k=',k,' idxbgn=',idxbgn,' idxend=',idxend,' nelemreduce:',nelemreduce
!     FLUSH(6)
!   END IF
    CALL MPI_REDUCE(array(idxbgn),recv_buf(idxbgn),nelemreduce, &
                    MPI_REAL,reduce_op,root,comm_channel,ierr)
    !CALL MPI_BARRIER(comm_channel,ierr)
  END DO
  RETURN
END SUBROUTINE LGARRAY_REDUCER

SUBROUTINE vmaxmin(varray,nx,ny,ibgn,iend,jbgn,jend,vmin,vmax)
  IMPLICIT NONE
  INTEGER :: nx,ny
  REAL :: varray(nx,ny)
  INTEGER :: ibgn,iend
  INTEGER :: jbgn,jend
  REAL :: vmin,vmax

  INTEGER :: i,j

  vmin=varray(ibgn,jbgn)
  vmax=varray(ibgn,jbgn)
  DO j=jbgn,jend
    DO i=ibgn,iend
      vmin=min(varray(i,j),vmin) 
      vmax=max(varray(i,j),vmax) 
!     IF(varray(i,j) > 1000.) THEN
!       WRITE(6,'(a,f12.2,a,2i5)') ' Huge value of: ',varray(i,j),' at i,j: ',i,j
!     END IF
    END DO
  END DO
  RETURN
END SUBROUTINE vmaxmin
