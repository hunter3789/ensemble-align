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
!    - Can Read GEFS and HREF data
!    - Write shifted means (also include PM and LPM) and shift vectors
!      into NETCDF file
!

  USE eccodes
  USE netcdf

  IMPLICIT NONE
  include 'netcdf.inc'
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
  REAL, ALLOCATABLE :: ensfcst_shf(:,:,:,:)
  REAL, ALLOCATABLE :: ensmean(:,:)
  REAL, ALLOCATABLE :: ensmax(:,:)
  REAL, ALLOCATABLE :: ensshfmn(:,:,:)
  REAL, ALLOCATABLE :: ensshfmax(:,:,:)
  REAL, ALLOCATABLE :: enslpm(:,:), enspm(:,:)
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
! Dimenstions 
!
!-----------------------------------------------------------------------
!
  INTEGER :: nvar = 1
  INTEGER :: nx,ny,nxny,nelements,nx_grib,ny_grib,cx,cy
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

  INTEGER :: i,j,k,k1,k2,ii,jj
  INTEGER :: ictr,jctr,realdata,missing
  INTEGER :: lvldbg,ipass,mxzone,applyshft
  INTEGER :: itime,ifhr,ntstep,itstep,rtime
  INTEGER :: narg,membknt,nelem2d,nelem3d
  INTEGER :: lstring,rsize,bufelems,mpibufsize
  INTEGER :: istat,istatus,inqstat,ierr,ireturn
  INTEGER :: year,month,day,hour,minute
  INTEGER :: ilap,jlap,istep,jstep
  INTEGER :: odate(4), pdate(4), odate2(4)
  REAL :: dx,dy
  REAL :: ctrx0,ctry0,ctrx,ctry
  REAL :: delx,dely
  REAL :: refmax,radx,rady
  REAL :: piov2,radx2inv,rady2inv,radius
  REAL :: rninv,fnorm,fcstmin,fcstmax,favgmin,favgmax
  REAL :: time,vmin,vmax
  REAL :: cput0,cput1,cput2,cput3,cput4,cput5,cput6,cput7
  REAL, PARAMETER :: rmisg_data = 9.0E36

  integer             :: ifile, idx, igrib
  integer             :: iret, ncid, ndim, count
  integer,allocatable :: dimids(:),vardims(:),chunks(:),varids(:)
  integer,allocatable :: vardims2(:),chunks2(:),vardims3(:)
  REAL, ALLOCATABLE   :: values(:),projx(:),projy(:)

  CHARACTER(LEN=2)  :: acchrstr
  CHARACTER(LEN=4)  :: fhrstr
  CHARACTER(LEN=8)  :: datestr
  CHARACTER(LEN=6)  :: varid0
  CHARACTER(LEN=6)  :: varid
  CHARACTER(LEN=8)  :: varunits
  CHARACTER(LEN=8)  :: is_dir
  CHARACTER(LEN=20) :: validtime
  CHARACTER(LEN=20) :: varout
  CHARACTER(LEN=20) :: parm
  CHARACTER(LEN=40) :: fullvarname
  CHARACTER(LEN=20) :: runname
  CHARACTER(LEN=512):: fcstname
  CHARACTER(LEN=512):: outname
  CHARACTER(LEN=512):: fname
  CHARACTER(LEN=512):: grb2fname
  CHARACTER(LEN=256):: charg
  CHARACTER(LEN=256):: fname_input
  CHARACTER(LEN=256):: levsfc
  CHARACTER(LEN=256):: metadata
  CHARACTER(LEN=256):: command_string
  CHARACTER(LEN=20), allocatable :: fname_head(:)

  LOGICAL :: file_exists,dir_exists,lposdef

  INTEGER, EXTERNAL :: grb2_wrt
  
!
!-----------------------------------------------------------------------
!
! These flags indicate whether the quantities are calculated
!
!-------------------------------------------------------------------------
!
  INTEGER, PARAMETER :: isam = 1
  INTEGER, PARAMETER :: isamlpm = 1
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mpif.h'
  INCLUDE 'align.inc'
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

  fname_input = ' '
  IF(myproc == root) THEN
    narg=COMMAND_ARGUMENT_COUNT()
    IF( narg > 1 ) THEN
      CALL GET_COMMAND_ARGUMENT(1, charg, istat, istatus )
      IF( charg(1:2) == '-i' .OR. charg(1:6) == '-input') THEN
        CALL GET_COMMAND_ARGUMENT(2, charg, istat, istatus )
        fname_input = charg
      END IF
      WRITE(stdout,'(1x,2a)') 'Namelist file from command line = ',TRIM(fname_input)
    END IF
  END IF
  
  piov2 = 2.0*atan(1.0)
  time = 0.
  varid0 = 'tp'
  varout = 'P03M'
  varunits = 'mm'
  runname = 'enstest'
  realdata = 3
  applyshft = 0

  call getarg(1,datehrstr)
  call getarg(2,fhrstr)
  read(fhrstr,*) ifhr

  READ(datehrstr,'(i4,i2,i2,i2)') odate(1), odate(2), odate(3), odate(4)
  CALL GETH_NEWDATE(pdate, odate, ifhr*3600)
  call calc_hour(pdate, itime)
  call calc_hour(odate, rtime)
!
!-------------------------------------------------------------------------
!
! Get NAMELIST variables
!
!-----------------------------------------------------------------------
!
  !CALL initalign(nx,ny,dx,dy,lvldbg,fname_input,istatus)
  !WRITE(6,'(i0,i0,f0.0,f0.0)') nx, ny, dx, dy  

  IF ( realdata > 0 ) THEN
    IF ( realdata == 1) THEN
      WRITE(fname_input,'(a)') '/home/clee/METplus/Tutorial/data/ecen_t063_kod1_h024.202209020000.gb1'

      call codes_open_file(ifile,fname_input,'r')
      call codes_index_create(idx,fname_input,'shortName,level,number')    

      call codes_index_select(idx,'shortName',varid0)
      call codes_index_select(idx,'level',0)
      call codes_index_select(idx,'number',0)

      call codes_new_from_index(idx,igrib,iret)
      call codes_get(igrib,'Ni',nx)
      call codes_get(igrib,'Nj',ny)
      call codes_get(igrib,'numberOfForecastsInEnsemble',nmembers)
      call codes_release(igrib)
    ELSE IF (realdata == 2) THEN
      nx_grib = 0
      ny_grib = 0
      nx = 154
      ny = 82
      nmembers = 30
    ELSE IF (realdata == 3) THEN
      !WRITE(fname_input,'(a,i4.4,3(i2.2),a,i2.2,a,i2.2,a,i3.3)') '/non/clee/gefs/', & 
      !        odate(1), odate(2), odate(3), odate(4), '/m', k, '/gep', k, '.t00z.pgrb2a.0p50.f', ifhr

      !call codes_open_file(ifile,fname_input,'r')
      !call codes_index_create(idx,fname_input,'shortName,level,lengthOfTimeRange')

      !call codes_index_select(idx,'shortName',varid0)
      !call codes_index_select(idx,'level',0)
      !call codes_index_select(idx,'lengthOfTimeRange',1)

      !call codes_new_from_index(idx,igrib,iret)
      !call codes_get(igrib,'Ni',nx)
      !call codes_get(igrib,'Nj',ny)
      !call codes_get(igrib,'numberOfForecastsInEnsemble',nmembers)
      !call codes_release(igrib)

      nx_grib = 0
      ny_grib = 0
      nx = 1799
      ny = 1059
      nmembers = 10
    END IF
  ELSE
    nx=43
    ny=43
    nmembers = 5
  END IF

  dx = 3000.!50000.0
  dy = 3000.!50000.0
  
  alignopt = 1
  calcmean = 1
  iposdef = 1
  ifmtout = 7
  minkdratio = 0.8
  slnratio0h = 0.3
  slnratio48h = 0.5

  DO i=1,mx_shfpass
    nizone(i)=10
    njzone(i)=10
  END DO

  DO i=1,nvarshf
    wgtvar(i)=1.0
  END DO

  lpmopt = 1
  patch_nx = 10
  patch_ny = 10
  ovx = 30
  ovy = 30
  gauss_sigma = 1
  filt_min = 0.

  shfoutdir='./shfmn'
  lpmoutdir='./lpm'
  shflpmoutdir='./shflpm'
  wrtorgmemb=0
  wrtshfmemb=0
  wrtshfmean=1


  nshfpass = 2
  nbaksmth = 0
  nshfsmth = 0
  noutsmth = 1
  minkdat = 3
  !slnratio = 0.5
  hrzlap = 0.5

  ibgn = 1
  jbgn = 1
  iend = nx
  jend = ny

  lvldbg = 0

  lposdef=(iposdef > 0)
  nxny=nx*ny

!
!#######################################################################
!
! Shift Zones
!
! nizone: number of zones in x direction for each shift iteration
! njzone: number of zones in y direction for each shift iteration
!
!#######################################################################


  !nizone(1)= 16;nizone(2)= 32;nizone(3)=64;nizone(4)=25;nizone(5)=10;
  !njzone(1)= 8; njzone(2)= 16;njzone(3)=32;njzone(4)=25;njzone(5)=10;

  iborder(1)=0;iborder(2)=0;iborder(3)=0;iborder(4)=5;
  izsize(1)=200;izsize(2)=75;izsize(3)=50;izsize(4)=10;
  jborder(1)=0;jborder(2)=0;jborder(3)=0;jborder(4)=5;
  jzsize(1)=200;jzsize(2)=75;jzsize(3)=50;jzsize(4)=10;

  mxzone=0
  DO ipass=1,nshfpass
    ilap=MAX(IFIX((izsize(ipass)*hrzlap)+0.5),1)
    jlap=MAX(IFIX((jzsize(ipass)*hrzlap)+0.5),1)
    istep=izsize(ipass)-ilap
    jstep=jzsize(ipass)-jlap
    nizone(ipass)=(nx-(2*iborder(ipass)))/istep
    njzone(ipass)=(ny-(2*jborder(ipass)))/jstep
    mxzone=MAX(mxzone,(nizone(ipass)*njzone(ipass)))
    IF (myproc == root) THEN
      WRITE(6,'(4(a,i0))') 'nizone(', ipass, ') : ', nizone(ipass), &
          ', njzone(', ipass, '): ', njzone(ipass)
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
  ALLOCATE(posdef(nvar))
  ALLOCATE(xs(nx))
  ALLOCATE(ys(ny))
  ALLOCATE(ensfname(mx_members))
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
  ALLOCATE(lats(nx,ny),lons(nx,ny))
  ALLOCATE(ensshfmn(nx,ny,nshfpass),ensshfmax(nx,ny,nshfpass))
  ALLOCATE(enslpm(nx,ny),enspm(nx,ny))
  ALLOCATE(ensshflpm(nx,ny,nshfpass),ensshfpm(nx,ny,nshfpass))
  ALLOCATE(gs_weight(nx*ny))

  ALLOCATE(ensfcst(nx,ny,nmembers))
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
!
!-----------------------------------------------------------------------
!
!  Initialize Arrays
!
!-----------------------------------------------------------------------
!
  IF(myproc == root) WRITE(6,'(a)') 'Initialize Arrays' 
  CALL CPU_TIME(cput0)
  !datestr=datehrstr(1:8)
  !READ(datehrstr,'(i4,i2,i2,i2)') year,month,day,hour
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
  runname='ar'//datehrstr
  levsfc='surface'
  !alignopt = 1
  !calcmean = 1 

! bufelems=nx*ny
! CALL MPI_PACK_SIZE(bufelems,MPI_REAL,MPI_COMM_WORLD,rsize,ierr)
! print *, ' MPI PACK REAL size: ',rsize
! mpibufsize=nprocs*(rsize+MPI_BSEND_OVERHEAD)
! print *, ' mpibufsize: ',mpibufsize
! ALLOCATE(mpi_work_buf(mpibufsize),STAT=ierr)
! print *, ' mpi_work_buf status: ',ierr
  CALL gaussian_weight(nx,ny,gauss_sigma,gs_weight)
!
!-----------------------------------------------------------------------
!
!  Get Ensemble Forecast Data
!
!-----------------------------------------------------------------------
!
  IF(myproc == root) WRITE(6,'(a)') 'Get Ensemble Forecast Data' 
  DO i=1,nx
    xs(i)=(float(i-2)+0.5)*dx
  END DO
  DO j=1,ny
    ys(j)=(float(j-2)+0.5)*dy
  END DO

  CALL setdxdy(nx,ny,1,nx-1,1,ny-1,                                    &
           xs,ys,dxfld,dyfld,rdxfld,rdyfld)
  IF(myproc == root) WRITE(6,'(a)') 'Set dx dy' 

  IF ( realdata > 0 ) THEN
    !ntstep=1+((tend-tbegin)/tintvl)
  ELSE
    membknt=min(nmembers,5)
    !tbegin=0
    !tend=0
    !tintvl=10800
    !ntstep=1
  END IF

  tbegin=0
  tend=0
  tintvl=10800
  ntstep=1

  IF( myproc == root ) THEN
    WRITE(6,'(a,i5,a,i5)') ' MPI nprocs:',nprocs,'  myproc:',myproc
    CALL CPU_TIME(cput1)
    WRITE(6,'(a,f10.2,a)') ' Start-up CPU time: ',(cput1-cput0),' seconds'
  END IF

  DO itstep=1,ntstep

    IF( myproc == root) CALL CPU_TIME(cput1)
    !itime=tbegin+(itstep-1)*tintvl
    !time=float(itime)
    !ifhr = itime/3600
    !IF( myproc == root) WRITE(6,'(//a,i9,a,i3)') &
    !  'Processing forecast time: ',itime,'  Hour:',ifhr
    ensfcst(:,:,:) = 0.
    !WRITE(hh,'(i2.2)') ifhr
    !WRITE(hhm3,'(i2.2)') (ifhr-3)
    membknt=nmembers

    IF( myproc == root) THEN
      IF ( realdata > 0 ) THEN
        IF (realdata == 1) THEN
          DO k=1,nmembers
            call codes_index_select(idx,'shortName',varid0)
            call codes_index_select(idx,'level',0)
            call codes_index_select(idx,'number',k-1)
            call codes_new_from_index(idx,igrib,iret)

            allocate(values(Nx*Ny))
            call codes_get(igrib,'values',values)
 
            do i=1,Nx
             do j=1,Ny
              ensfcst(i,j,k) = values((j-1)*Nx+i)
             enddo
            enddo

            deallocate(values)
            call codes_release(igrib)
          END DO

          call codes_close_file(ifile)
        ELSE IF (realdata == 2) THEN
          DO k=1,nmembers
            WRITE(fname_input,'(a,i4.4,3(i2.2),a,i2.2,a,i2.2,a,i3.3)') '/non/clee/gefs/', &
              odate(1), odate(2), odate(3), odate(4), '/m', k, '/gep', k, '.t00z.pgrb2a.0p50.f', ifhr 
            !print *, fname_input

            call codes_open_file(ifile,fname_input,'r')
            call codes_index_create(idx,fname_input,'shortName,level')

            call codes_index_select(idx,'shortName',varid0)
            call codes_index_select(idx,'level',0)
            call codes_new_from_index(idx,igrib,iret)

            IF (nx_grib .eq. 0 .OR. ny_grib .eq. 0) THEN
              call codes_get(igrib,'Ni',nx_grib)
              call codes_get(igrib,'Nj',ny_grib)

              do i=1,Nx
                do j=1,Ny
                  lons(i,j) = float((452+i))*0.5
                  lats(i,j) = 90.0 - float((j-1+60))*0.5
                enddo
              enddo
            END IF
            allocate(values(nx_grib*ny_grib))
            call codes_get(igrib,'values',values)

            do i=1,Nx
              do j=1,Ny
                ensfcst(i,j,k) = values((j-1+60)*nx_grib+i+452)
              enddo
            enddo

            deallocate(values)
            call codes_release(igrib)
            call codes_close_file(ifile)
          END DO
        ELSE IF (realdata == 3) THEN
          allocate(fname_head(nmembers))
          write(fname_head(1),'(a)') 'hrrr_ncep_'
          write(fname_head(2),'(a)') 'hiresw_conusnssl_'
          write(fname_head(3),'(a)') 'hiresw_conusarw_'
          write(fname_head(4),'(a)') 'hiresw_conusfv3_'
          write(fname_head(5),'(a)') 'nam_conusnest_'
          write(fname_head(6),'(a)') 'hrrr_ncep_'
          write(fname_head(7),'(a)') 'hiresw_conusnssl_'
          write(fname_head(8),'(a)') 'hiresw_conusarw_'
          write(fname_head(9),'(a)') 'hiresw_conusfv3_'
          write(fname_head(10),'(a)') 'nam_conusnest_'
       
          DO k=1,nmembers
            WRITE(fname_input,'(a,i4.4,2(i2.2),a,a,i4.4,3(i2.2),a,i3.3,a)') '/archive2/2022_HMT_summer/href/', &
              odate(1), odate(2), odate(3), '/', trim(fname_head(k)), odate(1), odate(2), odate(3), odate(4), &
              'f', ifhr, '.grib2'
            print *, trim(fname_input)

            call codes_open_file(ifile,fname_input,'r')
            IF (k.eq.5 .or. k.eq.10) THEN
              call codes_index_create(idx,fname_input,'shortName,level')
            ELSE
              call codes_index_create(idx,fname_input,'shortName,level,lengthOfTimeRange')
            END IF

            call codes_index_select(idx,'shortName',varid0)
            call codes_index_select(idx,'level',0)
            IF (k.eq.5 .or. k.eq.10) THEN
            ELSE
              call codes_index_select(idx,'lengthOfTimeRange',1)
            END IF
            call codes_new_from_index(idx,igrib,iret)

            IF (nx_grib .eq. 0 .OR. ny_grib .eq. 0) THEN
              call codes_get(igrib,'Nx',nx_grib)
              call codes_get(igrib,'Ny',ny_grib)
              call codes_get(igrib,'missingValue',missing)
              !print *, nx_grib, ny_grib

              allocate(values(nx_grib*ny_grib))
              call codes_get(igrib,'latitudes',values)

              do i=1,Nx
               do j=1,Ny
                lats(i,j) = values((j-1)*Nx+i)
               enddo
              enddo

              call codes_get(igrib,'longitudes',values)

              do i=1,Nx
               do j=1,Ny
                lons(i,j) = values((j-1)*Nx+i)
               enddo
              enddo

              deallocate(values)
            END IF

            allocate(values(nx_grib*ny_grib))
            call codes_get(igrib,'values',values)

            do i=1,Nx
             do j=1,Ny
              ensfcst(i,j,k) = values((j-1)*Nx+i)
             enddo
            enddo

            deallocate(values)
            call codes_release(igrib)
            call codes_close_file(ifile)

            IF (k.eq.5 .or. k.eq.10) THEN
              !IF (MOD(ifhr,3) .ne. 1) THEN

              !  WRITE(fname_input,'(a,i4.4,2(i2.2),a,a,i4.4,3(i2.2),a,i3.3,a)') '/archive2/2022_HMT_summer/href/', &
              !    odate(1), odate(2), odate(3), '/', trim(fname_head(k)), odate(1), odate(2), odate(3), odate(4), &
              !    'f', ifhr-1, '.grib2'
              !  print *, trim(fname_input)

              !  call codes_open_file(ifile,fname_input,'r')
              !  call codes_index_create(idx,fname_input,'shortName,level')

              !  call codes_index_select(idx,'shortName',varid0)
              !  call codes_index_select(idx,'level',0)
              !  call codes_new_from_index(idx,igrib,iret)

              !  allocate(values(nx_grib*ny_grib))
              !  call codes_get(igrib,'values',values)

                do i=1,Nx
                 do j=1,Ny
                  if (ensfcst(i,j,k).eq.missing) then
                    ensfcst(i,j,k) = 0.
                  end if
                  !ensfcst(i,j,k) = ensfcst(i,j,k) - values((j-1)*Nx+i)
                 enddo
                enddo

              !  deallocate(values)
              !  call codes_release(igrib)
              !  call codes_close_file(ifile)

              !END IF
            ELSE
                WRITE(fname_input,'(a,i4.4,2(i2.2),a,a,i4.4,3(i2.2),a,i3.3,a)') '/archive2/2022_HMT_summer/href/', &
                  odate(1), odate(2), odate(3), '/', trim(fname_head(k)), odate(1), odate(2), odate(3), odate(4), &
                  'f', ifhr-1, '.grib2'
                print *, trim(fname_input)

                call codes_open_file(ifile,fname_input,'r')
                call codes_index_create(idx,fname_input,'shortName,level,lengthOfTimeRange')

                call codes_index_select(idx,'shortName',varid0)
                call codes_index_select(idx,'level',0)
                call codes_index_select(idx,'lengthOfTimeRange',1)
                call codes_new_from_index(idx,igrib,iret)

                allocate(values(nx_grib*ny_grib))
                call codes_get(igrib,'values',values)

                do i=1,Nx
                 do j=1,Ny
                  ensfcst(i,j,k) = ensfcst(i,j,k) + values((j-1)*Nx+i)
                 enddo
                enddo

                deallocate(values)
                call codes_release(igrib)
                call codes_close_file(ifile)


                WRITE(fname_input,'(a,i4.4,2(i2.2),a,a,i4.4,3(i2.2),a,i3.3,a)') '/archive2/2022_HMT_summer/href/', &
                  odate(1), odate(2), odate(3), '/', trim(fname_head(k)), odate(1), odate(2), odate(3), odate(4), &
                  'f', ifhr-2, '.grib2'
                print *, trim(fname_input)

                call codes_open_file(ifile,fname_input,'r')
                call codes_index_create(idx,fname_input,'shortName,level,lengthOfTimeRange')

                call codes_index_select(idx,'shortName',varid0)
                call codes_index_select(idx,'level',0)
                call codes_index_select(idx,'lengthOfTimeRange',1)
                call codes_new_from_index(idx,igrib,iret)

                allocate(values(nx_grib*ny_grib))
                call codes_get(igrib,'values',values)

                do i=1,Nx
                 do j=1,Ny
                  ensfcst(i,j,k) = ensfcst(i,j,k) + values((j-1)*Nx+i)
                 enddo
                enddo

                deallocate(values)
                call codes_release(igrib)
                call codes_close_file(ifile)
            END IF

            IF (k.eq.nmembers/2) THEN
              ifhr = ifhr + 12
              CALL GETH_NEWDATE(odate, odate, -12*3600)
            END IF
          END DO
        END IF
      ELSE
        refmax=50.
        ctrx0=60.0E03
        ctry0=60.0E03
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
            ctrx=ctrx0
            ctry=ctry0-12.0E03
          ELSE IF(k==5) THEN
            ctrx=ctrx0+12.0E03
            ctry=ctry0+12.0E03
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
      END IF
    END IF

    IF ( nprocs > 1 ) THEN
      !CALL MPI_BCAST(membknt,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      !lstring=LEN(ensfname(1))
      !CALL MPI_BCAST(ensfname,(lstring*nmembers),MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
      !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      nelem3d=nx*ny*nmembers
      nelem2d=nx*ny
      IF( nelem3d < max_elem_send) THEN
        CALL MPI_BCAST(ensfcst,nelem3d,MPI_REAL,root,MPI_COMM_WORLD,ierr)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ELSE IF( nelem2d < max_elem_send) THEN
        DO k=1,nmembers
          CALL MPI_BCAST(ensfcst(1,1,k),nelem2d,MPI_REAL,root, &
                     MPI_COMM_WORLD,ierr)
          !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        END DO
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ELSE
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        DO k=1,nmembers
          CALL LGARRAY_BCASTR(ensfcst(1,1,k),nelem2d,max_elem_send, &
                          myproc,root,MPI_COMM_WORLD,ierr)
        END DO
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      END IF
      IF( myproc == root) THEN
        WRITE(6,'(a)') 'ensfcst broadcast complete'
      ELSE
        IF (realdata == 3) THEN
          ifhr = ifhr + 12
        END IF
      END IF
    END IF

    !IF ( realdata == 2) THEN
    !  do i=1,Nx
    !    do j=1,Ny
    !      lons(i,j) = float((452+i))*0.5
    !      lats(i,j) = 90.0 - float((j-1+60))*0.5
    !    enddo
    !  enddo
    !END IF

!
!-----------------------------------------------------------------------
!
!  Calculate simple mean
!
!-----------------------------------------------------------------------
!
    IF( myproc == root ) THEN
      CALL CPU_TIME(cput2)
      WRITE(6,'(a,f10.2,a)') ' Set-up CPU time: ',(cput2-cput1),' seconds'
      WRITE(6,'(a,i5)') ' Align Option: ',alignopt
    END IF
    IF( alignopt == 2 .OR. calcmean == 1 .OR. lpmopt == 1) THEN
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
        IF( wrtorgmemb == 1 .AND. myproc == root ) THEN
          WRITE(6,'(a,i2,a,f6.1,a,f6.1)') ' Ensemble index: ',k,       &
              '   fcst min: ',fcstmin,'  fcst max: ',fcstmax
          WRITE(varid,'(a2,a2,i2.2)') varid0(1:2),acchrstr,k
          WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
          varunits='mm'
          WRITE(outname,'(a,a,i2.2)') TRIM(runname),'_memb',k
          !CALL wrtvar2(nx,ny,1,ensfcst(1,1,k),varid,varout,varunits,&
          !                  ifmtout,time,outname,shfoutdir,istatus)
        END IF
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
      IF( myproc == root ) THEN
        WRITE(6,'(/a,f6.1,a,f6.1)') ' Ensemble simple mean: fcst min:',&
              favgmin,'  fcst max: ',favgmax
        WRITE(varid,'(a2,a)') acchrstr,'havg'
        WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
        WRITE(outname,'(a,a)') TRIM(runname),'_ensmn'
        !CALL wrtvar2(nx,ny,1,ensmean,varid,varout,varunits,         &
        !                  ifmtout,time,outname,shfoutdir,istatus)
        CALL CPU_TIME(cput3)
        WRITE(6,'(a,f10.2,a)') ' Simple mean CPU time: ',              &
             (cput3-cput2),' seconds'
      END IF
    ELSE
      IF( myproc == root) CALL CPU_TIME(cput3)
    END IF
!
!-----------------------------------------------------------------------
!
!   Calculate shift vectors
!   alignopt = 1   Align to all individual members
!   alignopt = 2   Align to ensemble mean location
!
!-----------------------------------------------------------------------
!
    IF( alignopt == 1) THEN
      wgtvar(:) = 1.
      ensfcst_shf(:,:,:,:) = 0.
      ibkshift(:,:,:,:) = 0
      jbkshift(:,:,:,:) = 0
      !fnorm=-1.0/float(membknt)
      fnorm=1.0/float(membknt)

      DO ipass=1,nshfpass
        DO k1=1,membknt
          !IF( myproc == root ) WRITE(6,'(a,i4)') 'Processing ensemble menber: ',k1
          xshiftmn(:,:,k1,ipass) = 0.
          yshiftmn(:,:,k1,ipass) = 0.

          time = real(ifhr)

          IF (realdata == 3) THEN
            IF (k1.le.membknt/2) THEN
              time = real(ifhr) - 12.
              IF( myproc == root ) WRITE(6,'(a,i4,a,i4)') 'Processing ensemble member: ',k1, ', leadtime: ', nint(time)
            ELSE
              IF( myproc == root ) WRITE(6,'(a,i4,a,i4,a)') &
                                   'Processing ensemble member: ',k1, ', leadtime: ', nint(time), ' (time-lagged)'
            END IF
          ELSE
            IF( myproc == root ) WRITE(6,'(a,i4,a,i4)') 'Processing ensemble member: ',k1, ', leadtime: ', nint(time)
          END IF

          DO k2=1,membknt
            IF( k2 /= k1 ) THEN
              IF ( myproc == root ) THEN
                print *, ''
                print *, k1, k2          
              END IF
              xshift(:,:) = 0.
              yshift(:,:) = 0.
              CALL rshift2dgrd(nx,ny,nvar,ipass,mx_shfpass,mxzone,       &
                 max_elem_send,applyshft,posdef,time,xs,ys,              &
                 ensfcst(:,:,k1),ensfcst(:,:,k2),tem2d1,tem2d2,          &
                 ensfcst_shf(:,:,k1,ipass),                              &
                 istart,jstart,ifinish,jfinish,                          &
                 ibkshift(:,:,k1,ipass),jbkshift(:,:,k1,ipass),          &
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

          END IF
!
!-----------------------------------------------------------------------
!
! Apply shift to ensemble using mean shift and sum to calculate
! ensemble mean
!
!-----------------------------------------------------------------------
!
          IF( myproc == root ) THEN
            print *, ''
            print *, ''   
            !print *, 'member', k1     
            IF (lvldbg > 10) THEN
              CALL a2dmax0(ensfcst(:,:,k1),1,nx,ibgn,iend,           &
                           1,ny,jbgn,jend,vmin,vmax)
              WRITE(6,'(1x,2(a,f13.4))')                                   &
                    'Pre-shift min = ', vmin,',  max =',vmax
            END IF
            CALL movegr(nx,ny,ensfcst(:,:,k1),tem2d3,                    &
                      xshiftmn(:,:,k1,ipass),yshiftmn(:,:,k1,ipass),     &
                      ensfcst_shf(:,:,k1,ipass),                         &
                      ibgn,iend,jbgn,jend,                               &
                      dxfld,dyfld,rdxfld,rdyfld,                         &
                      slopey,alphay,betay)

            DO i=1,noutsmth
              CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                          ensfcst_shf(:,:,k1,ipass),tem2d3,ensfcst_shf(:,:,k1,ipass))
            END DO

            IF (lvldbg > 10) THEN
              CALL a2dmax0(ensfcst_shf(:,:,k1,ipass),1,nx,ibgn,iend,       &
                         1,ny,jbgn,jend,vmin,vmax)
              WRITE(6,'(1x,2(a,f13.4))')                                   &
                    'Post-shift min = ', vmin,',  max =',vmax
            END IF

            IF( wrtshfmemb == 1 ) THEN
              WRITE(varid,'(a,i3.3)') 'shf',k1
              WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
              varunits='mm'
              WRITE(outname,'(a,a,i2.2)') TRIM(runname),'_shf',k
              !CALL wrtvar2(nx,ny,1,ensfcst_shf(1,1,k1),varid,varout,varunits,&
              !                     ifmtout,time,outname,shfoutdir,istatus)
            END IF
            ensshfmn(:,:,ipass)=ensshfmn(:,:,ipass)+ensfcst_shf(:,:,k1,ipass)
            !IF (ipass.lt.nshfpass) THEN
            !  ensfcst(:,:,k1,ipass+1) = ensfcst_shf(:,:,k1,ipass)
            !END IF

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

              IF( myproc == root) WRITE(6,'(a)') 'back-shift broadcast complete'
            END IF

          END IF
        END DO  ! k1

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
        ENDIF

      END DO  ! ipass

      IF( myproc == root ) THEN
        IF( wrtshfmean == 1 ) THEN
          WRITE(varid,'(a3,a2,a)') varid0(1:3),acchrstr,'_'
          WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
          varunits='mm'
          WRITE(outname,'(a,a,i1)') TRIM(runname),'_sam',alignopt
          !CALL wrtvar2(nx,ny,1,ensshfmn,varid,varout,varunits,              &
          !                   ifmtout,time,outname,shfoutdir,istatus)
        END IF
        CALL CPU_TIME(cput4)
        WRITE(6,'(a,f10.2,a)') ' Ensemble aligned mean CPU time: ',(cput4-cput3),' seconds'
      END IF
    ELSE
      wgtvar(:) = 1.
      ensfcst_shf(:,:,:,:) = 0.
      ibkshift(:,:,:,:) = 0
      jbkshift(:,:,:,:) = 0.
!
!     Align each ensemble member to the mean
!
      DO ipass=1,nshfpass
      DO k1=1,membknt
        IF( myproc == root ) WRITE(6,'(a,i4)') 'Processing ensemble menber: ',k1
        xshift(:,:) = 0.
        yshift(:,:) = 0.
        IF (lvldbg > 10) THEN
          CALL a2dmax0(ensfcst(:,:,k1),1,nx,ibgn,iend,             &
                       1,ny,jbgn,jend,vmin,vmax)
          WRITE(6,'(1x,2(a,f13.4))')                &
                'Pre-shift pre-call min = ', vmin,',  max =',vmax
        END IF
        CALL rshift2dgrd(nx,ny,nvar,ipass,mx_shfpass,mxzone,           &
               max_elem_send,applyshft,posdef,time,xs,ys,              &
               ensfcst(:,:,k1),ensmean,tem2d1,tem2d2,                  &
               ensfcst_shf(:,:,k1,ipass),                              &
               istart,jstart,ifinish,jfinish,                          &
               ibkshift(:,:,k1,ipass),jbkshift(:,:,k1,ipass),          &
               xshift,yshift,xsum,ysum,wgtsum,                         &
               dxfld,dyfld,rdxfld,rdyfld,                              &
               slopey,alphay,betay,                                    &
               tem2d3,recv_buf,lvldbg)

        IF (ipass .gt. 1) THEN
          DO j=1,ny
            DO i=1,nx
              !xshift(i,j)=xshift(i,j)+FLOAT(ibkshift(i,j,k1,ipass))
              !yshift(i,j)=yshift(i,j)+FLOAT(jbkshift(i,j,k1,ipass))
            END DO
          END DO
        END IF

        IF (ipass .ne. nshfpass) THEN
          DO j=1,ny
            DO i=1,nx
              !ibkshift(i,j,k1,ipass+1)=NINT(xshift(i,j))
              !jbkshift(i,j,k1,ipass+1)=NINT(yshift(i,j))
            END DO
          END DO
        END IF
!
        IF( myproc == root ) THEN
          IF (lvldbg > 10) THEN
            CALL a2dmax0(ensfcst(:,:,k1),1,nx,ibgn,iend,             &
                       1,ny,ibgn,iend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))')                                   &
                  'Pre-shift min = ', vmin,',  max =',vmax
          END IF
          CALL movegr(nx,ny,ensfcst(:,:,k1),tem2d3,              &
                    xshift,yshift,ensfcst_shf(:,:,k1,ipass),           &
                    ibgn,iend,jbgn,jend,                               &
                    dxfld,dyfld,rdxfld,rdyfld,                         &
                    slopey,alphay,betay)
          IF( myproc == root ) THEN
            CALL a2dmax0(ensfcst_shf(:,:,k1,ipass),1,nx,ibgn,iend,       &
                         1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))') 'Post-shift min = ', vmin,        &
                  ',  max =',vmax
          END IF
          IF( wrtshfmemb == 1 ) THEN
            WRITE(varid,'(a,i3.3)') 'shf',k1
            WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
            varunits='mm'
            WRITE(outname,'(a,a,i2.2)') TRIM(runname),'_shf',k1
            !CALL wrtvar2(nx,ny,1,ensfcst_shf(1,1,k1),varid,varout,varunits,&
            !                   ifmtout,time,outname,shfoutdir,istatus)
          END IF
        END IF
        ensshfmn(:,:,ipass)=ensshfmn(:,:,ipass)+ensfcst_shf(:,:,k1,ipass)
        !IF (ipass.lt.nshfpass) THEN
        !  ensfcst(:,:,k1,ipass+1) = ensfcst_shf(:,:,k1,ipass)
        !END IF
        !ensshfmn(:,:)=ensshfmn(:,:)+ensfcst_shf(:,:,k1)

        IF (ipass.eq.nshfpass) THEN
          IF( nprocs > 1 ) THEN
            nelem2d=nx*ny
            IF( nelem2d < max_elem_send) THEN
              CALL MPI_BCAST(ensshfmn,nelem2d,MPI_REAL,root,MPI_COMM_WORLD,ierr)
              !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            ELSE
              CALL LGARRAY_BCASTR(ensshfmn,nelem2d,max_elem_send, &
                            myproc,root,MPI_COMM_WORLD,ierr)
            END IF
            !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            IF( myproc == root) WRITE(6,'(a)') 'ensfcst broadcast complete'
          END IF
        END IF
      END DO  ! k1
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
      END DO  ! ipass

      IF( myproc == root ) THEN
        IF( wrtshfmean == 1 ) THEN
          WRITE(varid,'(a3,a2,a)') varid0(1:3),acchrstr,'_'
          WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
          varunits='mm'
          !CALL wrtvar2(nx,ny,1,ensshfmn,varid,varout,varunits,              &
          !                     ifmtout,time,runname,shfoutdir,istatus)
        END IF
        CALL CPU_TIME(cput4)
        WRITE(6,'(a,f10.2,a)') ' Ensemble aligned to ensemble mean CPU time: ',(cput4-cput3),' seconds'
      END IF
    END IF


    IF( myproc == root ) THEN
      IF( lpmopt > 0 ) THEN
        CALL CPU_TIME(cput5)
        !enslpm(:,:) = 0.
        CALL lpm_mean_la(nx,ny,ibgn,iend,jbgn,jend, &
                  membknt,patch_nx,patch_ny,ovx,ovy, &
                  filt_min,gauss_sigma,gs_weight,        &
                  ensfcst(:,:,1:membknt),ensmean,ensmax,     &
                  enslpm)
        IF (lvldbg > 10) THEN
          CALL a2dmax0(enslpm,1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
          WRITE(6,'(1x,2(a,f13.4))') 'LPM output min = ', vmin,',  max =',vmax
        END IF
    !    WRITE(varid,'(a3,a2,a)') varid0(1:3),acchrstr,'_'
    !    WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
    !    varunits='mm'
    !    !CALL wrtvar2(nx,ny,1,enslpm,varid,varout,varunits,              &
    !    !                       ifmtout,time,runname,lpmoutdir,istatus)
        CALL CPU_TIME(cput6)
        WRITE(6,'(a,f10.2,a)') ' Ensemble LPM CPU time: ',(cput6-cput5),' seconds'
        CALL CPU_TIME(cput5)
        !ensshflpm(:,:,:) = 0.
        DO ipass=1,nshfpass
          CALL lpm_mean_la(nx,ny,ibgn,iend,jbgn,jend, &
                  membknt,patch_nx,patch_ny,ovx,ovy, &
                  filt_min,gauss_sigma,gs_weight,        &
                  ensfcst_shf(:,:,1:membknt,ipass),ensshfmn(:,:,ipass),ensshfmax(:,:,ipass), &
                  ensshflpm(:,:,ipass))
          IF (lvldbg > 10) THEN
            CALL a2dmax0(ensshflpm(:,:,ipass),1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))') 'Shift LPM output min = ', vmin,',  max=',vmax
          END IF
        ENDDO
    !    WRITE(varid,'(a3,a2,a)') varid0(1:3),acchrstr,'_'
    !    WRITE(varout,'(a,a2,a)') 'P',acchrstr,'M'
    !    varunits='mm'
    !    !CALL wrtvar2(nx,ny,1,ensshflpm,varid,varout,varunits,              &
    !    !                         ifmtout,time,runname,shflpmoutdir,istatus)
        CALL CPU_TIME(cput6)
        WRITE(6,'(a,f10.2,a)') ' Ensemble aligned LPM CPU time:',(cput6-cput5),' seconds'


        CALL pm_mean(nx,ny,membknt,RESHAPE(ensfcst(:,:,1:membknt),(/nx*ny*membknt/)), &
                  ensmean,enspm)
        IF (lvldbg > 10) THEN
          CALL a2dmax0(enspm,1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
          WRITE(6,'(1x,2(a,f13.4))') 'PM output min = ', vmin,',  max =',vmax
        END IF

        DO ipass=1,nshfpass
          CALL pm_mean(nx,ny,membknt,RESHAPE(ensfcst(:,:,1:membknt),(/nx*ny*membknt/)), &
                  ensshfmn(:,:,ipass),ensshfpm(:,:,ipass))
          IF (lvldbg > 10) THEN
            CALL a2dmax0(ensshfpm(:,:,ipass),1,nx,ibgn,iend,1,ny,jbgn,jend,vmin,vmax)
            WRITE(6,'(1x,2(a,f13.4))') 'Shift PM output min = ', vmin,',  max=',vmax
          END IF
        ENDDO

      END IF
    END IF
    
  END DO ! itime

!================================
!write netcdf file
!================================
  IF( myproc == root ) THEN

! Write file(netcdf)
    if (realdata == 3) then
      ifhr = ifhr - 12
    end if

    WRITE(fcstname,'(5a,i3.3,a)') '/non/clee/ens/data/', datehrstr(1:8), '/href_ens_', datehrstr, '_f', ifhr, '.nc'
    iret = nf_create(fcstname, NF_NETCDF4, ncid)
    print *, trim(fcstname)

! Define dimensions
    ndim = 5
    if (realdata == 3) then
      nvar = 17
    else
      nvar = 14
    end if
    allocate(dimids(ndim), vardims(ndim), chunks(ndim), varids(nvar))
    allocate(vardims2(ndim), chunks2(ndim), vardims3(1))

    if (realdata == 3) then
      cx = int(real(Nx)/2)
      cy = int(real(Ny)/2)
      !print *, cx, cy
      allocate(projx(nx), projy(ny))
      do i=1,Nx
        projx(i) = 3.*(i-1-cx)
      enddo
      do j=1,Ny
        projy(j) = 3.*(j-1-cy)
      enddo
    end if

    if (realdata == 3) then
      iret = nf_def_dim(ncid, 'x', nx, dimids(1))
      iret = nf_def_dim(ncid, 'y', ny, dimids(2))
    else
      iret = nf_def_dim(ncid, 'lon', nx, dimids(1))
      iret = nf_def_dim(ncid, 'lat', ny, dimids(2))
    end if
    iret = nf_def_dim(ncid, 'nmembers', nmembers, dimids(3))
    iret = nf_def_dim(ncid, 'nshfpass', nshfpass, dimids(4))
    iret = nf_def_dim(ncid, 'time', 1, dimids(5))  

    vardims(1) = dimids(1)
    vardims(2) = dimids(2)
    vardims(3) = dimids(3)
    vardims(4) = dimids(4)

    vardims2(1) = dimids(1)
    vardims2(2) = dimids(2)
    vardims2(3) = dimids(4)
    vardims2(4) = dimids(3)

    vardims3(1) = dimids(5)

    chunks(1) = nx
    chunks(2) = ny
    chunks(3) = nmembers
    chunks(4) = nshfpass

    chunks2(1) = nx
    chunks2(2) = ny
    chunks2(3) = nshfpass
    chunks2(4) = nmembers

! Define variables
    iret = nf_def_var(ncid, 'lat', nf_float, 2, vardims, varids(1))
    iret = nf_def_var(ncid, 'lon', nf_float, 2, vardims, varids(2))
    iret = nf_def_var(ncid, 'precipitation', nf_float, 3, vardims, varids(3))
    iret = nf_def_var(ncid, 'shiftvec_u', nf_float, 4, vardims, varids(4))
    iret = nf_def_var(ncid, 'shiftvec_v', nf_float, 4, vardims, varids(5))
    iret = nf_def_var(ncid, 'shift_prcp', nf_float, 4, vardims, varids(6))
    iret = nf_def_var(ncid, 'ensmean', nf_float, 2, vardims, varids(7))
    iret = nf_def_var(ncid, 'ensshfmean', nf_float, 3, vardims2, varids(8))
    iret = nf_def_var(ncid, 'enslpm', nf_float, 2, vardims, varids(9))
    iret = nf_def_var(ncid, 'ensshflpm', nf_float, 3, vardims2, varids(10))
    iret = nf_def_var(ncid, 'enspm', nf_float, 2, vardims, varids(11))
    iret = nf_def_var(ncid, 'ensshfpm', nf_float, 3, vardims2, varids(12))
    iret = nf_def_var(ncid, 'time', nf_int, 1, vardims3, varids(13))
    iret = nf_def_var(ncid, 'forecast_reference_time', nf_int, 1, vardims3, varids(14))

    if (realdata == 3) then
      iret = nf_def_var(ncid, 'grid_mapping', nf_int, 0, vardims3, varids(15))
      iret = nf_def_var(ncid, 'x', nf_float, 1, dimids(1), varids(16))
      iret = nf_def_var(ncid, 'y', nf_float, 1, dimids(2), varids(17))
    end if

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
    do k = 1, nvar
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
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(1), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if

    iret = nf_put_att_text(ncid, varids(2), 'long_name', len_trim('longitude'), &
                         'longitude')
    iret = nf_put_att_text(ncid, varids(2), 'standard_name', len_trim('longitude'), &
                         'longitude')
    iret = nf_put_att_text(ncid, varids(2), 'units', len_trim('degrees_east'), 'degrees_east')
    iret = nf_put_att_real(ncid, varids(2), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(2), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(2), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  
    iret = nf_put_att_text(ncid, varids(3), 'long_name', len_trim('precipitation'), &
                         'precipitation')
    iret = nf_put_att_text(ncid, varids(3), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(3), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(3), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(3), 'grid_mapping',len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(3), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(4), 'long_name', len_trim('shiftvec_u'), &
                         'shiftvec_u')
    iret = nf_put_att_text(ncid, varids(4), 'units', len_trim('grid'), 'grid')
    iret = nf_put_att_real(ncid, varids(4), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(4), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(4), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(4), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(5), 'long_name', len_trim('shiftvec_v'), &
                         'shiftvec_v')
    iret = nf_put_att_text(ncid, varids(5), 'units', len_trim('grid'), 'grid')
    iret = nf_put_att_real(ncid, varids(5), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(5), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(5), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(5), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(6), 'long_name', len_trim('shifted precipitation'), &
                         'shifted precipitation')
    iret = nf_put_att_text(ncid, varids(6), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(6), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(6), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(6), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(6), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(7), 'long_name', len_trim('ensmean'), &
                         'ensmean')
    iret = nf_put_att_text(ncid, varids(7), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(7), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(7), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(7), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(7), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(8), 'long_name', len_trim('ensshfmean'), &
                         'ensshfmean')
    iret = nf_put_att_text(ncid, varids(8), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(8), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(8), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(8), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(8), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(9), 'long_name', len_trim('enslpm'), &
                         'enslpm')
    iret = nf_put_att_text(ncid, varids(9), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(9), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(9), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(9), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(9), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(10), 'long_name', len_trim('ensshflpm'), &
                         'ensshflpm')
    iret = nf_put_att_text(ncid, varids(10), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(10), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(10), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(10), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(10), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(11), 'long_name', len_trim('enspm'), &
                         'enspm')
    iret = nf_put_att_text(ncid, varids(11), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(11), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(11), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(11), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(11), 'valid_time', len_trim(validtime), validtime)

    iret = nf_put_att_text(ncid, varids(12), 'long_name', len_trim('ensshfpm'), &
                         'ensshfpm')
    iret = nf_put_att_text(ncid, varids(12), 'units', len_trim('mm'), 'mm')
    iret = nf_put_att_real(ncid, varids(12), '_FillValue', nf_float, 1, -999.)
    iret = nf_put_att_text(ncid, varids(12), 'coordinates', len_trim('lat lon'), 'lat lon')
    if (realdata == 3) then
      iret = nf_put_att_text(ncid, varids(12), 'grid_mapping', len_trim('grid_mapping'), &
                         'grid_mapping')
    end if
  !iret = nf_put_att_text(ncid, varids(12), 'valid_time', len_trim(validtime), validtime)

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


    if (realdata == 3) then
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
    end if



    iret = nf_put_att_text(ncid, NCGLOBAL, 'Conventions', len_trim('CF-1.6'), 'CF-1.6')

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
    if (realdata == 3) then
      iret = nf_put_var_int(ncid, varids(15), 0)
      iret = nf_put_var_real(ncid, varids(16), projx)
      iret = nf_put_var_real(ncid, varids(17), projy)
    end if

! Close Netcdf File
    iret = nf_close(ncid)
    deallocate(dimids, vardims, chunks, varids, chunks2)
    if (realdata == 3) then
      deallocate(projx, projy)
    endif

  END IF

  CALL MPI_FINALIZE(ierr)
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

SUBROUTINE readvar2net(nx,ny,nz,fname,varid,var,varname,varunits,istatus)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz
  CHARACTER(LEN=*), INTENT(IN) :: fname
  CHARACTER(LEN=*), INTENT(IN) :: varid
  REAL, INTENT(OUT) :: var(nx,ny,nz)
  CHARACTER(LEN=*), INTENT(OUT) :: varname
  CHARACTER(LEN=*), INTENT(OUT) :: varunits
  INTEGER, INTENT(OUT) :: istatus
!
  INTEGER :: imode,nxid,nyid,nzid,ncvarid,ncid,istat
  INTEGER :: nxin,nyin,nzin
  INCLUDE 'netcdf.inc'
  istatus = 0
  var = -999.
  varname = 'NULL'
  varunits = 'NULL'
  imode = NF_NOWRITE
  istat = NF_OPEN(fname,imode,ncid)
  IF( istat == NF_NOERR) THEN
    nxin = -1
    nyin = -1
    nzin = -1
    istat = NF_INQ_DIMID(ncid,'x',nxid)
    istat = NF_INQ_DIMLEN(ncid,nxid,nxin)
    istat = NF_INQ_DIMID(ncid,'y',nyid)
    istat = NF_INQ_DIMLEN(ncid,nyid,nyin)
    istat = NF_INQ_DIMID(ncid,'z',nzid)
    istat = NF_INQ_DIMLEN(ncid,nzid,nzin)
    IF( nxin == nx .AND. nyin == ny .AND. nzin == nz) THEN
      istat = NF_INQ_VARID(ncid,varid,ncvarid)
      IF( istat == NF_NOERR ) THEN
        istat = NF_GET_VAR_REAL(ncid,ncvarid,var)
        IF( istat == NF_NOERR ) THEN
          istat = NF_GET_ATT_TEXT(ncid,ncvarid,'long_name',varname)
          istat = NF_GET_ATT_TEXT(ncid,ncvarid,'units',varunits)
        ELSE
          WRITE(6,'(3a,i4)') 'Error reading variable ',TRIM(varid),' istat=',istat
          istatus = -4
          RETURN
        END IF
      ELSE
        WRITE(6,'(3a,i4)') 'Error finding variable ',TRIM(varid),' istat=',istat
        istatus = -3
        RETURN
      END IF
    ELSE
      WRITE(6,'(a,3i6)') 'Dimensions in file: ',nxin,nyin,nzin
      WRITE(6,'(a,3i6)') 'Do not match dims: ',nx,ny,nz
      istatus=-2
      RETURN
    END IF
  ELSE
    WRITE(6,'(3a,i4)') 'Error opening file ',TRIM(fname),' istat=',istat
    istatus = -1
  END IF
  RETURN 
END SUBROUTINE readvar2net

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
