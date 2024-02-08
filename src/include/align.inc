!
!#######################################################################
!
!     Alignment shift control parameters
!     Input via namelist in INITALIGN 
!
!#######################################################################
!
  INTEGER :: nprocs
  INTEGER :: myproc

  COMMON /mpi_vars/ nprocs,myproc

  INTEGER, PARAMETER :: nvarshf=1
  INTEGER, PARAMETER :: mx_shfpass=5
  INTEGER, PARAMETER :: mx_members=25

  CHARACTER(LEN=256) :: basedir
  CHARACTER(LEN=256) :: shfoutdir
  CHARACTER(LEN=256) :: lpmoutdir
  CHARACTER(LEN=256) :: shflpmoutdir
  CHARACTER(LEN=256) :: grb2_template
  CHARACTER(LEN=10) :: datehrstr
  CHARACTER(LEN=20) :: varfname
  CHARACTER(LEN=20) :: varname
  CHARACTER(LEN=20) :: membname(mx_members)
  INTEGER :: ibgn,iend,jbgn,jend
  INTEGER :: tintvl,tbegin,tend,nmembers

  COMMON /ens_info/ ibgn,iend,jbgn,jend,tintvl,tbegin,tend,nmembers
  COMMON /ens_name/ basedir,shfoutdir,lpmoutdir,shflpmoutdir, &
                    datehrstr,varfname,varname,membname

  INTEGER :: alignopt,calcmean,iposdef,nshfpass,minkdat

  INTEGER :: nbaksmth(mx_shfpass),nshfsmth,noutsmth,ifmtout
  INTEGER :: wrtorgmemb,wrtshfmemb,wrtshfmean,grb2out
  REAL :: hrzlap,minkdratio,slnratio0h,slnratio48h,slen,reflmin

  COMMON /align_parm/ alignopt,calcmean,iposdef,nshfpass,minkdat,      &
                      nbaksmth,nshfsmth,noutsmth,ifmtout,hrzlap,       &
                      minkdratio,slnratio0h,slnratio48h,slen,reflmin,  &
                      wrtorgmemb,wrtshfmemb,wrtshfmean,grb2out,        &
                      grb2_template

  INTEGER :: iborder(mx_shfpass),jborder(mx_shfpass), &
             izsize(mx_shfpass),jzsize(mx_shfpass),  &
             nizone(mx_shfpass),njzone(mx_shfpass), &
             loopstep(mx_shfpass),procspg(mx_shfpass)
  COMMON /align_zones/ iborder,jborder,izsize,jzsize,nizone,njzone,    &
                       loopstep,procspg

  REAL :: wgtvar(nvarshf),threshvar(nvarshf)
  INTEGER :: thresh_flag(nvarshf)
  COMMON /shift_wgts/ wgtvar,thresh_flag,threshvar
!
!-----------------------------------------------------------------------
!
!   Localized Probability-Matched Mean Parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: patch_nx,patch_ny,ovx,ovy,gauss_sigma
  REAL :: filt_min
  COMMON /lpm_parm/ patch_nx,patch_ny,ovx,ovy,gauss_sigma,filt_min
