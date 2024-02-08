!**************************************************************************
!  This subroutine support various utility module, such as 
!    o get unit
!                                    coded by Jun-Tae Choi Jan. 10, 2005
!==========================================================================
!  Purpose :
!  Find a free FORTRAN I/O unit from a list and return that unit.
!  which is one of 11 to 49
!                                    coded by Jun-Tae Choi Jan. 15, 2005
!*************************************************************************
!-------------------------------------------------------------------------
  SUBROUTINE GET_UNIT(unit_file)
!-------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------
! Declarations of Subroutine Argument

  integer            :: unit_file

!-------------------------------------------------------------------------
! Variables Declarations

  integer, PARAMETER :: unit_first=11, unit_last=49
  logical            :: used
 
!-------------------------------------------------------------------------
! Control Variables Declaration

  integer            :: loop

!-------------------------------------------------------------------------


  do loop = unit_first, unit_last
     INQUIRE(unit=loop, opened=used)
     if (.not.used) then
         unit_file = loop
!        write(*,'(a,i2,a)') "Fortran I/O unit ", unit_file, &
!                   " picked from the free list."
         return
     endif
  enddo

  write(*,*) "Failing to get Fortran I/O unit and Program stop"
!
  END SUBROUTINE GET_UNIT
!
!========================================================================
  SUBROUTINE CAL_LOCALTIME(date,ndate,IDLT)
!========================================================================
!
  IMPLICIT NONE
  INTEGER, DIMENSION(12) :: DAY = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER                :: IZMN,IZDY,IZHR,IDLT, IZYR, IADY
  INTEGER                :: IYR,IMN,IDY,IHR
  INTEGER, DIMENSION(4)  :: ndate
  INTEGER, DIMENSION(4)  :: date 
!
      IYR = date(1) ; IMN = date(2) ; IDY = date(3) ; IHR = date(4)
!
      IF((MOD(IYR,4).EQ.0.AND.MOD(IYR,100).NE.0).OR. &
                              MOD(IYR,400).EQ.0) DAY(2)=29
!
      IZYR=IYR
      IZMN=IMN
      IZHR=MOD(IHR+IDLT,24)
      IADY=(IHR+IDLT)/24  
      IZDY=IDY+IADY

!---CHECK DAY,MONTH,YEAR
      IF(IZDY.GT.DAY(IZMN))THEN
         IF(IZMN.GE.1.AND.IZMN.LE.11)THEN
           IZDY=IZDY-DAY(IZMN)
           IZMN=IZMN+1
         ELSEIF(IZMN.GE.12)THEN
           IZDY=IZDY-DAY(12)
           IZMN=1
           IZYR=IZYR+1
         ENDIF
      ELSEIF(IZDY.LE.0)THEN
         IF(IZMN.GE.2.AND.IZMN.LE.12)THEN
           IZMN=IZMN-1
           IZDY=IZDY+DAY(IZMN)
         ELSEIF(IZMN.LE.1)THEN
           IZMN=12
           IZDY=IZDY+DAY(IZMN)
           IZYR=IZYR-1
         ENDIF
      ENDIF
!
      ndate(1) = IZYR ; ndate(2) = IZMN ; ndate(3) = IZDY ; ndate(4) = IZHR  
!
  END SUBROUTINE CAL_LOCALTIME
!
!========================================================================
  SUBROUTINE CAL_DATE(date,ndate,IDLD,IDLT)
!========================================================================
!
  IMPLICIT NONE
  INTEGER, DIMENSION(12) :: DAY = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER                :: IZMN,IZDY,IZHR,IDLT, IZYR, IADY, IDLD
  INTEGER                :: IYR,IMN,IDY,IHR
  INTEGER, DIMENSION(4)  :: ndate
  INTEGER, DIMENSION(4)  :: date
!
      IYR = date(1) ; IMN = date(2) ; IDY = date(3) ; IHR = date(4)
!
      IF((MOD(IYR,4).EQ.0.AND.MOD(IYR,100).NE.0).OR. &
                              MOD(IYR,400).EQ.0) DAY(2)=29
!
      IZYR=IYR
      IZMN=IMN
      IZHR=MOD(IHR+IDLT,24)
      IADY=(IHR+IDLT)/24+IDLD
      IZDY=IDY+IADY

!---CHECK DAY,MONTH,YEAR
      IF(IZDY.GT.DAY(IZMN))THEN
         IF(IZMN.GE.1.AND.IZMN.LE.11)THEN
           IZDY=IZDY-DAY(IZMN)
           IZMN=IZMN+1
         ELSEIF(IZMN.GE.12)THEN
           IZDY=IZDY-DAY(12)
           IZMN=1
           IZYR=IZYR+1
         ENDIF
      ELSEIF(IZDY.LE.0)THEN
         IF(IZMN.GE.2.AND.IZMN.LE.12)THEN
           IZMN=IZMN-1
           IZDY=IZDY+DAY(IZMN)
         ELSEIF(IZMN.LE.1)THEN
           IZMN=12
           IZDY=IZDY+DAY(IZMN)
           IZYR=IZYR-1
         ENDIF
      ENDIF
!
      ndate(1) = IZYR ; ndate(2) = IZMN ; ndate(3) = IZDY ; ndate(4) = IZHR
!
  END SUBROUTINE CAL_DATE
