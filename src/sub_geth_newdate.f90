SUBROUTINE GETH_NEWDATE(ndate, odate, idt)
!-------------------------------------------------------------------------
!  Purpose
!    To get new date with old date and offset idt
!
!  Variables
!    ndate(4) : new date (output)
!    odate(4) : old date (input)
!    idt      : offset in second
!               3hour after : idt=3*3600
!               2day and 3hour before : idt=(2*24+3)*3600
!
!  History 
!    Aug. 24, 2008 JunTae Choi
!       Coded to provide standard module
!-------------------------------------------------------------------------
  IMPLICIT NONE

! Declarations of subroutine argument

  INTEGER, DIMENSION(4) :: ndate, odate
  INTEGER               :: idt

! Declaration of local varibles

  INTEGER      :: nlen, olen
  INTEGER      :: yrnew, monew, dynew, hrnew, minew, scnew, frnew
  INTEGER      :: yrold, moold, dyold, hrold, miold, scold, frold
  INTEGER      :: mday(12), nday, nhour, nmin, nsec, nfrac, i, ifrc
  LOGICAL      :: opass
  CHARACTER*10 :: hfrc
  CHARACTER*1  ::  sp
       
      !  Local Variables
       
      !  yrold    -  indicates the year associated with "odate"
      !  moold    -  indicates the month associated with "odate"
      !  dyold    -  indicates the day associated with "odate"
      !  hrold    -  indicates the hour associated with "odate"
      !  miold    -  indicates the minute associated with "odate"
      !  scold    -  indicates the second associated with "odate"
       
      !  yrnew    -  indicates the year associated with "ndate"
      !  monew    -  indicates the month associated with "ndate"
      !  dynew    -  indicates the day associated with "ndate"
      !  hrnew    -  indicates the hour associated with "ndate"
      !  minew    -  indicates the minute associated with "ndate"
      !  scnew    -  indicates the second associated with "ndate"
       
      !  mday     -  a list assigning the number of days in each month
      
      !  i        -  loop counter
      !  nday     -  the integer number of days represented by "idt"
      !  nhour    -  the integer number of hours in "idt" after taking out
      !              all the whole days
      !  nmin     -  the integer number of minutes in "idt" after taking out
      !              all the whole days and whole hours.
      !  nsec     -  the integer number of minutes in "idt" after taking out
      !              all the whole days, whole hours, and whole minutes.
       
!-------------------------------------------------------------------------
      
      mday( 1) = 31
      mday( 2) = 28
      mday( 3) = 31
      mday( 4) = 30
      mday( 5) = 31
      mday( 6) = 30
      mday( 7) = 31
      mday( 8) = 31
      mday( 9) = 30
      mday(10) = 31
      mday(11) = 30
      mday(12) = 31
      
      !  Break down old hdate into parts
      
      hrold = 0
      miold = 0
      scold = 0
      frold = 0
      
      !  Use internal READ statements to convert the CHARACTER string
      !  date into INTEGER components.
      yrold = odate(1)
      moold = odate(2)
      dyold = odate(3)
      hrold = odate(4)
      miold = 0
      scold = 0 
      frold = 0 
   
      !  Set the number of days in February for that year.
      
!     mday(2) = nfeb(yrold)
      CALL      nfeb(yrold, mday(2))
      
      !  Check that ODATE makes sense.
      
      opass = .TRUE.
      
      !  Check that the month of ODATE makes sense.
      
      IF ((moold.GT.12).or.(moold.LT.1)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Month of ODATE = ', moold
         opass = .FALSE.
      END IF
      
      !  Check that the day of ODATE makes sense.
      
      IF ((dyold.GT.mday(moold)).or.(dyold.LT.1)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Day of ODATE = ', dyold
         opass = .FALSE.
      END IF
      
      !  Check that the hour of ODATE makes sense.
      
      IF ((hrold.GT.23).or.(hrold.LT.0)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Hour of ODATE = ', hrold
         opass = .FALSE.
      END IF
      
      !  Check that the minute of ODATE makes sense.
      
      IF ((miold.GT.59).or.(miold.LT.0)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Minute of ODATE = ', miold
         opass = .FALSE.
      END IF
      
      !  Check that the second of ODATE makes sense.
      
      IF ((scold.GT.59).or.(scold.LT.0)) THEN
         WRITE(*,*) 'GETH_NEWDATE:  Second of ODATE = ', scold
         opass = .FALSE.
      END IF
      
      !  Check that the fractional part  of ODATE makes sense.
      
      !KWM      IF ((scold.GT.59).or.(scold.LT.0)) THEN
      !KWM         WRITE(*,*) 'GETH_NEWDATE:  Second of ODATE = ', scold
      !KWM         opass = .FALSE.
      !KWM      END IF
      
      
      !  Date Checks are completed.  Continue.
      
      
      !  Compute the number of days, hours, minutes, and seconds in idt
      
         ifrc = 1
         nday   = ABS(idt)/86400 ! Integer number of days in delta-time
         nhour  = MOD(ABS(idt),86400)/3600
         nmin   = MOD(ABS(idt),3600)/60
         nsec   = MOD(ABS(idt),60)
         nfrac  = 0
      
      IF (idt.GE.0) THEN
      
         frnew = frold + nfrac
         IF (frnew.GE.ifrc) THEN
            frnew = frnew - ifrc
            nsec = nsec + 1
         END IF
      
         scnew = scold + nsec
         IF (scnew .GE. 60) THEN
            scnew = scnew - 60
            nmin  = nmin + 1
         END IF
      
         minew = miold + nmin
         IF (minew .GE. 60) THEN
            minew = minew - 60
            nhour  = nhour + 1
         END IF
      
         hrnew = hrold + nhour
         IF (hrnew .GE. 24) THEN
            hrnew = hrnew - 24
            nday  = nday + 1
         END IF
      
         dynew = dyold
         monew = moold
         yrnew = yrold
         DO i = 1, nday
            dynew = dynew + 1
            IF (dynew.GT.mday(monew)) THEN
               dynew = dynew - mday(monew)
               monew = monew + 1
               IF (monew .GT. 12) THEN
                  monew = 1
                  yrnew = yrnew + 1
                  ! If the year changes, recompute the number of days in February
                  CALL nfeb(yrnew, mday(2))
               END IF
            END IF
         END DO
      
      ELSE IF (idt.LT.0) THEN
      
         frnew = frold - nfrac
         IF (frnew .LT. 0) THEN
            frnew = frnew + ifrc
            nsec = nsec - 1
         END IF
      
         scnew = scold - nsec
         IF (scnew .LT. 00) THEN
            scnew = scnew + 60
            nmin  = nmin + 1
         END IF
      
         minew = miold - nmin
         IF (minew .LT. 00) THEN
            minew = minew + 60
            nhour  = nhour + 1
         END IF
      
         hrnew = hrold - nhour
         IF (hrnew .LT. 00) THEN
            hrnew = hrnew + 24
            nday  = nday + 1
         END IF
      
         dynew = dyold
         monew = moold
         yrnew = yrold
         DO i = 1, nday
            dynew = dynew - 1
            IF (dynew.eq.0) THEN
               monew = monew - 1
               IF (monew.eq.0) THEN
                  monew = 12
                  yrnew = yrnew - 1
                  ! If the year changes, recompute the number of days in February
!                 mday(2) = nfeb(yrnew)
                  CALL nfeb(yrnew, mday(2))
               END IF
               dynew = mday(monew)
            END IF
         END DO
      END IF
      
      !  Now construct the new mdate
      
      ndate(1) = yrnew
      ndate(2) = monew
      ndate(3) = dynew
      ndate(4) = hrnew
      
END SUBROUTINE GETH_NEWDATE

!-------------------------------------------------------------------------

      SUBROUTINE nfeb ( year, num_days )
   
      ! Compute the number of days in February for the given year
   
      IMPLICIT NONE
   
      INTEGER  year
      INTEGER  num_days
   
      num_days = 28 ! By default, February has 28 days ...
      IF (MOD(year,4).eq.0) THEN  
         num_days = 29  ! But every four years, it has 29 days ...
         IF (MOD(year,100).eq.0) THEN
            num_days = 28  ! Except every 100 years, when it has 28 days ...
            IF (MOD(year,400).eq.0) THEN
               num_days = 29  ! Except every 400 years, when it has 29 days.
            END IF
         END IF
      END IF
   
      END SUBROUTINE nfeb

!-------------------------------------------------------------------------

      SUBROUTINE mk_date(date,date1)
      INTEGER day, date(4), date1(4), nnn
      INTEGER no_day(12)
      DATA no_day/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      nnn = date(3)
      do i = 1,date(2)-1
         nnn = nnn + no_day(i)
      enddo

      if ((date(1).eq.2000).and.(date(2).gt.2)) nnn = nnn + 1
      day = nnn

      open(unit=3, file='./jtime')
      write(3,'(i4,i3.3,i2.2,i4,3i2.2)') &
            date(1), day, date(4), (date1(i),i=1,4)
      close(3)

      return
      END SUBROUTINE mk_date


!-------------------------------------------------------------------------

      SUBROUTINE calc_hour(date,totalhour)
      INTEGER :: date(4)
      INTEGER :: totaldate, totalhour
      INTEGER :: yyyy, mm, dd, hh
      integer :: jday1 = 2440588                    ! = 1970/01/01

      yyyy = date(1)
      mm = date(2)
      dd = date(3)
      hh = date(4)

      totaldate = dd-32075+1461*(yyyy+4800+(mm-14)/12)/4   &
                 +367*(mm-2-(mm-14)/12*12)/12              &
                 -3*((yyyy+4900+(mm-14)/12)/100)/4  -jday1

      totalhour = totaldate*24 + hh

      END SUBROUTINE calc_hour

