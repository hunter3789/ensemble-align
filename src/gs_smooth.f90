SUBROUTINE gaussian_weight(nx,ny,gs_scale,gs_weight)
  IMPLICIT NONE
  INTEGER :: nx,ny
  INTEGER :: gs_scale
  REAL :: gs_weight(nx*ny)
!
! Misc Local Variables
!
  INTEGER :: ng,nmax,nn,i0,j0
  REAL :: pi
  REAL :: sigma,sigmasq,sqng,distsq
!
! Beginning of executable code
!
!
! Initializations
!
  pi=4.0*atan(1.0)
  nmax=nx*ny
  gs_weight(:) = 0.
!
  IF(gs_scale > 0) THEN
    sigma = float(gs_scale)
    sigmasq = sigma*sigma
    ng = 2*gs_scale
    sqng = float(ng*ng)
    nn=0
    DO i0=-ng,ng
      DO j0=-ng,ng
        nn = nn+1
        IF( nn > nmax ) EXIT
        distsq = float(i0*i0+j0*j0)
        IF(distsq <= sqng) THEN
          gs_weight(nn) = exp(-0.5*distsq/sigmasq)/(2.*pi*sigmasq)
        END IF
      END DO
    END DO
  END IF
END SUBROUTINE gaussian_weight
SUBROUTINE gaussian_smooth(nx,ny,var,sigma,partweight,smprob)
!
! Author: Fanyou Kong
! 4/5/2010
! Rewrite based on the raw code (makePP) from Patrick Marsh
!
! PURPOSE: Apply a 2D Gaussian kernel smoothing to neighborhood
!          probability to generate 'Practically Perfect' (PP)
!          forecast [Silverman 1986; Brooks et al. ??]
!
! INPUTS
!   var is an array of unsmoothed probability of exceedance
!   sigma is the smoothing parameter (in terms of gridpoints)
!   partweight is the weighting array
!
! OUTPUT
!   smprob is an array of smoothed probability of exceedance

  INTEGER :: nx,ny
  INTEGER :: i,j,kk,i1,j1,ng,nn
  REAL :: sigma
  REAL :: var(nx,ny),smprob(nx,ny)
  REAL :: partweight(nx*ny)

  IF(int(sigma) == 0) THEN  ! No gaussian smoothing is applied
    smprob = var
    print *,'No gaussian smoothing!'
    RETURN
  END IF

  ng = int(2.*sigma)   ! 2 times sigma
!  ng = int(5.*sigma)   ! 5 times sigma

  smprob = 0.0

! Now go through the Practically Perfect Gaussian code...
!
! Now apply the weights at all points within 5*sigma of a "hit"...
!
  do i=1,nx
    do j=1,ny
      if(var(i,j) > 0.0) then
        nn = 0
!if(i==700 .and. j==380) print *,'i-ng,i+ng,j-ng,j+ng:',i-ng,i+ng,j-ng,j+ng
if(i==1 .and. j==1) print *,''
        do i1=i-ng,i+ng
          do j1=j-ng,j+ng
            nn = nn+1
            if(i1 < 1 .or. i1 > nx .or. j1 < 1 .or. j1 > ny) cycle

            smprob(i1,j1) = smprob(i1,j1)+partweight(nn)*var(i,j)
          enddo
        enddo
      endif
    enddo
  enddo
!print *,'smprob(700,380),var(700,380):',smprob(700,380),var(700,380)

  RETURN
END SUBROUTINE gaussian_smooth
