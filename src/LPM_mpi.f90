
SUBROUTINE lpm_mean_la(nx,ny,ibgn,iend,jbgn,jend, &
                    n_ens,patch_nx,patch_ny,ovx, ovy, &
                    filt_min,gauss_sigma,partweight,        &
                    var2d_ens,var2d_enmean,var2d_enmax,     &
                    var2d_lpm)
!
!Author: Nate Snook (28 Apr. 2017)
!
!Purpose: Calculate localized PM mean of a 2D field.  Localized PM mean
!         uses all data in a non-exclusive calculation area to obtain
!         PM mean for a smaller patch within that region of dimensions
!         (patch_nx, patch_ny).  The size of the calculation area is
!         controled by the overlap variables ovx and ovy.  
!
!         Here is a diagram of the design:
!
!         +------------------------------------------------------------+
!         |  full domain                                               |
!         |        +---------------------------+                       |
!         |        | calculation ^             |                       |
!         |        | area        |             |                       |
!         |        |            ovy            |                       |
!         |        |             |             |                       |
!         |        |             v             |                       |
!         |        |         +-------+         |                       |
!         |        |         | patch |         |                       |
!         |        |<--ovx-->|       |<--ovx-->|                       |
!         |        |         +-------+         |                       |
!         |        |             ^             |                       |
!         |        |             |             |                       |
!         |        |            ovy            |                       |
!         |        |             |             |                       |
!         |        |             v             |                       |
!         |        +---------------------------+                       |
!         |                                                            |
!         |                                                            |
!         |                                                            |
!         |                                                            |
!         +------------------------------------------------------------+
!
!         This calculation is repeated for patches in the x- and y-direction
!         such that the entire domain is covered, and then a smoother is run
!         on the resulting field to ensure a smooth 2D output.  For best 
!         results, patch_nx and patch_ny of no larger than 8 are recommended.
!
! History
!
! 6/2/2017 Fanyou Kong
!   - Code revision, bug fix, and clean up
! 6/7/2017 Fanyou Kong
!   - Add a constraint with ensemble max (var2d_enmax)
! 6/18/2021 Keith Brewster
!   - _la (limited area) version with custom ranges for i,j, not necessarily CAPS scalar
!   ranges and able to handle off-size edge patches
! 12/2/2022 ChangJae Lee
!   - Fix the offset error in the edge case
! 12/15/2023 ChangJae Lee
!   - revise the program to applicable to use MPI
!

IMPLICIT NONE

INCLUDE 'mpif.h'

!Inputs
INTEGER :: nx, ny, n_ens    ! number of x, y gridpoints, ensemble members
INTEGER :: ibgn,iend        ! i index range of valid physical (not BC) points (ARPS: 2,nx-2)
INTEGER :: jbgn,jend        ! j index range of valid physical (not BC) points (ARPS: 2,ny-2)
INTEGER :: patch_nx, patch_ny ! size of actual patches in gridpoints (recommended value <= 8)
                            !       For HMT QPF, patch_nx = patch_ny = 6 works rather well.
INTEGER :: ovx, ovy         ! how many points of patch "overlap" to include in calculation area.
                            !       Recommended value: 3 to 10 times patch_nx and patch_ny
                            !       For HMT QPF, a value of ovx = ovy = 30 works quite well.
REAL :: filt_min            ! LPM is set to zero where enmean < filt_min  -- this helps to remove
                            !       areas of noise that can occur due to patch overlap.  For QPF,
                            !       filt_min = 0.1 mm for 3hr accumulated precip works well.
REAL :: gauss_sigma         ! Parameter sigma (standard deviation) to be passed to gauss_rad
                            !       Recommended value = 2 (2 gridpoints) -- this needs some 
                            !       testing
REAL :: partweight(nx*ny)   ! Needed input parameter for gaussian_smooth
REAL :: var2d_ens(nx, ny, n_ens)  ! 2D (nx, ny) field for each member (n_ens)
REAL :: var2d_enmean(nx, ny)      ! 2D (nx, ny) ensemble mean field
REAL :: var2d_enmax(nx, ny)       ! 2D (nx, ny) ensemble max field
!Outputs
REAL :: var2d_lpm(nx, ny)         ! 2D (nx, ny) LPM field

!Calculated and used in this code
INTEGER :: ipatches, jpatches  ! number of patches in i, j directions
INTEGER :: ipatch, jpatch      ! iterator for patches in i, j directions
INTEGER :: pe_west, pe_east, pe_north, pe_south ! patch edge locations (i or j point)
INTEGER :: ce_west, ce_east, ce_north, ce_south ! calculation area edge locations
INTEGER :: calc_nx, calc_ny    ! The number of x, y gridpoints within the calculation area
INTEGER :: offset_i, offset_j  ! The distance of the patch edge from the W, S edge of calc area

REAL, ALLOCATABLE :: lpm_smoothed(:,:) ! 2D (nx, ny) LPM field after gaussian smoothing
REAL, ALLOCATABLE :: lpm_calc(:,:)  !The LPM over the calculation area.  Must be allocated since
                                    !we do not know ahead of time how big the calc area will be.
REAL, ALLOCATABLE :: var0(:)   ! The values from ALL ensemble members at all points on the 2D
                               ! slice of interest crammed into a 1D array (needed for call to
                               ! pm_mean)
                               ! Since this is limited to calc area, size is unknown ahead of 
                               ! time.
INTEGER :: ix,jy
INTEGER :: ibgnp1,iendm1,jbgnp1,jendm1
INTEGER :: iendpatch,jendpatch

!! ChangJae Lee
INTEGER :: add_nx, add_ny                         ! The distance of patch edge
INTEGER :: mpatch, buf_size, ierr, nprocs, myproc ! MPI varibales
INTEGER, PARAMETER :: root = 0                    ! root processor
REAL, ALLOCATABLE  :: recv_buf(:,:)               ! buffer array
!! I used common variables, but you can get those variables as arguments
COMMON /mpi_vars/ nprocs, myproc                  ! MPI varibales
!!!!!!!!!

!---------------------------------------------------------------!
!Initialize output array as zero as points outside valid area not assigned otherwise.
var2d_lpm = 0.

!---------------------------------------------------------------!
!Determine whether specified patch_nx, patch_ny divide evenly into the domain
!and if they do not, return an error.
!!! Kong comment: This might not be necessary - need to loosen
iendpatch = 0
IF (MOD(((iend-ibgn)+1), patch_ny) /= 0) THEN
!   CALL arpsstop('LPM MEAN ERROR: patches do not divide evenly into domain (y dir)', 1)
    iendpatch = 1
END IF

jendpatch = 0
IF (MOD(((jend-jbgn)+1), patch_nx) /= 0) THEN
!   CALL arpsstop('LPM MEAN ERROR: patches do not divide evenly into domain (x dir)', 1)
    jendpatch = 1
END IF

IF( myproc == root) THEN
    ALLOCATE(lpm_smoothed(nx,ny))
END IF

IF ( nprocs > 1 ) THEN
  buf_size=nx*ny
  ALLOCATE(recv_buf(nx,ny))
END IF

ipatches=iendpatch + (((iend-ibgn)+1)/patch_nx)
jpatches=jendpatch + (((jend-jbgn)+1)/patch_ny)

!print *,'LPM -- ipatches,jpatches:',ipatches,jpatches
!print *,'LPM -- ibgn,iend:',ibgn,iend
!print *,'LPM -- jbgn,jend:',jbgn,jend

ibgnp1 = ibgn + 1
iendm1 = iend - 1
jbgnp1 = jbgn + 1
jendm1 = jbgn + 1

!Now loop over each patch and perform the PM mean in that patch's calculation area
!   NOTE: where calculation areas extend outside the domain boundary, they are clipped to fit
DO ipatch=1, ipatches
    DO jpatch=1, jpatches

        !! ChangJae LEE
        !! Assign each patch's job to each processor
        mpatch = (ipatch - 1) * jpatches + jpatch
        IF( mod(mpatch,nprocs) == myproc ) THEN

            !Calculate patch edge locations (i for pe_west, pe_east; j for pe_south, pe_north)
            !   NOTE: the "2 + " is to deal with the non-physical gridpoint on S, W arps domain edges
            !         (it excludes the westmost/southmost point (which has index "1" in FORTRAN)
            pe_west = ibgn + ((ipatch-1) * patch_nx)
            pe_east = pe_west + patch_nx
            pe_east = min(pe_east,iend)
            pe_south= jbgn + ((jpatch-1) * patch_ny)
            pe_north= pe_south + patch_ny
            pe_north= min(pe_north,jend)

            !Calculate nominal edges of calculation area (before clipping, if any is needed)
            ce_west = pe_west - ovx
            ce_east = pe_east + ovx
            ce_south = pe_south - ovx
            ce_north = pe_north + ovx

            !! ChangJae Lee
            !! Catch edge cases:
            offset_i = ovx + 1    !If there is no clipping, offset_i is the same as ovx
            offset_j = ovy + 1    !If there is no clipping, offset_j is the same as ovy
            add_nx = patch_nx
            add_ny = patch_ny
            !!!!!!!!!

            IF (ce_west < ibgnp1) THEN    !West edge case
                offset_i = (pe_west - ibgn) + 1
                ce_west = ibgn
            END IF
       
            IF (ce_east > iend) THEN  !East edge case
                ce_east = iend   !No need to change offset in this case (only cares about offset
                                   !from west side of calcuation area, not east side)
            END IF
       
            IF (ce_south < jbgnp1) THEN    !South edge case
                offset_j = (pe_south - jbgn) + 1
                ce_south = jbgn
            END IF
       
            IF (ce_north > jend) THEN    !North edge case
                ce_north = jend   !No need to change offset in this case (only cares about offset
                                    !from south side of calcuation area, not north side)
            END IF

            !! ChangJae Lee
            !! Catch edge cases:
            IF (ce_west + offset_i + patch_nx - 1 > ce_east) THEN
                add_nx = ce_east - (ce_west + offset_i - 1)
            END IF

            IF (ce_south + offset_j + patch_ny - 1 > ce_north) THEN
                add_ny = ce_north - (ce_south + offset_j - 1)
            END IF   

            !! No overlapping:
            if (pe_east .ne. iend) then
                pe_east = pe_east - 1
                add_nx = add_nx - 1
            end if

            if (pe_north .ne. jend) then
                pe_north = pe_north - 1
                add_ny = add_ny - 1
            end if
            !!!!!!!!!
 
            !Define calc_nx, calc_ny:
            calc_nx = 1 + (ce_east - ce_west)
            calc_ny = 1 + (ce_north - ce_south)

            !Allocate LPM storage for calculation area and var0 using calc_nx, calc_ny:
            ALLOCATE(lpm_calc(calc_nx, calc_ny))
            ALLOCATE(var0(calc_nx*calc_ny*n_ens))

            !Calculate PM mean over calculation area:
            !First, if all member values in calculation area are zero, then LPM is zero
            IF(maxval(var2d_ens(ce_west:ce_east, ce_south:ce_north, :)) <= 0 ) THEN
                lpm_calc = 0.0
            ELSE
                var0 = RESHAPE(var2d_ens(ce_west:ce_east,ce_south:ce_north,:), &
                               (/calc_nx*calc_ny*n_ens/))
                !This call calculates PM mean over the calc area and stores it in lpm_calc:
                CALL pm_mean(calc_nx, calc_ny, n_ens, var0, &
                   var2d_enmean(ce_west:ce_east, ce_south:ce_north), lpm_calc)

                !Take the result and store the patch into the array for var2d_lpm:
                !! ChangJae Lee
                !! Catch edge cases:
                var2d_lpm(pe_west:pe_east, pe_south:pe_north) =  &
                    lpm_calc(offset_i:offset_i + add_nx, offset_j:offset_j + add_ny)  
                !!!!!!!!!
            END IF
 
            !Deallocate values for this patch; will re-allocate again for the next one.
            DEALLOCATE(lpm_calc)
            DEALLOCATE(var0)
        END IF
    END DO
END DO

!! ChangJae Lee
!! Congregate results to root processor
IF ( nprocs > 1 ) THEN
    recv_buf=0.
    CALL MPI_REDUCE(var2d_lpm,recv_buf,buf_size, &
                    MPI_REAL,MPI_SUM,root,MPI_COMM_WORLD,ierr)

    IF( myproc == root ) THEN
        var2d_lpm(:,:) = recv_buf(:,:)
    END IF
END IF


IF( myproc == root) THEN
    IF( filt_min > 0.0 ) THEN
    !Filter the LPM so that it is used only above a threshold value
      DO ix=1, nx
        DO jy=1, ny
            !IF(var2d_ens(ix, jy) < filt_min) THEN
            IF(var2d_enmean(ix, jy) < filt_min) THEN
                var2d_lpm(ix, jy) = 0.0
            END IF  
        END DO
      END DO
    END IF 

    ! Apply ensemble max as a ceiling (6/7/2017 - Fanyou Kong)
    var2d_lpm=min(var2d_lpm,var2d_enmax)

    IF ( gauss_sigma > 0 ) THEN
    !Apply a gaussian smoother to the filtered var2d_lpm to get rid of near-grid-scale noise
        CALL gaussian_smooth(nx,ny,var2d_lpm,gauss_sigma,partweight, &
                             lpm_smoothed)
        var2d_lpm=lpm_smoothed
    END IF

    DEALLOCATE(lpm_smoothed)
END IF

END SUBROUTINE lpm_mean_la
!---------------------------------------------------------------!
SUBROUTINE pm_mean(nx,ny,nmember,var0,var2dm,var2d_pm)

use m_mrgrnk

INTEGER :: nx,ny,nmember
INTEGER :: i,j,ij,i0,j0,ii0,ii1,ii2,ii
REAL :: pm0
REAL :: var2dm(nx,ny),var2d_pm(nx,ny)
REAL :: var0(nx*ny*nmember)

REAL, ALLOCATABLE :: var_2d(:)
INTEGER, ALLOCATABLE :: idx_mn(:),idx_2d(:)

ALLOCATE(var_2d(nx*ny))
ALLOCATE(idx_mn(nx*ny*nmember),idx_2d(nx*ny))

call mrgrnk(var0,idx_mn)
do ii=1,nx*ny*nmember
i=nx*ny*nmember - ii + 1
!print *,ii,idx_mn(i),var_2d(idx_mn(i))
ii2=ii-1
if(var0(idx_mn(i))<= 0.0) exit
enddo

do i=1,nx
do j=1,ny
ij=i+nx*(j-1)
var_2d(ij) = var2dm(i,j)
enddo
enddo
call mrgrnk(var_2d,idx_2d)

var2d_pm = 0.0
do ii=1,nx*ny
  i=nx*ny - ii + 1
  j0=(idx_2d(i)-1)/nx+1
  i0=idx_2d(i)-nx*(j0-1)
  ii0 = nx*ny*nmember - nmember*(ii-1)

! skip every mn points
!pm(i0,j0) = var(idx_mn(ii0))

! average over number of members (mn)
  pm0=0.0
  do ll=ii0,ii0-nmember+1,-1
  pm0=pm0+var0(idx_mn(ll))
  enddo
  var2d_pm(i0,j0) = pm0/float(nmember)

  ii1 = ii-1
  if(var_2d(idx_2d(i))<= 0.0) exit
enddo
!print *,'ii1,ii2,ii2/mn:',ii1,ii2,ii2/nmember
!print *,'Total rain points of all members, obs: ',itotle,itobs

  DEALLOCATE(var_2d, idx_mn, idx_2d)

  RETURN
END SUBROUTINE pm_mean
