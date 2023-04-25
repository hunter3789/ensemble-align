!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE MOVEGR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE movegr(nx,ny, var,wrk, xshf,yshf, varout,                    &
           ibgn,iend,jbgn,jend,                                         &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Using the shift vectors, xshf and yshf, translate the variables
!  in array var horizontally.
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!    nx,ny,nz Array dimensions for forecast field.
!
!  OUTPUT :
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny
  REAL, INTENT(IN) :: var(nx,ny)
  REAL, INTENT(OUT) :: wrk(nx,ny)
  REAL, INTENT(IN) :: xshf(nx,ny)
  REAL, INTENT(IN) :: yshf(nx,ny)
  REAL, INTENT(OUT) :: varout(nx,ny)
  INTEGER, INTENT(IN) :: ibgn,iend
  INTEGER, INTENT(IN) :: jbgn,jend
  REAL, INTENT(IN) :: dxfld(nx)
  REAL, INTENT(IN) :: dyfld(ny)
  REAL, INTENT(IN) :: rdxfld(nx)
  REAL, INTENT(IN) :: rdyfld(ny)
  REAL, INTENT(OUT) :: slopey(nx,ny)
  REAL, INTENT(OUT) :: alphay(nx,ny)
  REAL, INTENT(OUT) :: betay(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,idev,jdev,ii,jj
  REAL :: xdev,ydev,delx,dely
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Compute y-derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(nx,ny,1,                                                 &
               ibgn,iend,jbgn,jend,1,1,                                 &
               dyfld,rdyfld,var,                                        &
               slopey,alphay,betay)

  wrk(:,:) = var(:,:)
  DO j=jbgn,jend
    DO i=ibgn,iend
      idev=nint(xshf(i,j)-0.5)
      jdev=nint(yshf(i,j)-0.5)
      xdev=xshf(i,j)-FLOAT(idev)
      ydev=yshf(i,j)-FLOAT(jdev)
!    print *, '  i,j = ',i,j
!    print *, '  xshift,idev,xdev= ',u(i,j),idev,xdev
      jj=j+jdev
      ii=i+idev
!
!  Periodic x boundary conditions
!
!    IF(ii.ge.nx) ii=ii-ixlen
!    IF(ii.lt.1)  ii=ii+ixlen
!
!  Mirror boundary conditions   Mirror at j=3 and j=ny-2
!
!    IF(jj.ge.ny) THEN
!      jj=twnym2-jj
!      ydev=1.-ydev
!    ELSE IF(jj.lt.1) THEN
!      jj=5-jj
!      ydev=1.-ydev
!    END IF
!
!  For now assume zero gradiant boundaries
!
      ii = min(ii,iend-1)
      ii = max(ii,ibgn)
      jj = min(jj,jend)
      jj = max(jj,jbgn)

!    c1=xdev
!    c2=ydev
!    c3=1.-xdev
!    c4=1.-ydev
!    wrk(i,j)=
!    +        c3*(c4*var(  ii,jj)+c2*var(  ii,jj+1))+
!    +        c1*(c4*var(ii+1,jj)+c2*var(ii+1,jj+1))
!
      delx=xdev*dxfld(ii)
      dely=ydev*dyfld(jj)
!     print *, ' i,j,ii,jj:',i,j,ii,jj
      wrk(i,j)=(1.-delx*rdxfld(ii))*                                    &
               (var(ii  ,jj)+slopey(ii  ,jj)*dely)+                     &
               (delx*rdxfld(ii))*                                       &
               (var(ii+1,jj)+slopey(ii+1,jj)*dely)

    END DO
  END DO
!
!  Transfer shifted and original array into output array
!
  varout(:,:) = wrk(:,:)
  RETURN
END SUBROUTINE movegr
!

SUBROUTINE smooth2d(nx,ny,ibgn,iend,jbgn,jend,s,zin,zwork,zout)
!
!  Performs symmetrical two-dimensional smoothing of input field
!  zin which is output as zout, the smoothed field.  A work array
!  zwork, dimension (nx,ny) is required.
!
!  K. Brewster, October, 1991
!
  IMPLICIT NONE
!
!  Arguments
!
  INTEGER, INTENT(IN) :: nx,ny
  INTEGER, INTENT(IN) :: ibgn,iend,jbgn,jend
  REAL, INTENT(IN)    :: s
  REAL, INTENT(IN)    :: zin(nx,ny)
  REAL, INTENT(OUT)   :: zwork(nx,ny)
  REAL, INTENT(OUT)   :: zout(nx,ny)
!
!  Misc internal variables
!
  INTEGER :: i,j
  REAL :: wcen,wsid
!
  zwork(:,:)=zin(:,:)
  wcen=1.-s
  wsid=s*0.5
  DO j=jbgn,jend
    DO i=ibgn+1,iend-1
      zwork(i,j)=zin(i  ,j)*wcen +                                      &
                 zin(i+1,j)*wsid +                                      &
                 zin(i-1,j)*wsid
    END DO
  END DO
!
!
!
  zout(:,:)=zwork(:,:)
  DO j=jbgn+1, jend-1
    DO i=ibgn, iend
      zout(i,j)=zwork(i  ,j)*wcen +                                     &
                zwork(i,j+1)*wsid +                                     &
                zwork(i,j-1)*wsid
    END DO
  END DO
  RETURN
END SUBROUTINE smooth2d

SUBROUTINE smth3dl(nx,ny,nz, ibgn,iend,jbgn,jend,kbgn,kend, zin,zwork,zout,s)
!
!  Performs symmetrical three-dimensional smoothing of input field
!  zin which is output as zout, the smoothed field.  A work array
!  zwork, dimension (nx,ny,nz) is required.
!
!  K. Brewster, October, 1991
!
  IMPLICIT NONE
!
!  Arguments
!
  INTEGER :: nx,ny,nz
  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
  REAL :: zin(nx,ny,nz),zwork(nx,ny,nz),zout(nx,ny,nz)
  REAL :: s
!
!  Misc internal variables
!
  INTEGER :: i,j,k
  REAL :: wcen,wsid
!
  wcen=1.-s
  wsid=s*0.5

  zout = zin
!
  DO k=kbgn,kend
!
!  horizontal part, x direction
!
    DO j=jbgn,jend
      DO i=ibgn+1, iend-1
        zout(i,j,k)=zin(i  ,j,k)*wcen +                                &
                    zin(i+1,j,k)*wsid +                                &
                    zin(i-1,j,k)*wsid
      END DO
    END DO
!
!  horizontal part, y direction
!
    zwork(:,:,k) = zout(:,:,k)
    DO j=jbgn+1, jend-1
      DO i=ibgn, iend
        zwork(i,j,k)=zout(i  ,j,k)*wcen +                              &
                     zout(i,j+1,k)*wsid +                              &
                     zout(i,j-1,k)*wsid
      END DO
    END DO
  END DO
!
! vertical part
!
  zout = zwork
  DO k=kbgn+1, kend-1
    DO j=jbgn, jend
      DO i=ibgn, iend
        zout(i,j,k)=zwork(i  ,j,k)*wcen +                               &
                    zwork(i,j,k+1)*wsid +                               &
                    zwork(i,j,k-1)*wsid
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE smth3dl
