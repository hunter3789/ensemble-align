!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE INTONEF                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE intonef(nx,ny,nz,vnx,vny,vnz,                                &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                         &
           iorder,                                                      &
           x2d,y2d,zp, vx,vy,vzp,                                       &
           vprt, vbar, aprt, abar,                                      &
           iloc,jloc,zpver,dxfld,dyfld,rdxfld,rdyfld,                   &
           slopey,alphay,betay,                                         &
           ireturn )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Intfield interpolates scalars from a single field (the verification
!    fields, "verif") having Cartesian coordinates described by vx,vy,vzp
!    to a second set of fields described by cartesion coordinates x,y,zp.
!    It is assumed that x,y,zp and vx,vy,vzp are monotonically increasing
!    with increasing index.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster  OU School of Meteorology.  Feb 1992
!
!  MODIFICATION HISTORY:
!    12 Aug 1992  (KB) changed from arps2.5 to arps3.0
!    19 May 1993  (KB) changed from arps3.1 to arps3.2
!    24 May 1993  (KB) changed to special version for scalars only.
!
!     9 Sep 1995  (KB) added processing of sfc (soil) fields
!    07 Nov 1996  (KB) changed interpolation scheme
!                      Reordered sequence of variables in call.
!    15 Oct 1999  (KB via Eric Kemp) corrected some dimension
!                      statements (no effect on results)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number grid pts in the x-direction (east/west)
!    ny       Number grid pts in the y-direction (north/south)
!    nz       Number grid pts in the vertical
!
!    vnx      Number grid pts in the x-direction (east/west)
!    vny      Number grid pts in the y-direction (north/south)
!    vnz      Number grid pts in the vertical
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!    kbeg,kend   Range of z index to do interpolation
!
!    ivbeg,ivend   Range of x index to use in field to be interpolated
!    jvbeg,jvend   Range of y index to use in field to be interpolated
!    kvbeg,kvend   Range of z index to use in field to be interpolated
!
!    iorder   Interpolation parameter.
!             iorder specifies the order of interpolation
!             1 = bi-linear
!             2 = bi-quadratic
!
!    x2d      x coordinate (m) of interpolation points
!    y2d      y coordinate (m) of interpolation points
!    z        z coordinate (m) of interpolation points
!
!    vx       x coordinate (m) of field to be interpolated
!    vy       y coordinate (m) of field to be interpolated
!    zpver    z coordinate (m) of field to be interpolated
!
!    vprt     perturbation variable to be interpolated
!    vbar     mean variable to be interpolated
!
!  WORK ARRAYS:
!    iloc     I-index of interpolation points in field to be interpolated
!    jloc     J-index of interpolation points in field to be interpolated
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!  OUTPUT:
!    aprt     perturbation variable interpolated to ARPS grid
!    abar     mean variable interpolated to ARPS grid
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: zp(nx,ny,nz)
!
  INTEGER :: vnx,vny,vnz
  INTEGER :: ivbeg,ivend,jvbeg,jvend,kvbeg,kvend
  REAL :: vx(vnx)
  REAL :: vy(vny)
  REAL :: vzp(vnx,vny,vnz)
!
  INTEGER :: iorder
!
  REAL :: vprt(vnx,vny,vnz)
  REAL :: vbar(vnx,vny,vnz)
  REAL :: aprt(nx,ny,nz)
  REAL :: abar(nx,ny,nz)
!
  INTEGER :: iloc(nx,ny)
  INTEGER :: jloc(nx,ny)
  REAL :: zpver(nx,ny,vnz)
!
  REAL :: dxfld(vnx)
  REAL :: dyfld(vny)
  REAL :: rdxfld(vnx)
  REAL :: rdyfld(vny)
  REAL :: slopey(vnx,vny,vnz)
  REAL :: alphay(vnx,vny,vnz)
  REAL :: betay(vnx,vny,vnz)
!
  INTEGER :: ireturn
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,korder
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Find i,j locations
!
!-----------------------------------------------------------------------
!
  CALL setijloc(nx,ny,vnx,vny,x2d,y2d,vx,vy,iloc,jloc)
!
!-----------------------------------------------------------------------
!
!  Create array of verification heights at
!  forecast x,y locations
!
!-----------------------------------------------------------------------
!
  korder=MIN(iorder,2)
  DO k=1,vnz-1
    CALL fldint2d(nx,ny,vnx,vny,                                        &
                  ibeg,iend,jbeg,jend,                                  &
                  ivbeg,ivend,jvbeg,jvend,                              &
                  korder,x2d,y2d,vzp(1,1,k),vx,vy,iloc,jloc,            &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  slopey(1,1,k),alphay(1,1,k),betay(1,1,k),             &
                  zpver(1,1,k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate 3d fields
!
!-----------------------------------------------------------------------
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,x2d,y2d,zp,vprt,vx,vy,zpver,iloc,jloc,           &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aprt)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,x2d,y2d,zp,vbar,vx,vy,zpver,iloc,jloc,           &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                abar)
  RETURN
END SUBROUTINE intonef
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE FLDINT3D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fldint3d(nx,ny,nz,vnx,vny,vnz,                               &
           ibeg,iend,jbeg,jend,kstart,kend,                             &
           ivstart,ivend,jvstart,jvend,kvstart,kvend,                   &
           iorder,x2d,y2d,z,var,vx,vy,zpver,iloc,jloc,                  &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay,                                         &
           varint)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Intfield interpolates scalars from a set of fields (the verification
!    fields, "verif") having Cartesian coordinates described by vx,vy,vzp
!    to a second set of fields described by cartesion coordinates x,y,zp.
!    It is assumed that x,y,zp and vx,vy,vzp are monotonically increasing
!    with increasing index.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster  OU School of Meteorology.  Feb 1992
!
!  MODIFICATION HISTORY:
!    12 Aug 1992  (KB) changed from arps2.5 to arps3.0
!    19 May 1993  (KB) changed from arps3.1 to arps3.2
!    24 May 1993  (KB) changed to special version for scalars only.
!
!     9 Sep 1995  (KB) added processing of sfc (soil) fields
!    26 Apr 1996  (KB) Version 2.0 -- Uses Gauss Forward routines for
!                      interpolation rather than piecewise linear.
!    07 Nov 1996  (KB) Changed interpolation scheme
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!    nz       Number of model grid points in the vertical
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!    kbeg,kend   Range of z index to do interpolation
!
!    iorder   Interpolation parameter.
!             iorder specifies the order of interpolation
!             1 = bi-linear
!             2 = bi-quadratic
!
!    x2d      x coordinate (m) of interpolation points
!    y2d      y coordinate (m) of interpolation points
!    z        z coordinate (m) of interpolation points
!
!    vx       x coordinate (m) of field to be interpolated
!    vy       y coordinate (m) of field to be interpolated
!    zpver    z coordinate (m) of field to be interpolated
!
!    iloc     I-index of interpolation points in field to be interpolated
!    jloc     J-index of interpolation points in field to be interpolated
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!  WORK ARRAYS:
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!  OUTPUT:
!    varint   3-d array of interpolated values
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: vnx,vny,vnz
  INTEGER :: ibeg,iend,jbeg,jend,kstart,kend
  INTEGER :: ivstart,ivend,jvstart,jvend,kvstart,kvend
  INTEGER :: iorder
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: z(nx,ny,nz)
  REAL :: var(vnx,vny,vnz)
  REAL :: vx(vnx)
  REAL :: vy(vny)
  REAL :: zpver(nx,ny,vnz)
  INTEGER :: iloc(nx,ny)
  INTEGER :: jloc(nx,ny)
  REAL :: dxfld(vnx)
  REAL :: dyfld(vny)
  REAL :: rdxfld(vnx)
  REAL :: rdyfld(vny)
  REAL :: slopey(vnx,vny,vnz)
  REAL :: alphay(vnx,vny,vnz)
  REAL :: betay(vnx,vny,vnz)
  REAL :: varint(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kv,kbot,ktop
  REAL :: wtop,varbot,vartop
  REAL :: pntint2d
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Compute derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(vnx,vny,vnz,                                             &
               ivstart,ivend,jvstart,jvend,kvstart,kvend,               &
               dyfld,rdyfld,var,                                        &
               slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  Compute interpolated values
!
!-----------------------------------------------------------------------
!
  DO k=kstart,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
!
!-----------------------------------------------------------------------
!
!  Find location in second height array
!  Assign linear weight
!
!-----------------------------------------------------------------------
!
        DO kv=kvstart+1,kvend-1
          IF(zpver(i,j,kv) > z(i,j,k)) EXIT
        END DO
!        51    CONTINUE
        ktop=kv
        kbot=kv-1
        wtop=(z(i,j,k)-zpver(i,j,kbot))/                                &
             (zpver(i,j,ktop)-zpver(i,j,kbot))
        varbot=pntint2d(vnx,vny,                                        &
                 ivstart,ivend,jvstart,jvend,                           &
                 iorder,vx,vy,x2d(i,j),y2d(i,j),                        &
                 iloc(i,j),jloc(i,j),var(1,1,kbot),                     &
                 dxfld,dyfld,rdxfld,rdyfld,                             &
                 slopey(1,1,kbot),alphay(1,1,kbot),betay(1,1,kbot))
        vartop=pntint2d(vnx,vny,                                        &
                 ivstart,ivend,jvstart,jvend,                           &
                 iorder,vx,vy,x2d(i,j),y2d(i,j),                        &
                 iloc(i,j),jloc(i,j),var(1,1,ktop),                     &
                 dxfld,dyfld,rdxfld,rdyfld,                             &
                 slopey(1,1,ktop),alphay(1,1,ktop),betay(1,1,ktop))
        varint(i,j,k)=((1.-wtop)*varbot) + (wtop*vartop)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE fldint3d

!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE FLDINT2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fldint2d(nx,ny,vnx,vny,                                      &
           ibeg,iend,jbeg,jend,                                         &
           ivstart,ivend,jvstart,jvend,                                 &
           iorder,x2d,y2d,var,vx,vy,iloc,jloc,                          &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay,                                         &
           varint)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Interpolate a 2d field for several points 2-dimensions.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster  CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!
!    vnx      Number of verif grid points in the x-direction (east/west)
!    vny      Number of verif grid points in the y-direction (north/south)
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!
!    ivbeg,ivend   Range of x index to use in verification array
!    jvbeg,jvend   Range of y index to use in verification array
!
!    iorder   Interpolation parameter.
!             iorder specifies the order of interpolation
!             1 = bi-linear
!             2 = bi-quadratic
!
!    x2d      x coordinate (m) of interpolation points
!    y2d      y coordinate (m) of interpolation points
!
!    var      variable to be interpolated
!
!    vx       x coordinate (m) of field to be interpolated
!    vy       y coordinate (m) of field to be interpolated
!
!  WORK ARRAYS:
!    iloc     I-index of interpolation points in field to be interpolated
!    jloc     J-index of interpolation points in field to be interpolated
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!  OUTPUT
!    varint   Interpolated data array
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny
  INTEGER :: vnx,vny
  INTEGER :: ibeg,iend,jbeg,jend
  INTEGER :: ivstart,ivend,jvstart,jvend
  INTEGER :: iorder
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: var(vnx,vny)
  REAL :: vx(vnx)
  REAL :: vy(vny)
  INTEGER :: iloc(nx,ny)
  INTEGER :: jloc(nx,ny)
  REAL :: dxfld(vnx)
  REAL :: dyfld(vny)
  REAL :: rdxfld(vnx)
  REAL :: rdyfld(vny)
  REAL :: slopey(vnx,vny)
  REAL :: alphay(vnx,vny)
  REAL :: betay(vnx,vny)
  REAL :: varint(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,ii,jj
  REAL :: delx,dely
  REAL :: alpha,beta,rtwodx
  REAL :: varm1,var00,varp1
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
  CALL setdrvy(vnx,vny,1,                                               &
               ivstart,ivend,jvstart,jvend,1,1,                         &
               dyfld,rdyfld,var,                                        &
               slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  Compute bilinear interpolated value
!
!-----------------------------------------------------------------------
!
  IF(iorder == 1) THEN
    DO j=jbeg,jend
      DO i=ibeg,iend
        ii=MIN(MAX(iloc(i,j),ivstart),(ivend-1))
        jj=MIN(MAX(jloc(i,j),jvstart),(jvend-1))
        delx=(x2d(i,j)-vx(ii))
        dely=(y2d(i,j)-vy(jj))
        varint(i,j)=(1.-delx*rdxfld(ii))*                               &
                    (var(ii  ,jj)+slopey(ii  ,jj)*dely)+                &
                    (delx*rdxfld(ii))*                                  &
                    (var(ii+1,jj)+slopey(ii+1,jj)*dely)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Compute biquadratic
!
!-----------------------------------------------------------------------
!
  ELSE
    DO j=jbeg,jend
      DO i=ibeg,iend
        ii=MIN(MAX(iloc(i,j),(ivstart+1)),(ivend-1))
        jj=MIN(MAX(jloc(i,j),(jvstart+1)),(jvend-1))
        delx=(x2d(i,j)-vx(ii))
        dely=(y2d(i,j)-vy(jj))
!
!-----------------------------------------------------------------------
!
!    Stencil is ii-1 to ii+1 and jj-1 to jj + 1
!
!    Interpolate in y.
!
!-----------------------------------------------------------------------
!
        varm1=(alphay(ii-1,jj)*dely+betay(ii-1,jj))*dely+var(ii-1,jj)
        var00=(alphay(ii  ,jj)*dely+betay(ii  ,jj))*dely+var(ii  ,jj)
        varp1=(alphay(ii+1,jj)*dely+betay(ii+1,jj))*dely+var(ii+1,jj)
!
!-----------------------------------------------------------------------
!
!    Interpolate intermediate results in x.
!
!-----------------------------------------------------------------------
!
        rtwodx=1./(dxfld(ii-1)+dxfld(ii))
        alpha=((varp1-var00)*rdxfld(ii  ) +                             &
               (varm1-var00)*rdxfld(ii-1))*rtwodx
        beta=(varp1-var00)*rdxfld(ii) -                                 &
                 dxfld(ii)*alpha
        varint(i,j)=(alpha*delx+beta)*delx+var00
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE fldint2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION PNTINT2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION pntint2d(vnx,vny,                                            &
             ivbeg,ivend,jvbeg,jvend,                                   &
             iorder,vx,vy,xpnt,ypnt,iloc,jloc,var,                      &
             dxfld,dyfld,rdxfld,rdyfld,                                 &
             slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Interpolate a 2-d field for a single point on that plane.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    vnx       Number of model grid points in the x-direction (east/west)
!    vny       Number of model grid points in the y-direction (north/south)
!
!    ivbeg,ivend   Range of x index to use in verification array
!    jvbeg,jvend   Range of y index to use in verification array
!
!    iorder   Interpolation parameter.
!             iorder specifies the order of interpolation
!             1 = bi-linear
!             2 = bi-quadratic
!
!    vx       x coordinate of verif scalar grid points in physical space (m)
!    vy       y coordinate of verif scalar grid points in physical space (m)
!
!    xpnt     x coordinate (m) of interpolation point
!    ypnt     y coordinate (m) of interpolation point
!
!    iloc     I-index of interpolation point in field to be interpolated
!    jloc     J-index of interpolation point in field to be interpolated
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: pntint2d
  INTEGER :: vnx,vny
  INTEGER :: ivbeg,ivend,jvbeg,jvend
  INTEGER :: iorder
  REAL :: vx(vnx)
  REAL :: vy(vny)
  REAL :: xpnt
  REAL :: ypnt
  INTEGER :: iloc
  INTEGER :: jloc
  REAL :: var(vnx,vny)
  REAL :: dxfld(vnx)
  REAL :: dyfld(vny)
  REAL :: rdxfld(vnx)
  REAL :: rdyfld(vny)
  REAL :: slopey(vnx,vny)
  REAL :: alphay(vnx,vny)
  REAL :: betay(vnx,vny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ii,jj
  REAL :: delx,dely
  REAL :: alpha,beta,rtwodx
  REAL :: varm1,var00,varp1
  REAL :: varint
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Compute bilinear interpolated value
!
!-----------------------------------------------------------------------
!
  IF(iorder == 1) THEN
    ii=MIN(MAX(iloc,ivbeg),(ivend-1))
    jj=MIN(MAX(jloc,jvbeg),(jvend-1))
    delx=(xpnt-vx(ii))
    dely=(ypnt-vy(jj))
    varint=(1.-delx*rdxfld(ii))*                                        &
               (var(ii  ,jj)+slopey(ii  ,jj)*dely)+                     &
               (delx*rdxfld(ii))*                                       &
               (var(ii+1,jj)+slopey(ii+1,jj)*dely)
!
!-----------------------------------------------------------------------
!
!  Compute biquadratic
!
!-----------------------------------------------------------------------
!
  ELSE
    ii=MIN(MAX(iloc,(ivbeg+1)),(ivend-1))
    jj=MIN(MAX(jloc,(jvbeg+1)),(jvend-1))
    delx=(xpnt-vx(ii))
    dely=(ypnt-vy(jj))
!
!-----------------------------------------------------------------------
!
!    Stencil is ii-1 to ii+1 and jj-1 to jj + 1
!
!    Interpolate in y.
!
!-----------------------------------------------------------------------
!
    varm1=(alphay(ii-1,jj)*dely+betay(ii-1,jj))*dely+var(ii-1,jj)
    var00=(alphay(ii  ,jj)*dely+betay(ii  ,jj))*dely+var(ii  ,jj)
    varp1=(alphay(ii+1,jj)*dely+betay(ii+1,jj))*dely+var(ii+1,jj)
!
!-----------------------------------------------------------------------
!
!    Interpolate intermediate results in x.
!
!-----------------------------------------------------------------------
!
    rtwodx=1./(dxfld(ii-1)+dxfld(ii))
    alpha=((varp1-var00)*rdxfld(ii  ) +                                 &
           (varm1-var00)*rdxfld(ii-1))*rtwodx
    beta=(varp1-var00)*rdxfld(ii) -                                     &
             dxfld(ii)*alpha
    varint=(alpha*delx+beta)*delx+var00
  END IF
  pntint2d=varint
  RETURN
END FUNCTION pntint2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETIJLOC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setijloc(nx,ny,nx_ext,ny_ext,                                &
           x2d,y2d,x_ext,y_ext,iloc,jloc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Find i,j indices in verfication grid of each forecast point
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster  OU School of Meteorology.  Feb 1992
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!
!    nx_ext   Number of external grid pts in the x-direction (east/west)
!    ny_ext   Number of external grid pts in the y-direction (north/south)
!
!    x2d      x coordinate (m) of interpolation points
!    y2d      y coordinate (m) of interpolation points
!
!    x_ext    x-coordinate (m) of external grid pts
!    y_ext    y-coordinate (m) of external grid pts
!
!  OUTPUT:
!
!    iloc     i-index of interpolation points in field to be interpolated
!    jloc     j-index of interpolation points in field to be interpolated
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nx_ext,ny_ext
  REAL :: x2d(nx,ny)
  REAL :: y2d(nx,ny)
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)
  INTEGER :: iloc(nx,ny)
  INTEGER :: jloc(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,i_ext,j_ext
  INTEGER :: imid,jmid
  REAL :: xmid,ymid
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  imid=nx_ext/2
  xmid=x_ext(imid)
  jmid=ny_ext/2
  ymid=y_ext(jmid)
!
  DO j=1,ny
    DO i=1,nx
      IF(x2d(i,j) < xmid) THEN
        DO i_ext=imid,2,-1
          IF(x_ext(i_ext) <= x2d(i,j)) EXIT
        END DO
        iloc(i,j)=i_ext
      ELSE
        DO i_ext=imid,nx_ext-1
          IF(x_ext(i_ext) >= x2d(i,j)) EXIT
        END DO
        iloc(i,j)=i_ext-1
      END IF
!
      IF(y2d(i,j) < ymid) THEN
        DO j_ext=jmid,2,-1
          IF(y_ext(j_ext) <= y2d(i,j)) EXIT
        END DO
        jloc(i,j)=j_ext
      ELSE
        DO j_ext=jmid,ny_ext-1
          IF(y_ext(j_ext) >= y2d(i,j)) EXIT
        END DO
        jloc(i,j)=j_ext-1
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE setijloc
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETDXDY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setdxdy(nx,ny,                                               &
           ibeg,iend,jbeg,jend,                                         &
           x1d,y1d,dxfld,dyfld,rdxfld,rdyfld)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the local delta-x, delta-y and their inverses.
!    Precalculating these variables speeds up later calculations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!
!    x1d     Array of x-coordinate grid locations (m)
!    y1d     Array of y-coordinate grid locations (m)
!
!  OUTPUT:
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny
  INTEGER :: ibeg,iend
  INTEGER :: jbeg,jend
  REAL :: x1d(nx)
  REAL :: y1d(ny)
  REAL :: dxfld(nx)
  REAL :: dyfld(ny)
  REAL :: rdxfld(nx)
  REAL :: rdyfld(ny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,istop,jstop
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istop=MIN((iend-1),(nx-1))
  DO i=ibeg,istop
    dxfld(i)=(x1d(i+1)-x1d(i))
    rdxfld(i)=1./(x1d(i+1)-x1d(i))
  END DO
  jstop=MIN((jend-1),(ny-1))
  DO j=jbeg,jstop
    dyfld(j)=(y1d(j+1)-y1d(j))
    rdyfld(j)=1./(y1d(j+1)-y1d(j))
  END DO
  RETURN
END SUBROUTINE setdxdy
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETDRVY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setdrvy(nx,ny,nz,                                            &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           dyfld,rdyfld,var,                                            &
           slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the coefficients of interpolating polynomials
!    in the y-direction.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!    nz       Number of model grid points in the vertical
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!    kbeg,kend   Range of z index to do interpolation
!
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    var      variable to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL :: dyfld(ny)
  REAL :: rdyfld(ny)
  REAL :: var(nx,ny,nz)
  REAL :: slopey(nx,ny,nz)
  REAL :: alphay(nx,ny,nz)
  REAL :: betay(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: jstart,jstop
  REAL :: rtwody
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  jstart=MAX(jbeg,2)
  jstop=MIN((jend-1),(ny-2))
  DO k=kbeg,kend
    DO j=jstart,jstop
      DO i=ibeg,iend
        slopey(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j)
        rtwody=1./(dyfld(j-1)+dyfld(j))
        alphay(i,j,k)=((var(i,j+1,k)-var(i,j,k))*rdyfld(j) +            &
                 (var(i,j-1,k)-var(i,j,k))*rdyfld(j-1))*rtwody
        betay(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j) -              &
                   dyfld(j)*alphay(i,j,k)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE setdrvy
