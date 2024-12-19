!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
!
!  This is the "interp.F" script (based on the interpolation implemented in the
!  Roms Rutgers version)
!------------------------------------------------------------------------------
!  This file contains the subfunctions enabling the online interpolation of
!  forcing datasets on the simulation domain using linear or cubic approach.
!  These functions applied for all discretisations of the domain, MPI or OPENMP.
!------------------------------------------------------------------------------
#include "cppdefs.h"

      SUBROUTINE linterp2d (LBx, UBx, LBy, UBy,
     &                      Xinp, Yinp, Finp,
     &                      Istr, Iend, Jstr, Jend,
     &                      Xout, Yout, Fout)
!
!=======================================================================
!                                                                      !
!  Given any gridded 2D field, Finp, this routine linearly interpolate !
!  to locations (Xout,Yout).  To facilitate the  interpolation  within !
!  any irregularly gridded 2D field,  the fractional grid cell indices !
!  (Iout,Jout) with respect Finp are needed at input.  Notice that the !
!  routine "hindices" can be used to compute these indices.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     LBx        I-dimension lower bound of gridded field, Finp.       !
!     UBx        I-dimension upper bound of gridded field, Finp.       !
!     LBy        J-dimension lower bound of gridded field, Finp.       !
!     UBy        J-dimension upper bound of gridded field, Finp.       !
!     Xinp       X-locations of gridded field, Finp.                   !
!     Yinp       Y-locations of gridded field, Finp.                   !
!     Finp       2D field to interpolate from.                         !
!     Istr       Starting data I-index to interpolate, Fout.           !
!     Iend       Ending   data I-index to interpolate, Fout.           !
!     Jstr       Starting data J-index to interpolate, Fout.           !
!     Jend       Ending   data J-index to interpolate, Fout.           !
!     Xout       X-locations to interpolate, Fout.                     !
!     Yout       Y-locations to interpolate, Fout.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated 2D field.                                !
!                                                                      !
!=======================================================================
!
!
      implicit none
#include "param.h"
!
!  Imported variable declarations.
!
      integer, intent(in) :: LBx, UBx, LBy, UBy
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(kind=8), intent(in) :: Xinp(LBx:UBx)
      real(kind=8), intent(in) :: Yinp(LBy:UBy)
      real(kind=8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
!
      real(kind=8), intent(in) :: Xout(GLOBAL_2D_ARRAY)
      real(kind=8), intent(in) :: Yout(GLOBAL_2D_ARRAY)
!
      real(kind=8), intent(out) :: Fout(GLOBAL_2D_ARRAY)
!
!  Local variable declarations.
!
      integer i, i1, i2, j, j1, j2, ii, jj
      real(kind=8) cff, x, x1, x2, y, y1, y2
!
!-----------------------------------------------------------------------
!  Linearly interpolate requested field
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!          i1=INT(Iout(i,j))
!          i2=i1+1
!          j1=INT(Jout(i,j))
!          j2=j1+1
           DO ii=LBx,(UBx-1)
             if ((Xinp(ii).le.Xout(i,j)).and.
     &           (Xinp(ii+1).gt.Xout(i,j))) then
               i1=ii
               i2=ii+1
               goto 10
             endif
           enddo
           print*, 'Did not find i1 and i2'
           goto 100
 10        continue
           DO jj=LBy,(UBy-1)
             if ((Yinp(jj).le.Yout(i,j)).and.
     &           (Yinp(jj+1).gt.Yout(i,j))) then
               j1=jj
               j2=jj+1
               goto 20
             endif
           enddo
           print*, 'Did not find j1 and j2'
           goto 100
 20        continue

          IF (((LBx.le.i1).and.(i1.le.UBx)).and.
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
            x1=Xinp(i1)
            x2=Xinp(i2)
            y1=Yinp(j1)
            y2=Yinp(j2)
            x=Xout(i,j)
            y=Yout(i,j)

            cff= Finp(i1,j1)*(x2-x )*(y2-y )
     &          +Finp(i2,j1)*(x -x1)*(y2-y )
     &          +Finp(i1,j2)*(x2-x )*(y -y1)
     &          +Finp(i2,j2)*(x -x1)*(y -y1)

            Fout(i,j)=cff/((x2-x1)*(y2-y1))
          END IF
        END DO
      END DO
      RETURN
 100  continue
      print*, 'error in linterp2d'
      END SUBROUTINE linterp2d

! ******************************************************************************
      SUBROUTINE cinterp2d (LBx, UBx, LBy, UBy,
     &                      Xinp, Yinp, Finp,
     &                      Istr, Iend, Jstr, Jend,
     &                      Xout, Yout, Fout)
!
!=======================================================================
!                                                                      !
!  Given any gridded 2D field,  Finp, at locations (Xinp,Yinp) this    !
!  routine performs bicubic interpolation at locations (Xout,Yout).    !
!  To facilitate the interpolation within any  irregularly  gridded    !
!  field, the fractional grid cell indices (Iout,Jout) with respect    !
!  Finp are needed at input. Notice that the routine "hindices" can    !
!  be used to compute these indices.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     LBx        I-dimension lower bound of gridded field, Finp.       !
!     UBx        I-dimension upper bound of gridded field, Finp.       !
!     LBy        J-dimension lower bound of gridded field, Finp.       !
!     UBy        J-dimension upper bound of gridded field, Finp.       !
!     Xinp       X-locations of gridded field, Finp.                   !
!     Yinp       Y-locations of gridded field, Finp.                   !
!     Finp       2D field to interpolate from.                         !
!     Istr       Starting data I-index to interpolate, Fout.           !
!     Iend       Ending   data I-index to interpolate, Fout.           !
!     Jstr       Starting data J-index to interpolate, Fout.           !
!     Jend       Ending   data J-index to interpolate, Fout.           !
!     Xout       X-locations to interpolate, Fout.                     !
!     Yout       Y-locations to interpolate, Fout.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated 2D field.                                !
!                                                                      !
!=======================================================================
!
      implicit none
#include "param.h"
#include "grid.h"
!
!  Imported variable declarations.
!
      integer, intent(in) :: LBx, UBx, LBy, UBy
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(kind=8), intent(in) :: Xinp(LBx:UBx)
      real(kind=8), intent(in) :: Yinp(LBy:UBy)
      real(kind=8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
!
      real(kind=8), intent(in) :: Xout(GLOBAL_2D_ARRAY)
      real(kind=8), intent(in) :: Yout(GLOBAL_2D_ARRAY)
!
      real(kind=8), intent(out) :: Fout(GLOBAL_2D_ARRAY)
!
!  Local variable declarations.
!
      integer i, ic, iter, i1, i2, j, jc, j1, j2, ii, jj

      real(kind=8) :: a11, a12, a21, a22
      real(kind=8) :: e11, e12, e21, e22
      real(kind=8) :: cff, d1, d2, dfc, dx, dy, eta, xi, xy, yx
      real(kind=8) :: f0, fx, fxx, fxxx, fxxy, fxy, fxyy, fy, fyy, fyyy

      real(kind=8), parameter :: C01 = 1.0/48.0
      real(kind=8), parameter :: C02 = 1.0/32.0
      real(kind=8), parameter :: C03 = 0.0625                  ! 1/16
      real(kind=8), parameter :: C04 = 1.0/6.0
      real(kind=8), parameter :: C05 = 0.25
      real(kind=8), parameter :: C06 = 0.5
      real(kind=8), parameter :: C07 = 0.3125                  ! 5/16
      real(kind=8), parameter :: C08 = 0.625                   ! 5/8
      real(kind=8), parameter :: C09 = 1.5
      real(kind=8), parameter :: C10 = 13.0/24.0

      real(kind=8), parameter :: LIMTR = 3.0
      real(kind=8), parameter :: spv = 0.0            ! HGA need work

      real(kind=8), dimension(-1:2,-1:2) :: dfx, dfy, ff
!
!-----------------------------------------------------------------------
!  Interpolates requested field locations (Xout,Yout).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!          i1=INT(Iout(i,j))
!          i2=i1+1
!          j1=INT(Jout(i,j))
!          j2=j1+1
           DO ii=LBx,(UBx-1)
             if ((Xinp(ii).le.Xout(i,j)).and.
     &           (Xinp(ii+1).gt.Xout(i,j))) then
               i1=ii
               i2=ii+1
               goto 10
             endif
           enddo
           print*, 'Did not find i1 and i2',
     &           Istr,Iend,Jstr,Jend,i,j,Xout(i,j),Xout(i-1,j)
           goto 100
 10        continue
           DO jj=LBy,UBy-1
             if ((Yinp(jj).le.Yout(i,j)).and.
     &           (Yinp(jj+1).gt.Yout(i,j))) then
               j1=jj
               j2=jj+1
               goto 20
             endif
           enddo
           print*, 'Did not find j1 and j2'
           goto 100
 20        continue

          IF (((LBx.le.i1).and.(i1.le.UBx)).and.
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
!
!  Determine local fractional coordinates (xi,eta) corresponding to
!  the target point (Xout,Yout) on the grid (Xinp,Yinp). Here, "xi"
!  and "eta" are defined, in such a way, that xi=eta=0 corresponds
!  to the middle of the cell (i1:i1+1,j1:j1+1), while xi=+/-1/2 and
!  eta=+/-1/2 (any combination +/- signs) corresponds to the four
!  corner points of the cell. Inside the cell it is assumed that
!  (Xout,Yout) are expressed via bi-linear functions of (xi,eta),
!  where term proportional to xi*eta does not vanish because
!  coordinate transformation may be at least weakly non-orthogonal
!  due to discretization errors. The associated non-linear system
!  is solved by iterative method of Newton.
!
            xy=Xinp(i2)-Xinp(i1)-Xinp(i2)+Xinp(i1)
            yx=Yinp(j2)-Yinp(j2)-Yinp(j1)+Yinp(j1)
            dx=Xout(i,j)-0.25*(Xinp(i2)+Xinp(i1)+
     &                         Xinp(i2)+Xinp(i1))
            dy=Yout(i,j)-0.25*(Yinp(j2)+Yinp(j2)+
     &                         Yinp(j1)+Yinp(j1))
!
!  The coordinate transformation matrix:
!
!           e11 e12
!           e21 e22
!
!  contains derivatives of (Xinp,Yinp) with respect to (xi,eta). Because
!  the coordinates may be non-orthogonal (at least due to discretization
!  errors), the nonlinear system
!
!           e11*xi+e12*eta+xy*xi*eta=dx
!           e21*xi+e22*eta+yx*xi*eta=dy
!
!  needs to be solved in order to retain symmetry.
!
            e11=0.5*(Xinp(i2)-Xinp(i1)+Xinp(i2)-Xinp(i1))
            e12=0.5*(Xinp(i2)+Xinp(i1)-Xinp(i2)-Xinp(i1))
            e21=0.5*(Yinp(j2)-Yinp(j2)+Yinp(j1)-Yinp(j1))
            e22=0.5*(Yinp(j2)+Yinp(j2)-Yinp(j1)-Yinp(j1))
!
            cff=1.0/(e11*e22-e12*e21)
            xi=cff*(e22*dx-e12*dy)
            eta=cff*(e11*dy-e21*dx)
!
            DO iter=1,4
              d1=dx-e11*xi-e12*eta-xy*xi*eta
              d2=dy-e21*xi-e22*eta-yx*xi*eta
              a11=e11+xy*eta
              a12=e12+xy*xi
              a21=e21+yx*eta
              a22=e22+yx*xi
              cff=1.0/(a11*a22-a12*a21)
              xi =xi +cff*(a22*d1-a12*d2)
              eta=eta+cff*(a11*d2-a21*d1)
            END DO

#ifndef CUBIC_MASKED
!
!  Genuinely two-dimensional, isotropic cubic interpolation scheme
!  using 12-point stencil.  In the code below the interpolated field,
!  Fout, is expanded into two-dimensional Taylor series of local
!  fractional coordinates "xi" and "eta", retaining all terms of
!  combined power up to third order (that is, xi, eta, xi^2, eta^2,
!  xi*eta, xi^3, eta^3, xi^2*eta, and xi*eta^2), with all
!  coefficients (i.e, derivatives) computed via           x  x
!  two-dimensional finite difference expressions          |  |
!  of "natural" order of accuracy: 4th-order for       x--x--x--x
!  the field itself and its first derivatives in          |  |
!  both directions; and 2nd-order for all higher-      x--x--x--x
!  order derivatives. The permissible range of            |  |
!  of coordinates is -1/2 < xi,eta < +1/2, which          x--x
!  covers the central cell on the stencil, while
!  xi=eta=0 corresponds to its center. This interpolation scheme has
!  the property that if xi,eta=+/-1/2 (any combination of +/- signs)
!  it reproduces exactly value of the function at the corresponding
!  corner of the central "working" cell. However, it does not pass
!  exactly through the  extreme points of the stencil, where either
!  xi=+/-3/2 or eta+/-3/2. And, unlike a split-directional scheme,
!  when interpolating along the line eta=+/-1/2 (similarly xi=+/-1/2),
!  it has non-zero contribution from points on the side from the line,
!  except if xi=-1/2; 0; +1/2 (similarly eta=-1/2; 0; +1/2).
!
            DO jc=-1,2
              DO ic=-1,2
                ff(ic,jc)=Finp(MAX(1,MIN(UBx,i1+ic)),
     &                         MAX(1,MIN(UBy,j1+jc)))
              END DO
            END DO

            f0=C07*(ff(1,1)+ff(1,0)+ff(0,1)+ff(0,0))-
     &         C02*(ff(2,0)+ff(2,1)+ff(1,2)+ff(0,2)+
     &              ff(-1,1)+ff(-1,0)+ff(0,-1)+ff(1,-1))

            fx=C08*(ff(1,1)+ff(1,0)-ff(0,1)-ff(0,0))-
     &         C01*(ff(2,1)+ff(2,0)-ff(-1,1)-ff(-1,0))-
     &         C03*(ff(1,2)-ff(0,2)+ff(1,-1)-ff(0,-1))

            fy=C08*(ff(1,1)-ff(1,0)+ff(0,1)-ff(0,0))-
     &         C01*(ff(1,2)+ff(0,2)-ff(1,-1)-ff(0,-1))-
     &         C03*(ff(2,1)-ff(2,0)+ff(-1,1)-ff(-1,0))

            fxy=ff(1,1)-ff(1,0)-ff(0,1)+ff(0,0)

            fxx=C05*(ff(2,1)-ff(1,1)-ff(0,1)+ff(-1,1)+
     &               ff(2,0)-ff(1,0)-ff(0,0)+ff(-1,0))

            fyy=C05*(ff(1,2)-ff(1,1)-ff(1,0)+ff(1,-1)+
     &               ff(0,2)-ff(0,1)-ff(0,0)+ff(0,-1))

            fxxx=C06*(ff(2,1)+ff(2,0)-ff(-1,1)-ff(-1,0))-
     &           C09*(ff(1,1)+ff(1,0)-ff(0,1)-ff(0,0))

            fyyy=C06*(ff(1,2)+ff(0,2)-ff(1,-1)-ff(0,-1))-
     &           C09*(ff(1,1)-ff(1,0)+ff(0,1)-ff(0,0))

            fxxy=C06*(ff(2,1)-ff(1,1)-ff(0,1)+ff(-1,1)-
     &                ff(2,0)+ff(1,0)+ff(0,0)-ff(-1,0))

            fxyy=C06*(ff(1,2)-ff(1,1)-ff(1,0)+ff(1,-1)-
     &                ff(0,2)+ff(0,1)+ff(0,0)-ff(0,-1))
#else
!
!  Algorithm below is equivalent to the one above, except that special
!  care is taken to avoid interpolation accross land. This is achieved
!  by shortening the stencil and reducing order of polynomial, if
!  extreme points of the stencil touch land. This is achieved by
!  expressing all f0,fx,fy,...,fxyy in terms of values of interpolated
!  field at the four corners of central cell (which already checked to
!  stay away from land), and eight one-sided differences dfx,dfy (see
!  below) in such a way that field values at the extreme points of the
!  12-point stencil do not participate directly into f0,fx,...,fxyy.
!  Should an extreme point of the stencil touch land, thus making it
!  impossible to compute the corresponding one-sided difference, this
!  difference is retracted toward the center of the stencil.
!
!  Optionally, a slope-limiting algorithm may be employed to prevent
!  spurious oscillations of the interpolant. This is a valuable property,
!  if dealing with rough data, however, as a side effect, it turns off
!  high-order interpolation in the vicinity of extrema.
!
!  The slope-limiting algorithm employed here checks that two consecutive
!  elementary differences, "dfx" and "dfc" have the same sign and differ
!  in magnitude by no more than factor of 3.
!
            ff(0,0)=Finp(i1,j1)
            ff(1,0)=Finp(i2,j1)
            ff(0,1)=Finp(i1,j2)
            ff(1,1)=Finp(i2,j2)
!
            dfc=ff(1,1)-ff(0,1)
            IF (i1+2.gt.UBx) THEN
              dfx(1,1)=dfc
            ELSE IF (Finp(i1+2,j2).eq.spv) THEN
              dfx(1,1)=dfc
            ELSE
              dfx(1,1)=Finp(i1+2,j2)-ff(1,1)
# ifdef LIMTR
              IF ((dfx(1,1)*dfc).lt.0.0) THEN
                dfx(1,1)=0.0
              ELSE IF (ABS(dfx(1,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(1,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=ff(1,0)-ff(0,0)
            IF ((i1+2).gt.UBx) THEN
              dfx(1,0)=dfc
            ELSE IF (Finp(i1+2,j1).eq.spv) THEN
              dfx(1,0)=dfc
            ELSE
              dfx(1,0)=Finp(i1+2,j1)-ff(1,0)
# ifdef LIMTR
              IF ((dfx(1,0)*dfc).lt.0.0) THEN
                dfx(1,0)=0.0
              ELSE IF (ABS(dfx(1,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(1,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=ff(1,1)-ff(0,1)
            IF (i1-1.lt.1) THEN
              dfx(0,1)=dfc
            ELSE IF (Finp(i1-1,j2).eq.spv) THEN
              dfx(0,1)=dfc
            ELSE
              dfx(0,1)=ff(0,1)-Finp(i1-1,j2)
# ifdef LIMTR
              IF ((dfx(0,1)*dfc).lt.0.0) THEN
                dfx(0,1)=0.0
              ELSE IF (ABS(dfx(0,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(0,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=ff(1,0)-ff(0,0)
            IF (i1-1.lt.1) THEN
              dfx(0,0)=dfc
            ELSE IF (Finp(i1-1,j1).eq.spv) THEN
              dfx(0,0)=dfc
            ELSE
              dfx(0,0)=ff(0,0)-Finp(i1-1,j1)
# ifdef LIMTR
              IF ((dfx(0,0)*dfc).lt.0.0) THEN
                dfx(0,0)=0.0
              ELSE IF (ABS(dfx(0,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(0,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=ff(1,1)-ff(1,0)
            IF (j1+2.gt.UBy) THEN
              dfy(1,1)=dfc
            ELSE IF (Finp(i2,j1+2).eq.spv) THEN
              dfy(1,1)=dfc
            ELSE
              dfy(1,1)=Finp(i2,j1+2)-ff(1,1)
# ifdef LIMTR
              IF ((dfy(1,1)*dfc).lt.0.0) THEN
                dfy(1,1)=0.0
              ELSEIF (ABS(dfy(1,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(1,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=ff(0,1)-ff(0,0)
            IF (j1+2.gt.UBy) THEN
              dfy(0,1)=dfc
            ELSE IF (Finp(i1,j1+2).eq.spv) THEN
              dfy(0,1)=dfc
            ELSE
              dfy(0,1)=Finp(i1,j1+2)-ff(0,1)
# ifdef LIMTR
              IF ((dfy(0,1)*dfc).lt.0.0) THEN
                dfy(0,1)=0.0
              ELSE IF (ABS(dfy(0,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(0,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=ff(1,1)-ff(1,0)
            IF (j1-1.lt.1) THEN
              dfy(1,0)=dfc
            ELSE IF (Finp(i2,j1-1).eq.spv) THEN
              dfy(1,0)=dfc
            ELSE
              dfy(1,0)=ff(1,0)-Finp(i2,j1-1)
# ifdef LIMTR
              IF ((dfy(1,0)*dfc).lt.0.0) THEN
                dfy(1,0)=0.0
              ELSE IF (ABS(dfy(1,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(1,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=ff(0,1)-ff(0,0)
            IF (j1-1.lt.1) THEN
              dfy(0,0)=dfc
            ELSE IF (Finp(i1,j1-1).eq.spv) THEN
              dfy(0,0)=dfc
            ELSE
              dfy(0,0)=f(0,0)-Finp(i1,j1-1)
# ifdef LIMTR
              IF ((dfy(0,0)*dfc).lt.0.0) THEN
                dfy(0,0)=0.0
              ELSEIF (ABS(dfy(0,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(0,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            f0=C05*(ff(1,1)+ff(1,0)+ff(0,1)+ff(0,0))-
     &         C02*(dfx(1,1)+dfx(1,0)-dfx(0,1)-dfx(0,0)+
     &              dfy(1,1)-dfy(1,0)+dfy(0,1)-dfy(0,0))

            fx=C10*(ff(1,1)-ff(0,1)+ff(1,0)-ff(0,0))-
     &         C01*(dfx(1,1)+dfx(1,0)+dfx(0,1)+dfx(0,0))-
     &         C03*(dfy(1,1)-dfy(0,1)-dfy(1,0)+dfy(0,0))

            fy=C10*(ff(1,1)-ff(1,0)+ff(0,1)-ff(0,0))-
     &         C01*(dfy(1,1)+dfy(0,1)+dfy(1,0)+dfy(0,0))-
     &         C03*(dfx(1,1)-dfx(1,0)-dfx(0,1)+dfx(0,0))

            fxy=ff(1,1)-ff(1,0)-ff(0,1)+ff(0,0)

            fxx=C05*(dfx(1,1)-dfx(0,1)+dfx(1,0)-dfx(0,0))

            fyy=C05*(dfy(1,1)-dfy(1,0)+dfy(0,1)-dfy(0,0))

            fxxx=C06*(dfx(1,1)+dfx(1,0)+dfx(0,1)+dfx(0,0))-
     &           ff(1,1)+ff(0,1)-ff(1,0)+ff(0,0)

            fyyy=C06*(dfy(1,1)+dfy(0,1)+dfy(1,0)+dfy(0,0))-
     &           ff(1,1)+ff(1,0)-ff(0,1)+ff(0,0)

            fxxy=C06*(dfx(1,1)-dfx(0,1)-dfx(1,0)+dfx(0,0))

            fxyy=C06*(dfy(1,1)-dfy(1,0)-dfy(0,1)+dfy(0,0))
#endif
            Fout(i,j)=f0+
     &                fx*xi+
     &                fy*eta+
     &                C06*fxx*xi*xi+
     &                fxy*xi*eta+
     &                C06*fyy*eta*eta+
     &                C04*fxxx*xi*xi*xi+
     &                C06*fxxy*xi*xi*eta+
     &                C04*fyyy*eta*eta*eta+
     &                C06*fxyy*xi*eta*eta
          END IF
        END DO
      END DO

      RETURN

 100  continue
      print*, 'error in cinterp2d'
      END SUBROUTINE cinterp2d

! ******************************************************************************
! FUNCTION TO INTERPOLATE IRREGULAR GRID ON CROCO_GRID
!
      SUBROUTINE irrinterp2d (LBx, UBx, LBy, UBy,
     &                      Xinp, Yinp, Finp,
     &                      Istr, Iend, Jstr, Jend,
     &                      Xout, Yout, Fout)

      implicit none
#include "param.h"
!
!  Imported variable declarations.
!
      integer, intent(in) :: LBx, UBx, LBy, UBy
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(kind=8), intent(in) :: Xinp(LBx:UBx,LBy:UBy)
      real(kind=8), intent(in) :: Yinp(LBx:UBx,LBy:UBy)
      real(kind=8), intent(in) :: Finp(LBx:UBx,LBy:UBy)
!
      real(kind=8), intent(in) :: Xout(GLOBAL_2D_ARRAY)
      real(kind=8), intent(in) :: Yout(GLOBAL_2D_ARRAY)
!
      real(kind=8), intent(out) :: Fout(GLOBAL_2D_ARRAY)

      real(kind=8), dimension(Istr:Iend,Jstr:Jend) :: xi,eta
      integer(kind=4), dimension(Istr:Iend,Jstr:Jend) :: ip,jp
      real(kind=8), dimension(LBx:UBx,LBy:UBy) :: srX,srY,sXY,sYX
!
!!!!!!!!!!!!
! ORIGINAL !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call r2r_interp_init_thread( nx,ny,x,y, ncx,ncy, xc,yc,
!     &                                                ip,jp, xi,eta)
!
!      call spln2d(nx,ny, src, srX,srY,sXY,sYX)
!
!      call spln2d_interp(nx,ny, src, srX, srY, sXY,
!     &                           ncx,ncy, ip,jp, xi,eta, targ) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call spln2d(UBx,UBy, Finp, srX,srY,sXY,sYX)


      call r2r_interp_init_thread(UBx,UBy,Xinp,Yinp,Istr,Iend,
     &                             Jstr,Jend, Xout,Yout ,ip,jp, xi,eta)



      call spln2d_interp(UBx,UBy, Finp, srX, srY, sXY,
     &                      Istr, Iend, Jstr, Jend, ip,jp, xi,eta, Fout)

            

      RETURN
      END SUBROUTINE irrinterp2d


      SUBROUTINE r2r_interp_init_thread( nx,ny,x,y,Istr, ncx,Jstr,ncy,
     &                                            xc,yc, ip,jp, xi,eta)

! Initialize interpolation between two arbitrary oriented curvilinear
! grids: given arrays of coordinates [x(nx,ny),y(nx,ny)] for points of
! a non-Cartesian "source" grid, and similar [xc(ncx,ncy),yc(ncx,ncy)]
! for the "target" grid (hereafter "c" stands for "child") find arrays
! of indices [ip(ncx,ncy),jp(ncx,ncy)] and arrays fractional distances
! [xi(ncx,ncy),eta(ncx,ncy)] such that each point of the target grid
! (xc,yc) is surrounded by 4 points of the source,
!
!              [x,y](ip,jp+1) --__
!                     /             --__ [x,y](ip+1,jp+1)
!                    /   [xc,yc]            /
!                   /                      /
!              [x,y](ip,jp) -- __         /
!                               --__ [x,y](ip+1,jp)
!
! such that bi-linear interpolation of [x,y] into location of [xc,yc]
! yields [xc,yc] themselves,
!
!             xc = (1-xi)*(1-eta)*x(i,j) + xi*(1-eta)*x(i+1,j)
!                 +(1-xi)* eta *x(i,j+1) + xi* eta *x(i+1,j+1)
!
!             yc = (1-xi)*(1-eta)*y(i,j) + xi*(1-eta)*y(i+1,j)
!                 +(1-xi)* eta *y(i,j+1) + xi* eta *y(i+1,j+1)
!
! in other words, (ip+xi,jp+eta) are coordinates of point (xc,yc) in
! "continuous index" space (i,j) of the source grid.

! Normally the target grid is expected to be entirely within the area
! covered by the source, however it is not strictly required by the
! algorithm below: if some portions of the target grid are outside the
! source (hence interpolation is impossible for these locations) then
! ip=-1 and jp=-1 are set as special values.

      implicit none
      integer nx,ny, ncx,ncy,Istr,Jstr
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer icmin,icmax,jcmin,jcmax
C$    integer numthreads, trd, chunk_size
C$    integer omp_get_num_threads, omp_get_thread_num
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$    chunk_size=(ncy+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx ; jcmin=1 ; jcmax=ncy
C$    jcmin=1+trd*chunk_size ; jcmax=min(jcmin+chunk_size-1,ncy)

      call search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     &                                       ncx,ncy, xc,yc, ip,jp)

      call comp_offsets_tile(  icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     &                               ncx,ncy, xc,yc, ip,jp, xi,eta)
      END SUBROUTINE


      SUBROUTINE search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny,
     &                                  x,y, ncx,ncy, xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc         !   4 <---- 3
      integer(kind=4), dimension(ncx,ncy) :: ip,jp      !   !       ^ 
      integer i,j, ic,jc, inew,jnew, iinc               !   !       !
      real(kind=8) dx4,dy4, dx3,dy3, r41,r34,           !   v       !
     &             dx1,dy1, dx2,dy2, r12,r23            !   1 ----> 2

! The search algorithm is organized as follows:
!
!    initialization: select starting point ic,jc
!    and compute cyclic vector products r12...r41
!
!    do while(.true.)
!      do while( r12...r41 >= 0 -- meaning that point
!                 ic,jc is inside of [i:i+1]x[j:j+1] )
!
!        record ip,jp for this point, and proceed to the next point
!        of target grid using by incrementing only ONE index, either
!        "ic" or "jc" -- if "ic" already reached its bound, while
!        simultaneously reversing the direction of "ic" sweep -- the
!        reversing is done to make sure that the next target point is
!        always nearby, to maximize the probability that it is still
!        inside of [i:i+1]x[j:j+1] cell, so the inner while-loop would
!        not break off too often. Recompute r12...r41 for the next
!        ic,jc point to check logical condition during next iteration
!        of while loop; [Note that there is no attempt to change
!        source-grid indices i,j inside this loop.]
!
!      enddo
!
!      once the while(r12...r41 >= 0) loop breaks off, it means that
!      the point ic,jc is no longer inside [i:i+1]x[j:j+1], so i and/or
!      j index of source grid must be incremented, depending which of
!      the r12...r41 is negative.  The indices i,j are incremented by
!      one and are restricted to be within the range of source grid.
!      This leads to two possibilities:
!      either
!             (i) at least one of them, i,j is incremented, then keep
!                 ic,jc the same and proceed with recomputing r12...r41
!                 to resume while(r12...r41 >= 0) loop;
!      or
!            (ii) both increments of i and j are canceled by the
!                 restriction, so i,j get "stuck"  which is detected
!                 by inew,jnew both having the same values. This means
!                 that this ic,jc point cannot be bounded properly
!                 because its location is outside the source grid, so
!                 mark it by special values of ip,jp (essentially skip
!                 this point) and proceed to the next one by
!                 incrementing  ip,jp;
!
!    enddo !<-- while(.true.)
!
! Note (a) that not having special care for (i) leads to an infinite
!          loop in the algorithm if some of the target points cannot be
!          bounded; and
!      (b) assuming that the target grid has finer resolution than the
!          source, it is expected that most of the computing time is
!          spent inside while(r12...r41 >= 0) loop, hence its body is
!          minimized at the expense of somewhat awkwardness of what is
!          around it.
      i=nx/2 ; j=ny/2 ; ic=icmin ; jc=jcmin ; iinc=+1

      dx1=x(i+1,j  )-xc(ic,jc) ; dy1=y(i+1,j  )-yc(ic,jc)   !  3 <-- 2
      dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)   !  !     !
      dx3=x(i  ,j+1)-xc(ic,jc) ; dy3=y(i  ,j+1)-yc(ic,jc)   !  !     !
      dx4=x(i  ,j  )-xc(ic,jc) ; dy4=y(i  ,j  )-yc(ic,jc)   !  4 --> 1

      r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
      r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      do while(.true.)
        do while( r12 >= 0.D0 .and. r23 >= 0.D0
     &      .and. r34 >= 0.D0 .and. r41 >= 0.D0 )

          ip(ic,jc)=i ; jp(ic,jc)=j        !<-- record bounding indices

          if (iinc > 0 .and. ic < icmax) then    !--> proceed to the
            ic=ic+1                              !    next target point
          elseif (iinc < 0 .and. ic > icmin) then
            ic=ic-1
          else
            jc=jc+1 ; iinc=-iinc ; if (jc > jcmax) return
          endif

          dx1=x(i+1,j  )-xc(ic,jc) ; dy1=y(i+1,j  )-yc(ic,jc)
          dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
          dx3=x(i  ,j+1)-xc(ic,jc) ; dy3=y(i  ,j+1)-yc(ic,jc)
          dx4=x(i  ,j  )-xc(ic,jc) ; dy4=y(i  ,j  )-yc(ic,jc)

          r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
          r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4
        enddo
        inew=i ; jnew=j
        if (r12 < 0.D0) inew=min(inew+1, nx-1)
        if (r23 < 0.D0) jnew=min(jnew+1, ny-1)
        if (r34 < 0.D0) inew=max(inew-1,    1)
        if (r41 < 0.D0) jnew=max(jnew-1,    1)

        if (inew == i .and. jnew == j) then
          ip(ic,jc)=-1 ; jp(ic,jc)=-1       !<-- cannot be bounded
          if (iinc > 0 .and. ic < icmax) then  !--> proceed to the
            ic=ic+1                            !    next target point
          elseif (iinc < 0 .and. ic > icmin) then
            ic=ic-1
          else
            jc=jc+1 ; iinc=-iinc ; if (jc > jcmax) return
          endif
        else
          i=inew ; j=jnew       !<-- accept move of source location.
        endif

        dx1=x(i+1,j  )-xc(ic,jc) ; dy1=y(i+1,j  )-yc(ic,jc)
        dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
        dx3=x(i  ,j+1)-xc(ic,jc) ; dy3=y(i  ,j+1)-yc(ic,jc)
        dx4=x(i  ,j  )-xc(ic,jc) ; dy4=y(i  ,j  )-yc(ic,jc)

        r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
        r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      enddo
      END SUBROUTINE


      SUBROUTINE comp_offsets_tile(icmin,icmax,jcmin,jcmax, nx,ny,
     &                         x,y, ncx,ncy, xc,yc, ip,jp, xi,eta)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer ic,jc, i,j, iter
      real(kind=8) a11,a12,a21,a22, det, dX,dY, p,p1,q,q1
      do jc=jcmin,jcmax
        do ic=icmin,icmax
          i=ip(ic,jc) ; j=jp(ic,jc)
          if (i > 0 .and. j > 0) then
            a11=0.5D0*(x(i+1,j+1)-x(i,j+1) +x(i+1,j)-x(i,j))
            a12=0.5D0*(x(i+1,j+1)-x(i+1,j) +x(i,j+1)-x(i,j))
            a21=0.5D0*(y(i+1,j+1)-y(i,j+1) +y(i+1,j)-y(i,j))
            a22=0.5D0*(y(i+1,j+1)-y(i+1,j) +y(i,j+1)-y(i,j))

            dX=xc(ic,jc)-0.25D0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))
            dY=yc(ic,jc)-0.25D0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))

            det=1.D0/(a11*a22-a12*a21)
            xi(ic,jc) =0.5D0 + det*(a22*dX-a12*dY)
            eta(ic,jc)=0.5D0 + det*(a11*dY-a21*dX)
          else
            xi(ic,jc)=-1.D0 ; eta(ic,jc)=-1.D0  !<-- special values
          endif
        enddo
        do iter=1,10
          do ic=icmin,icmax
            i=ip(ic,jc) ; j=jp(ic,jc)
            if (i > 0 .and. j > 0) then
              p=xi(ic,jc)  ; p1=1.D0-p
              q=eta(ic,jc) ; q1=1.D0-q

              a11=q*(x(i+1,j+1)-x(i,j+1)) +q1*(x(i+1,j)-x(i,j))
              a12=p*(x(i+1,j+1)-x(i+1,j)) +p1*(x(i,j+1)-x(i,j))
              a21=q*(y(i+1,j+1)-y(i,j+1)) +q1*(y(i+1,j)-y(i,j))
              a22=p*(y(i+1,j+1)-y(i+1,j)) +p1*(y(i,j+1)-y(i,j))

              dX=xc(ic,jc) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     &                      -p*q1*x(i+1,j)  -p1*q1*x(i,j)
              dY=yc(ic,jc) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     &                      -p*q1*y(i+1,j)  -p1*q1*y(i,j)

              det=1.D0/(a11*a22-a12*a21)

              xi(ic,jc) =p + det*(a22*dX-a12*dY)
              eta(ic,jc)=q + det*(a11*dY-a21*dX)
!            else
!            print *,'xi, eta = -1 pour ic,jc = ',ic,jc
            endif
          enddo
        enddo
      enddo
      END SUBROUTINE






      SUBROUTINE spln2d(nx,ny, src, srX,srY,sXY,sYX)
! The following routine takes a two-dimensional array "src" and
! computes its derivatives
!
!        d src           d src           d^2 src           d^2 src
! srX = -------   srY = -------   sXY = ---------   sYX = ---------
!          dx             dy              dx dy             dy dx
!
! to construct  bi-cubic spline interpolation.  All the derivatives
! are defined at the same grid locations as the field itself and are
! obtained from the condition of continuity of the second derivative
! in the same direction.  The spline algorithm is one-dimensional
! Gaussian elimination procedure applied individually in each
! direction, so the mixed derivatives are computed as derivative of
! previously computed derivative in transversal direction.
! It is expected that sXY == sYX, however this routine computes both
! of them (just a matter of which direction goes first and which
! follows) and compares the outcome.
!
! All the derivatives are defined as "with respect to array index",
! i.e., dx=dy=1 in both directions.

! "Natural" boundary conditions are assumed on the sides (the second
! derivative vanishes at the ends).

! Array "src" which can be be either single or double precision,
! matching TYPE in the name of the subroutine, however all computed
! derivatives srX,srY,sXY,sYX are stored in double precision arrays
! regardless of type of "src".


      implicit none
      integer nx,ny, i,j
      real(kind=8), dimension(nx,ny) :: src
      real(kind=8), dimension(nx,ny) :: srX,srY,sXY,sYX
      real(kind=8) cff, error
      real(kind=8), dimension(nx) :: CF
      real(kind=8), dimension(ny) :: FC
!      allocate(CF(1:nx))
!      allocate(FC(1:ny))
      CF(1)=0.5D0
      do j=1,ny
        srX(1,j)=1.5D0*(src(2,j)-src(1,j))
      enddo
      do i=2,nx-1,+1
        CF(i)=1.D0/(4.D0-CF(i-1))
        cff=CF(i)
        do j=1,ny
          srX(i,j)=cff*(3.D0*(src(i+1,j)-src(i-1,j)) -srX(i-1,j))
        enddo
      enddo
      cff=1.D0/(2.D0-CF(nx-1))
      do j=1,ny
        srX(nx,j)=cff*(3.D0*(src(nx,j)-src(nx-1,j))-srX(nx-1,j))
      enddo
      do i=nx-1,1,-1
        cff=CF(i)
        do j=1,ny
          srX(i,j)=srX(i,j)-cff*srX(i+1,j)
        enddo
      enddo
      FC(1)=0.5D0
      do i=1,nx
        srY(i,1)=1.5D0*(src(i,2)-src(i,1))
        sXY(i,1)=1.5D0*(srX(i,2)-srX(i,1))
      enddo
      do j=2,ny-1,+1
        FC(j)=1.D0/(4.D0-FC(j-1))
        cff=FC(j)
        do i=1,nx
          srY(i,j)=cff*( 3.D0*(src(i,j+1)-src(i,j-1)) -srY(i,j-1))
          sXY(i,j)=cff*( 3.D0*(srX(i,j+1)-srX(i,j-1)) -sXY(i,j-1))
        enddo
      enddo
      cff=1.D0/(2.D0-FC(ny-1))
      do i=1,nx
        srY(i,ny)=cff*( 3.D0*(src(i,ny)-src(i,ny-1)) -srY(i,ny-1))
        sXY(i,ny)=cff*( 3.D0*(srX(i,ny)-srX(i,ny-1)) -sXY(i,ny-1))
      enddo
      do j=ny-1,1,-1
        cff=FC(j)
        do i=1,nx
          srY(i,j)=srY(i,j)-cff*srY(i,j+1)
          sXY(i,j)=sXY(i,j)-cff*sXY(i,j+1)
        enddo
      enddo
      CF(1)=0.5D0
      do j=1,ny
        sYX(1,j)=1.5D0*(srY(2,j)-srY(1,j))
      enddo
      do i=2,nx-1,+1
        CF(i)=1.D0/(4.D0-CF(i-1))
        cff=CF(i)
        do j=1,ny
          sYX(i,j)=cff*(3.D0*(srY(i+1,j)-srY(i-1,j)) -sYX(i-1,j))
        enddo
      enddo
      cff=1.D0/(2.D0-CF(nx-1))
      do j=1,ny
        sYX(nx,j)=cff*(3.D0*(srY(nx,j)-srY(nx-1,j))-sYX(nx-1,j))
      enddo
      do i=nx-1,1,-1
        cff=CF(i)
        do j=1,ny
          sYX(i,j)=sYX(i,j)-cff*sYX(i+1,j)
        enddo
      enddo
      error=0.D0
      do j=1,ny
        do i=1,nx
          sXY(1,j)=0.5D0*(sXY(i,j)+sYX(i,j))
        enddo
      enddo
      END SUBROUTINE spln2d



      SUBROUTINE spln2d_interp(nx,ny, src, srX, srY, sXY,
     &                          Istr, ncx,Jstr,ncy, ip,jp, xi,eta, targ)
      implicit none
      integer nx,ny, ncx,ncy, Istr, Jstr
      real(kind=8), dimension(nx,ny) :: src
      real(kind=8),    dimension(nx,ny) :: srX,srY,sXY
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      real(kind=8),    dimension(ncx,ncy) :: xi,eta
      real(kind=8), dimension(ncx,ncy) :: targ
      integer icmin,icmax,jcmin,jcmax, chunk_size
      icmin=1 ; icmax=ncx ; jcmin=1 ; jcmax=ncy ; chunk_size=0

      call spln2d_interp_tile(nx,ny, src, srX, srY, sXY,
     &                           ncx,ncy, ip,jp, xi,eta, targ,
     &                               icmin,icmax,jcmin,jcmax)
      END SUBROUTINE spln2d_interp

      SUBROUTINE spln2d_interp_tile(nx,ny, src, srX, srY, sXY,
     &                                 ncx,ncy, ip,jp, xi,eta, targ,
     &                                     icmin,icmax,jcmin,jcmax)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=8), dimension(nx,ny) :: src
      real(kind=8),    dimension(nx,ny) :: srX,srY,sXY
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      real(kind=8),    dimension(ncx,ncy) :: xi,eta
      real(kind=8), dimension(ncx,ncy) :: targ
!-->  internal variables
      integer ic,jc, iL,iR, jL,jR
      real(kind=8) px,qx,pqx, HxR,HxL,GxR,GxL,
     &             py,qy,pqy, HyR,HyL,GyR,GyL
     
      do jc=jcmin,jcmax
        do ic=icmin,icmax
          if (ip(ic,jc) > 0 .and. jp(ic,jc) > 0) then
            iL=ip(ic,jc)  ;  px=xi(ic,jc)
            iR=iL+1       ;  qx=1.D0-px ; pqx=px*qx
            GxR=-px*pqx   ;  HxR=px*(px+2.D0*pqx)
            GxL= qx*pqx   ;  HxL=qx*(qx+2.D0*pqx)

            jL=jp(ic,jc)  ;  py=eta(ic,jc)
            jR=jL+1       ;  qy=1.D0-py ; pqy=py*qy
            GyR=-py*pqy   ;  HyR=py*(py+2.D0*pqy)
            GyL= qy*pqy   ;  HyL=qy*(qy+2.D0*pqy)

            targ(ic,jc) = HxL*HyL*src(iL,jL) + HxR*HyL*src(iR,jL)
     &                  + HxL*HyR*src(iL,jR) + HxR*HyR*src(iR,jR)

     &                  + GxL*HyL*srX(iL,jL) + GxR*HyL*srX(iR,jL)
     &                  + GxL*HyR*srX(iL,jR) + GxR*HyR*srX(iR,jR)

     &                  + HxL*GyL*srY(iL,jL) + HxR*GyL*srY(iR,jL)
     &                  + HxL*GyR*srY(iL,jR) + HxR*GyR*srY(iR,jR)

     &                  + GxL*GyL*sXY(iL,jL) + GxR*GyL*sXY(iR,jL)
     &                  + GxL*GyR*sXY(iL,jR) + GxR*GyR*sXY(iR,jR)
          else
            targ(ic,jc) = 0.D0
          endif
        enddo
      enddo

      END SUBROUTINE spln2d_interp_tile

