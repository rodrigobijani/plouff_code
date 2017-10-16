MODULE gravmag_mod
! Module to calculate gravity and magnetics quantities of a polygonal prism

USE geometric_mod, ONLY: polyprism_cls
USE kernels_mod

IMPLICIT NONE
! Variables:
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(p=4, r=4)    ! Single precision 
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(p=8, r=8)    ! Double precision
REAL(KIND=DP), PARAMETER:: CM = 10d-7, T2NT = 10d9 
! CM   = constant used in the magnetic method in henry/m (SI)
! T2NT = constant to convert from Tesla to nano Tesla (SI)                    
PUBLIC

CONTAINS
! CREATE THE PLOUFF CALCULATIONS OF THE POTENCIAL OF 1/R FOR A PRISM WITH POLIGONAL CROSS-SECTION

SUBROUTINE total_field(xp,yp,zp,pp,inc,dec,azi,tf) 
! Computes the total-field anomaly of a polygonal prism using Plouff 1976 paper (implementation based on Uieda 2013)
IMPLICIT NONE
TYPE(polyprism_cls), INTENT(IN):: pp ! prism
REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp ! obs points
REAL(KIND=DP), INTENT(IN):: inc,dec,azi ! inclination, declination and azimuth of the regional field (degrees)
REAL(KIND,DP), DIMENSION(:), INTENT(OUT):: tf
! Local Variables:
REAL(KIND=DP):: fx,fy,fz ! cartesian components of the regional field
REAL(KIND=DP):: mx,my,mz ! cartesian components of the magnetization of the polyprism
REAL(KIND=DP), DIMENSION( SIZE(tf,1) ) :: v1,v2,v3,v4,v5,v6,bx,by,bz ! values from kernels_mod and components of regional fiel

! Check dimensions of the obs points:
IF( SIZE(xp,1) /= SIZE(yp,1) .AND. SIZE(xp,1) /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1) )THEN
  PRINT*, 'gravmag_mod: subroutine Dimension of input arrays should be the same.'
  STOP
ENDIF

IF( SIZE(tf,1) /= SIZE(xp,1) ) THEN
  PRINT*, 'gravmag_mod: Dimension of tf array should be the same as data points (xp,yp,zp)'

! Compute the cartesian components of the regional field:
 CALL dircos(inc,dec,azi,fx,fy,fz)

! Compute the total field anomaly:
tf = 0.0

! get magnetization of the polyprism:
 mx = pp%mx
 my = pp%my
 mz = pp%mz
! call kernels:
 CALL kernelxx(xp,yp,zp,v1)
 CALL kernelxy(xp,yp,zp,v2)
 CALL kernelxz(xp,yp,zp,v3)
 CALL kernelyy(xp,yp,zp,v4)
 CALL kernelyz(xp,yp,zp,v5)
 CALL kernelzz(xp,yp,zp,v6)
! compute the components of the magnetic field of the polyprism:
 bx = v1*mx + v2*my + v3*mz
 by = v2*mx + v4*my + v5*mz
 bz = v3*mx + v5*my + v6*mz
 tf = fx*bx + fy*by + fz*bz
 tf = tf * CM * T2NT     ! tf is computed as if the regional and the source field have the same directions

END SUBROUTINE total_field

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE bx(xp,yp,zp,pp,b)
IMPLICIT NONE

! Computes the x component of the magnetic field of a polygonal prism using Plouff 1976 paper (implementation based on Uieda 2013)
IMPLICIT NONE
TYPE(polyprism_cls), INTENT(IN):: pp ! prism
REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp ! obs points
REAL(KIND,DP), DIMENSION(:), INTENT(OUT):: b ! bx component 
! Local Variables:
REAL(KIND=DP):: mx,my,mz ! cartesian components of the magnetization of the polyprism
REAL(KIND=DP), DIMENSION( SIZE(tf,1) ) :: v1,v2,v3 ! values from kernels_mod and components of regional fiel

! Check dimensions of the obs points:
IF(SIZE(xp,1) /= SIZE(yp,1) .AND. SIZE(xp,1) /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1))THEN
  PRINT*, 'gravmag_mod: subroutine bx: Dimension of input arrays should be the same.'
  STOP
ENDIF

IF(SIZE(b,1) /= SIZE(xp,1) ) THEN
  PRINT*, 'gravmag_mod: Dimension of tf array should be the same as data points (xp,yp,zp)'

! compute bx:
 b = 0.0 ! cleans the array
! get magnetization of the polyprism:
 mx = pp%mx
 my = pp%my
 mz = pp%mz
! call kernels:
 CALL kernelxx(xp,yp,zp,v1)
 CALL kernelxy(xp,yp,zp,v2)
 CALL kernelxz(xp,yp,zp,v3)
 b = mx*v1 + my*v2 + mz*v3
 ! convert to nT:
 b = b * CM * T2NT

END SUBROUTINE bx

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE by(xp,yp,zp,pp,b)
IMPLICIT NONE

! Computes the y component of the magnetic field of a polygonal prism using Plouff 1976 paper (implementation based on Uieda 2013)
IMPLICIT NONE
TYPE(polyprism_cls), INTENT(IN):: pp ! prism
REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp ! obs points
REAL(KIND,DP), DIMENSION(:), INTENT(OUT):: b ! bx component 
! Local Variables:
REAL(KIND=DP):: mx,my,mz ! cartesian components of the magnetization of the polyprism
REAL(KIND=DP), DIMENSION( SIZE(tf,1) ) :: v2,v4,v5 ! values from kernels_mod and components of regional fiel

! Check dimensions of the obs points:
IF(SIZE(xp,1) /= SIZE(yp,1) .AND. SIZE(xp,1) /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1))THEN
  PRINT*, 'gravmag_mod: subroutine bx: Dimension of input arrays should be the same.'
  STOP
ENDIF

IF(SIZE(b,1) /= SIZE(xp,1) ) THEN
  PRINT*, 'gravmag_mod: Dimension of tf array should be the same as data points (xp,yp,zp)'

! compute bx:
 b = 0.0 ! cleans the array
! get magnetization of the polyprism:
 mx = pp%mx
 my = pp%my
 mz = pp%mz
! call kernels:
 CALL kernelyy(xp,yp,zp,v4)
 CALL kernelxy(xp,yp,zp,v2)
 CALL kernelyz(xp,yp,zp,v5)
 b = mx*v2 + my*v4 + mz*v5
 ! convert to nT:
 b = b * CM * T2NT

END SUBROUTINE bx

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE bz(xp,yp,zp,pp,b)
IMPLICIT NONE

! Computes the z component of the magnetic field of a polygonal prism using Plouff 1976 paper (implementation based on Uieda 2013)
IMPLICIT NONE
TYPE(polyprism_cls), INTENT(IN):: pp ! prism
REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp ! obs points
REAL(KIND,DP), DIMENSION(:), INTENT(OUT):: b ! bx component 
! Local Variables:
REAL(KIND=DP):: mx,my,mz ! cartesian components of the magnetization of the polyprism
REAL(KIND=DP), DIMENSION( SIZE(tf,1) ) :: v3,v5,v6 ! values from kernels_mod and components of regional fiel

! Check dimensions of the obs points:
IF(SIZE(xp,1) /= SIZE(yp,1) .AND. SIZE(xp,1) /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1))THEN
  PRINT*, 'gravmag_mod: subroutine bx: Dimension of input arrays should be the same.'
  STOP
ENDIF

IF(SIZE(b,1) /= SIZE(xp,1) ) THEN
  PRINT*, 'gravmag_mod: Dimension of tf array should be the same as data points (xp,yp,zp)'

! compute bx:
 b = 0.0 ! cleans the array
! get magnetization of the polyprism:
 mx = pp%mx
 my = pp%my
 mz = pp%mz
! call kernels:
 CALL kernelxz(xp,yp,zp,v3)
 CALL kernelyz(xp,yp,zp,v5)
 CALL kernelzz(xp,yp,zp,v6)
 b = mx*v3 + my*v5 + mz*v6
 ! convert to nT:
 b = b * CM * T2NT

END SUBROUTINE bz

! TODO: implement the gravity outputs (gx,gy,gz,gxx,gxy,gxz....)!
  


END MODULE gravmag_mod
