MODULE kernels

! Module to calculate geometric sources
USE geometric_mod, ONLY: polyprism_cls


IMPLICIT NONE
! Variables:
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(p=4, r=4)    ! Single precision 
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(p=8, r=8)    ! Double precision                   
PUBLIC

CONTAINS
! CREATE THE PLOUFF CALCULATIONS OF THE POTENCIAL OF 1/R FOR A PRISM WITH POLIGONAL CROSS-SECTION

SUBROUTINE kernelxx(xp,yp,zp,res) ! res é um vetor cuja a dimensao esta associada ao numero de pontos de obs
 IMPLICIT NONE
 TYPE(polyprism_cls), INTENT(IN) :: pp
 REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: res ! result

 REAL(KIND=DP):: dummy,z1,z2,nverts ! variables comming from the polyprism_class
 !REAL(KIND=DP):: z1e,z2e
 INTEGER(KIND=DP):: nx ,nv ! number of vertices in the pp class
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) x, y, z1e, z2e, z1_sqr, z2_sqr, x1, x2, y1, y2, dx, dy, n, g, dist, cross, p, d1, d2 &
 v1_sqr, v2_sqr, r11, r12, r21, r22 , atan_diff_d2, atan_diff_d1

! Check dimensionality:
nx = SIZE(xp,1)  ! number of data points
 IF(nx /= SIZE(yp,1) .AND. nx /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1) )THEN
   PRINT*,'Dimensionality error in subroutine kernelxx: xp,yp,zp must have the same dimension'
   STOP
 ENDIF  

! Dummy variable for singularity points:
 dummy = 1e-10
! Get polygon info:
 z1 = pp%z1 
 z2 = pp%z2
 nv = pp%nv
! Allocating local arrays:
 ALLOCATE(x(nv), y(nv), z1e(nx), z2e(nx), z1_sqr(nx), z2_sqr(nx), x1(nx), x2(nx), dx(nx), dy(nx), n(nx), g(nx), dist(nx), &
 cross(nx), p(nx), d1(nx), d2(nx), v1_sqr(nx), v2_sqr(nx), r11(nx), r12(nx), r21(nx), r22(nx), atan_diff_d2(nx), atan_diff_d1(nx),&
 tmp(nx) ) 
x  = pp%xv ! pp = object (polyprism_class) and xv is an array inside pp
y  = pp%yv
! Calculate the effect of the prism
! dealing with z coordinates of the polygon:
 z1e = z1 - zp
 z2e = z2 - zp
 ! Square:
 z1_sqr = z1e * z1e
 z2_sqr = z2e * z2e

 ! Cleaning up the res array:
 res = 0d0
 DO i=0,nv
   x1    = x(i) - xp          
   y1    = y(i) - yp
   x2    = x( MOD(i+1), nv ) - xp ! mod serve para criar um laço interno (gerenciando as arestas do poligono) 
   y2    = y( MOD(i+1), nv ) - yp
   dx    = x2 - x1 + dummy
   dy    = y2 - y1 + dummy
   n     = dx/dy
   g     = x1 - y1 * n
   dist  = DSQRT(dx*dx + dy*dy)
   cross = x1*y2 - x2*y1
   p     = cross/dist + dummy
   d1    = (dx*x1 + dy*y1) / dist + dummy
   d2    = (dx*x2 + dy*y2) / dist + dummy 
   
   v1_sqr = x1*x1 + y1*y1
   v2_sqr = x2*x2 + y2*y2
  
   r11 = DSQRT(v1_sqr + z1_sqr)
   r12 = DSQRT(v1_sqr + z2_sqr)
   r21 = DSQRT(v2_sqr + z1_sqr)
   r22 = DSQRT(v2_sqr + z2_sqr)

   atan_diff_d2 = DATAN2(z2*d2, p*r22) - DATAN2(z1*d2, p*r21)
   atan_diff_d1 = DATAN2(z2*d1, p*r12) - DATAN2(z1*d1, p*r11)
   
   tmp = g*y2*atan_diff_d2 /(p*d2) + n * p * atan_diff_d2 / (d2)
   tmp =  tmp - ( g * y1 * atan_diff_d1 /(p*d1) + n * p * atan_diff_d1/(d1) )
   tmp = tmp + n*DLOG( (z2 + r12)*(z1 + r21) / ( (z1 + r11)*(z2 + r22) + dummy ) + dummy ) 
   tmp = tmp * (-1.0 / (1.0 + n * n) )
   res = res + tmp
 ENDDO
END SUBROUTINE kernelxx

!##############################################################
!##############################################################
!##############################################################

SUBROUTINE kernelxy(xp,yp,zp,res) ! res é um vetor cuja a dimensao esta associada ao numero de pontos de obs
 IMPLICIT NONE
 TYPE(polyprism_cls), INTENT(IN) :: pp
 REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: res ! result

 REAL(KIND=DP):: dummy,z1,z2,nverts ! variables comming from the polyprism_class
 INTEGER(KIND=DP):: nx ,nv ! number of vertices in the pp class
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) x, y, z1e, z2e, z1_sqr, z2_sqr, x1, x2, y1, y2, dx, dy, n,g, gg dist, cross, p, d1, d2 &
 v1_sqr, v2_sqr, r11, r12, r21, r22 , atan_diff_d2, atan_diff_d1, g_sqr

! Check dimensionality:
nx = SIZE(xp,1)  ! number of data points
 IF(nx /= SIZE(yp,1) .AND. nx /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1) )THEN
   PRINT*,'Dimensionality error in subroutine kernelxx: xp,yp,zp must have the same dimension'
   STOP
 ENDIF  

! Dummy variable for avoiding singularity points:
 dummy = 1e-10
 ! Get the number of vertices fo allocation of working arrays:
 nv = pp%nv  ! number of vertices

! Allocating local arrays:
 ALLOCATE(x(nv), y(nv), z1e(nx), z2e(nx), z1_sqr(nx), z2_sqr(nx), x1(nx), x2(nx), dx(nx), dy(nx), n(nx), g(nx), gg(nx), dist(nx), &
 cross(nx), p(nx), d1(nx), d2(nx), v1_sqr(nx), v2_sqr(nx), r11(nx), r12(nx), r21(nx), r22(nx), atan_diff_d2(nx), atan_diff_d1(nx),&
 tmp(nx) ) 

! Get polygon info:
 z1 = pp%z1  ! top of prism
 z2 = pp%z2  ! bottom of prism
 x  = pp%xv ! pp = object (polyprism_class) and xv is an array inside pp
 y  = pp%yv

! Calculate the potencial of the prism:
 z1e = z1 - zp
 z2e = z2 - zp
 ! Square:
 z1_sqr = z1e * z1e
 z2_sqr = z2e * z2e

 ! Cleaning up the res array:
 res = 0.0
 DO i=0,nv
   x1    = x(i) - xp          
   y1    = y(i) - yp
   x2    = x( MOD(i+1), nv ) - xp ! mod serve para criar um laço interno (gerenciando as arestas do poligono) 
   y2    = y( MOD(i+1), nv ) - yp
   dx    = x2 - x1 + dummy
   dy    = y2 - y1 + dummy
   n     = dx / dy
   g     = x1 - y1 * n
   gg    = g*g
   dist  = DSQRT(dx*dx + dy*dy)
   cross = x1*y2 - x2*y1
   p     = cross/dist + dummy
   d1    = (dx*x1 + dy*y1) / dist + dummy
   d2    = (dx*x2 + dy*y2) / dist + dummy 
   
   v1_sqr = x1*x1 + y1*y1    ! square of vertices values
   v2_sqr = x2*x2 + y2*y2
  
   r11 = DSQRT(v1_sqr + z1_sqr)
   r12 = DSQRT(v1_sqr + z2_sqr)
   r21 = DSQRT(v2_sqr + z1_sqr)
   r22 = DSQRT(v2_sqr + z2_sqr)

   atan_diff_d2 = DATAN2(z2*d2, p*r22) - DATAN2(z1*d2, p*r21)
   atan_diff_d1 = DATAN2(z2*d1, p*r12) - DATAN2(z1*d1, p*r11)
   
   tmp = ( gg + g*n*y2) * atan_diff_d2/(p*d2) - p*atan_diff_d2/d2
   tmp = tmp - ( (gg + g*n*y1) * atan_diff_d1/(p*d1) - p*atan_diff_d1/d1 ) 
   tmp = tmp + DLOG( (z2+r22)*(z1+r11) / ( (z1+r21)*(z2+r12) + dummy ) + dummy )
   tmp = tmp * ( 1.0/ (1.0 + n*n) )
   res = res + tmp
 ENDDO
END SUBROUTINE kernelxy

!##############################################################
!##############################################################
!##############################################################

SUBROUTINE kernelyy(xp,yp,zp,res) ! res é um vetor cuja a dimensao esta associada ao numero de pontos de obs
 IMPLICIT NONE
 TYPE(polyprism_cls), INTENT(IN) :: pp
 REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: res ! result

 REAL(KIND=DP):: dummy,z1,z2,nverts ! variables comming from the polyprism_class
 INTEGER(KIND=DP):: nx ,nv ! number of vertices in the pp class
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) x, y, z1e, z2e, z1_sqr, z2_sqr, x1, x2, y1, y2, dx, dy, m, c, dist, cross, p, d1, d2 &
 v1_sqr, v2_sqr, r11, r12, r21, r22 , atan_diff_d2, atan_diff_d1, g_sqr

! Check dimensionality:
nx = SIZE(xp,1)  ! number of data points
 IF(nx /= SIZE(yp,1) .AND. nx /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1) )THEN
   PRINT*,'Dimensionality error in subroutine kernelxx: xp,yp,zp must have the same dimension'
   STOP
 ENDIF  

! Dummy variable for avoiding singularity points:
 dummy = 1e-10
 ! Get the number of vertices fo allocation of working arrays:
 nv = pp%nv  ! number of vertices

! Allocating local arrays:
 ALLOCATE(x(nv), y(nv), z1e(nx), z2e(nx), z1_sqr(nx), z2_sqr(nx), x1(nx), x2(nx), dx(nx), dy(nx), m(nx), c(nx), dist(nx), &
 cross(nx), p(nx), d1(nx), d2(nx), v1_sqr(nx), v2_sqr(nx), r11(nx), r12(nx), r21(nx), r22(nx), atan_diff_d2(nx), atan_diff_d1(nx),&
 tmp(nx) ) 

! Get polygon info:
 z1 = pp%z1  ! top of prism
 z2 = pp%z2  ! bottom of prism
 x  = pp%xv ! pp = object (polyprism_class) and xv is an array inside pp
 y  = pp%yv

! Calculate the potencial of the prism:
 z1e = z1 - zp
 z2e = z2 - zp
 ! Square:
 z1_sqr = z1e * z1e
 z2_sqr = z2e * z2e

 ! Cleaning up the res array:
 res = 0.0
 DO i=0,nv
   x1    = x(i) - xp          
   y1    = y(i) - yp
   x2    = x( MOD(i+1), nv ) - xp ! mod serve para criar um laço interno (gerenciando as arestas do poligono) 
   y2    = y( MOD(i+1), nv ) - yp
   dx    = x2 - x1 + dummy
   dy    = y2 - y1 + dummy
   m     = dy / dx
   c     = y1 - x1 * m
   dist  = DSQRT(dx*dx + dy*dy)
   cross = x1*y2 - x2*y1
   p     = cross/dist + dummy
   d1    = (dx*x1 + dy*y1) / dist + dummy
   d2    = (dx*x2 + dy*y2) / dist + dummy 
   
   v1_sqr = x1*x1 + y1*y1    ! square of vertices values
   v2_sqr = x2*x2 + y2*y2
  
   r11 = DSQRT(v1_sqr + z1_sqr)
   r12 = DSQRT(v1_sqr + z2_sqr)
   r21 = DSQRT(v2_sqr + z1_sqr)
   r22 = DSQRT(v2_sqr + z2_sqr)

   atan_diff_d2 = DATAN2(z2*d2, p*r22) - DATAN2(z1*d2, p*r21)
   atan_diff_d1 = DATAN2(z2*d1, p*r12) - DATAN2(z1*d1, p*r11)
   
   tmp = c*x2*atan_diff_d2 / (p*d2) - m*p*atan_diff_d2 / d2
   tmp = tmp - c*x1*atan_diff_d1 / (p*d1) - m*p*atan_diff_d1 / d1
   tmp = tmp + m*DLOG( (z2+r12) * (z1+r21) / ( (z2+r22)*(z1+r11) ) + dummy )
   tmp = tmp * 1.0 / (1.0 + m*m)   
   res = res + tmp
 ENDDO
END SUBROUTINE kernelyy

!##############################################################
!##############################################################
!##############################################################

SUBROUTINE kernelxz(xp,yp,zp,res) ! res é um vetor cuja a dimensao esta associada ao numero de pontos de obs
 IMPLICIT NONE
 TYPE(polyprism_cls), INTENT(IN) :: pp
 REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: res ! result

 REAL(KIND=DP):: dummy,z1,z2,nverts ! variables comming from the polyprism_class
 INTEGER(KIND=DP):: nx ,nv ! number of vertices in the pp class
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) x, y, z1e, z2e, z1e2, z2e2, x1, x2, y1, y2, dx, dy, n, g, dist, p, d1, d2 &
 v1_sqr, v2_sqr, r11, r12, r21, r22 , log_diff_d2, log_diff_d1, ng, nn_p1

! Check dimensionality:
nx = SIZE(xp,1)  ! number of data points
 IF(nx /= SIZE(yp,1) .AND. nx /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1) )THEN
   PRINT*,'Dimensionality error in subroutine kernelxx: xp,yp,zp must have the same dimension'
   STOP
 ENDIF  

! Dummy variable for avoiding singularity points:
 dummy = 1e-10
! Get polygon info:
 z1 = pp%z1  ! top of prism
 z2 = pp%z2  ! bottom of prism
 nv = pp%nv  ! number of vertices

! Allocating local arrays:
 ALLOCATE(x(nv), y(nv), z1e(nx), z2e(nx), z1_sqr(nx), z2_sqr(nx), x1(nx), x2(nx), dx(nx), dy(nx), n(nx), g(nx), dist(nx), &
 p(nx), d1(nx), d2(nx), v1_sqr(nx), v2_sqr(nx), r11(nx), r12(nx), r21(nx), r22(nx), log_diff_d2(nx), log_diff_d1(nx),&
 tmp(nx), ng(nx), nn_p1(nx) )
 
 ! Coordinates of the prism:
 x  = pp%xv ! pp = object (polyprism_class) and xv is an array inside pp
 y  = pp%yv

! Calculate the potencial of the prism:
 z1e = z1 - zp
 z2e = z2 - zp
 ! Square:
 z1_sqr = z1e * z1e
 z2_sqr = z2e * z2e

 ! Cleaning up the res array:
 res = 0.0
 DO i=0,nv ! check if we do need to start from zero as well
   x1    = x(i) - xp          
   y1    = y(i) - yp
   x2    = x( MOD(i+1), nv ) - xp ! mod serve para criar um laço interno (gerenciando as arestas do poligono) 
   y2    = y( MOD(i+1), nv ) - yp
   dx    = x2 - x1 + dummy
   dy    = y2 - y1 + dummy
   n     = dx / dy
   nn_p1 = n*n + 1.0
   g     = x1 - y1 * n
   ng    = n * g
   dist  = DSQRT(dx*dx + dy*dy)
   p     = cross/dist + dummy
   d1    = (dx*x1 + dy*y1) / dist + dummy
   d2    = (dx*x2 + dy*y2) / dist + dummy 
   
   v1_sqr = x1*x1 + y1*y1    ! square of vertices values
   v2_sqr = x2*x2 + y2*y2
  
   r11 = DSQRT(v1_sqr + z1_sqr)
   r12 = DSQRT(v1_sqr + z2_sqr)
   r21 = DSQRT(v2_sqr + z1_sqr)
   r22 = DSQRT(v2_sqr + z2_sqr)
 
   log_r22 = DLOG( (r22 - d2)/ (r22 + d2) + dummy )
   log_r21 = DLOG( (r21 - d2)/ (r21 + d2) + dummy )
   log_r12 = DLOG( (r12 - d1)/ (r12 + d1) + dummy )
   log_r11 = DLOG( (r11 - d1)/ (r11 + d1) + dummy ) 

   log_diff_d1 = (0.5/d1) * (log_r12 - log_r11) 
   log_diff_d2 = (0.5/d2) * (log_r22 - log_r21)
   
   tmp = (y2 * nn_p1 + ng) * log_diff_d2
   tmp = tmp - ( (y1 * nn_p1 + ng) * log_diff_d1 )
   tmp = tmp * ( -1.0 / (nn_p1) )
   res = res + tmp
   
 ENDDO
END SUBROUTINE kernelxz

!##############################################################
!##############################################################
!##############################################################

SUBROUTINE kernelyz(xp,yp,zp,res) ! res é um vetor cuja a dimensao esta associada ao numero de pontos de obs
 IMPLICIT NONE
 TYPE(polyprism_cls), INTENT(IN) :: pp
 REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: res ! result

 REAL(KIND=DP):: dummy,z1,z2,nverts ! variables comming from the polyprism_class
 INTEGER(KIND=DP):: nx ,nv ! number of vertices in the pp class
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) x, y, z1e, z2e, z1e2, z2e2, x1, x2, y1, y2, dx, dy, m, c, dist, p, d1, d2 &
 v1_sqr, v2_sqr, r11, r12, r21, r22 , log_diff_d2, log_diff_d1, mm_p1, cm

! Check dimensionality:
nx = SIZE(xp,1)  ! number of data points
 IF(nx /= SIZE(yp,1) .AND. nx /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1) )THEN
   PRINT*,'Dimensionality error in subroutine kernelxx: xp,yp,zp must have the same dimension'
   STOP
 ENDIF  

! Dummy variable for avoiding singularity points:
 dummy = 1e-10

 ! Get the number of vertices for allocation issues:
 nv = pp%nv  ! number of vertices

! Allocating local arrays:
 ALLOCATE(x(nv), y(nv), z1e(nx), z2e(nx), z1_sqr(nx), z2_sqr(nx), x1(nx), x2(nx), dx(nx), dy(nx), m(nx), c(nx), dist(nx), &
 p(nx), d1(nx), d2(nx), v1_sqr(nx), v2_sqr(nx), r11(nx), r12(nx), r21(nx), r22(nx), log_diff_d2(nx), log_diff_d1(nx),&
 tmp(nx), mm_p1(nx), cm(nx) )
 
 ! Get polygon info:
 z1 = pp%z1  ! top of prism
 z2 = pp%z2  ! bottom of prism
 x  = pp%xv ! pp = object (polyprism_class) and xv is an array inside pp
 y  = pp%yv

! Calculate the potencial of the prism:
 z1e = z1 - zp
 z2e = z2 - zp
 ! Square:
 z1_sqr = z1e * z1e
 z2_sqr = z2e * z2e

 ! Cleaning up the res array:
 res = 0.0
 DO i=0,nv ! check if we do need to start from zero as well
   x1    = x(i) - xp          
   y1    = y(i) - yp
   x2    = x( MOD(i+1), nv ) - xp ! mod serve para criar um laço interno (gerenciando as arestas do poligono) 
   y2    = y( MOD(i+1), nv ) - yp
   dx    = x2 - x1 + dummy
   dy    = y2 - y1 + dummy
   m     = dy/dx
   mm_p1 = m*m + 1.0
   c     = y1 - x1 * m
   cm    = c * m
   dist  = DSQRT(dx*dx + dy*dy)
   p     = cross/dist + dummy
   d1    = (dx*x1 + dy*y1) / dist + dummy
   d2    = (dx*x2 + dy*y2) / dist + dummy 
   
   v1_sqr = x1*x1 + y1*y1    ! square of vertices values
   v2_sqr = x2*x2 + y2*y2
  
   r11 = DSQRT(v1_sqr + z1_sqr)
   r12 = DSQRT(v1_sqr + z2_sqr)
   r21 = DSQRT(v2_sqr + z1_sqr)
   r22 = DSQRT(v2_sqr + z2_sqr)
 
   log_r22 = DLOG( (r22 - d2)/ (r22 + d2) + dummy )
   log_r21 = DLOG( (r21 - d2)/ (r21 + d2) + dummy )
   log_r12 = DLOG( (r12 - d1)/ (r12 + d1) + dummy )
   log_r11 = DLOG( (r11 - d1)/ (r11 + d1) + dummy ) 

   log_diff_d1 = (0.5/d1) * (log_r12 - log_r11) 
   log_diff_d2 = (0.5/d2) * (log_r22 - log_r21)
   
   tmp = (x2*mm_p1+cm) * (0.5/d2) * (log_r22 - log_r21)
   tmp = tmp - ( (x1*mm_p1+cm) * (0.5/d1) * (log_r12 - log_r11) )   
   tmp = tmp * (1.0/mm_p1)
   res = res + tmp
   
 ENDDO
END SUBROUTINE kernelyz

!##############################################################
!##############################################################
!##############################################################

SUBROUTINE kernelzz(xp,yp,zp,res) ! res é um vetor cuja a dimensao esta associada ao numero de pontos de obs
 IMPLICIT NONE
 TYPE(polyprism_class), INTENT(IN) :: pp
 REAL(KIND=DP), DIMENSION(:), INTENT(IN):: xp,yp,zp
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: res ! result

 REAL(KIND=DP):: dummy,z1,z2,nverts ! variables comming from the polyprism_class
 INTEGER(KIND=DP):: nx ,nv ! number of vertices in the pp class
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) x, y, z1e, z2e, z1e2, z2e2, x1, x2, y1, y2, dx, dy, dist, p, d1, d2 &
 v1_sqr, v2_sqr, r11, r12, r21, r22

! Check dimensionality:
nx = SIZE(xp,1)  ! number of data points
 IF(nx /= SIZE(yp,1) .AND. nx /= SIZE(zp,1) .AND. SIZE(yp,1) /= SIZE(zp,1) )THEN
   PRINT*,'Dimensionality error in subroutine kernelxx: xp,yp,zp must have the same dimension'
   STOP
 ENDIF  

! Dummy variable for avoiding singularity points:
 dummy = 1e-10

 ! Get the number of vertices for allocation issues:
 nv = pp%nv  ! number of vertices

! Allocating local arrays:
 ALLOCATE(x(nv), y(nv), z1e(nx), z2e(nx), z1_sqr(nx), z2_sqr(nx), x1(nx), x2(nx), dx(nx), dy(nx), dist(nx), &
 p(nx), d1(nx), d2(nx), v1_sqr(nx), v2_sqr(nx), r11(nx), r12(nx), r21(nx), r22(nx) )
 
 ! Get polygon info:
 z1 = pp%z1  ! top of prism
 z2 = pp%z2  ! bottom of prism
 x  = pp%xv ! pp = object (polyprism_class) and xv is an array inside pp
 y  = pp%yv

! Calculate the potencial of the prism:
 z1e = z1 - zp
 z2e = z2 - zp
 ! Square:
 z1_sqr = z1e * z1e
 z2_sqr = z2e * z2e

 ! Cleaning up the res array:
 res = 0.0
 DO i=0,nv ! check if we do need to start from zero as well
   x1    = x(i) - xp          
   y1    = y(i) - yp
   x2    = x( MOD(i+1), nv ) - xp ! mod serve para criar um laço interno (gerenciando as arestas do poligono) 
   y2    = y( MOD(i+1), nv ) - yp
   dx    = x2 - x1 ! + dummy
   dy    = y2 - y1 ! + dummy
   
   dist  = DSQRT(dx*dx + dy*dy) + dummy ! dist is only used in divisions, so we need to control division by zero with dummy factor
   cross = x1*y2 - x2*y1 
   p     = cross/dist  !+ dummy (no need in this case)
   d1    = (dx*x1 + dy*y1) / dist !+ dummy
   d2    = (dx*x2 + dy*y2) / dist !+ dummy 
   
   v1_sqr = x1*x1 + y1*y1    ! square of vertices values
   v2_sqr = x2*x2 + y2*y2
  
   r11 = DSQRT(v1_sqr + z1_sqr)
   r12 = DSQRT(v1_sqr + z2_sqr)
   r21 = DSQRT(v2_sqr + z1_sqr)
   r22 = DSQRT(v2_sqr + z2_sqr)
   res = res + (DATAN2(z2*d2,p*r22) - DATAN2(z1*d2,p*r21) - DATAN2(z2*d1,p*r12) + DATAN2(z1*d1,p*r11) )
 ENDDO
END SUBROUTINE kernelzz


END MODULE kernels
