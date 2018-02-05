MODULE obs_mod
! Module to create observation points:
! subroutine's names are inspired in python
IMPLICIT NONE
PUBLIC :: SP, DP

! Variables:
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(p=4, r=4)    ! Single precision 
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(p=15, r=10)    ! Double precision                   

CONTAINS

 SUBROUTINE arange(xi,xf,dx,x) 
 ! create an array with even equally-space data points:
 IMPLICIT NONE 
 REAL(KIND=DP), INTENT(IN):: xi,xf,dx
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: x ! array with data points
 ! Local variables:
 !REAL(KIND=DP):: nonint  ! to check if dx is a non-integer variable
 INTEGER:: i !, np
 
 !nonint = MOD(dx, 2d0)
 ! for consistency, dx should be a integer. otherwise, use linspace
 !IF(nonint /= 0 )THEN
 !  PRINT*, 'step is a non-integer number. You should try the subroutine linspace for gain of consistency.'
   !RETURN
 !ENDIF

 ! check if start is higher than stop:
 IF(xi<xf)THEN  ! start is lower than stop:
   x(1) = xi
   DO i=2, SIZE(x,1)
     x(i) = x(i-1) + dx
   ENDDO
 ELSE  ! start is higher than stop:
   x(1) = xi
   DO i=2, SIZE(x,1)
     x(i) = x(i-1) - dx
   ENDDO
 ENDIF   

 END SUBROUTINE arange
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

 SUBROUTINE linspace(xi, xf, np, x, endpoint)
 ! create an array x of even values ranging from xi to xf. inspired in python np.linspace script
 IMPLICIT NONE
 REAL(KIND=DP), INTENT(IN):: xi,xf  ! xi = start, xf = stop
 INTEGER, INTENT(IN):: np ! np = number of points
 REAL(KIND=DP), DIMENSION(:), INTENT(OUT):: x ! array with data points
 LOGICAL, INTENT(IN):: endpoint ! true will include xf into the x array; false, will not.

 REAL(KIND=DP):: dx ! step
 INTEGER:: i
 ! check array specifications:
 IF(np /= SIZE(x,1) )THEN
   PRINT*,'dimension of x array badly set. Verify for consistency.'
   RETURN
 ENDIF

 ! Only one point in the array:
 IF(np==1)THEN
   x = xf
 ENDIF
 ! check the endpoint condition:
 IF(endpoint)THEN ! include the xf value into the x array:
   dx = (xf - xi) / DFLOAT(np - 1)
 ELSE ! exclude the xf value from x array:
   dx = (xf - xi) / DFLOAT(np)
 ENDIF  
 ! generating the x array:
 IF(xi<xf)THEN
   DO i=1,np
     x(i) = xi + dx*DFLOAT(i-1)
   ENDDO
 ELSE
   DO i=np,1,-1
     x(i) = xi + dx*DFLOAT(i-1)
   ENDDO
 ENDIF 
END SUBROUTINE linspace
  
!--!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
 SUBROUTINE meshgrid(xi, xf, nx, yi, yf, ny, X, Y)
 ! create a 2D equally-spaced grid of points:
 IMPLICIT NONE
 REAL(KIND=DP), INTENT(IN):: xi,xf,yi,yf   ! x and y limits
 INTEGER, INTENT(IN):: nx, ny ! number of points in x and y directions
 REAL(KIND=DP), DIMENSION(nx,ny), INTENT(OUT) :: X,Y ! matrices 
 
 ! Local variables:
 REAL(KIND=DP), DIMENSION(nx):: wx ! one-dimensional array with x values
 REAL(KIND=DP), DIMENSION(ny):: wy ! one-dimensional array with y values
 INTEGER:: i,j
 
 ! Use subroutine linspace to create working arrays:
 CALL linspace(xi, xf, nx, wx, .TRUE.)
 CALL linspace(yi, yf, ny, wy, .TRUE.)
 
 ! create the grid:
 DO i=1,nx
   DO j=1,ny
     X(i,j) = wx(i)
     Y(i,j) = wy(j)
   ENDDO
 ENDDO

 END SUBROUTINE meshgrid
!--!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
END MODULE obs_mod
