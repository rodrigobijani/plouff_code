PROGRAM test_obs_mod
! test if all subroutines inside the module obs_mod.f95 are working properly
 USE obs_mod, ONLY: DP, arange, meshgrid, linspace
 IMPLICIT NONE
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:):: x 
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:):: Xi, Yi
 INTEGER:: i, j, np

 ! Set the allocation of x array in terms of number of points:
 np = 1 + CEILING( (10d0 - 0d0)/ 1d0 )
 ! Allocate array:
 ALLOCATE( x(np) )
 ! Cleanning the array:
 x = 0.0
! Test #01 : subroutine arrange for xi > xf :
 CALL arange (0d0, 10d0, 1d0, x)
 PRINT*, 'Test #01' 
 DO i=1,SIZE(x,1)
   PRINT*,x(i)
 ENDDO
 ! check if the ranges of x array are out of bounds:
 IF(ANY(x < 0.0) .OR. ANY(x > 10.0) )THEN; PRINT*,' Test #01 of subroutine arrange FAILED: (Array out of ranges) ' ; STOP ; ENDIF
 !DEALLOCATE(x)
 !-----------------------------------------------------------------------------------------------------------------------------------
 
 ! Test #02: subroutine arange for xi < xf :
 CALL arange(10d0, 0d0, 1d0, x)
 PRINT*, 'Test #02'
 DO i=1,SIZE(x,1)
   PRINT*,x(i)
 ENDDO
 IF(ANY(x < 0d0) .OR. ANY(x > 10d0) )THEN; PRINT*,' Test #02 of subroutine arrange FAILED: (Array out of ranges) ' ; STOP ; ENDIF
 DEALLOCATE(x)
 !-----------------------------------------------------------------------------------------------------------------------------------

 ! Test #03: subroutine arange dx < 1 :
 np = 1 + CEILING( (10d0 - 0d0)/ 5d-1 )
 ! Allocate array:
 ALLOCATE( x(np) )
 ! Cleanning the array:
 x = 0.0
 CALL arange(0d0, 10d0, 5d-1, x)
 PRINT*, 'Test #03'
 DO i=1,SIZE(x,1)
   PRINT*,x(i)
 ENDDO
 !IF(ANY(x/= 0d0) )THEN; PRINT*,' Test #03 of subroutine arrange FAILED: (dx < 1.0) ' ; STOP ; ENDIF
 DEALLOCATE(x)
 !-----------------------------------------------------------------------------------------------------------------------------------

  ! Test #04: subroutine arange for xi = xf :
  np = 1 + CEILING( (10d0 - 10d0)/ 1d0 )
  ! Allocate array:
  ALLOCATE( x(np) )
  ! Cleanning the array:
  x = 0.0
  CALL arange(10d0, 10d0, 1d0, x)
  IF( ALL(x/=10d0) )THEN; PRINT*,' Test #04 of subroutine arrange FAILED: (wrong value) ' ; STOP ; ENDIF
  IF(SIZE(x,1) /= 1 )THEN; PRINT*, ' Test #04 of subroutine arrange FAILED: (Allocation mismatch) ' ; STOP ; ENDIF
 ! DEALLOCATE(x) ! finishing tests with arange subroutine, so need to deallocate array
  PRINT*,'(1) - Tests with subroutine arange PASSED'
  DEALLOCATE(x) 
 !-----------------------------------------------------------------------------------------------------------------------------------

  ! Test #05: test subroutine linspace including the xf point into the array:
  ALLOCATE( x(20) )
  ! cleanning array:
  x = 0d0
  CALL linspace(0d0, 10d0, 20, x, .TRUE.)
  PRINT*,'Test #05'
  DO i=1,SIZE(x,1)
   PRINT*,x(i)
  ENDDO
  IF( x( SIZE(x,1) ) /= 10d0 )THEN; PRINT*, &
  ' Test #05 of subroutine linspace FAILED: (no extreme value in the array) ' ; STOP ; ENDIF
 !-----------------------------------------------------------------------------------------------------------------------------------
 
 ! Test #06 : subroutine linspace excluding xf of the array:
  CALL linspace(0d0, 10d0, 20, x, .FALSE.) 
  PRINT*,'Test #06'
  DO i=1,SIZE(x,1)
   PRINT*,x(i)
  ENDDO
  IF( x(SIZE(x,1) ) == 10d0 )THEN; PRINT*,& 
  ' Test #06 of subroutine linspace FAILED: (inclusion of xf value in the array) ' ; STOP ; ENDIF

!----------------------------------------------------------------------------------------------------------------------------------- 
 ! Test #07 : subroutine linspace with xf < xi:
  CALL linspace(10d0, 0d0, 20, x, .TRUE.) 
  PRINT*,'Test #07'
  DO i=1,SIZE(x,1)
   PRINT*,x(i)
  ENDDO
  IF( x(SIZE(x,1)) /= 0d0 )THEN; PRINT*,& 
  ' Test #07 of subroutine linspace FAILED: (inclusion of xi value in the array) ' ; STOP ; ENDIF
  DEALLOCATE(x) ! finishing tests with arange subroutine, so need to deallocate array
  PRINT*,'(2) - Tests with subroutine linspace PASSED' 
 !-----------------------------------------------------------------------------------------------------------------------------------
 
 ! Test #08: subroutine meshgrid for a square mesh:
  ALLOCATE(Xi(5,5), Yi(5,5) )
  CALL meshgrid(0d0, 10d0, 5, 0d0, 10d0, 5, Xi, Yi)
  PRINT*,'Xi'
  DO i=1,5
    WRITE(*,17) ( Xi(i,j), j=1,5)
  ENDDO

  PRINT*,'Yi'
  DO i=1,5
    WRITE(*,17) ( Yi(i,j) , j=1,5)
  ENDDO

  17 FORMAT(5(ES12.4E2,2x))		! formatação para matrizes em notação cientifica   

  IF(SIZE(Xi,1) /= SIZE(Xi,2) .OR. SIZE(Yi,1)/= SIZE(Yi,2) )THEN; PRINT*,& 
  ' Test #08 of subroutine meshgrid FAILED: (not a square mesh) ' ; STOP ; ENDIF
 PRINT*,'(3) - Tests with subroutine meshgrid PASSED' 
 !-----------------------------------------------------------------------------------------------------------------------------------
  
END PROGRAM test_obs_mod
