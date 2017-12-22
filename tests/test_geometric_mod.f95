PROGRAM test_geometric_mod
 ! test program to verify if all subroutines of geometric_mod is working properly:
 ! authors: Bijani and Carreira 
 USE geometric_mod
 IMPLICIT NONE
 TYPE(polyprism_type):: pp ! the object
 INTEGER:: n , fid
 CHARACTER(LEN=25)::inputfile

 ! 1) test #01 of geometric clear: see if all values are cleared by the subroutine:
 pp%z1    = 11.0
 pp%z2    = 11.0
 pp%nv    = 1000
 pp%ppval = 20.0
 pp%props = 'dens'
 !pp%inc    = 12.0
 !pp%dec    = 44.0
 CALL geometric_clear(pp)
 IF(pp%z1 /= 0.0 .OR. pp%z1 /= 0.0 )THEN; PRINT*,'Test 01 for z1 and z2 of geometric_clear FAILED' ; STOP ; ENDIF
 IF(pp%nv /= 0 )THEN; PRINT*,'Test 01 for nv of geometric_clear FAILED' ; STOP ; ENDIF
 IF(pp%props /= 'none')THEN; PRINT*,'Test 01 for props name of geometric_clear FAILED' ; STOP ; ENDIF
 IF(pp%ppval /= 0.0)THEN; PRINT*,'Test 01 for ppval of geometric_clear FAILED' ; STOP ; ENDIF
 IF(pp%inc /= 0.0)THEN; PRINT*,'Test 01 for inc of geometric_clear FAILED' ; STOP ; ENDIF
 IF(pp%dec /= 0.0)THEN; PRINT*,'Test 01 for dec of geometric_clear FAILED' ; STOP ; ENDIF

 PRINT*,'(1) - Tests on geometric_clear PASSED'
 PRINT *, ''

 ! 2) test #02: check if geometric_allocate does really allocate arrays:
 pp%nv = 5000
 n     = pp%nv
 CALL geometric_allocate(pp)
 IF(SIZE(pp%xv,1) /= n .OR. SIZE(pp%yv,1) /= n )THEN; PRINT*,'Test 01 of geometric_allocate FAILED' ; STOP ; ENDIF   
 PRINT*,'(2) - Tests on geometric_allocate PASSED'
 PRINT *, ''  

 ! 3) test for input file:
 inputfile = '../tests/test_dens.inp'
 fid = 1
 CALL geometric_read_inputs(pp,fid,inputfile) 
 IF(pp%z1 == 0.0 .OR. pp%z1 == 0.0 )THEN; PRINT*,'Test 01 for z1 and z2 of geometric_read_inputs FAILED' ; STOP ; ENDIF
 IF(pp%nv == 0 )THEN; PRINT*,'Test 01 for nv of  geometric_read_inputs FAILED' ; STOP ; ENDIF
 IF(pp%props == 'none')THEN; PRINT*,'Test 01 for props name of geometric_read_inputs FAILED' ; STOP ; ENDIF
 IF(pp%ppval == 0.0)THEN; PRINT*,'Test 01 for ppval of geometric_read_inputs FAILED' ; STOP ; ENDIF
 IF(pp%inc /= 0.0)THEN; PRINT*,'Test 01 for inc of  geometric_read_inputs FAILED' ; STOP ; ENDIF
 IF(pp%dec /= 0.0)THEN; PRINT*,'Test 01 for dec of  geometric_read_inputs FAILED' ; STOP ; ENDIF

 PRINT*,'(3) - Tests on geometric_read_inputs PASSED'
 PRINT *, '' 

END PROGRAM test_geometric_mod
