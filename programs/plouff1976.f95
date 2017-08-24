PROGRAM plouff1976

  IMPLICIT NONE
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !Projeto de inversão magnética
  !Criador: Bijani,R.S.
  !Colaborador: Carreira,V.R.
  !Este programa tem o reproduzir a rotina publicada em 1976 por Plouff para
  !inversão de dados magnéticos
  !Categoria: Inversão de dados magnéticos
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			!####################### LISTA DE VARIÁVEIS ###################!
      !H - Magnitude do campo magnético da Terra                              !
      !T - Vetor campo magnético total                                        !
      !x - componente horizontal do campo magnético                           !
			!y - componente horizontal do campo magnético   	                  !
			!z - componente vertical do campo magnético         	               !
			!r - distância do centro de coordenadas até um infinitésimo da fonte !
      !Rk -                                                          !
      !dx -                                                          !
      !dy -                                                          !
      !dz -                                                          !
      !J -                                                           !
      !V -                                                           !
      !F -                                                           !
      !P -                                                           !
      !W -                                                           !
      !l -                                                           !
      !m -                                                           !
      !n -                                                           !
			!##############################################################!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$$$$$$$$$$$$$$$$$$$$$ DECLARAÇÃO DAS VARIÁVEIS GLOBAIS $$$$$$$$$$$$$$$$$$$$$!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
  INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=8, r=8)
  INTEGER, PARAMETER:: DBL = SELECTED_REAL_KIND(p=4, r=4)
  INTEGER(KIND=DBL):: i, j, k
  REAL(KIND=SGL):: H, J, V, F, P, W, l, m, n, r, Rk, dx, dy, dz
  REAL(KIND=SGL):: inicial, final, custocomputacional
  REAL(KIND=SGL), ALLOCATABLE, DIMENSION(:):: x, y, z
  REAL(KIND=SGL), ALLOCATABLE, DIMENSION(:,:,:):: T


    !ALLOCATE()


	CALL cpu_time(inicial)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$ INÍCIO DO PROGRAMA $$$$$$$$$$$$$$$$$$$$$$$$$$$$!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  CALL cpu_time(final)
  custocomputacional=final-inicial
  PRINT*, 'Custo Computacional=',custocomputacional, 'segundos'
  PRINT*,'*********************** FIM ***************************'
END PROGRAM plouff1976
