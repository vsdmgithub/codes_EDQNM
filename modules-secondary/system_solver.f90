! <f
! --------------------------------------------------------------
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
! CODE BY:
! --------   |         |   ---------        /\        |\      |
! |          |         |  |                /  \       | \     |
! |          |         |  |               /    \      |  \    |
! --------   |         |  |   ------|    /------\     |   \   |
!         |  |         |  |         |   /        \    |    \  |
!         |  |         |  |         |  /          \   |     \ |
! ---------   ----------  ----------  /            \  |      \|
! --------------------------------------------------------------
! #########################
! MODULE: system_solver
! LAST MODIFIED: 21 June 2022
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER FOR  EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_solver
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and updates it by a step, using the subroutines
! 1. rk4_algorithm
! 2. time_derivative
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_basicfunctions

  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT NONE
  ! _________________________
  ! SOLVER ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::d_spec1, d_spec2, d_spec3, d_spec4
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec_temp,integrating_factor
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  CONTAINS
! </f>

	SUBROUTINE allocate_solver_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       allocate all solver arrays
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y     A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( d_spec1( N ), d_spec2( N ), d_spec3( N ), d_spec4( N ))
    ALLOCATE( spec_temp( N ) )

    IF ( viscosity .GT. tol_float ) THEN

	    ALLOCATE( integrating_factor( N ) )
	    integrating_factor = DEXP( - viscosity * laplacian_k * dt )
			! INTEGRATION FACTOR TO INCLUDE VISCOSITY EFFECTS

		END IF

  END
! </f>

	SUBROUTINE rk4_algorithm
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'spec(k,t)-> spec(k,t+1)'
  ! Alg: - Runga kutta 4th order
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'spec(k)''

    spec_temp = spec
    CALL time_derivative(d_spec1) ! This call provides the time derivative for the spectral

  	spec      = spec_temp + hf * d_spec1
    CALL time_derivative(d_spec2)

    spec      = spec_temp + hf * d_spec2
    CALL time_derivative(d_spec3)

		spec      = spec_temp + d_spec3
    CALL time_derivative(d_spec4)

    ! Final increment for 'v(k)'
    IF ( viscosity .GT. tol_float ) THEN
			spec      = ( spec_temp + ( d_spec1 + two * d_spec2 + two * d_spec3 + d_spec4 ) / six ) * integrating_factor
		ELSE
			spec      =   spec_temp + ( d_spec1 + two * d_spec2 + two * d_spec3 + d_spec4 ) / six
		END IF

  END
! </f>

  SUBROUTINE time_derivative(d_spec)
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the time derivative matrix for matrix 'spec(k_ind)'
  ! This is the EDQNM EQUATION implemented for numerical computation
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! TRANSFER  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION ,DIMENSION( N ),INTENT(OUT)::d_spec

    CALL transfer_term
    ! REF-> <<< system_basicfunctions >>>

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   E   D   Q   N   M          E   Q   N.
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    IF ( ( viscosity .GT. tol_float ) .AND. ( forcing_status .EQ. 1 ) ) THEN

	    d_spec  =   dt * ( transfer_spec + forcer )
	    ! The transfer term and the forcing term

		ELSE

	    d_spec  =   dt * transfer_spec
	    ! The transfer term

		END IF

	END
! </f>

	SUBROUTINE deallocate_solver_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       deallocate all solver arrays
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y     D  E   A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(d_spec1, d_spec2, d_spec3, d_spec4)
    DEALLOCATE(spec_temp)

    IF ( viscosity .GT. tol_float ) THEN
	    DEALLOCATE(integrating_factor)
		END IF

  END
! </f>

END MODULE system_solver
