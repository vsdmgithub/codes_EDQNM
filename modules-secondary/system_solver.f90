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
! MODULE NAME: system_solver
! LAST MODIFIED: 15 NOV 2022
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
! *. rk4_algorithm
! *. time_derivative
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_solver_eqns

	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT NONE
	! _________________________
	! SOLVER ARRAYS
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dV_1, dV_2, dV_3, dV_4, dV
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dB_1, dB_2, dB_3, dB_4, dB
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
		ALLOCATE( dV_1( N ), dV_2( N ), dV_3( N ), dV_4( N ))
		ALLOCATE( dB_1( N ), dB_2( N ), dB_3( N ), dB_4( N ))
		ALLOCATE( dV( N ), dB( N ) )
		ALLOCATE( spec_temp_V( N ) )
		ALLOCATE( spec_temp_B( N ) )

		IF ( visc_status .EQ. 1 ) THEN

			ALLOCATE( integ_factor_V( N ) )
			integ_factor_V = DEXP( - two * visc * laplacian_k * dt_cur )
			! INTEGRATION FACTOR TO INCLUDE VISCOSITY EFFECTS

		END IF
		IF ( diff_status .EQ. 1 ) THEN

			ALLOCATE( integ_factor_B( N ) )
			integ_factor_B = DEXP( - two * diff * laplacian_k * dt_cur )
			! INTEGRATION FACTOR TO INCLUDE DIFFUSIVITY EFFECTS

		END IF

	END
	! </f>

	SUBROUTINE rk4_algorithm_V
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'en_spec(k,t)-> en_spec(k,t+1)'
	! Alg: - Runga kutta 4th order
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'en_spec(k)''

		spec_temp_V = en_spec_V

		CALL time_derivative_V(dV_1) ! This call provides the time derivative for the spectral
		en_spec_V      = spec_temp_V + hf * dV_1

		CALL time_derivative_V(dV_2)
		en_spec_V      = spec_temp_V + hf * dV_2

		CALL time_derivative_V(dV_3)
		en_spec_V      = spec_temp_V + dV_3

		CALL time_derivative_V(dV_4)

  ! Final increment for 'v(k)'
    IF ( visc_status .EQ. 1 ) THEN
			en_spec_V      = ( spec_temp_V + ( dV_1 + two * dV_2 + two * dV_3 + dV_4 ) / six ) * integ_factor_V
		ELSE
			en_spec_V      =   spec_temp_V + ( dV_1 + two * dV_2 + two * dV_3 + dV_4 ) / six
		END IF

		DO dum_ind = 0, N
			IF ( en_spec_V( dum_ind ) .LT. zero ) THEN
				print*,"Kinetic spectrum got negative at k=",dum_ind
				nan_status = 1
				EXIT
			END IF
		END DO

	END
	! </f>

	SUBROUTINE rk4_algorithm_B
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'en_spec(k,t)-> en_spec(k,t+1)'
	! Alg: - Runga kutta 4th order
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'en_spec(k)''

		spec_temp_B = en_spec_B

		CALL time_derivative_B(dB_1) ! This call provides the time derivative for the spectral
		en_spec_B   = spec_temp_B + hf * dB_1

		CALL time_derivative_B(dB_2)
		en_spec_B   = spec_temp_B + hf * dB_2

		CALL time_derivative_B(dB_3)
		en_spec_B   = spec_temp_B + dV_3

		CALL time_derivative_B(dB_4)

    ! Final increment for 'B(k)'
    IF ( diff_status .EQ. 1 ) THEN
			en_spec_B = ( spec_temp_B + ( dB_1 + two * dB_2 + two * dB_3 + dB_4 ) / six ) * integ_factor_B
		ELSE
			en_spec_B = spec_temp_B + ( dB_1 + two * dB_2 + two * dB_3 + dB_4 ) / six
		END IF

		DO dum_ind = 0, N
			IF ( en_spec_B( dum_ind ) .LT. zero ) THEN
				print*,"Magnetic spectrum got negative at k=",dum_ind
				nan_status = 1
				EXIT
			END IF
		END DO

	END
	! </f>

	SUBROUTINE rk4_algorithm
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'en_spec(k,t)-> en_spec(k,t+1)'
	! Alg: - Runga kutta 4th order
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'en_spec(k)''


		spec_temp_V = en_spec_V
		spec_temp_B = en_spec_B

		CALL time_derivative_B(dB_1)
		CALL time_derivative_V(dV_1)
		en_spec_V   = spec_temp_V + hf * dV_1
		en_spec_B   = spec_temp_B + hf * dB_1

		CALL time_derivative_B(dB_2)
		CALL time_derivative_V(dV_2)
		en_spec_V   = spec_temp_V + hf * dV_2
		en_spec_B   = spec_temp_B + hf * dB_2

		CALL time_derivative_B(dB_3)
		CALL time_derivative_V(dV_3)
		en_spec_V   = spec_temp_V + dV_3
		en_spec_B   = spec_temp_B + dV_3

		CALL time_derivative_B(dB_4)
		CALL time_derivative_V(dV_4)

    ! Final increment for 'v(k)'
    IF ( diff_status .EQ. 1 ) THEN
			en_spec_B = ( spec_temp_B + ( dB_1 + two * dB_2 + two * dB_3 + dB_4 ) / six ) * integ_factor_B
		ELSE
			en_spec_B = spec_temp_B + ( dB_1 + two * dB_2 + two * dB_3 + dB_4 ) / six
		END IF

    ! Final increment for 'v(k)'
    IF ( visc_status .EQ. 1 ) THEN
			en_spec_V = ( spec_temp_V + ( dV_1 + two * dV_2 + two * dV_3 + dV_4 ) / six ) * integ_factor_V
		ELSE
			en_spec_V = spec_temp_V + ( dV_1 + two * dV_2 + two * dV_3 + dV_4 ) / six
		END IF

		DO dum_ind = 0, N
			IF ( en_spec_V( dum_ind ) .LT. zero ) THEN
				print*,"Kinetic spectrum got negative at k=",dum_ind
				en_spec_V( dum_ind ) = spec_temp_V( dum_ind )
				nan_status = 1
				EXIT
			END IF
			IF ( en_spec_B( dum_ind ) .LT. zero ) THEN
				print*,"Magnetic spectrum got negative at k=",dum_ind
				en_spec_B( dum_ind ) = spec_temp_B( dum_ind )
				nan_status = 1
				EXIT
			END IF
		END DO

	END
	! </f>

	SUBROUTINE ab4_algorithm
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE Adam Bashforth (predictor corrector) algorithm to move one step forward in time for the matrix 'en_spec(k,t)-> en_spec(k,t+1)'
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'en_spec(k)''
		IF ( t_step .EQ. 0 ) THEN
			CALL time_derivative_B(dB)
			CALL time_derivative_V(dV)
			CALL rk4_algorithm
			dV_4 = dV
			dB_4 = dB
		ELSE IF ( t_step .EQ. 1 ) THEN
			CALL time_derivative_B(dB)
			CALL time_derivative_V(dV)
			CALL rk4_algorithm
			dV_3 = dV
			dB_3 = dB
		ELSE IF ( t_step .EQ. 2 ) THEN
			CALL time_derivative_B(dB)
			CALL time_derivative_V(dV)
			CALL rk4_algorithm
			dV_2 = dV
			dB_2 = dB
		ELSE

			spec_temp_V = en_spec_V
			spec_temp_B = en_spec_B

			CALL time_derivative_B(dB_1)
			CALL time_derivative_V(dV_1)
	    IF ( visc_status .EQ. 1 ) THEN
				en_spec_V = integ_factor_V * ( spec_temp_V + ( - 9.0D0 * dV_4 + 37.0D0 * dV_3 - 59.0D0 * dV_2 + 55.0D0 * dV_1 ) / 24.0D0 )
			ELSE
				en_spec_V = spec_temp_V + ( - 9.0D0 * dV_4 + 37.0D0 * dV_3 - 59.0D0 * dV_2 + 55.0D0 * dV_1 ) / 24.0D0
			END IF
	    IF ( diff_status .EQ. 1 ) THEN
				en_spec_B = integ_factor_B * ( spec_temp_B + ( - 9.0D0 * dB_4 + 37.0D0 * dB_3 - 59.0D0 * dB_2 + 55.0D0 * dB_1 ) / 24.0D0 )
			ELSE
				en_spec_B = spec_temp_B + ( - 9.0D0 * dB_4 + 37.0D0 * dB_3 - 59.0D0 * dB_2 + 55.0D0 * dB_1 ) / 24.0D0
			END IF

			CALL time_derivative_V(dV_4)
			CALL time_derivative_B(dB_4)

	    IF ( visc_status .EQ. 1 ) THEN
				en_spec_V =  integ_factor_V * ( spec_temp_V + ( dV_3 - 5.0D0 * dV_2 + 19.0D0 * dV_1 + 9.0D0 * dV_4 ) / 24.0D0 )
			ELSE
				en_spec_V = spec_temp_V + ( dV_3 - 5.0D0 * dV_2 + 19.0D0 * dV_1 + 9.0D0 * dV_4 ) / 24.0D0
			END IF
	    IF ( diff_status .EQ. 1 ) THEN
				en_spec_B = integ_factor_B * ( spec_temp_B + ( dB_3 - 5.0D0 * dB_2 + 19.0D0 * dB_1 + 9.0D0 * dB_4 ) / 24.0D0 )
			ELSE
				en_spec_B = spec_temp_B + ( dB_3 - 5.0D0 * dB_2 + 19.0D0 * dB_1 + 9.0D0 * dB_4 ) / 24.0D0
			END IF

			dV_4 = dV_3
			dV_3 = dV_2
			dV_2 = dV_1

			dB_4 = dB_3
			dB_3 = dB_2
			dB_2 = dB_1

		END IF

		DO dum_ind = 0, kD_V_ind
			IF ( en_spec_V( dum_ind ) .LT. zero ) THEN
				print*,"Kinetic spectrum got negative at k=",dum_ind
				en_spec_V( dum_ind ) = spec_temp_V( dum_ind )
				nan_status = 1
				EXIT
			END IF
		END DO
		DO dum_ind = kD_V_ind, N
			IF ( en_spec_V( dum_ind ) .LT. zero ) THEN
				en_spec_V( dum_ind ) = spec_temp_V( dum_ind )
			END IF
		END DO
		DO dum_ind = 0, kD_B_ind
			IF ( en_spec_B( dum_ind ) .LT. zero ) THEN
				print*,"Magnetic spectrum got negative at k=",dum_ind
				en_spec_B( dum_ind ) = spec_temp_B( dum_ind )
				nan_status = 1
				EXIT
			END IF
		END DO
		DO dum_ind = kD_B_ind, N
			IF ( en_spec_B( dum_ind ) .LT. zero ) THEN
				en_spec_B( dum_ind ) = spec_temp_B( dum_ind )
			END IF
		END DO

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
		DEALLOCATE(dV_1, dV_2, dV_3, dV_4)
		DEALLOCATE(dB_1, dB_2, dB_3, dB_4)
		DEALLOCATE(dB, dV)
		DEALLOCATE(spec_temp_V,spec_temp_B)
		IF ( visc_status .EQ. 1 ) THEN
			DEALLOCATE(integ_factor_V)
		END IF
		IF ( diff_status .EQ. 1 ) THEN
			DEALLOCATE(integ_factor_B)
		END IF

	END
! </f>

END MODULE system_solver
