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
	USE system_basicfunctions

	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT NONE
	! _________________________
	! SOLVER ARRAYS
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::d_spec1, d_spec2, d_spec3, d_spec4
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec_temp,integ_factor
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

		IF ( visc_status .EQ. 1 ) THEN

			ALLOCATE( integ_factor( N ) )
			integ_factor = DEXP( - two * visc * laplacian_k * dt )
			! INTEGRATION FACTOR TO INCLUDE VISCOSITY EFFECTS

		END IF

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

		spec_temp = en_spec
		CALL time_derivative(d_spec1) ! This call provides the time derivative for the spectral

		en_spec      = spec_temp + hf * d_spec1
		CALL time_derivative(d_spec2)

		en_spec      = spec_temp + hf * d_spec2
		CALL time_derivative(d_spec3)

		en_spec      = spec_temp + d_spec3
		CALL time_derivative(d_spec4)

    ! Final increment for 'v(k)'
    IF ( visc_status .EQ. 1 ) THEN
			en_spec      = ( spec_temp + ( d_spec1 + two * d_spec2 + two * d_spec3 + d_spec4 ) / six ) * integ_factor
		ELSE
			en_spec      =   spec_temp + ( d_spec1 + two * d_spec2 + two * d_spec3 + d_spec4 ) / six
		END IF

	END
	! </f>

	SUBROUTINE time_derivative(d_spec)
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'en_spec(k_ind)'
	! This is the EDQNM EQUATION implemented for numerical computation
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION ,DIMENSION( N ),INTENT(OUT)::d_spec

		CALL compute_transfer_term
		! REF-> <<< system_basicfunctions >>>

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   E   D   Q   N   M          E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		IF ( forc_status .EQ. 1 ) THEN

			d_spec = dt * ( tr_spec + fr_spec )
			! The transfer term and the forcing term

		ELSE

			d_spec = dt * tr_spec
			! The transfer term

		END IF

	END
	! </f>

	SUBROUTINE compute_transfer_term
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the transfer term for every k
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   T   R   A   N   S   F   E   R       T   E   R   M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		tr_spec   =   zero
		! Reseting the transfer term

		DO k_ind = 1 , N
		DO q_ind = 1 , N
		DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

			CALL transfer_term_integrand

			tr_spec( k_ind ) = tr_spec( k_ind ) + integrand
			! Summation terms over all possible q,p for a given k.

	END DO
	END DO
	END DO
	END
! </f>

	SUBROUTINE transfer_term_integrand
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the integrand for a given k,q,p
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION:: p_d,k_d,E_k,E_p,E_q

		CALL eddy_damping_factor
		! Eddy damping term for the triad

		p_d = wno( p_ind )**( dim - one )
		k_d = wno( k_ind )**( dim - one )
		E_k = en_spec( k_ind )
		E_p = en_spec( p_ind )
		E_q = en_spec( q_ind )
		integrand  =  ( k_d * E_p - p_d * E_k ) * E_q
		integrand  =  integrand * geom_B( k_ind, q_ind, p_ind ) * eddy_damping
		integrand  =  integrand * weightage(k_ind, q_ind, p_ind)
		IF (integrand .GT. 1E5) THEN
			print*,"NaN here ",integrand,weightage(k_ind,q_ind,p_ind),k_ind,q_ind,p_ind
			nan_status=1
		END IF

	END
! </f>

	SUBROUTINE eddy_damping_factor
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the damping factor for the third wnoent
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   E  D  D  Y            F  R  E  Q  U  E  N  C  Y
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		eddy_array   = laplacian_k * en_spec * wno_band
		eddy_k       = SUM( eddy_array( : k_ind) ) ** hf
		eddy_q       = SUM( eddy_array( : q_ind) ) ** hf
		eddy_p       = SUM( eddy_array( : p_ind) ) ** hf
		eddy_freq    = eddy_const * ( eddy_k + eddy_q + eddy_p )

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! D A M P I N G       F A C T O R
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		visc_freq    = laplacian_k( k_ind ) + laplacian_k( q_ind ) + laplacian_k( p_ind )
		visc_freq    = visc * visc_freq
		eddy_damping = one - DEXP( -( eddy_freq + visc_freq ) * time_now )
		eddy_damping = eddy_damping / ( eddy_freq + visc_freq )

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
		IF ( visc_status .EQ. 1 ) THEN
			DEALLOCATE(integ_factor)
		END IF

	END
! </f>

END MODULE system_solver
