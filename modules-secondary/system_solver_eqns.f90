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

MODULE system_solver_eqns
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Gives the time-derivative at a given time.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_basicfunctions

	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT NONE
	! _________________________
	! LOCAL  VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION:: p_d,k_d
	DOUBLE PRECISION:: k_pq,p_kq
	DOUBLE PRECISION:: E_V_k,E_V_p,E_V_q
	DOUBLE PRECISION:: E_B_k,E_B_p,E_B_q
	DOUBLE PRECISION:: eddy_freq
	DOUBLE PRECISION:: eddy_damping_V,eddy_damping_Bk
	DOUBLE PRECISION:: eddy_damping_Bq,eddy_damping_Bp
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec_temp_V,spec_temp_B
	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	CONTAINS
	! </f>
	SUBROUTINE time_derivative_V(d_spec)
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

		CALL compute_transfer_term_V
		! REF-> <<< system_basicfunctions >>>

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   E   D   Q   N   M          E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		IF ( forc_status .EQ. 1 ) THEN

			d_spec = dt_cur * ( tr_spec_V + fr_spec )
			! The transfer term and the forcing term

		ELSE

			d_spec = dt_cur * tr_spec_V
			! The transfer term

		END IF

	END
	! </f>

	SUBROUTINE time_derivative_B(d_spec)
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

		CALL compute_transfer_term_B
		! REF-> <<< system_basicfunctions >>>

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  M  H  D          E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
			d_spec = dt_cur * tr_spec_B
			! The transfer term

	END
	! </f>

	SUBROUTINE compute_transfer_term_V
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

		tr_spec_V_self   =   zero
		tr_spec_V_intr   =   zero
		! Reseting the transfer term

		DO k_ind = 1 , N
		DO q_ind = 1, N
		DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

			CALL transfer_term_integrand_V

			IF (nan_status .EQ. 1) THEN
				! EXIT
				GOTO 911
			END IF

			tr_spec_V_self( k_ind ) = tr_spec_V_self( k_ind ) + integrand_V_self
			IF ( coupling_status .GE. 1 ) THEN
				tr_spec_V_intr( k_ind ) = tr_spec_V_intr( k_ind ) + integrand_V_intr
				! Summation terms over all possible q,p for a given k.
			END IF

		END DO
		END DO
		END DO

		tr_spec_V = tr_spec_V_self + tr_spec_V_intr
		911 CONTINUE

	END
! </f>

	SUBROUTINE compute_transfer_term_B
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

		tr_spec_B_self   =   zero
		tr_spec_B_intr   =   zero
		! Reseting the transfer term

		DO k_ind = 1 , N
		DO q_ind = 1 , N
		DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

			CALL transfer_term_integrand_B

			IF (nan_status .EQ. 1) THEN
				! EXIT
				GOTO 911
			END IF

			tr_spec_B_self( k_ind ) = tr_spec_B_self( k_ind ) + integrand_B_self
			tr_spec_B_intr( k_ind ) = tr_spec_B_intr( k_ind ) + integrand_B_intr
			! Summation terms over all possible q,p for a given k.

		END DO
		END DO
		END DO

		tr_spec_B = tr_spec_B_self + tr_spec_B_intr
		911 CONTINUE

	END
! </f>

	SUBROUTINE transfer_term_integrand_V
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the integrand for a given k,q,p
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		eddy_freq       = eddy_V( k_ind ) + eddy_V( q_ind ) + eddy_V( p_ind )
		eddy_damping_V  = ( one - DEXP( - time_now * eddy_freq ) ) / eddy_freq
		eddy_freq       = eddy_V( k_ind ) + eddy_B( q_ind ) + eddy_B( p_ind )
		eddy_damping_Bk = ( one - DEXP( - time_now * eddy_freq ) ) / eddy_freq

		p_d                = wno( p_ind )**( dim - one )
		k_d                = wno( k_ind )**( dim - one )
		k_pq               = wno( k_ind ) / ( wno( p_ind ) * wno( q_ind ) )
		E_V_k              = en_spec_V( k_ind )
		E_V_p              = en_spec_V( p_ind )
		E_V_q              = en_spec_V( q_ind )
		E_B_p              = en_spec_B( p_ind )
		E_B_q              = en_spec_B( q_ind )

		integrand_V_self   = k_pq * geom_b( k_ind, q_ind, p_ind ) * ( k_d * E_V_p - p_d * E_V_k ) * E_V_q
		integrand_V_self   = eddy_damping_V * triad_weightage( k_ind, q_ind, p_ind ) * integrand_V_self
		integrand_V_self   = integrand_V_self * wno_band( q_ind ) * wno_band( p_ind )

		integrand_V_intr   = zero

		IF ( coupling_status .EQ. 1 ) THEN
			! integrand_V_intr = + k_pq * geom_b( k_ind, q_ind, p_ind ) * k_d * E_B_p * E_B_q ! Nonlinear term
			! integrand_V_intr = - k_pq * geom_c( k_ind, q_ind, p_ind ) * p_d * E_V_k * E_B_q ! Linear term
			integrand_V_intr = k_pq * ( geom_b( k_ind, q_ind, p_ind ) * k_d * E_B_p - geom_c( k_ind, q_ind, p_ind ) * p_d * E_V_k ) * E_B_q ! Full
			integrand_V_intr = eddy_damping_Bk * triad_weightage( k_ind, q_ind, p_ind ) *  integrand_V_intr
			integrand_V_intr = integrand_V_intr * wno_band( q_ind ) * wno_band( p_ind )
		END IF

		IF ( integrand_V_self .NE. integrand_V_self ) THEN
			print*,"NaN here at V-transfer term calc at ",k_ind,q_ind,p_ind,integrand_V_self,integrand_V_intr
			nan_status=1
		END IF

	END
! </f>

	SUBROUTINE transfer_term_integrand_B
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the integrand for a given k,q,p
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		eddy_freq       = eddy_V( q_ind ) + eddy_B( k_ind ) + eddy_B( p_ind )
		eddy_damping_Bq = ( one - DEXP( - time_now * eddy_freq ) ) / eddy_freq
		eddy_freq       = eddy_V( p_ind ) + eddy_B( k_ind ) + eddy_B( q_ind )
		eddy_damping_Bp = ( one - DEXP( - time_now * eddy_freq ) ) / eddy_freq

		p_d              = wno( p_ind )**( dim - one )
		k_d              = wno( k_ind )**( dim - one )
		k_pq             = wno( k_ind ) / ( wno( p_ind ) * wno( q_ind ) )
		p_kq             = wno( p_ind ) / ( wno( k_ind ) * wno( q_ind ) )
		E_V_p            = en_spec_V( p_ind )
		E_V_q            = en_spec_V( q_ind )
		E_B_p            = en_spec_B( p_ind )
		E_B_q            = en_spec_B( q_ind )
		E_B_k            = en_spec_B( k_ind )

		! integrand_B_self = zero
		integrand_B_self = k_pq * geom_h( k_ind, q_ind, p_ind ) * ( k_d * E_B_p - p_d * E_B_k ) * E_V_q
		integrand_B_self = eddy_damping_Bq * triad_weightage( k_ind, q_ind, p_ind ) * integrand_B_self
		! integrand_B_self = integrand_B_self * wno_band( q_ind ) * wno_band( p_ind )
		integrand_B_self = integrand_B_self * wno_band( q_ind ) * wno_band( p_ind )
		! integrand_B_intr = zero

		! integrand_B_intr = - p_kq * geom_c( p_ind, q_ind, k_ind ) * p_d * E_B_k * E_B_q  ! Nonlinear term
		! integrand_B_intr = + p_kq * geom_c( p_ind, q_ind, k_ind ) * k_d * E_V_p * E_B_q  !  Linear terms
		integrand_B_intr = p_kq * geom_c( p_ind, q_ind, k_ind ) * ( k_d * E_V_p - p_d * E_B_k ) * E_B_q ! Full
		integrand_B_intr = eddy_damping_Bp * triad_weightage( k_ind, q_ind, p_ind ) * integrand_B_intr
		integrand_B_intr = integrand_B_intr * wno_band( q_ind ) * wno_band( p_ind )

		IF (integrand_B_intr .NE. integrand_B_intr) THEN
			print*,"NaN here at B-transfer term calc at ",k_ind,q_ind,p_ind,&
			triad_weightage(k_ind,q_ind, p_ind),integrand_B_self,integrand_B_intr
			nan_status=1
		END IF

	END
! </f>

	SUBROUTINE compute_adaptive_time_step
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO compute new time-step
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT  NONE
		DOUBLE PRECISION::dt_cur_tem

		dt_cur     = dt
		dt_cur_tem = dt
		DO dum_ind = 1, N
			IF ( DABS(tr_spec_B( dum_ind )) .GT. zero ) THEN
				cfl_cur = en_spec_B( dum_ind ) / ( dt_cur * DABS( tr_spec_B( dum_ind ) ) )
				IF ( cfl_cur .LE. cfl_sys ) THEN
					dt_cur_tem = en_spec_B( dum_ind ) / ( cfl_sys * DABS( tr_spec_B( dum_ind ) ) )
					IF ( dt_cur_tem .LT. dt_cur ) THEN
						dt_cur = dt_cur_tem
					END IF
				END IF
			END IF
		END DO

		IF ( dt .NE. dt_cur ) THEN
			dt_cur_tem = dt_cur
			CALL find_time_step( dt_cur_tem, dt_cur )
			PRINT*,'dt changed to ',dt_cur,'at t=',time_now
		END IF

	! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	END
! </f>

	SUBROUTINE dynamo_rate_calc
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the dynamo rate for each wno0
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  D Y N A M O   T E R M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		dyn_rate_spec   =   zero
		! Reseting the  term

		DO q_ind = 1 , N
		DO k_ind = 1 , N
		DO p_ind = p_ind_min( q_ind, k_ind ), p_ind_max( k_ind, q_ind )

		eddy_freq              = eddy_V( p_ind ) + eddy_B( k_ind ) + eddy_B( q_ind )
		eddy_damping_Bp        = ( one - DEXP( - time_now * eddy_freq ) ) / eddy_freq

		k_d                    = wno( k_ind )**( dim - one )
		p_kq                   = wno( p_ind ) / ( wno( k_ind ) * wno( q_ind ) )
		E_V_p                  = en_spec_V( p_ind )

		dyn_rate_integrand     = + p_kq * geom_c( p_ind, q_ind, k_ind ) * k_d * E_V_p
		dyn_rate_integrand     = eddy_damping_Bp * triad_weightage( k_ind, q_ind, p_ind ) * dyn_rate_integrand
		dyn_rate_integrand     = dyn_rate_integrand * wno_band( k_ind ) * wno_band( p_ind )

		dyn_rate_spec( q_ind ) = dyn_rate_spec( q_ind ) + dyn_rate_integrand

		END DO
		END DO
		END DO

	END
! </f>

END MODULE system_solver_eqns
