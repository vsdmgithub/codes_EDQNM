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
	! LOCAL  VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION:: p_d,k_d
	DOUBLE PRECISION:: k_pq,p_kq
	DOUBLE PRECISION:: E_V_k,E_V_p,E_V_q
	DOUBLE PRECISION:: E_B_k,E_B_p,E_B_q
	DOUBLE PRECISION:: visc_freq,diff_freq,eddy_freq,alfven_freq
	DOUBLE PRECISION:: alfven_k,alfven_p,alfven_q
	DOUBLE PRECISION:: eddy_V_k,eddy_V_q,eddy_V_p
	DOUBLE PRECISION:: eddy_B_k,eddy_B_q,eddy_B_p
	DOUBLE PRECISION:: eddy_damping_V,eddy_damping_Bk
	DOUBLE PRECISION:: eddy_damping_Bq,eddy_damping_Bp
	! _________________________
	! SOLVER ARRAYS
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dV_spec1, dV_spec2, dV_spec3, dV_spec4
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dB_spec1, dB_spec2, dB_spec3, dB_spec4
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec_temp_V,spec_temp_B
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
		ALLOCATE( dV_spec1( N ), dV_spec2( N ), dV_spec3( N ), dV_spec4( N ))
		ALLOCATE( dB_spec1( N ), dB_spec2( N ), dB_spec3( N ), dB_spec4( N ))
		ALLOCATE( spec_temp_V( N ) )
		ALLOCATE( spec_temp_B( N ) )

		IF ( visc_status .EQ. 1 ) THEN

			ALLOCATE( integ_factor_V( N ) )
			integ_factor_V = DEXP( - two * visc * laplacian_k * dt )
			! INTEGRATION FACTOR TO INCLUDE VISCOSITY EFFECTS

		END IF
		IF ( diff_status .EQ. 1 ) THEN

			ALLOCATE( integ_factor_B( N ) )
			integ_factor_B = DEXP( - two * diff * laplacian_k * dt )
			! INTEGRATION FACTOR TO INCLUDE DIFFUSIVITY EFFECTS

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

		spec_temp_V = en_spec_V
		spec_temp_B = en_spec_B

		CALL time_derivative_V(dV_spec1) ! This call provides the time derivative for the spectral
		en_spec_V      = spec_temp_V + hf * dV_spec1
		CALL time_derivative_B(dB_spec1) ! This call provides the time derivative for the spectral
		en_spec_B   = spec_temp_B + hf * dB_spec1

		CALL time_derivative_V(dV_spec2)
		en_spec_V      = spec_temp_V + hf * dV_spec2
		CALL time_derivative_B(dB_spec2)
		en_spec_B   = spec_temp_B + hf * dB_spec2

		CALL time_derivative_V(dV_spec3)
		en_spec_V      = spec_temp_V + dV_spec3
		CALL time_derivative_B(dB_spec3)
		en_spec_B   = spec_temp_B + dV_spec3

		CALL time_derivative_V(dV_spec4)
		CALL time_derivative_B(dB_spec4)

    ! Final increment for 'v(k)'
    IF ( diff_status .EQ. 1 ) THEN
			en_spec_B = ( spec_temp_B + ( dB_spec1 + two * dB_spec2 + two * dB_spec3 + dB_spec4 ) / six ) * integ_factor_B
		ELSE
			en_spec_B = spec_temp_B + ( dB_spec1 + two * dB_spec2 + two * dB_spec3 + dB_spec4 ) / six
		END IF

    ! Final increment for 'v(k)'
    IF ( visc_status .EQ. 1 ) THEN
			en_spec_V      = ( spec_temp_V + ( dV_spec1 + two * dV_spec2 + two * dV_spec3 + dV_spec4 ) / six ) * integ_factor_V
		ELSE
			en_spec_V      =   spec_temp_V + ( dV_spec1 + two * dV_spec2 + two * dV_spec3 + dV_spec4 ) / six
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

		CALL time_derivative_V(dV_spec1) ! This call provides the time derivative for the spectral
		en_spec_V      = spec_temp_V + hf * dV_spec1

		CALL time_derivative_V(dV_spec2)
		en_spec_V      = spec_temp_V + hf * dV_spec2

		CALL time_derivative_V(dV_spec3)
		en_spec_V      = spec_temp_V + dV_spec3

		CALL time_derivative_V(dV_spec4)

  ! Final increment for 'v(k)'
    IF ( visc_status .EQ. 1 ) THEN
			en_spec_V      = ( spec_temp_V + ( dV_spec1 + two * dV_spec2 + two * dV_spec3 + dV_spec4 ) / six ) * integ_factor_V
		ELSE
			en_spec_V      =   spec_temp_V + ( dV_spec1 + two * dV_spec2 + two * dV_spec3 + dV_spec4 ) / six
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

			d_spec = dt * ( tr_spec_V + fr_spec )
			! The transfer term and the forcing term

		ELSE

			d_spec = dt * tr_spec_V
			! The transfer term

		END IF

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

		CALL time_derivative_B(dB_spec1) ! This call provides the time derivative for the spectral
		en_spec_B   = spec_temp_B + hf * dB_spec1

		CALL time_derivative_B(dB_spec2)
		en_spec_B   = spec_temp_B + hf * dB_spec2

		CALL time_derivative_B(dB_spec3)
		en_spec_B   = spec_temp_B + dV_spec3

		CALL time_derivative_B(dB_spec4)

    ! Final increment for 'B(k)'
    IF ( diff_status .EQ. 1 ) THEN
			en_spec_B = ( spec_temp_B + ( dB_spec1 + two * dB_spec2 + two * dB_spec3 + dB_spec4 ) / six ) * integ_factor_B
		ELSE
			en_spec_B = spec_temp_B + ( dB_spec1 + two * dB_spec2 + two * dB_spec3 + dB_spec4 ) / six
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
			d_spec = dt * tr_spec_B
			! The transfer term

		! DO dum_ind = 0, N
		! 	IF ( spec_temp_B( dum_ind ) .LT. 6.0D0 * d_spec( dum_ind ) ) THEN
		! 		! print*,"Magnetic spectrum got negative at k=",dum_ind
		! 		d_spec( dum_ind ) = zero
		! 	END IF
		! END DO

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
		DOUBLE PRECISION::dumI,dumS
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

		DO dum_ind = 0, N
			dumS = dt * tr_spec_B_self( dum_ind )
			dumI = dt * tr_spec_B_intr( dum_ind )
			IF ( spec_temp_B( dum_ind ) + dumS + dumI .LT. zero ) then
				print*,"Magnetic spectrum will become negative at k=",dum_ind,dumS,dumI,spec_temp_B(dum_ind)
				print*,growth_factor
		! 		d_spec( dum_ind ) = zero
			END IF
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

		CALL eddy_damping_factor
		! Eddy damping term for the triad

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
			print*,"NaN here at V-transfer term calc at ",k_ind,q_ind,p_ind,&
			triad_weightage(k_ind,q_ind, p_ind),eddy_damping,integrand_V_self,integrand_V_intr
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

		CALL eddy_damping_factor
		! Eddy damping term for the triad

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

		CALL eddy_damping_factor
		! Eddy damping term for the triad

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
		eddy_V_k          = DSQRT( DABS( SUM( en_spec_V( :k_ind ) * laplacian_k( :k_ind ) * wno_band( :k_ind ) ) ) )
		eddy_V_p          = DSQRT( DABS( SUM( en_spec_V( :p_ind ) * laplacian_k( :p_ind ) * wno_band( :p_ind ) ) ) )
		eddy_V_q          = DSQRT( DABS( SUM( en_spec_V( :q_ind ) * laplacian_k( :q_ind ) * wno_band( :q_ind ) ) ) )
		eddy_freq         = eddy_const * ( eddy_V_k + eddy_V_q + eddy_V_p )
		visc_freq         = visc * ( laplacian_k( k_ind ) + laplacian_k( q_ind ) + laplacian_k( p_ind ) )
		eddy_damping_V    = ( one - DEXP( - time_now * ( visc_freq + eddy_freq ) ) ) / ( visc_freq + eddy_freq )

		IF ( coupling_status .GT. 0 ) THEN

			eddy_B_k        = DSQRT( DABS( SUM( en_spec_B( :k_ind ) * laplacian_k( :k_ind ) * wno_band( :k_ind ) ) ) )
			eddy_B_p        = DSQRT( DABS( SUM( en_spec_B( :p_ind ) * laplacian_k( :p_ind ) * wno_band( :p_ind ) ) ) )
			eddy_B_q        = DSQRT( DABS( SUM( en_spec_B( :q_ind ) * laplacian_k( :q_ind ) * wno_band( :q_ind ) ) ) )
			alfven_k        = wno( k_ind ) * DSQRT( DABS( SUM( en_spec_B( :k_ind ) * wno_band( :k_ind ) ) ) )
			alfven_p        = wno( p_ind ) * DSQRT( DABS( SUM( en_spec_B( :p_ind ) * wno_band( :p_ind ) ) ) )
			alfven_q        = wno( q_ind ) * DSQRT( DABS( SUM( en_spec_B( :q_ind ) * wno_band( :q_ind ) ) ) )

			visc_freq       = visc * laplacian_k( k_ind )
			diff_freq       = diff * ( laplacian_k( q_ind ) + laplacian_k( p_ind ) )
			eddy_freq       = eddy_const * ( eddy_V_k + eddy_B_p + eddy_B_q )
			alfven_freq     = alfven_const * ( alfven_p + alfven_q )
			eddy_damping_Bk = ( one - DEXP( - time_now * ( visc_freq + diff_freq + eddy_freq + alfven_freq ) ) ) &
			                   / ( visc_freq + diff_freq + eddy_freq + alfven_freq )

			visc_freq       = visc * laplacian_k( p_ind )
			diff_freq       = diff * ( laplacian_k( k_ind ) + laplacian_k( q_ind ) )
			eddy_freq       = eddy_const * ( eddy_V_p + eddy_B_k + eddy_B_q )
			alfven_freq     = alfven_const * ( alfven_k + alfven_q )
			eddy_damping_Bp = ( one - DEXP( - time_now * ( visc_freq + diff_freq + eddy_freq + alfven_freq ) ) ) &
			                   / ( visc_freq + diff_freq + eddy_freq + alfven_freq )

			visc_freq       = visc * laplacian_k( q_ind )
			diff_freq       = diff * ( laplacian_k( k_ind ) + laplacian_k( p_ind ) )
			eddy_freq       = eddy_const * ( eddy_V_q + eddy_B_k + eddy_B_p )
			alfven_freq     = alfven_const * ( alfven_k + alfven_p )
			eddy_damping_Bq = ( one - DEXP( - time_now * ( visc_freq + diff_freq + eddy_freq + alfven_freq ) ) ) &
			                   / ( visc_freq + diff_freq + eddy_freq + alfven_freq )

		END IF
		eddy_damping      = 4.0D0 /( ( one / eddy_damping_V ) + ( one / eddy_damping_Bk ) + ( one / eddy_damping_Bp ) &
		 												 + ( one / eddy_damping_Bq ) )
		! eddy_damping_Bk   = eddy_damping
		! eddy_damping_Bp   = eddy_damping
		! eddy_damping_Bq   = eddy_damping
		! eddy_damping_V    = eddy_damping

		IF ( eddy_damping .NE. eddy_damping ) THEN
			print*,"NaN here at eddy damping calc ",k_ind,q_ind,p_ind,eddy_damping
			nan_status=1
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
		DEALLOCATE(dV_spec1, dV_spec2, dV_spec3, dV_spec4)
		DEALLOCATE(dB_spec1, dB_spec2, dB_spec3, dB_spec4)
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
