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
! MODULE NAME  : system_basicfunctions
! LAST MODIFIED: 09 OCT 2023
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MODULE CONTAINING SUBROUTINES, FUNCTIONS, PARAMETERS FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>
MODULE system_basicfunctions
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * All the system parameters, corresponding to the different variant of EDQNM equation we are studying
! is defined here. Note global variables is common amongst these different variants.
! * The EDQNM arrays are allocated and declared here.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_initialcondition
	USE system_basicoutput
	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT  NONE
	! _________________________
	! ARRAYS
	! !!!!!!!!!!!!!!!!!!!!!!!!!

	CONTAINS
! </f>

	SUBROUTINE allocate_edqnm_arrays
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	!       Allocate all arrays related to EDQNM algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  R  R  A  Y     A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE( p_ind_min( N, N ), p_ind_max( N, N ) )
		ALLOCATE( geom_b( N, N, N), geom_c( N, N, N ), geom_h( N, N, N ) )
		ALLOCATE( kqp_status( N, N, N ) )
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	END
! </f>

	SUBROUTINE init_edqnm_arrays
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to initialize all in-built arrays that need not change during the evolution of system.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION:: wno_p_min,wno_p_max
		DOUBLE PRECISION:: z_f, x_f, y_f
		DOUBLE PRECISION:: geo,rat


		triad_count   = 0
		kqp_status    = 0
		! Reseting to zero for safety

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  L I M I T S   O F    I N T E G R A T I O N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		DO k_ind = 1 , N         ! For each 'k'
		DO q_ind = 1 , N         ! First loop for 'q' from 0 till 'k'

			wno_p_min                           = DABS( wno( k_ind ) - wno( q_ind ) )
			wno_p_max                           = wno( k_ind ) + wno( q_ind )
			p_ind_min( k_ind, q_ind )           = upperindex( wno_p_min )
			p_ind_max( k_ind, q_ind )           = lowerindex( wno_p_max )
			! FINDING THE MIN/MAXIMUM OF 'p_ind' FOR EVERY PAIR OF 'k_ind,q_ind'

			DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

				IF ( triangle_compatibility( k_ind, q_ind, p_ind ) .EQ. 1 ) THEN
				! THIS HAS TO BE TRUE, OTHERWISE THERE IS AN AN ERROR IN THE UPPER AND LOWER INDEX DETERMINATION

					kqp_status( k_ind, q_ind, p_ind ) = 1
					x_f                               = + cosine( q_ind, p_ind, k_ind )
					y_f                               = + cosine( k_ind, q_ind, p_ind )
					z_f                               = + cosine( k_ind, p_ind, q_ind )
					! FINDING COSINES FOR ALL THREE SIDES ONCE IT IS APPROVED TO PARTICIPATE IN THE TRIAD INTERACTION

					! GEOMETRIC FACTORS DEFINED THROUGH X,Y,Z
					geo                               = wno( p_ind ) / wno( k_ind )
					geom_b( k_ind, q_ind, p_ind )     = geo * ( ( z_f ** thr ) + x_f * y_f )
					geom_c( k_ind, q_ind, p_ind )     = geo * ( z_f * ( one -  y_f * y_f ) )
					geom_h( k_ind, q_ind, p_ind )     = geo * ( z_f + x_f * y_f )

					triad_count                            = triad_count + 1
					! COUNTING THE TRIAD

				ELSE

					PRINT*,"TRIANGLE INCOMPATIBLE FOR ",wno(k_ind),wno(q_ind),wno(p_ind)

				END IF

			END DO

		END DO
		END DO

	END
! </f>

	SUBROUTINE triad_debug
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to check all the triads are initialized properly.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		er_V_self = zero
		er_B_self = zero
		er_VB     = zero
		!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   T  R  I  A  D      D  E  B  U  G
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k_ind = 1, N
		DO q_ind = 1, N

			IF (p_ind_max(k_ind,q_ind) .GT. N) THEN
				PRINT*,'ERROR IN UPPER BOUNDARY',p_ind_max(k_ind,q_ind)
				sys_status = 0
			END IF
			IF (p_ind_min(k_ind,q_ind) .LT. 1) THEN
				PRINT*,'ERROR IN LOWER BOUNDARY',p_ind_min(k_ind,q_ind),k_ind,q_ind
				sys_status = 0
			END IF

		DO p_ind = p_ind_min(k_ind,q_ind),p_ind_max(k_ind,q_ind)

			IF ( kqp_status( k_ind, q_ind, p_ind ) .NE. kqp_status( k_ind, p_ind, q_ind) ) THEN
				PRINT*,'ERROR IN TRIAD : kqp is ',kqp_status( k_ind, q_ind, p_ind ),' whereas kpq is ',kqp_status( k_ind, p_ind, q_ind )
				PRINT*,'at k=',wno( k_ind ),' q=',wno( q_ind ),' p=',wno( p_ind )
				sys_status = 0
			END IF
			IF ( kqp_status( k_ind, q_ind, p_ind ) .NE. kqp_status( q_ind, k_ind, p_ind) ) THEN
				PRINT*,'ERROR IN TRIAD : kqp is ',kqp_status( k_ind, q_ind, p_ind ),' whereas qkp is ',kqp_status( q_ind, k_ind, p_ind )
				PRINT*,'at k=',wno( k_ind ),' q=',wno( q_ind ),' p=',wno( p_ind )
				sys_status = 0
			END IF

			IF ( kqp_status( k_ind, q_ind, p_ind ) .EQ. 1 ) THEN
				er_V_self = er_V_self+ DABS( laplacian_k( k_ind ) * geom_b( k_ind, q_ind, p_ind ) &
				 													 - laplacian_k( p_ind ) * geom_b( p_ind, q_ind, k_ind ) )
				er_B_self = er_B_self+ DABS( laplacian_k( k_ind ) * geom_h( k_ind, q_ind, p_ind ) &
				 													 - laplacian_k( p_ind ) * geom_h( p_ind, q_ind, k_ind ) )
				er_VB     = er_VB + DABS( geom_b( k_ind, q_ind, p_ind ) + geom_b( k_ind, p_ind, q_ind ) &
				 										- geom_c( k_ind, q_ind, p_ind ) - geom_c( k_ind, p_ind, q_ind ) )
			END IF

		END DO
		END DO
		END DO

	END
! </f>

	SUBROUTINE compute_eddy_damping
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the damping factor for the third wnoent
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		DOUBLE PRECISION::eddy_0

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   E  D  D  Y            F  R  E  Q  U  E  N  C  Y
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		wno_diss_V = ( ds_rate_visc_V / visc ** thr ) **qtr
		IF ( wno_diss_V .GT. wno_forc ) THEN
			kD_V_ind   = FLOOR( DLOG( wno_diss_V / wno_base ) / wno_scale_log ) + 1
		END IF
		IF ( wno_diss_V .GT. wno( N ) ) THEN
			kD_V_ind   = N
		END IF

		!kD_V_ind = 37
		! FOR CUSTOM USE

		eddy_0 = DSQRT( SUM( ( en_spec_V( :kD_V_ind ) + en_spec_B( :kD_V_ind ) ) * laplacian_k( :kD_V_ind ) * wno_band( :kD_V_ind ) ) )
		! Eddy time-scale at the dissipative wavenumber

		DO k_ind = 1, N

			eddy_V( k_ind ) = DSQRT( SUM( ( en_spec_V( :k_ind ) + en_spec_B( :k_ind ) ) * laplacian_k( :k_ind ) * wno_band( :k_ind ) ) )
			! Eddy time-scale at the 'i'th wavenumber

			IF ( coupling_status .NE. 0 ) THEN

				! eddy_B( k_ind ) = eddy_const * ( eddy_0 ** eddy_exp_C ) * ( eddy_V( k_ind ) ** eddy_exp ) + diff * laplacian_k( k_ind )
				! When '\alpha' effect is included for 'B'.

				eddy_B( k_ind ) = eddy_const * eddy_V( k_ind ) + diff * laplacian_k( k_ind )
				! When '\alpha' effect is not included for 'B'.

				eddy_B( k_ind ) = eddy_B( k_ind ) + alfven_const * wno( k_ind ) * DSQRT( DABS(SUM( en_spec_B( :k_ind ) * wno_band( :k_ind ))))
				! Additional Alfven velocity time-scale

			END IF

			! eddy_V( k_ind ) = eddy_const * ( eddy_0 ** eddy_exp_C ) * ( eddy_V( k_ind ) ** eddy_exp ) + visc * laplacian_k( k_ind )
			! When '\alpha' effect is included for 'U'.

			eddy_V( k_ind ) = eddy_const * eddy_V( k_ind )  + visc * laplacian_k( k_ind )
			! When '\alpha' effect is not included for 'U'.

		END DO

	END
! </f>

	SUBROUTINE compute_kinetic_spectral_data
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO compute all spectral data and write them
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! S P E C T R A L    D A T A
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k_ind = 1, N
			fl_spec_V( k_ind ) 		=  SUM( tr_spec_V( k_ind : ) * wno_band( k_ind : ) )
			fl_spec_V_self( k_ind )	=  SUM( tr_spec_V_self( k_ind : ) * wno_band( k_ind : ) )
		END DO

		CALL write_kinetic_energy_spectrum()
		! REF-> <<< system_basicoutput >>>

	END
! </f>

	SUBROUTINE compute_magnetic_spectral_data
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO compute all spectral data and write them
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! S P E C T R A L    D A T A
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k_ind = 1, N
			fl_spec_B( k_ind ) 		= SUM( tr_spec_B( k_ind : ) * wno_band( k_ind : ) )
			fl_spec_B_intr( k_ind )	= SUM( tr_spec_B_intr( k_ind : ) * wno_band( k_ind : ) )
		END DO

		CALL write_magnetic_energy_spectrum()
		! REF-> <<< system_basicoutput >>>

	END
! </f>

	SUBROUTINE compute_kinetic_temporal_data
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO compute net energy, enstrophy, dissipation and skewness for every time step
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! T E M P O R A L    D A T A
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		energy_V     = SUM( en_spec_V * wno_band )
		enstrophy_V  = SUM( laplacian_k * en_spec_V * wno_band )

		IF ( visc_status .EQ. 1 ) THEN
			ds_rate_visc_V = two * visc * enstrophy_V
		END IF

		ds_rate_net_V  = ( energy_V - energy_V_prev ) / dt
		ds_rate_intr_V = SUM( tr_spec_V_intr * wno_band )
		ds_rate_self_V = SUM( tr_spec_V_self * wno_band )

		IF ( forc_status .EQ. 1 ) THEN
			ds_rate_forc   = SUM( fr_spec * wno_band )
		END IF

		IF ( forc_status .EQ. 1 ) THEN
			skewness = SUM( ( tr_spec_V + fr_spec ) * laplacian_k * wno_band )
		ELSE
			skewness = SUM( tr_spec_V * laplacian_k * wno_band )
		END IF

		skewness   = skewness * ( enstrophy_V ** ( -1.5D0 )) * skewness_const

		CALL write_kinetic_temporal_data
		! REF-> <<< system_basicoutput >>>

		energy_V_prev = energy_V

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE compute_magnetic_temporal_data
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO compute net energy, enstrophy, dissipation and skewness for every time step
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! T E M P O R A L    D A T A
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		energy_B     = SUM( en_spec_B * wno_band )
		enstrophy_B  = SUM( laplacian_k * en_spec_B * wno_band )
		potential_B  = SUM( ( en_spec_B / laplacian_k ) * wno_band )

		IF ( diff_status .EQ. 1 ) THEN
			ds_rate_diff_B  = two * diff * enstrophy_B
		END IF

		ds_rate_net_B   = ( energy_B - energy_B_prev ) / dt
		ds_rate_intr_B  = SUM( tr_spec_B_intr * wno_band )
		ds_rate_self_B  = SUM( tr_spec_B_self * wno_band )
		dynamo_exp      = DLOG( energy_B / energy_B_prev ) / dt
		! dynamo_exp_calc = SUM( dyn_rate_spec * en_spec_B * wno_band ) / energy_B

		CALL write_magnetic_temporal_data
		! REF-> <<< system_basicoutput >>>

		energy_B_prev = energy_B

		energy_tot = energy_V + energy_B

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE compute_forcing_spectrum
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO compute how to initiate forcing spectrum
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  F  O  R  C  I  N  G
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		fr_spec = zero

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! TYPE OF FORCING
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! 1. FORCED AT LARGE SCALE WITH A FIXED FORM, AND COEFFICIENT 
		!  	ADJUSTED TO MATCH THE NET DISSIPATION
		! --------------------------------------------------------------------------
		forcing_factor = ( ds_rate_visc_V + ds_rate_diff_B )
		DO k_ind = 1, N
			fr_spec( k_ind ) = forcing_factor * spec0( k_ind )
		END DO

		! 2. INJECTED LINEAR TO ENERGY SPECTRUM USING INTEGRATING FACTOR
		! TURN ON THE FBACK_COEF TO KEEP CONSTANT ENERGY
		! --------------------------------------------------------------------------
		! forcing_factor = ( ds_rate_visc_V + ds_rate_diff_B )
		! forcing_factor = ds_rate_ref_V / SUM( en_spec_V(:kF_ind) * wno_band(:kF_ind) )
		! DO k_ind = 1, kF_ind
		! 	fr_spec( k_ind ) = forcing_factor * en_spec_V( k_ind )
		! END DO
		! integ_factor_V = DEXP( (- two * visc * laplacian_k + fr_spec ) * dt )

		! 3. CONSTANT ENERGY FOR CERTAIN MODES
		! --------------------------------------------------------------------------
		! DO k_ind = 1, kF_ind
		! 	en_spec_V( k_ind ) = spec0( k_ind )
		! END DO

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


	END
! </f>

	SUBROUTINE prepare_initial_condition
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO initiate kinetic energy spectrum
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE

		! CALL IC_V_large_eddies
		CALL IC_V_power_law
		! REF-> <<< system_initialcondition >>>
		
		CALL compute_eddy_damping
		! REF-> <<< system_basicfunctions >>>

		CALL prepare_output
		! REF-> <<< system_basicoutput >>>

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	END
! </f>

	SUBROUTINE prepare_perturbation_dynamo
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! TO compute how to initiate perturbation for the dynamo testing
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  Y  N  A  M  O
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		energy_B_prev   = energy_B_0

		! CALL IC_B_small_scale_dynamo 
		CALL IC_B_large_eddies
		! This is for perturbing the field at small scales
		! REF-> <<< system_initialcondition >>>

		CALL compute_eddy_damping
		! REF-> <<< system_basicfunctions >>>

		CALL prepare_output_dynamo
		! REF-> <<< system_basicoutput >>>
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	END
! </f>

	SUBROUTINE perform_debug
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to check Nan in data
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	IMPLICIT NONE

	DO k_ind = 1, N
		IF ( en_spec_V ( k_ind) .NE. en_spec_V ( k_ind) ) THEN
			nan_count =  nan_count + 1
		END IF
		IF ( en_spec_B ( k_ind) .NE. en_spec_B ( k_ind) ) THEN
			nan_count =  nan_count + 1
		END IF
	END DO

	IF ( nan_count .GT. 0 ) THEN
		PRINT*,nan_count," NaN encountered around t=",time_now
		! IF any NaN is encountered, the loop is exited without any further continuation.
		nan_status = 1
	END IF

  END
! </f>

	SUBROUTINE deallocate_edqnm_arrays
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	!       deallocate all edqnm arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  R  R  A  Y     D  E   A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE( p_ind_min, p_ind_max )
		DEALLOCATE( geom_b, geom_c, geom_h )
		DEALLOCATE( kqp_status )
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	END
! </f>

	DOUBLE PRECISION FUNCTION cosine( i1, i2, i3 )
! <f
	! ------------
	! FUNCTION TO: Calculate the cosine of angle opposite to 'i3'
	! -------------
		IMPLICIT NONE
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		INTEGER(KIND=4),INTENT(IN)::i1, i2, i3

		cosine  = ( wno( i1 ) ** two ) + ( wno( i2 ) ** two  ) - ( wno( i3 ) ** two )
		cosine  = cosine / ( two * wno( i1 ) * wno( i2 ) )
		RETURN

	END
! </f>

  INTEGER FUNCTION lowerindex( wno0 )
! <f
	! ------------
	! FUNCTION TO: Calculate the largest index whose momentum is smaller than the given momentum 'wno0'
	! -------------
		DOUBLE PRECISION::wno0
		IF ( wno0 .GT. wno_max - tol_float ) THEN
			lowerindex = N
		! print*,'lowerindex',wno0, 'tolerance',lowerindex
		ELSE
			lowerindex = FLOOR ( ( DLOG( wno0 / wno_base) / wno_scale_log ) + tol_float ) + 1
		! print*,'lowerindex',wno0, 1+ DLOG( wno0/wno_base) / wno_scale_log,lowerindex
		END IF
		RETURN

  END
! </f>

  INTEGER FUNCTION upperindex( wno0 )
! <f
	! ------------
	! FUNCTION TO: Calculate the smallest index whose momentum is larger than the given momentum 'wno0'
	! -------------
		DOUBLE PRECISION::wno0
		IF ( wno0 .LT. wno_min + tol_float ) THEN
			upperindex = 1
		! print*,'upperindex',wno0, 'tolerance',upperindex
		ELSE
			upperindex  =  CEILING ( ( DLOG( wno0 / wno_base) / wno_scale_log )  ) + 1
		! print*,'upperindex',wno0, 1+ DLOG( wno0/wno_base) / wno_scale_log,upperindex
		END IF
		RETURN

  END
! </f>

	INTEGER FUNCTION triangle_compatibility( i1, i2, i3 )
! <f
	! ------------
	! FUNCTION TO: Calculate the compatibility that three momentum can form a triangle, 1 is yes, 0 means no.
	! -------------
		IMPLICIT NONE
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		INTEGER(KIND=4),INTENT(IN)::i1, i2, i3

		triangle_compatibility =  0

		IF ( wno(i1) + wno(i2) .GE. wno(i3) ) THEN
		IF ( wno(i3) + wno(i2) .GE. wno(i1) ) THEN
		IF ( wno(i1) + wno(i3) .GE. wno(i2) ) THEN

			triangle_compatibility =  1

		END IF
		END IF
		END IF

		RETURN

		END
! </f>

END MODULE system_basicfunctions
