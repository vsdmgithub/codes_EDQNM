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
! LAST MODIFIED: 15 NOV 2022
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
	INTEGER(KIND=4),DIMENSION(:,:,:),ALLOCATABLE  ::kqp_status

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
		ALLOCATE( triad_weightage( N, N, N ) )
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
		DOUBLE PRECISION:: wp,wq,wk
		DOUBLE PRECISION:: geom_factor
		DOUBLE PRECISION:: geo,triad_factor

		ALLOCATE( kqp_status( N, N, N ) )

		triad_count   = 0
		triad_deleted = 0
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
					geo                               = ( ( z_f ** thr ) + x_f * y_f + dim_min_3 * hf * ( z_f + x_f * y_f ) )
					geom_b( k_ind, q_ind, p_ind )     = geo * wno( p_ind ) / wno( k_ind )
					geo                               = ( z_f * ( one -  y_f * y_f ) + dim_min_3 * hf * ( z_f + x_f * y_f ) )
					geom_c( k_ind, q_ind, p_ind )     = geo * wno( p_ind ) / wno( k_ind )
					geo                               = ( dim - one ) * hf * ( z_f + x_f * y_f )
					geom_h( k_ind, q_ind, p_ind )     = geo * wno( p_ind ) / wno( k_ind )

					! GEOMETRIC FACTORS DEFINED THROUGH K,P,Q
					! wk                            = wno( k_ind )
					! wp                            = wno( p_ind )
					! wq                            = wno( q_ind )
					! geom_factor                   = (wk+wp+wq)*(wk+wp-wq)*(wk-wp+wq)*(-wk+wp+wq)/(8.0D0*wk*wk*wq*wq)
					! geo                           = geom_factor * ( (dim-one)*wk*wk*wp*wp - ( wk*wk + wp*wp )*wq*wq + wq**4.0d0 )
					! geom_b( k_ind, q_ind, p_ind ) = geo/(wk*wk*wp*wp)
					! geo                           = geom_factor * ( (dim-two)*wk*wk + wp*wp - wq*wq )
					! geom_c( k_ind, q_ind, p_ind ) = geo/(wk*wk)
					! geo                           = geom_factor * (dim - one)
					! geom_h( k_ind, q_ind, p_ind ) = geo

					IF ( dim_min_3 .LT. -0.01D0 ) THEN
					  IF ( DABS( x_f ** two - one ) .LT. tol_float ) THEN
							triad_factor  = zero
							triad_deleted = triad_deleted + 1
							! TO AVOID NAN
						ELSE
							triad_factor  = DSQRT( DABS( one - ( x_f ** two ) ) / ( wno( k_ind ) ** two ) ) ** dim_min_3
						END IF
					ELSE
						triad_factor    = DSQRT( DABS( one - ( x_f ** two ) ) / ( wno( k_ind ) ** two ) ) ** dim_min_3
					END IF
					! FACTOR of ( (1-x^2)/k^2 ) ** (d-3)/2

					triad_weightage( k_ind, q_ind, p_ind ) = dim_const * triad_factor
					triad_count                            = triad_count + 1
					! COUNTING THE TRIAD

					IF (triad_weightage( k_ind, q_ind, p_ind ) .NE. triad_weightage( k_ind, q_ind, p_ind ) ) THEN
						PRINT*,"NaN HERE AT TRIAD WEIGHTAGE",k_ind,q_ind,p_ind,triad_weightage(k_ind,q_ind,p_ind)
						nan_status=1
					END IF

				ELSE

					PRINT*,"TRIANGLE INCOMPATIBLE FOR ",wno(k_ind),wno(q_ind),wno(p_ind)

				END IF

			END DO

		END DO
		END DO

		triad_count = triad_count - triad_deleted
		! Subracting the triads that are deleted

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
		! DOUBLE PRECISION:: dum,wp,wq,wk
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

		DO p_ind = 1, 5

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

		DEALLOCATE( kqp_status )
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
			fl_spec_V( k_ind ) = - SUM( tr_spec_V( : k_ind) * wno_band( : k_ind) )
	  END DO

		! UNCOMMENT TO WRITE IN SEPERATE FILES
		! CALL write_spectrum('energy_V',en_spec_V)
		! CALL write_spectrum('transfer_V',tr_spec_V)
		! CALL write_spectrum('flux_V',fl_spec_V)
		! ! REF-> <<< system_basicoutput >>>

		! IF ( visc_status .EQ. 1 ) THEN
		! 	CALL write_spectrum('dissipation',two * visc * laplacian_k * en_spec_V) !  DISSIPATION FILE
		! END IF

		! UNCOMMENT TO WRITE IN SINGLE FILE
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
			fl_spec_B( k_ind ) = - SUM( tr_spec_B( : k_ind) * wno_band( : k_ind) )
	  END DO

		! UNCOMMENT TO WRITE IN SEPERATE FILES
		! CALL write_spectrum('energy_B',en_spec_B)
		! CALL write_spectrum('transfer_B',tr_spec_B)
		! CALL write_spectrum('flux_B',fl_spec_B)
		! ! REF-> <<< system_basicoutput >>>

		! IF ( diff_status .EQ. 1 ) THEN
		! 	CALL write_spectrum('diffusion',  two * diff * laplacian_k * en_spec_B) !  DISSIPATION FILE
		! END IF

		! UNCOMMENT TO WRITE IN SINGLE FILE
		CALL write_magnetic_energy_spectrum()
		! REF-> <<< system_basicoutput >>>

	END
! </f>

	SUBROUTINE compute_temporal_data
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

		IF ( forc_status .EQ. 1 ) THEN
			ds_rate_forc   = SUM( fr_spec * wno_band )
		END IF

		IF ( forc_status .EQ. 1 ) THEN
			skewness = SUM( ( tr_spec_V + fr_spec ) * laplacian_k * wno_band )
		ELSE
			skewness = SUM( tr_spec_V * laplacian_k * wno_band )
		END IF

		skewness   = skewness * ( enstrophy_V ** ( -1.5D0 )) * skewness_const

		IF ( coupling_status .NE. 0 ) THEN

			energy_B     = SUM( en_spec_B * wno_band )
			enstrophy_B  = SUM( laplacian_k * en_spec_B * wno_band )

			IF ( diff_status .EQ. 1 ) THEN
				ds_rate_diff_B  = two * diff * enstrophy_B
			END IF

			ds_rate_net_B   = ( energy_B - energy_B_prev ) / dt
			ds_rate_intr_B  = SUM( tr_spec_B_intr * wno_band )
			dynamo_exp      = ds_rate_net_B / energy_B
			dynamo_exp_calc = SUM( dyn_rate_spec * en_spec_B * wno_band ) / energy_B

			CALL write_magnetic_temporal_data
			! REF-> <<< system_basicoutput >>>

		END IF

		energy_tot   =  energy_V + energy_B
		CALL write_kinetic_temporal_data
		! REF-> <<< system_basicoutput >>>

		energy_V_prev = energy_V
		energy_B_prev = energy_B

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

		! IF ( MOD(t_step,t_step_forcing) .EQ. 0 ) THEN
		! END IF
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! TYPE OF FORCING
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! 1. Injected as constant
		! F(k)= f * E(k) for 0<k<kF_ind; f = ds_rate / 0_to_kF ind E(k)dk
		! --------------------------------------------------------------------------
		! forcing_factor = forcing_factor * ( ( 1 - ( energy_tot - energy_0 ) * fback_coef ) **0.25D0 )
		! fact = 0.2, fdback=15 working for d=2

		forcing_factor = ( ds_rate_visc_V + ds_rate_diff_B ) * ( ( 1 - ( energy_tot - energy_0 ) ) ** two )
		DO k_ind = 1, N
			fr_spec( k_ind ) = forcing_factor * spec0( k_ind )
		END DO

		! 2. Injected linear to energy spectrum using integrating factor
		! turn on the fback_coef to keep constant energy
		! --------------------------------------------------------------------------
		! ds_rate_ref_V = ds_rate_V - fback_coef * ( energy_V - energy_V_0 )
		! ds_rate_ref_V  = - ( energy_V - energy_V_0 ) / dt
		! ds_rate_ref_V  = ds_rate_visc_V - ds_rate_intr_V
		! forcing_factor = ds_rate_ref_V / SUM( en_spec_V(:kF_ind) * wno_band(:kF_ind) )
		! fr_spec        = zero
		! DO k_ind = 1, kF_ind
		! 	fr_spec( k_ind ) = forcing_factor * en_spec_V( k_ind )
		! END DO
		! integ_factor_V = DEXP( (- two * visc * laplacian_k + forcing_factor ) * dt )
		! energy_V_0     = energy_V

		! 3. Constant energy for certain modes
		! --------------------------------------------------------------------------
		! DO k_ind = 1, kF_ind
		! 	en_spec_V( k_ind ) = spec0( k_ind )
		! END DO

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
		coupling_status = 1
		energy_B_prev   = energy_B_0
		en_spec_V       = en_spec_V * ( energy_V_0 / SUM( en_spec_V * wno_band ) )

		! CALL IC_B_large_eddies_single_mode
		CALL IC_B_large_eddies
		! CALL IC_B_read_from_file
		! CALL IC_B_equipartition
		! REF-> <<< system_initialcondition >>>

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
		DEALLOCATE( triad_weightage )
		DEALLOCATE( geom_b, geom_c, geom_h )
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

! <f
	INTEGER FUNCTION triangle_compatibility( i1, i2, i3 )
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
