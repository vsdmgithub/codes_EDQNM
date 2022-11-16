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
		ALLOCATE( kqp_status( N, N, N ), geom_B( N, N, N ), weightage( N, N, N ) )
		ALLOCATE( eddy_array( N ) )
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
		DOUBLE PRECISION:: wei,geo

		kqp_status = 0
		geom_B     = zero
		triad_count= 0
		! Reseting arrays to zero for safety

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

			! print*,wno(k_ind),wno(q_ind),wno(p_ind_min(k_ind,q_ind)),wno(p_ind_max(k_ind,q_ind))

			DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

				IF ( triangle_compatibility( k_ind, q_ind, p_ind ) .EQ. 1 ) THEN

					kqp_status( k_ind, q_ind, p_ind ) = 1
					! MEANING, THE GIVEN THREE MOMENTUM 'p,k,q' CAN FORM A TRIANGLE.

					x_f                               = - cosine( q_ind, p_ind, k_ind )
					y_f                               = + cosine( k_ind, q_ind, p_ind )
					z_f                               = + cosine( k_ind, p_ind, q_ind )
					! FINDING COSINES FOR ALL THREE SIDES ONCE IT IS APPROVED TO PARTICIPATE IN THE TRIAD INTERACTION

					geo                               = ( ( z_f ** thr ) + x_f * y_f + dim_min_3 * hf * ( z_f + x_f * y_f ) )
					geom_B( k_ind, q_ind, p_ind )     = geo * wno( p_ind ) / wno( k_ind )
					! GEOMETRIC FACTOR 'B' IN THE E.D.Q.N.M

					wei                               = wno_band( q_ind ) * wno_band( p_ind )
					! FACTOR OF AREA OF THE STRIP IN 'QP' PLANE
					wei                               = wei * sym_const
					! FACTOR OF 8*S_(D-2)/S_(D-1)
					wei                               = wei * wno( k_ind ) / ( wno( q_ind ) * wno( p_ind ) )
					! FACTOR OF K/PQ
					weightage( k_ind, q_ind, p_ind )  = wei * ( DSQRT( ( one - ( x_f ** two ) ) / ( wno( k_ind ) ** two ) ) ** dim_min_3 )
					! FACTOR of ( (1-x^2)/k^2 ) ** (d-3)/2

					triad_count                       = triad_count + 1
					! COUNTING THE TRIAD

				END IF

			END DO

		END DO
		END DO
! print*,N*N*N,triad_count
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

		DO p_ind = 1, N

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

		END DO
		END DO
		END DO

	END
! </f>

	SUBROUTINE compute_spectral_data
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
			fl_spec( k_ind ) = - SUM( tr_spec( : k_ind) * wno_band( : k_ind) )
	  END DO

		! CALL write_spectrum('energy',en_spec) !  ENERGY FILE
		! CALL write_spectrum('transfer',tr_spec) !  ENERGY TRANSFER FILE
		! CALL write_spectrum('flux',fl_spec)     !  FLUX  FILE
		! ! REF-> <<< system_basicoutput >>>
		! IF ( visc_status .EQ. 1 ) THEN
		! 	CALL write_spectrum('dissipation',two * visc * laplacian_k * en_spec) !  DISSIPATION FILE
		! END IF

		CALL write_energy_spectrum()

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! T E M P O R A L    D A T A
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		energy     = SUM( en_spec * wno_band )
		enstrophy  = SUM( laplacian_k * en_spec * wno_band )

		IF ( visc_status .EQ. 1 ) THEN
			ds_rate  = two * visc * enstrophy
		END IF

		IF ( forc_status .EQ. 1 ) THEN
			skewness = SUM( ( tr_spec + fr_spec ) * laplacian_k * wno_band )
		ELSE
			skewness = SUM( tr_spec * laplacian_k * wno_band )
		END IF

		skewness   = skewness * ( enstrophy ** ( -1.5D0 )) * skewness_const

		CALL write_temporal_data
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
		IF ( en_spec ( k_ind) .NE. en_spec ( k_ind) ) THEN
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
		DEALLOCATE( kqp_status, geom_B, weightage )
		DEALLOCATE( eddy_array )
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

		IF ( wno(i1) + wno(i2) .GT. wno(i3) ) THEN
		IF ( wno(i3) + wno(i2) .GT. wno(i1) ) THEN
		IF ( wno(i1) + wno(i3) .GT. wno(i2) ) THEN

			triangle_compatibility =  1

		END IF
		END IF
		END IF

		RETURN

		END
! </f>

END MODULE system_basicfunctions
