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
! MODULE NAME: system_advfunctions
! LAST MODIFIED: 18 OCT 2023
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER FOR  EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_advfunctions
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * Using transfer spectrum - splits into local and nonlocal contributions
! * Using transfer spectrum - splits into positive and negative contributions
! * Mix localness and positiveness
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_solver

	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT NONE
	! _________________________
	! ADVANCED FUNCTION VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(3)::triad_sides
	! _________________________
	! ARRAYS
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::fl_pos,fl_neg
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::fl_loc,fl_non
	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	CONTAINS
! </f>

	SUBROUTINE allocate_flux_decomposition_arrays
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	!       allocate all arrays used in this module
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  R  R  A  Y     A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE( fl_pos( N ), fl_neg( N ) )
		ALLOCATE( fl_loc( N ), fl_non( N ) )

	END
! </f>

	SUBROUTINE compute_flux_decomposition_V
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to decompose the net flux at 'k' into positive and negative flux.
	! CALL this to decompose the net flux at 'k' into local and non-local flux.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		DOUBLE PRECISION::min_max

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   F L U X     D E C O M P O S I T I O N     
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		fl_pos          = zero
		fl_neg          = zero
		fl_loc          = zero
		fl_non          = zero

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   P  O  S  I  T  I  V  E     F  L  U  X
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k2_ind = 1 , N
		DO k_ind  = k2_ind , N
		DO q_ind  = 1 , k2_ind
		DO p_ind  = 1 , k2_ind

			IF ( kqp_status( k_ind, q_ind, p_ind )  .EQ. 1 ) THEN

				CALL transfer_term_integrand_V
				fl_pos( k2_ind ) = fl_pos( k2_ind ) + integrand_V_self * wno_band( k_ind )

			END IF

		END DO
		END DO
		END DO
		END DO

		fl_neg    =   fl_pos    -   fl_spec_V

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  L O C A L     F L U X
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k2_ind = 1, N
		DO k_ind = k2_ind, N
		DO q_ind = 1, N
		DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

			triad_sides =   (/ wno( k_ind ), wno( q_ind ), wno( p_ind ) /)
			! Array consisting the 3 sides of the triangle

			min_max   =   MINVAL( triad_sides ) / MAXVAL( triad_sides )
			! Finding the ratio of minimum to maximum sides of triangle

			CALL transfer_term_integrand_V

			IF ( min_max .GT. loc_par ) THEN

				fl_loc( k2_ind ) = fl_loc( k2_ind ) + integrand_V_self * wno_band( k_ind )
				! Summation terms over all possible k,q,p for a given k' in the region of local interactions.

			ELSE 

				fl_non( k2_ind ) = fl_non( k2_ind ) + integrand_V_self * wno_band( k_ind )

			END IF

		END DO
		END DO
		END DO
		END DO

		CALL write_flux_decomposition_V

	END
! </f>

	SUBROUTINE write_flux_decomposition_V
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Write all the flux decomposition data
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  P  R  I  N   T          O  U  T  P  U  T
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_sp )) &
		            // 'flux_V_t_'// TRIM( ADJUSTL( file_time ) ) // '.dat'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OPEN(UNIT = 282, FILE = file_name)

		DO k_ind = 1, N

			WRITE(282,f_d12p6,ADVANCE  = 'no')  wno( k_ind )
			WRITE(282,f_d32p17,ADVANCE = 'no')fl_spec_V_self( k_ind )
			WRITE(282,f_d32p17,ADVANCE = 'no')fl_pos( k_ind )
			WRITE(282,f_d32p17,ADVANCE = 'no')fl_neg( k_ind )
			WRITE(282,f_d32p17,ADVANCE = 'no')fl_loc( k_ind )
			WRITE(282,f_d32p17,ADVANCE = 'yes')fl_non( k_ind )

		END DO

		CLOSE(282)

  END
! </f>

	SUBROUTINE compute_flux_decomposition_B_self
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to decompose the net flux at 'k' into positive and negative flux for the self transfer term.
	! CALL this to decompose the net flux at 'k' into local and nonlocal flux for the self transfer term.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		DOUBLE PRECISION::min_max

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   F L U X     D E C O M P O S I T I O N     
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		fl_pos          = zero
		fl_neg          = zero
		fl_loc          = zero
		fl_non          = zero

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   P  O  S  I  T  I  V  E     F  L  U  X
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k2_ind = 1 , N
		DO k_ind  = k2_ind , N
		DO q_ind  = 1 , k2_ind
		DO p_ind  = 1 , k2_ind

			IF ( kqp_status( k_ind, q_ind, p_ind )  .EQ. 1 ) THEN

				CALL transfer_term_integrand_B
				fl_pos( k2_ind ) = fl_pos( k2_ind ) + integrand_B_self * wno_band( k_ind )

			END IF

		END DO
		END DO
		END DO
		END DO

		fl_neg    =   fl_pos    -   (fl_spec_B - fl_spec_B_intr)

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  L O C A L     F L U X
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k2_ind = 1, N
		DO k_ind = k2_ind, N
		DO q_ind = 1, N
		DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

			triad_sides =   (/ wno( k_ind ), wno( q_ind ), wno( p_ind ) /)
			! Array consisting the 3 sides of the triangle

			min_max   =   MINVAL( triad_sides ) / MAXVAL( triad_sides )
			! Finding the ratio of minimum to maximum sides of triangle

			CALL transfer_term_integrand_B

			IF ( min_max .GT. loc_par ) THEN

				fl_loc( k2_ind ) = fl_loc( k2_ind ) + integrand_B_self * wno_band( k_ind )
				! Summation terms over all possible k,q,p for a given k' in the region of local interactions.

			ELSE 

				fl_non( k2_ind ) = fl_non( k2_ind ) + integrand_B_self * wno_band( k_ind )

			END IF

		END DO
		END DO
		END DO
		END DO

		CALL write_flux_decomposition_B_self

	END
! </f>

	SUBROUTINE write_flux_decomposition_B_self
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Write all the flux decomposition data
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  P  R  I  N   T          O  U  T  P  U  T
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_sp )) &
		            // 'flux_B_self_t_'// TRIM( ADJUSTL( file_time ) ) // '.dat'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OPEN(UNIT = 382, FILE = file_name)

		DO k_ind = 1, N

			WRITE(382,f_d12p6,ADVANCE  = 'no')  wno( k_ind )
			WRITE(382,f_d32p17,ADVANCE = 'no')(fl_spec_B( k_ind ) - fl_spec_B_intr( k_ind ))
			WRITE(382,f_d32p17,ADVANCE = 'no')fl_pos( k_ind )
			WRITE(382,f_d32p17,ADVANCE = 'no')fl_neg( k_ind )
			WRITE(382,f_d32p17,ADVANCE = 'no')fl_loc( k_ind )
			WRITE(382,f_d32p17,ADVANCE = 'yes')fl_non( k_ind )

		END DO

		CLOSE(382)

  END
! </f>

	SUBROUTINE compute_flux_decomposition_B_intr
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to decompose the net flux at 'k' into local and nonlocal flux for the coupling transfer term.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		DOUBLE PRECISION::min_max

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   F L U X     D E C O M P O S I T I O N     
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		fl_loc          = zero
		fl_non          = zero

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  L O C A L     F L U X
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k2_ind = 1, N
		DO k_ind = k2_ind, N
		DO q_ind = 1, N
		DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

			triad_sides =   (/ wno( k_ind ), wno( q_ind ), wno( p_ind ) /)
			! Array consisting the 3 sides of the triangle

			min_max   =   MINVAL( triad_sides ) / MAXVAL( triad_sides )
			! Finding the ratio of minimum to maximum sides of triangle

			CALL transfer_term_integrand_B

			IF ( min_max .GT. loc_par ) THEN

				fl_loc( k2_ind ) = fl_loc( k2_ind ) + integrand_B_intr * wno_band( k_ind )
				! Summation terms over all possible k,q,p for a given k' in the region of local interactions.

			ELSE 

				fl_non( k2_ind ) = fl_non( k2_ind ) + integrand_B_intr * wno_band( k_ind )

			END IF

		END DO
		END DO
		END DO
		END DO

		CALL write_flux_decomposition_B_intr

	END
! </f>

	SUBROUTINE write_flux_decomposition_B_intr
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Write all the flux decomposition data
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  P  R  I  N   T          O  U  T  P  U  T
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_sp )) &
		            // 'flux_B_intr_t_'// TRIM( ADJUSTL( file_time ) ) // '.dat'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OPEN(UNIT = 482, FILE = file_name)

		DO k_ind = 1, N

			WRITE(482,f_d12p6,ADVANCE  = 'no')  wno( k_ind )
			WRITE(482,f_d32p17,ADVANCE = 'no')fl_spec_B_intr( k_ind )
			WRITE(482,f_d32p17,ADVANCE = 'no')fl_loc( k_ind )
			WRITE(482,f_d32p17,ADVANCE = 'yes')fl_non( k_ind )

		END DO

		CLOSE(482)

  END
! </f>

	SUBROUTINE deallocate_flux_decomposition_arrays
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! To deallocate all the flux decomposition arrays.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  R  R  A  Y        D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(fl_pos, fl_neg)
		DEALLOCATE(fl_loc, fl_non)

	END
! </f>

END MODULE system_advfunctions
