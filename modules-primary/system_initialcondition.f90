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
! MODULE NAME  : system_initialcondition
! LAST MODIFIED: 15 NOV 2022
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MODULE TO GET INITIAL CONDITION FOR THE SPECTRUM
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>
MODULE system_initialcondition
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * Either read from a file .
! * Or choose a type of IC.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT  NONE

  CONTAINS
! </f>

	SUBROUTINE IC_V_large_eddies
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with large eddies. Same as one of the forcing template
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		en_spec_V = energy_V_0 * spec0
		! en_spec_V =  spec0
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_V_kolmo
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with large eddies. Same as one of the forcing template
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		en_spec_V = energy_V_0 * specK
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_V_equipartition
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with large eddies. Same as one of the forcing template
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k_ind = 1, N
			en_spec_V( k_ind ) = wno( k_ind ) ** ( dim - one )
		END DO

		en_spec_V  = en_spec_V * ( energy_V_0 / SUM( en_spec_V * wno_band ) )
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_B_large_eddies
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with large eddies. Same as one of the forcing template
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		en_spec_B = energy_B_0 * spec0
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_B_small_eddies
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with small eddies.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION::s_exp,dum
		DOUBLE PRECISION,DIMENSION(N)::spec_B
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		s_exp     = 2.0D0 ! Integral scale spectrum exponent
		dum       = hf / ( wno_diss_B** two )
		spec_B    = ( wno ** s_exp ) * DEXP( - dum * laplacian_k )
		spec_B    = spec0 / SUM( spec_B * wno_band )
		en_spec_B = energy_B_0 * spec_B
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_B_large_eddies_single_mode
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with a single mode
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! en_spec_B           = tol_float
		! en_spec_B( N-1 ) = energy_B_0 / wno_band( N-1 )
		en_spec_B( kI_ind ) = energy_B_0 / wno_band( kI_ind )
		! en_spec_B( kD_ind ) = energy_B_0 / wno_band( kD_ind )
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_B_equipartition
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with large eddies. Same as one of the forcing template
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		DO k_ind = 1, N
			en_spec_B( k_ind ) = wno( k_ind ) ** ( dim - one )
		END DO

		en_spec_B  = en_spec_B * ( energy_B_0 / SUM( en_spec_B * wno_band ) )
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

  SUBROUTINE IC_V_read_from_file
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition from a given data from a file
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION::A0,dum
		CHARACTER(LEN=100)::spec0_address
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		spec0_address='IC_V_N37'
		! where initial condition is stored

		OPEN( UNIT = 2001 ,FILE = TRIM(ADJUSTL(spec0_address))//'.dat' )

		DO k_ind = 1, N

			READ( 2001, '(F12.6)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='NO') en_spec_V( k_ind )
			READ( 2001, '(F32.17)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='YES')  dum

		END DO

		CLOSE( 2001 )

		A0        = energy_V_0 / SUM( en_spec_V * wno_band )
		! Adjustng the normalization constant.

		en_spec_V = A0 * en_spec_V
		! UNCOMMENT TO NORMALIZE ENERGY

		PRINT*,"IC READ SUCCESFULLY"
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

END MODULE system_initialcondition
