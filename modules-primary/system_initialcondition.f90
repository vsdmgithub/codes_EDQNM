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
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_V_power_law
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with large eddies. Same as one of the forcing template
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		DOUBLE PRECISION::spec_exp
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! spec_exp = -5.0D0/3.0D0
		spec_exp = - ( 16.0D0 - 9.0D0 * eddy_damping_exp ) / ( 8.0D0 - 3.0D0 * eddy_damping_exp )
		en_spec_V = wno**spec_exp
		en_spec_V = energy_V_0 * en_spec_V / SUM( en_spec_V * wno_band )
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
		DOUBLE PRECISION::en_noise
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		en_noise = zero
		! Set it to zero, if no noise is required to avoid NaN in the magnetic spectrum

		en_spec_B = ( energy_B_0 - en_noise ) * spec0 +  en_noise / ( N * wno_band )
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE IC_B_small_scale_dynamo
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Initialize initial condition with a single mode
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION,DIMENSION(N)::spec_B
		DOUBLE PRECISION::en_noise
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		en_noise = 1E-8
		! Set it to zero, if no noise is required to avoid NaN in the magnetic spectrum

		en_spec_B    = DEXP( - ( ( wno-wno_diss_B ) / wno_base ) ** two ) 
		en_spec_B    = en_spec_B / SUM( en_spec_B * wno_band )
		en_spec_B    = ( energy_B_0 - en_noise ) * en_spec_B +  en_noise / ( N * wno_band )
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
		CHARACTER(LEN=150)::spec0_address
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		spec0_address='IC_V_N41'

		! where initial condition is stored

		OPEN( UNIT = 2001 ,FILE = TRIM(ADJUSTL(spec0_address))//'.dat' )

		DO k_ind = 1, N

			READ( 2001, '(F12.6)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='NO') en_spec_V( k_ind )
			READ( 2001, '(F32.17)',ADVANCE='NO')   dum
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

  SUBROUTINE IC_B_read_from_file
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
		DOUBLE PRECISION::dum
		CHARACTER(LEN=100)::spec0_address
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		spec0_address='IC_B_N33'
		! where initial condition is stored

		OPEN( UNIT = 2001 ,FILE = TRIM(ADJUSTL(spec0_address))//'.dat' )

		DO k_ind = 1, N

			READ( 2001, '(F12.6)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='NO') en_spec_B( k_ind )
			READ( 2001, '(F32.17)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='NO')   dum
			READ( 2001, '(F32.17)',ADVANCE='YES')  dum

		END DO

		CLOSE( 2001 )

		en_spec_B = energy_B_0 * en_spec_B / SUM( en_spec_B * wno_band )
		! Adjustng the normalization constant.

		PRINT*,"IC READ SUCCESFULLY"
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

END MODULE system_initialcondition
