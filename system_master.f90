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
! PROGRAM NAME : EDQNM d-dim
! LAST MODIFIED: 15 NOV 2022
! _____________________________________
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! PROGRAM FOR SOLVING 'D' DIMENSIONAL EDQNM EQUATION
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! </f>
PROGRAM EDQNM
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * This program solves the 'd' dimensional EDQNM equation, refer Orsag 1970.
! * It is equipped to solve in 'd' dimension, where 'd' is a parameter that can be specified.
! * Other inputs involve, total energy, initial condition and viscosity.
! * The system is solved in 'N' wavenumbers, spaced logarithmically. Not recommended to change the spacing.
! * Subroutines and functions are referenced to their corresponding modules below the calling statement, unless they are in the same module or program.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
!  MODULES
!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
USE system_main
USE system_timer

IMPLICIT NONE

CHARACTER(LEN=5)::run
! </f>

	CALL start_run_timer
	! REF-> <<< system_timer >>>

	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!  I  N  I  T  I  A  L  I  Z  A  T  I  O  N
	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	CALL read_input
	! REF-> <<< system_basicvarables >>>

	U_GRID=2
	W_GRID=2

	CALL init_global_variables
	! REF-> <<< system_basicvariables >>>

	CALL init_global_arrays
	! REF-> <<< system_basicvariables >>>

	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!  E  V  O  L  U  T  I  O  N
	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	CALL pre_analysis
	! REF-> <<< system_main >>>

	run = 'y' ! Easy way to stop the evolution with only initiation

	IF ( (run .EQ. 'y') .AND. ( sys_status .EQ. 1 ) ) THEN

		CALL  time_evolution ! Solve the EDQNM equation, in discrete time
		! REF-> <<< system_main >>>

	END IF

	IF ( sim_status .EQ. 1 ) THEN

		CALL post_analysis
		! REF-> <<< system_main >>>

		CALL end_run_timer
		! REF-> <<< system_timer >>>

	END IF

END PROGRAM EDQNM
