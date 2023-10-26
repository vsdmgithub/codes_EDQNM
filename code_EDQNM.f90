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
! PROGRAM NAME : EDQNM_MHD_alpha_D3
! _____________________________________
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! PROGRAM FOR SOLVING 3-DIMENSIONAL EDQNM-MHD ALPHA VARIANT EQUATION  
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! </f>
PROGRAM EDQNM_MHD_alpha_D3
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! - This code solves the 3-dimensional Eddy Damped Quasi-Normal Markovian Magneto Hydrodynamical (EDQNM-MHD) equation, refer to *Orsag 1970, Pouquet 1976* for more details.
! - In this model, a variation is implemented in the eddy damping timescale that is parameterised by $\alpha$, hence the name.
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

!==================================================================
	! FOR CHANGING KINETIC AND MAGNETIC REYNOLDS NUMBER IN BATCH SCRIPTS
	! U_GRID=XXX
	! W_GRID=YYY

	! USE THESE FOR A SINGLE KINETIC AND MAGNETIC REYNOLDS NUMBER
	U_GRID = 1
	W_GRID = 1
!==================================================================

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

		CALL  time_evolution ! Solve the equation, in discrete time steps
		! REF-> <<< system_main >>>

	END IF

	IF ( sim_status .EQ. 1 ) THEN

		CALL post_analysis
		! REF-> <<< system_main >>>

		CALL end_run_timer
		! REF-> <<< system_timer >>>

	END IF

END PROGRAM EDQNM_MHD_alpha_D3
