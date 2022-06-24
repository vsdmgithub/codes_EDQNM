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
! PROGRAM : EDQNM
! LAST MODIFIED: 21 JUNE 2022
! _____________________________________
! LIST OF MODULES USED :
! 1. main
! 2. solver
! 3. initialcondition
! 4. basicvariables
! 5. basicfunctions
! 6. basicoutput
! 7. advfunctions
! 8. timer
! 9. constants
! 10.auxilaries
! _____________________________________
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! PROGRAM FOR SOLVING EDQNM EQUATION
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! </f>
PROGRAM EDQNM
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This program solves the EDQNM equation
! All the work is done in the modules. Calling a few would finish the code.
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
	! REF-> <<< system_basicvariables >>>

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

  IF ( (run .EQ. 'y') .AND. ( all_set .EQ. 1 ) ) THEN

    PRINT*,'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
    PRINT*,'   S  I  M  U  L  A  T  I  O  N        S  T  A  R  T  S '
    PRINT*,'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'

    CALL  time_evolution ! Solve the EDQNM equation, in discrete time
		! REF-> <<< system_main >>>

    PRINT*,'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    PRINT*,'   S  I  M  U  L  A  T  I  O  N        E  N  D  S '
    PRINT*,'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'

  END IF

  IF ( state_sim .EQ. 1 ) THEN

	  CALL post_analysis ! Does the post-analysis, deallocating
		! REF-> <<< system_main >>>

	  CALL end_run_timer
		! REF-> <<< system_timer >>>

	END IF

END PROGRAM EDQNM
