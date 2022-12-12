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

! ##################
! MODULE NAME: system_main
! LAST MODIFIED: 15 NOV 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MAIN MODULE FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_main
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * This is the main module, containing the flow of simulation.
! * Simply divided into three parts
! * -> Pre-analysis
! * -> Evolution
! * -> Inter-analysis
! * -> Post-analysis
! MAIN MODULE
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_advfunctions

  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
  IMPLICIT  NONE

  CONTAINS
! </f>

  SUBROUTINE pre_analysis
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! Call this to check the time_step, initial condition and output folder creation.
	! Does time_step check, initial condition and writing details of simulation
	! Allocating the evolution arrays, if everything is set, 'sys_status' will be 1.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	IMPLICIT NONE

	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!       T    I    M     E              S    T    E    P              C   H    E   C   K
	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	IF ( dt .LT. MIN( time_visc, time_rms_V, time_diff ) ) THEN

		sys_status = 1

		CALL allocate_edqnm_arrays
		! REF-> <<< system_basicfunctions >>>

		CALL init_edqnm_arrays
		! REF-> <<< system_basicfunctions >>>

		CALL triad_debug
		! REF-> <<< system_basicfunctions >>>

		IF ( sys_status .EQ. 1 ) THEN ! Checked again in triad debug

			CALL IC_V_large_eddies
			! CALL IC_V_read_from_file
			! en_spec_V = zero
			! REF-> <<< system_initialcondition >>>

			IF ( coupling_status .NE. 0 ) THEN
				CALL IC_B_large_eddies_single_mode
				! CALL IC_B_large_eddies
				! CALL IC_B_small_eddies
				! REF-> <<< system_initialcondition >>>
			END IF

			CALL allocate_solver_arrays
			! REF-> <<< system_solver >>>

			! CALL allocate_flux_decomposition_arrays
			! REF-> <<< system_advfunctions >>>

			CALL prepare_output
			! REF-> <<< system_basicoutput >>>

		END IF

	ELSE

		CALL print_error_timestep
		! REF-> <<< system_basicoutput >>>

	END IF

	END
! </f>

	SUBROUTINE time_evolution
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! Loop of time steps, where at each step the spectral velocities
	! are updated through any of the algoritm. Meanwhile, inter_analysis and
	! outputs are printed respectively.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		! ================================================================
		!             S        T         A         R        T
		! 8888888888888888888888888888888888888888888888888888888888888888
		CALL write_sim_start
		! REF-> <<< system_basicoutput >>>

		GOTO 922

		DO t_step = 0, t_step_total

			CALL inter_analysis

			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			!  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			IF ( coupling_status .EQ. 1 ) THEN
				CALL rk4_algorithm
			ELSE IF ( coupling_status .EQ. 0 ) THEN
				CALL rk4_algorithm_V
			ELSE IF ( coupling_status .EQ. 2 ) THEN
				CALL rk4_algorithm_B
			END IF
			! Updates velocity and magnetic field spectrum as per EDQNM-MHD equation for next time step

			IF ( nan_status .EQ. 1 ) THEN

				CALL print_error_nan
				! REF-> <<< system_basicoutput >>>
				EXIT ! Meaning 'NaN' is encountered during the Debug

			END IF

		END DO

		922 CONTINUE
		CALL prepare_perturbation_dynamo
		! REF-> <<< system_basicfunctions >>>

		DO t_step = 0, t_step_total

			CALL inter_analysis

			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			!  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			CALL rk4_algorithm
			! Updates velocity and magnetic field spectrum as per EDQNM-MHD equation for next time step

			IF ( nan_status .EQ. 1 ) THEN

				CALL print_error_nan
				! REF-> <<< system_basicoutput >>>
				EXIT ! Meaning 'NaN' is encountered during the Debug

			END IF

		END DO
		! ================================================================
		!                    E     N     D
		! 8888888888888888888888888888888888888888888888888888888888888888
		CALL write_sim_end
		! REF-> <<< system_basicoutput >>>

		sim_status = 1
		! Stating that the simulation has ended.

	END
! </f>

  SUBROUTINE inter_analysis
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! This does all the inter_analysis, making calls to write output during the evolution, debug and statistics part.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		CALL step_to_time_convert(t_step, time_now, dt)
		! Converts the 't_step' to actual time 'time_now'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  S  A  V  I  N  G    D  A  T  A
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		IF (MOD(t_step,t_step_save) .EQ. 0) THEN

			! Writes 'time_now' as a CHARACTER
			WRITE (file_time,f_d8p4) time_now

			CALL compute_transfer_term_V
			! REF-> <<< system_solver >>>

			CALL compute_kinetic_spectral_data
			! REF-> <<< system_basicfunctions >>>

			! CALL flux_decomposition
			! REF-> <<< system_advfunctions >>>

			IF ( coupling_status .NE. 0 ) THEN

				CALL compute_transfer_term_B
				! REF-> <<< system_solver >>>

				CALL compute_magnetic_spectral_data
				! REF-> <<< system_basicfunctions >>>

			END IF

		END IF

		CALL compute_temporal_data
		! REF-> <<< system_basicfunctions >>>

		IF ( forc_status .EQ. 1 ) THEN

			CALL compute_forcing_spectrum
			! REF-> <<< system_basicfunctions >>>

		END IF

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  B  U  G         F  O  R          N  a   N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF (MOD(t_step,t_step_debug) .EQ. 0) then

			CALL perform_debug
			! REF-> <<< system_basicfunctions >>>

			CALL print_running_status
			! REF-> <<< system_basicoutput >>>

		END IF

	END
! </f>

	SUBROUTINE post_analysis
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! This does all the post analysis, making calls to write output after the evolution, debug and statistics part.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		CALL deallocate_edqnm_arrays
		! REF-> <<< system_basicfunctions >>>

		CALL deallocate_solver_arrays
		! REF-> <<< system_solver >>>

		! CALL deallocate_flux_decomposition_arrays
		! REF-> <<< system_advfunctions >>>

		CALL deallocate_system_arrays
		! REF-> <<< system_basicvariables >>>

	END
! </f>

END MODULE system_main
