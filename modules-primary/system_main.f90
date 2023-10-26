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
! LAST MODIFIED: 15 OCT 2023
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MAIN MODULE FOR EDQNM-MHD EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_main
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * This is the main module, containing the flow of the simulation.
! * Divided into four main parts (SUBROUTINES):
! * -> pre_analysis
! * -> time_evolution
! * -> inter_analysis
! * -> post_analysis

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
	! Call this to do the following 
	! * -> Check the value of time step
	! * -> Initializes the parameters in the EDQNM problem
	! If everything is set, 'sys_status' will be 1.
	! * -> Initialize the initial condition
	! * -> Output folder creation.
	! * -> Allocating the evolution arrays, 
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	IMPLICIT NONE

	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!       T    I    M     E              S    T    E    P              C   H    E   C   K
	!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	IF ( dt .LT. MIN( time_visc, time_rms_V, time_rms_B , time_diff ) ) THEN

		sys_status = 1

		CALL allocate_edqnm_arrays
		! REF-> <<< system_basicfunctions >>>

		CALL init_edqnm_arrays
		! REF-> <<< system_basicfunctions >>>

		CALL triad_debug
		! REF-> <<< system_basicfunctions >>>

		IF ( sys_status .EQ. 1 ) THEN ! Checked again in triad debug

			CALL prepare_initial_condition
			! REF-> <<< system_basicfunctions >>>

			CALL allocate_solver_arrays
			! REF-> <<< system_solver >>>

			CALL allocate_flux_decomposition_arrays
			! REF-> <<< system_advfunctions >>>

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

		coupled_code = 0

		! ================================================================
		!             S        T         A         R        T
		! 8888888888888888888888888888888888888888888888888888888888888888
		CALL write_sim_start
		! REF-> <<< system_basicoutput >>>

		!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  K I N E T I C    S P E C T R U M    E V O L U T I O N
		!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		! INITIALIZATION OF THE LOOP
		time_now  = -dt
		save_ind  = 0
		t_step    = 0

		DO WHILE ( time_now .LT. time_total )

			CALL inter_analysis

			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			!  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			! CALL rk4_algorithm_V
			CALL ab4_algorithm_V
			! REF-> <<< system_solver >>>
			! Updates kinetic energy spectrum as per EDQNM equation for next time step

			IF ( nan_status .EQ. 1 ) THEN

				CALL print_error_nan
				! REF-> <<< system_basicoutput >>>
				EXIT ! Meaning 'NaN' is encountered during the Debug

			END IF

			IF ( t_step .GT. t_step_jump ) THEN 
			! Helps get out the loop in the middle, set it in <<< system_basicvariables >>>
				EXIT
			END IF

			t_step = t_step + 1

		END DO

		CALL compute_transfer_term_V
		! REF-> <<< system_solver_eqns >>>

		IF ( coupling_status .NE. 0 ) THEN 

			CALL prepare_perturbation_dynamo
			! Gets the magnetic spectrum initialized and ready for evolution, either coupled or frozen
			! REF-> <<< system_basicfunctions >>>

			coupled_code = 1
			! Stating that the coupling is happening.

		ELSE 

			GOTO 923
			! Skip the MHD evolution

		END IF

		!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  K I N E T I C    A N D     M A G N E T I C    S P E C T R U M
		!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		! INITIALIZATION OF THE LOOP
		time_now  = -dt
		save_ind  = 0
		t_step    = 0

		DO WHILE ( time_now .LT. time_total )

			CALL inter_analysis

			! CALL compute_adaptive_time_step
			! Use it if the dt is chosen large, and energy spectrum becomes negative at some time.
		  ! REF-> <<< system_solver_eqns >>>

			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			!  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			IF ( coupling_status .EQ. 1 ) THEN 

				! CALL rk4_algorithm
				CALL ab4_algorithm
				! REF-> <<< system_solver >>>
				! Updates velocity and magnetic field spectrum as per EDQNM-MHD equation for next time step

			ELSE 

				! CALL rk4_algorithm_B
				CALL ab4_algorithm_B
				! REF-> <<< system_solver >>>
				! Updates only the magnetic field spectrum as per EDQNM-MHD equation for next time step

			END IF

			IF ( nan_status .EQ. 1 ) THEN

				CALL print_error_nan
				! REF-> <<< system_basicoutput >>>
				EXIT ! Meaning 'NaN' is encountered during the Debug

			END IF

			t_step = t_step + 1

		END DO
		! ================================================================
		!                    E     N     D
		! 8888888888888888888888888888888888888888888888888888888888888888

		923 CONTINUE 
		
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
	! This does the analysis of functions in between the time steps,
	! making calls to write output during the evolution, debug and statistics part.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		time_now = time_now + dt_cur
		! Converts the 't_step' to actual time 'time_now'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  S  A  V  I  N  G    D  A  T  A
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		save_pointer = CEILING( ( time_now + tol_double ) / time_save )
		! A pointer that shows when the data has to be saved 

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		! SAVING SPECTRUM AT SELECTED INTERVALS
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF ( save_pointer .GT. save_ind ) THEN

			WRITE (file_time,f_d8p4) save_ind * time_save
			! Writes 'time_now' as a CHARACTER

			save_ind = save_ind + 1

			IF ( coupled_code .EQ. 0 ) THEN

				CALL compute_transfer_term_V
					! REF-> <<< system_solver_eqns >>>

				CALL compute_kinetic_spectral_data
				! REF-> <<< system_basicfunctions >>>

				CALL compute_flux_decomposition_V
				! REF-> <<< system_advfunctions >>>

			ELSE IF ( coupling_status .EQ. 1 ) THEN

				CALL compute_transfer_term_V
					! REF-> <<< system_solver_eqns >>>

				CALL compute_transfer_term_B
				! REF-> <<< system_solver_eqns >>>

				CALL compute_kinetic_spectral_data
				! REF-> <<< system_basicfunctions >>>

				CALL compute_magnetic_spectral_data
				! REF-> <<< system_basicfunctions >>>

				CALL compute_flux_decomposition_B_self
				CALL compute_flux_decomposition_B_intr
				! REF-> <<< system_advfunctions >>>

				! CALL dynamo_rate_calc
				! REF-> <<< system_solver_eqns >>>

			ELSE IF ( coupling_status .EQ. 2 ) THEN

				CALL compute_transfer_term_B
				! REF-> <<< system_solver_eqns >>>

				CALL compute_magnetic_spectral_data
				! REF-> <<< system_basicfunctions >>>

				CALL compute_flux_decomposition_B_self
				CALL compute_flux_decomposition_B_intr
				! REF-> <<< system_advfunctions >>>

				! CALL dynamo_rate_calc
				! REF-> <<< system_solver_eqns >>>

			END IF

		END IF

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  SAVING TEMPORAL DATA AT EVERY TIMESTEP
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF ( coupled_code .EQ. 0 ) THEN

			CALL compute_kinetic_temporal_data
			! REF-> <<< system_basicfunctions >>>

		ELSE IF ( coupling_status .EQ. 1 ) THEN

			CALL compute_kinetic_temporal_data
			! REF-> <<< system_basicfunctions >>>

			CALL compute_magnetic_temporal_data
			! REF-> <<< system_basicfunctions >>>

		ELSE IF ( coupling_status .EQ. 2 ) THEN

			CALL compute_magnetic_temporal_data
			! REF-> <<< system_basicfunctions >>>

		END IF

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  CALCULATING FORCING
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF ( forc_status .EQ. 1 ) THEN
			IF ( coupled_code .EQ. 0 ) THEN
				CALL compute_forcing_spectrum
				! REF-> <<< system_basicfunctions >>>
			ELSE IF( coupling_status .EQ. 1) THEN
				CALL compute_forcing_spectrum
				! REF-> <<< system_basicfunctions >>>
			END IF
			! No need for the frozen model in the second half
		END IF

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  CALCULATING EDDY DAMPING TIMESCALES
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		CALL compute_eddy_damping
		! REF-> <<< system_basicfunctions >>>

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  B  U  G         F  O  R          N  a   N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF (MOD(t_step,t_step_debug) .EQ. 0) then

			CALL perform_debug
			! REF-> <<< system_basicfunctions >>>

			energy_V     = SUM( en_spec_V * wno_band )
			energy_B     = SUM( en_spec_B * wno_band )
			energy_tot   = energy_V + energy_B

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

		CALL deallocate_flux_decomposition_arrays
		! REF-> <<< system_advfunctions >>>

		CALL deallocate_system_arrays
		! REF-> <<< system_basicvariables >>>

	END
! </f>

END MODULE system_main
