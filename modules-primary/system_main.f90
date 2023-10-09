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
	INTEGER(KIND=4)::coupled_code 

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

			! CALL IC_V_large_eddies
			CALL IC_V_read_from_file
			! CALL IC_V_power_law
			! REF-> <<< system_initialcondition >>>

			CALL compute_eddy_damping
			! REF-> <<< system_basicfunctions >>>

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

		coupled_code = 0

		! ================================================================
		!             S        T         A         R        T
		! 8888888888888888888888888888888888888888888888888888888888888888
		CALL write_sim_start
		! REF-> <<< system_basicoutput >>>

		IF ( coupling_status .EQ. 2 ) THEN
			GOTO 922
			! Skip to the frozen model
		END IF

		time_now  = -dt
		save_ind  = 0
		t_step    = 0
		! INITIALIZATION OF THE LOOP

		DO WHILE ( time_now .LT. time_total )

			CALL inter_analysis

			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			!  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			CALL rk4_algorithm_V
			! REF-> <<< system_solver >>>
			! Updates kinetic energy spectrum as per EDQNM equation for next time step

			IF ( nan_status .EQ. 1 ) THEN

				CALL print_error_nan
				! REF-> <<< system_basicoutput >>>
				EXIT ! Meaning 'NaN' is encountered during the Debug

			END IF

			IF ( t_step .GT. t_step_jump ) THEN ! Helps get out the loop in the middle

				EXIT

			END IF

			t_step = t_step + 1

		END DO

		922 CONTINUE

		CALL compute_transfer_term_V
		! REF-> <<< system_solver_eqns >>>

		CALL compute_kinetic_temporal_data
		! REF-> <<< system_basicfunctions >>>

		IF ( coupling_status .NE. 0 ) THEN 

			CALL prepare_perturbation_dynamo
			! REF-> <<< system_basicfunctions >>>

		ELSE 
			GOTO 923
			! Skip the dynamo model
		END IF

		time_now = -dt
		save_ind = 0
		t_step   = 0
		! INITIALIZATION OF THE LOOP

		DO WHILE ( time_now .LT. time_total )

			coupled_code = 1

			CALL inter_analysis

			! CALL compute_adaptive_time_step
		  ! REF-> <<< system_solver_eqns >>>

			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			!  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
			!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				CALL ab4_algorithm
				! REF-> <<< system_solver >>>
				! Updates velocity and magnetic field spectrum as per EDQNM-MHD equation for next time step

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
	! This does all the inter_analysis, making calls to write output during the evolution, debug and statistics part.
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

		IF ( save_pointer .GT. save_ind ) THEN

			WRITE (file_time,f_d8p4) save_ind * time_save
			! Writes 'time_now' as a CHARACTER

			save_ind = save_ind + 1

			IF ( coupling_status .NE. 2 ) THEN

				CALL compute_transfer_term_V
					! REF-> <<< system_solver_eqns >>>

				CALL compute_kinetic_spectral_data
				! REF-> <<< system_basicfunctions >>>

				! CALL compute_flux_decomposition
				! REF-> <<< system_advfunctions >>>

			END IF

			IF ( ( coupling_status .NE. 0 ) .AND. ( coupled_code .EQ. 1 ) ) THEN

				CALL compute_transfer_term_B
				! REF-> <<< system_solver_eqns >>>

				! CALL dynamo_rate_calc
				! REF-> <<< system_solver_eqns >>>

				CALL compute_magnetic_spectral_data
				! REF-> <<< system_basicfunctions >>>

			END IF

		END IF

		IF ( coupling_status .NE. 2 ) THEN

			CALL compute_kinetic_temporal_data
			! REF-> <<< system_basicfunctions >>>

		END IF

		IF ( ( coupling_status .NE. 0 ) .AND. ( coupled_code .EQ. 1 ) ) THEN

			CALL compute_magnetic_temporal_data
			! REF-> <<< system_basicfunctions >>>

		END IF

		IF ( forc_status .EQ. 1 ) THEN

			CALL compute_forcing_spectrum
			! REF-> <<< system_basicfunctions >>>

		END IF

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
			energy_tot   =  energy_V + energy_B

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
