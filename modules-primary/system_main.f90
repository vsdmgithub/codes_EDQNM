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
! MODULE: system_main
! LAST MODIFIED: 21 JUNE 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MAIN MODULE FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_main
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This is the main module. All the other modules are sub-modules to this.
! Then nesting of all sub-modules is as follows
! MAIN MODULE
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_solver
  USE system_advfunctions

  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
	IMPLICIT  NONE
  ! _________________________
  !  VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ---------------------------------------------------------
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  CONTAINS
! </f>

  SUBROUTINE pre_analysis
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Call this to check the time_step, initial condition and output folder creation.
  ! Does time_step check, initial condition and writing details of simulation
  ! Allocating the evolution arrays, if everything is set, 'all_set' will be 1.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  IMPLICIT NONE

  !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !       T    I    M     E              S    T    E    P              C   H    E   C   K
  !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF ( dt  .LT. MIN( time_visc, time_rms ) ) THEN

    all_set =  1

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL make_initial_condition
    ! REF-> <<< system_initialcondition >>>

    ! CALL read_initial_condition
    ! REF-> <<< system_initialcondition >>>

    CALL prepare_output
    ! REF-> <<< system_basicoutput >>>

    CALL allocate_edqnm_arrays
    ! REF-> <<< system_basicfunctions >>>

    CALL init_edqnm_arrays
    ! REF-> <<< system_basicfunctions >>>

    CALL allocate_solver_arrays
    ! REF-> <<< system_solver >>>

    CALL allocate_flux_decomposition_arrays
    ! REF-> <<< system_advfunctions >>>

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
    DO t_step = 0, t_step_total

      CALL inter_analysis

      !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
      !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL rk4_algorithm
      ! Updates velocity field as per EDQNM equation for next time step

      IF ( NaN_count .GT. 0) THEN

        CALL print_error_nan
      ! REF-> <<< system_basicoutput >>>

        EXIT ! Meaning 'NaN' is encountered during the Debug

      END IF

    END DO
    ! ================================================================
    !                    E     N     D
    ! 8888888888888888888888888888888888888888888888888888888888888888

    state_sim = 1
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

      WRITE (file_time,f_d8p4) time_now
      ! Writes 'time_now' as a CHARACTER

      CALL compute_spectral_data
      ! REF-> <<< system_basicfunctions >>>

      CALL flux_decomposition
      ! REF-> <<< system_advfunctions >>>

    END IF

    CALL compute_temporal_data
    ! REF-> <<< system_basicfunctions >>>

    IF ( ( viscosity .GT. tol_float ) .AND. ( forcing_status .EQ. 1 ) ) THEN

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  F  O  R  C  I  N  G
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      net_dissipation_rate = dissipation_rate - ( energy - init_energy )/dt
      forcer               = net_dissipation_rate * forcer_template
      ! FORCING

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

    CALL write_temporal_data
    ! REF-> <<< system_basicoutput >>>

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
