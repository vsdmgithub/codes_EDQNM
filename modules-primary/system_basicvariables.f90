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
! MODULE: system_basicvariables
! LAST MODIFIED: 21 JUNE 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MODULE CONTAINING ALL PRIMARY VARIABLES FOR THE CODE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_basicvariables
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space are declared and given values
! here, wheras temporary (IF necessary) are declared within the subroutines.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_constants
  USE system_auxilaries
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT  NONE
  ! _________________________
  ! SPECTRAL SPACE VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::N
	INTEGER (KIND=4)::k_ind,q_ind,p_ind
	INTEGER (KIND=4)::k2_ind
  CHARACTER(LEN=60)::N_char
  ! ---------------------------------------------------------
  DOUBLE PRECISION::wno_base,max_wno,min_wno
  DOUBLE PRECISION::log_lambda,lambda
  ! _________________________
  ! TIME (SIMULATION) VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::t_step,t_step_total
	INTEGER (KIND=4)::t_step_save
	INTEGER (KIND=4)::t_step_debug
  INTEGER (KIND=4)::t_step_forcing
	INTEGER (KIND=4)::no_of_saves,no_of_debug
  ! ---------------------------------------------------------
  DOUBLE PRECISION::time_total,time_now
  DOUBLE PRECISION::time_save,time_forcing,time_min
  DOUBLE PRECISION::time_visc,time_rms
  DOUBLE PRECISION::dt,dt_max
  ! _________________________
  ! SYSTEM VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)::state_sim,all_set,NaN_count
  INTEGER(KIND=4)::forcing_status
  INTEGER(KIND=4)::ind_integral,ind_dissipation
  INTEGER(KIND=4)::no_of_triads
  INTEGER(KIND=4)::CFL_min,CFL
  ! ---------------------------------------------------------
  DOUBLE PRECISION::viscosity
  DOUBLE PRECISION::forcing
  DOUBLE PRECISION::init_energy
  DOUBLE PRECISION::energy,enstrophy
  DOUBLE PRECISION::net_dissipation_rate,dissipation_rate,skewness
  DOUBLE PRECISION::eddy_constant
  DOUBLE PRECISION::localness_cutoff_ratio
  ! _________________________
  ! SOLVER VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION::viscous_freq,eddy_freq
  DOUBLE PRECISION::eddy_k,eddy_q,eddy_p
  DOUBLE PRECISION::integrand
  DOUBLE PRECISION::damping
  DOUBLE PRECISION::min_max_ratio
  ! _________________________
  ! GLOBAL ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::wno, wno_band
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::time
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::laplacian_k
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::en_time,es_time,ds_time,sk_time
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec,transfer_spec,flux
  ! _________________________
  ! EDQNM ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! _________________________
  DOUBLE PRECISION,DIMENSION(3)::triad_sides
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::geom_fac
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::forcer,forcer_template
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::eddy_array
  INTEGER(KIND=4),DIMENSION(:,:,:),ALLOCATABLE::kqp_status
  INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE::p_ind_min,p_ind_max
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
! </f>

  SUBROUTINE read_input
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read simulation parameters from a file 'input_file'
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER( LEN = 60 )::input_file

    input_file  = 'parameters.dat'
    ! This file contains all major Input parameters to be fed from outside file

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    OPEN(unit=1001,file=TRIM(ADJUSTL(input_file)))

    READ(1001,f_d8p4,ADVANCE='yes')
    READ(1001,f_i6,ADVANCE='yes') N
    ! No of momentum shells

    READ(1001,f_d8p4,ADVANCE='yes')
    READ(1001,f_d12p6,ADVANCE='yes') dt
    ! ! Time step size

    READ(1001,f_d8p4,ADVANCE='yes')
    READ(1001,f_d8p4,ADVANCE='yes') time_total
    ! Total time to simulate

    READ(1001,f_d8p4,ADVANCE='yes')
    READ(1001,f_i4,ADVANCE='yes') no_of_saves
    ! No of saves for the velocity data.

    CLOSE(1001)
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE init_global_variables
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       Initialize all the global variables that are used through out the code.
  ! PREREQUISITE: SUBROUTINE 'read_input' has to be called before initializing  these variables.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! S  P  E  C  T  R  U  M        A  N  D         T  I  M  E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    lambda                 = two ** ( 0.25D0 )
    ! Ratio of consecutive shells

    log_lambda             = DLOG( lambda )
    ! Logarithm of Lambda

    wno_base               = two ** ( - thr )
    ! Base wavenumber

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! NOTES:
    ! 1. This is forced viscous EDQNM, with forcing given in the first
    ! few shells matching the dissipation rate with slight fluctuations
    ! to keep it random. The time averaged net energy remains constant.
    ! 2. Viscosity levels for resolutions
    ! N45 - Minimum of 0.0005.
    ! 3. Eddy constant is generally not changed.
    ! 4. Two timescales are derived, one from net energy, other from viscosity
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    init_energy            = one
    ! INitial energy

    ind_integral           = CEILING( DBLE( N ) / 10.0D0 )
    ! Index (position) of integral scale

    ind_dissipation        = N - FLOOR( DBLE( N ) / 4 )
    ! Index (position) of dissipation scale

    viscosity              = zero
    ! Viscosity

    min_wno                = wno_base * lambda
    ! Min wave number

    max_wno                = wno_base * ( lambda ** ( N - 1) )
    ! Max wave number

    time_rms               = one / DSQRT( init_energy * ( max_wno ** two ) )
    ! Time scale from energy and largest momentum

    time_visc              = one / ( viscosity * ( max_wno ** two ) + tol_float )
    ! Time scale from viscosity and largest momentum, tol_float is added to avoid NaN in case of inviscid

    CFL_min                = 10
    ! Minimum of CFL

    time_min               = MIN( time_rms, time_visc )

    time_forcing           = time_min
    ! Time scale of forcing

    dt_max                 = time_min / DBLE( CFL_min )
    ! Maximum time step to satisfy the CFL condition.

    CALL find_time_step( dt_max, dt )
    ! REF-> <<< system_auxilaries >>>

    CFL                    = time_min / dt
    ! Actual CFL of the system

    CALL time_to_step_convert( time_total, t_step_total, dt)
    ! REF-> <<< system_auxilaries >>>
    ! Converts time to steps

    CALL time_to_step_convert( time_forcing, t_step_forcing, dt)
    ! REF-> <<< system_auxilaries >>>

  	t_step_save            = t_step_total / no_of_saves
    ! Determines how many time steps after the save has to be made.

    no_of_debug            = 10
    ! No of times the data checked for any NaN

  	t_step_debug           = t_step_total / no_of_debug
    ! No of times, the data will be checked for 'NaN' during the simul

    CALL step_to_time_convert( t_step_save, time_save, dt)
    ! REF-> <<< system_auxilaries >>>
    ! Converts steps to time

    WRITE (N_char, f_i8) N
    ! converting resolution value to CHARACTER

    eddy_constant          = 0.54D0
    ! Eddy constant in its expression

    localness_cutoff_ratio = 0.4
    ! Ratio of min to max triad sides, to say it is a nonlocal triad interactions
    ! This is userdefined . Has to be << 1 is must.

    forcing_status         = 0
    ! '1' to activate forcing, '0' to deactivate forcing

    state_sim              = 0
    ! This being the first variable to start the simulation. At last, it will be set to '1'

    all_set                = 0
    ! Checking the readiness of system

    NaN_count              = 0
    ! No of NaN count

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

	SUBROUTINE init_global_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       Initialize all the global arrays that are declared/allotted here.
  ! PREREQUISITE: SUBROUTINE 'init_global_variables' has to be called before initializing  these arrays.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    CALL allocate_system_arrays

    ! Defining the band for each discrete k_i as k^{u}_i,k^{l}_i  for above and below
    ! EQN:- $\Delta k_i=k_i \ln(\lambda)$
    ! EQN:- $k_i^{l}=\frac{\Delta k_i}{\lambda -1}
    ! EQN:- $k_i^{u}=\frac{\lambda \Delta k_i}{\lambda -1}

    DO k_ind                = 1, N

      wno( k_ind )          = wno_base * ( lambda ** (k_ind - 1) )
      ! Value of wnoentum for the corresponding index.

      wno_band( k_ind )     = wno( k_ind ) * log_lambda
      ! Gap around the wnoentum

      laplacian_k( k_ind )  = wno( k_ind ) ** two
      ! Laplacian in spectral space

    END DO

  END
! </f>

	SUBROUTINE allocate_system_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       allocate all system arrays
  ! PREREQUISITE: SUBROUTINE 'init_global_variables' has to be called before initializing  these arrays.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y     A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( wno( N ) , wno_band ( N ) )
    ALLOCATE( laplacian_k( N ) )
    ALLOCATE( spec( N ), transfer_spec( N ), flux( N ) )

  END
! </f>

	SUBROUTINE deallocate_system_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       deallocate all system arrays
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y     D  E   A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( wno, wno_band )
    DEALLOCATE( laplacian_k, spec, transfer_spec, flux )

  END
! </f>

END MODULE system_basicvariables
