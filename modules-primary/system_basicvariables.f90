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
! MODULE NAME  : system_basicvariables
! LAST MODIFIED: 15 NOV 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MODULE CONTAINING ALL PRIMARY VARIABLES AND ARRAYS FOR THE CODE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_basicvariables
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * All the global variables for the simulation are declared and given values for the needed variables.
! * All the arrays are declared here. Allocation and Deallocation of system arrays is done here.
! * Other arrays are allocated/deallocated in "system_basicfunctions"
! * Other temporary (if necessary) are declared within the subroutines/functions.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
!  SUB-MODULES
!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
USE system_auxilaries
! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

IMPLICIT  NONE
	! _________________________
	! SPECTRAL SPACE VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::N
	INTEGER (KIND=4)::k_ind,q_ind,p_ind
	INTEGER (KIND=4)::k2_ind,dum_ind
	CHARACTER(LEN=60)::N_char,dim_char
	! ---------------------------------------------------------
	DOUBLE PRECISION::dim,dim_min_3
	DOUBLE PRECISION::wno_base,wno_max,wno_min
	DOUBLE PRECISION::wno_forc,wno_diss,wno_int
	DOUBLE PRECISION::wno_scale_log,wno_scale
	! _________________________
	! TIME (SIMULATION) VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::t_step,t_step_total
	INTEGER (KIND=4)::t_step_save
	INTEGER (KIND=4)::t_step_debug
	INTEGER (KIND=4)::no_of_saves,no_of_debug
	! ---------------------------------------------------------
	DOUBLE PRECISION::time_total,time_now
	DOUBLE PRECISION::time_save
	DOUBLE PRECISION::time_visc,time_rms_V
	DOUBLE PRECISION::time_diff,time_rms_B
	DOUBLE PRECISION::dt,dt_max
	! _________________________
	! SYSTEM VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER(KIND=4)::nan_count
	INTEGER(KIND=4)::sim_status,sys_status,nan_status
	INTEGER(KIND=4)::visc_status,diff_status,forc_status
	INTEGER(KIND=4)::coupling_status
	INTEGER(KIND=4)::kI_ind,kD_ind,kF_ind
	INTEGER(KIND=4)::triad_count
	INTEGER(KIND=4)::triad_deleted
	INTEGER(KIND=4)::cfl_sys
	! ---------------------------------------------------------
	DOUBLE PRECISION::visc,diff,prandl_no
	DOUBLE PRECISION::forcing_factor
	DOUBLE PRECISION::energy_V_0,energy_B_0
	DOUBLE PRECISION::energy_V,energy_B
	DOUBLE PRECISION::energy_tot
	DOUBLE PRECISION::energy_0
	DOUBLE PRECISION::enstrophy_V,enstrophy_B
	DOUBLE PRECISION::ds_rate_V,ds_rate_B
	DOUBLE PRECISION::ds_rate_ref_V
	DOUBLE PRECISION::skewness
	DOUBLE PRECISION::skewness_const
	DOUBLE PRECISION::er_V_self,er_B_self,er_VB
	DOUBLE PRECISION::dim_const
	DOUBLE PRECISION::kol_const
	DOUBLE PRECISION::fback_coef
	! _________________________
	! SOLVER VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION::eddy_const,magn_const
	DOUBLE PRECISION::integrand_V,integrand_B
	DOUBLE PRECISION::eddy_damping
	! _________________________
	! GLOBAL ARRAYS
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::wno, wno_band
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::wno_right, wno_left
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::laplacian_k
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::en_spec_V,tr_spec_V,fl_spec_V
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::en_spec_B,tr_spec_B,fl_spec_B
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::fr_spec,spec0
	! _________________________
	! EDQNM ARRAYS
	! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! _________________________
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::geom_b,geom_c,geom_h,triad_weightage
	INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE   ::p_ind_min,p_ind_max
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

		input_file  = 'input.dat'
		! This file contains all major Input parameters to be fed from outside file

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		OPEN(unit=1001,file=TRIM(ADJUSTL(input_file)))

		READ(1001,f_d8p4,ADVANCE='yes')
		READ(1001,f_i6,ADVANCE='yes') N
		! No of momentum shells

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
		! _________________________
		! LOCAL VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION::time_min,visc_ref,diff_ref
		INTEGER(KIND=4)::N_ref,cfl_ref

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! NOTES:
		! 1. This is forced viscous EDQNM, with forcing given in the first
		! few shells matching the dissipation rate.
		! 2. Viscosity levels for resolutions
		! N45 - Minimum of 0.0005.
		! 3. Eddy constant is generally not changed.
		! 4. Two timescales are derived, one from net energy, other from viscosity
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		dim                                   = 3.0D0
		dim_min_3                             = dim - thr
		! Dimension of the space in which EDQNM is computed

		visc_status                           = 1
		! '1' to include viscosity, '0' to do inviscid

		diff_status                           = 1
		! '1' to include diffusivity, '0' to do inviscid

		forc_status                           = 1
		! '1' to activate forcing, '0' to deactivate forcing (only for kinetic spectrum)

		coupling_status                       = 0
		! '1' for MHD EDQNM, '0' for only kinetic EDQNM, '2' for only MHD with fixed E(k)

		N_ref                                 = 45
		! Reference resolution

		! visc_ref                            = 5E-4
		visc_ref                              = 1E-3
		! Viscosity standard (minimum) for N  =45

		wno_scale                             = two ** ( 0.25D0 )

		visc                                  = visc_ref * ( wno_scale ** ( N_ref - N ) )
		! Adjusted minimum viscosity for the current N

		! visc                                = 0.001
		! UNCOMMENT FOR CUSTOM VISCOSITY

		diff_ref                              = 1E-3
		! Reference diffusivity

		diff                                  = diff_ref * ( wno_scale ** ( N_ref - N ) )
		! Adjusted minimum diffusivity for the current N

		! diff                                = 0.001
		! UNCOMMENT FOR CUSTOM VISCOSITY

		prandl_no															= visc / diff
		! Prandl number

		energy_0 															= one
		energy_V_0                            = energy_0
		! Initial kinetic energy

		energy_B_0                            = 1E-5
		! Initial magnetic energy

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! S P E C T R U M A N D T I M E
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		wno_scale_log                         = DLOG( wno_scale )
		! Ratio of consecutive shells

		wno_base                              = two ** ( - thr )
		! Base wavenumber

		wno_min                               = wno_base
		! Min wave number

		wno_max                               = wno_base * ( wno_scale ** ( N - 1) )
		! Max wave number

		kI_ind                                = 2
		! Index (position) of integral scale

		kF_ind                                = 4
		! Index (position) of forcing scale

		kD_ind                                = N - 10
		! Index (position) of dissipation scale

		ds_rate_ref_V                         = one
		! REF-> compute_forcing_spectrum in  <<< system_basicfunctions >>>

		fback_coef                            = 0.4
		! REF-> compute_forcing_spectrum in  <<< system_basicfunctions >>>
		! Feedback of current energy trend to force accordingly, '0' means no feedback

		time_rms_V                            = one / DSQRT( energy_V_0 * ( wno_max ** two ) )
		time_rms_B                            = one / DSQRT( energy_B_0 * ( wno_max ** two ) )
		! Time scale from energy and largest momentum

		time_visc                             = one / ( visc * ( wno_max ** two ) + tol_float )
		time_diff                             = one / ( diff * ( wno_max ** two ) + tol_float )
		! Time scales from viscosity and diffusivity

		cfl_ref                               = 20
		! Minimum of CFL

		time_min                              = MIN( time_rms_V, time_rms_B, time_visc, time_diff )

		dt_max                                = time_min / DBLE( cfl_ref )
		! Maximum time step to satisfy the CFL condition.

		CALL find_time_step( dt_max, dt )
		! REF-> <<< system_auxilaries >>>

		! dt                                  = 0.005
		! UNCOMMENT TO GIVE CUSTOM 'dt'

		cfl_sys                               = FLOOR( time_min / dt )
		! Actual CFL of the system

		CALL time_to_step_convert( time_total, t_step_total, dt)
		! REF-> <<< system_auxilaries >>>

		t_step_save                           = t_step_total / no_of_saves
		! Determines how many time steps after the save has to be made.

		no_of_debug                           = 5
		! No of times the data checked for any NaN

		t_step_debug                          = t_step_total / no_of_debug
		! No of times, the data will be checked for 'NaN' during the simul

		CALL step_to_time_convert( t_step_save, time_save, dt)
		! REF-> <<< system_auxilaries >>>

		WRITE (N_char, f_i8) N
		! converting resolution value to CHARACTER

		WRITE (dim_char, f_d5p2) dim
		! converting dimension to CHARACTER

		kol_const                             = 1.7D0
		eddy_const                            = 0.19D0 * DSQRT( ( kol_const ) ** thr )
		magn_const                            = DSQRT( two / thr )

		skewness_const                        = DSQRT(135.0D0/98.0D0)
		! Constant appearing in the calc. of skewness

		dim_const                             = 8.0D0 * solid_angle( dim - two ) / solid_angle( dim - one )
		! REF-> <<< system_auxilaries >>>
		! Const of integration in 'd' dimension for the transfer term

		sim_status                            = 0
		! This being the first variable to start the simulation. At last, it will be set to '1'

		sys_status                            = 0
		! Checking the parameters of system, after pre_analysis subroutine it will be set to '1'

		nan_count                             = 0
		nan_status                            = 0
		! No of NaN count and status of NaN (if any)

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
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION::s_exp,dum

		CALL allocate_system_arrays

		! NOTE:==================================>>>
		! Defining the band for each discrete k_i as k^{u}_i,k^{l}_i  for above and below
		! EQN:- $\Delta k_i=k_i \ln(\wno_scale)$
		! EQN:- $k_i^{l}=\frac{\Delta k_i}{\wno_scale -1}
		! EQN:- $k_i^{u}=\frac{\wno_scale \Delta k_i}{\wno_scale -1}
		! =======================================<<<
		DO k_ind                = 1, N

			wno( k_ind )          = wno_base * ( wno_scale ** DBLE(k_ind - 1) )
			! Value of momentum for the corresponding index.

			wno_band( k_ind )     = wno( k_ind ) * wno_scale_log
			! Gap around the momentum

			wno_left( k_ind )     = wno_band( k_ind ) / ( wno_scale - one )
			! Right edge around the momentum band

			wno_right( k_ind )    = wno_scale * wno_band( k_ind ) / ( wno_scale - one )
			! Left edge around the momentum band

			laplacian_k( k_ind )  = wno( k_ind ) ** two
			! Laplacian in spectral space

		END DO

		wno_forc                              = wno( kF_ind )
		! Wavenumber of forcing scale

		wno_int                               = wno( kI_ind )
		! Wavenumber of forcing scale

		wno_diss                              = wno( kD_ind )
		! Wavenumber of dissipation scale

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  F  O  R  C  I  N  G       T  E  M  P  L  A  T  E
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		s_exp = 4.0D0 ! Integral scale spectrum exponent
		dum   = hf / ( wno_int ** two )
		spec0 = ( (wno / wno_int) ** s_exp ) * DEXP( - dum * laplacian_k )
		spec0 = spec0 / SUM( spec0 * wno_band )
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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
		ALLOCATE( wno_right( N ) , wno_left ( N ) )
		ALLOCATE( laplacian_k( N ) )
		ALLOCATE( spec0( N ), en_spec_V( N ), en_spec_B( N ) )
		ALLOCATE( tr_spec_V( N ), tr_spec_B( N ) )
		ALLOCATE( fl_spec_V( N ), fl_spec_B( N ) )
		IF ( forc_status .EQ. 1 ) THEN
			ALLOCATE( fr_spec( N ) )
		END IF

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
		DEALLOCATE( wno , wno_band )
		DEALLOCATE( wno_right , wno_left  )
		DEALLOCATE( laplacian_k )
		DEALLOCATE( spec0, en_spec_V, en_spec_B )
		DEALLOCATE( tr_spec_V, tr_spec_B, fl_spec_V, fl_spec_B )
		IF ( forc_status .EQ. 1 ) THEN
			DEALLOCATE( fr_spec )
		END IF

	END
! </f>

END MODULE system_basicvariables
