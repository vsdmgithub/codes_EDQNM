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
! MODULE NAME  : system_basicoutput
! LAST MODIFIED: 15 NOV 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MODULE CONTAINING ALL OUTPUT FROM THE SIMULATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_basicoutput
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! * Directories are created respectively
! * Spectral outputs
! * Temporal outputs
! * Error messages(if any) and Running status of the simulation
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_basicvariables
	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT  NONE
	! _________________________
	! OUTPUT VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!
	CHARACTER(LEN =140)::file_name
	CHARACTER(LEN =40) ::file_time
	CHARACTER(LEN =40) ::path_dir
	CHARACTER(LEN =80) ::type_sim
	CHARACTER(LEN =40) ::name_sim
	CHARACTER(LEN =100)::file_address
	CHARACTER(LEN =40) ::sub_dir_sp
	CHARACTER(LEN =40) ::sub_dir

	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	CONTAINS
! </f>

	SUBROUTINE prepare_output
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Name the folders, create them, open files to write system_basicoutput .
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE

		CALL name_output_dir
		! Names all the directories where output is stored

		CALL create_output_dir
		! Creates the directories

		CALL write_simulation_details
		! Writes the parameters used in the simulation

	END
! </f>

	SUBROUTINE name_output_dir
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Name the folders, file address
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT  NONE

		path_dir    =   '../data/'
		! path of the main directory relative to this file.

		sub_dir_sp  =   'spec/'
		! Sub directory name to store spectral data

		type_sim    =  'N' // TRIM( ADJUSTL( N_char ) ) // '/'
		! type_sim    =  'D' // TRIM( ADJUSTL( dim_char ) ) // '_U' // TRIM( ADJUSTL( U_char ) ) // 'W' // TRIM( ADJUSTL( W_char ) ) // '/'
		! type of simulation, the data is storing

		name_sim    =  'D' // TRIM( ADJUSTL( dim_char ) )
		! CALL get_simulation_name(name_sim)
		! REF-> <<< system_auxilaries >>>
		! Creating dated and timed name for the simulation for this particular type

		! name_sim    =   'test_sim'
		! Use this to give CUSTOM SIMULATION NAME

		file_address =   TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) //  &
		                 TRIM( ADJUSTL( name_sim ) ) // '/kin/'
		! Address should be added to all file names, if needed sub-dir can be declared later and appended to
		! this address

	END
! </f>

	SUBROUTINE create_output_dir
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Create the directories.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT  NONE

		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( path_dir ) ) )
		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) )
		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) // TRIM( ADJUSTL( name_sim ) ) )
		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) )
		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) )
		! Command to create the main directory and sub directories (name_sim) in the desired path
		! If exists already, it won't be an error

	END
! </f>

	SUBROUTINE prepare_output_dynamo
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Create the directories for dynamo study after EDQNM steady state.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT  NONE

		file_address =   TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) //  &
		                 TRIM( ADJUSTL( name_sim ) ) // '/dyn/'
		! Address should be added to all file names, if needed sub-dir can be declared later and appended to
		! this address
		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) )
		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) )
		! Command to create the main directory and sub directories (name_sim) in the desired path
		! If exists already, it won't be an error

		! CALL write_simulation_details
		! Writes the parameters used in the simulation

		WRITE(*,'(A63)')'-----------------------------------------------------------'
		WRITE(*,'(A63)')' DYNAMO EFFECT TESTING: PERTURBATION SEEDED IN MAG. SPECTRUM'
		WRITE(*,'(A63)')'-----------------------------------------------------------'

	END
! </f>

	SUBROUTINE write_simulation_details
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	!   Write the details of the simulation
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT  NONE

		file_name = TRIM(ADJUSTL(file_address))//'system_details'
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!   S  I  M  U  L  A  T  I  O  N     D  E  T  A  I  L  S
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OPEN(UNIT =233,FILE=TRIM( ADJUSTL( file_name ) ) // '.dat')
		WRITE(233,"(A40)")TRIM(ADJUSTL('--------------------------------------------------------------------'))
		WRITE(233,"(A40)")TRIM(ADJUSTL('------  EDQNM  EQUATION----------------------'))
		WRITE(233,"(A40)")TRIM(ADJUSTL(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'))
		WRITE(233,"(A40)")TRIM(ADJUSTL('-----------PARAMETERS OF SIMULATION------------'))
		WRITE(233,"(A40)")TRIM(ADJUSTL('--------------------------------------------------------------------'))
		WRITE(233,"(A1,A20,A2,I8)")      '*',' No of modes    ',                '= ',N
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' Dimension of sys',               '= ',dim
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' Time step   ',                   '= ',dt
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' Total time ',                    '= ',time_total
		WRITE(233,"(A1,A20,A2,I8)")      '*',' Total time steps   ',            '= ',t_step_total
		WRITE(233,"(A1,A20,A2,I8)")      '*',' Forcing status  ',               '= ',forc_status
		WRITE(233,"(A1,A20,A2,I8)")      '*',' Viscosity status   ',            '= ',visc_status
		WRITE(233,"(A1,A20,A2,I8)")      '*',' Diffusivity status   ',          '= ',diff_status
		WRITE(233,"(A1,A20,A2,F8.6)")    '*',' Viscosity  ',                    '= ',visc
		WRITE(233,"(A1,A20,A2,F8.6)")    '*',' Diffusivity  ',                  '= ',diff
		WRITE(233,"(A1,A20,A2,F8.6)")    '*',' Prandl No   ',                   '= ',prandl_no
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' Viscous timescale  ',            '= ',time_visc
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' Diffusive timescale  ',          '= ',time_diff
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' K.Energy timescale  ',           '= ',time_rms_V
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' M.Energy timescale  ',           '= ',time_rms_B
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' K. Diff rate        ',           '= ',ds_rate_ref_V
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' B. Diff rate        ',           '= ',ds_rate_ref_B
		WRITE(233,"(A1,A20,A2,I8)")      '*',' No of saves   ',                 '= ',no_of_saves
		WRITE(233,"(A1,A20,A2,I8)")      '*',' CFL System   ',                  '= ',cfl_sys
		WRITE(233,"(A1,A20,A2,F8.3)")    '*',' Initial kin.energy ',            '= ',energy_V_0
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' Initial mag.energy ',            '= ',energy_B_0
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' Smallest wavenumber',            '= ',wno_min
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' Largest wavenumber ',            '= ',wno_max
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' K.Diss. wavenumber  ',           '= ',wno( kD_V_ind )
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' M.Diss. wavenumber  ',           '= ',wno( kD_B_ind )
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' Int. wavenumber  ',              '= ',wno( kI_ind )
		WRITE(233,"(A1,A20,A2,F8.2)")    '*',' Forc wavenumber  ',              '= ',wno( kF_ind )
		WRITE(233,"(A1,A20,A2,I8)")      '*',' Total Triad count ',             '= ',triad_count
		WRITE(233,"(A1,A20,A2,I8)")      '*',' Deleted triads ',                '= ',triad_deleted
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' Error in G.fac "b" ',            '= ',er_V_self
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' Error in G.fac "c" ',            '= ',er_VB
		WRITE(233,"(A1,A20,A2,ES8.2)")   '*',' Error in G.fac "h" ',            '= ',er_B_self
		WRITE(233,"(A1,A20,A2,F8.3)")    '*',' Eddy const',                     '= ',eddy_const
		WRITE(233,"(A1,A20,A2,F8.3)")    '*',' Dim const',                      '= ',dim_const

		CLOSE(233)
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		WRITE(*,"(A40)") TRIM(ADJUSTL('--------------------------------------------------------------------'))
		WRITE(*,"(A40)") TRIM(ADJUSTL('------  EDQNM  EQUATION----------------------'))
		WRITE(*,"(A40)") TRIM(ADJUSTL(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'))
		WRITE(*,"(A40)") TRIM(ADJUSTL('-----------PARAMETERS OF SIMULATION------------'))
		WRITE(*,"(A40)") TRIM(ADJUSTL('--------------------------------------------------------------------'))
		WRITE(*,"(A1,A20,A2,I8)")      '*',' No of modes    ',                '= ',N
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' Dimension of sys',               '= ',dim
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' Time step   ',                   '= ',dt
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' Total time ',                    '= ',time_total
		WRITE(*,"(A1,A20,A2,I8)")      '*',' Total time steps   ',            '= ',t_step_total
		WRITE(*,"(A1,A20,A2,I8)")      '*',' Forcing status   ',              '= ',forc_status
		WRITE(*,"(A1,A20,A2,I8)")      '*',' Viscosity status   ',            '= ',visc_status
		WRITE(*,"(A1,A20,A2,I8)")      '*',' Diffusivity status   ',          '= ',diff_status
		WRITE(*,"(A1,A20,A2,F8.6)")    '*',' Viscosity  ',                    '= ',visc
		WRITE(*,"(A1,A20,A2,F8.6)")    '*',' Diffusivity  ',                  '= ',diff
		WRITE(*,"(A1,A20,A2,F8.6)")    '*',' Prandl No   ',                   '= ',prandl_no
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' Viscous timescale  ',            '= ',time_visc
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' Diffusive timescale  ',          '= ',time_diff
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' K.Energy timescale  ',           '= ',time_rms_V
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' M.Energy timescale  ',           '= ',time_rms_B
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' K. Diff rate        ',           '= ',ds_rate_ref_V
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' B. Diff rate        ',           '= ',ds_rate_ref_B
		WRITE(*,"(A1,A20,A2,I8)")      '*',' No of saves   ',                 '= ',no_of_saves
		WRITE(*,"(A1,A20,A2,I8)")      '*',' CFL System   ',                  '= ',cfl_sys
		WRITE(*,"(A1,A20,A2,F8.3)")    '*',' Initial kin.energy ',            '= ',energy_V_0
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' Initial mag.energy ',            '= ',energy_B_0
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' Smallest wavenumber',            '= ',wno_min
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' Largest wavenumber ',            '= ',wno_max
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' K.Diss. wavenumber  ',           '= ',wno( kD_V_ind )
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' M.Diss. wavenumber  ',           '= ',wno( kD_B_ind )
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' Int. wavenumber  ',              '= ',wno( kI_ind )
		WRITE(*,"(A1,A20,A2,F8.2)")    '*',' Forc wavenumber  ',              '= ',wno( kF_ind )
		WRITE(*,"(A1,A20,A2,I8)")      '*',' Total Triad count ',             '= ',triad_count
		WRITE(*,"(A1,A20,A2,I8)")      '*',' Deleted triads ',                '= ',triad_deleted
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' Error in G.fac "b" ',            '= ',er_V_self
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' Error in G.fac "c" ',            '= ',er_VB
		WRITE(*,"(A1,A20,A2,ES8.2)")   '*',' Error in G.fac "h" ',            '= ',er_B_self
		WRITE(*,"(A1,A20,A2,F8.3)")    '*',' Eddy const',                     '= ',eddy_const
		WRITE(*,"(A1,A20,A2,F8.3)")    '*',' Dim const',                      '= ',dim_const

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		file_name = TRIM( ADJUSTL( file_address ) ) // 'wavenumbers.dat'
		OPEN(UNIT = 818, FILE = file_name)

		DO k_ind = 1, N
			WRITE(818,f_d16p8,ADVANCE = 'no')  wno_left( k_ind )
			WRITE(818,f_d16p8,ADVANCE = 'no')  wno( k_ind )
			WRITE(818,f_d16p8,ADVANCE = 'no')  wno_right( k_ind )
			WRITE(818,f_d16p8,ADVANCE = 'yes')  wno_band( k_ind )
		END DO

		CLOSE(818)
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // 'triads/' )
		DO k_ind = 1, N, 3
			CALL write_triad( k_ind )
		END DO
		! Writes all possible q,p for given k

	END
! </f>

	SUBROUTINE write_triad( k_ind0 )
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Write the list of q,p in the triad
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CHARACTER(LEN=10)::k_name
		INTEGER(KIND=4),INTENT(IN):: k_ind0

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  P  R  I  N   T          O  U  T  P  U  T
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		WRITE (k_name,f_i2) k_ind0
		! Writes 'k_ind0' as a CHARACTER

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		file_name = TRIM( ADJUSTL( file_address ) ) // 'triads/triad_' // TRIM( ADJUSTL ( k_name ) ) // '.dat'
		OPEN(UNIT = 828, FILE = file_name)

		DO q_ind = 1, N
		DO p_ind = p_ind_min( k_ind0, q_ind ), p_ind_max( k_ind0, q_ind )
				WRITE(828,f_i6,ADVANCE = 'no')  k_ind0
				WRITE(828,f_i6,ADVANCE = 'no')  q_ind
				WRITE(828,f_i6,ADVANCE = 'yes') p_ind
		END DO
		END DO

		CLOSE(828)
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END
! </f>

	SUBROUTINE write_spectrum( data_name, data_k )
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Write spectral data under the name of the file 'nam'
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CHARACTER(LEN=*),INTENT(IN)::data_name
		DOUBLE PRECISION,DIMENSION( N ),INTENT(IN):: data_k

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  P  R  I  N   T          O  U  T  P  U  T
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_sp )) &
		            // TRIM( ADJUSTL ( data_name ) ) // '_t_'// TRIM( ADJUSTL( file_time ) ) // '.dat'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OPEN(UNIT = 888, FILE = file_name)

		DO k_ind = 1, N
			WRITE(888,f_d12p6,ADVANCE = 'no')  wno( k_ind )
			WRITE(888,f_d32p17,ADVANCE = 'yes')data_k( k_ind )
		END DO

		CLOSE(888)

  END
! </f>

	SUBROUTINE write_kinetic_energy_spectrum()
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Write all spectral data in a single file (optional)
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  P  R  I  N   T          O  U  T  P  U  T
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_sp )) &
		            // 'spectral_data_V_t_'// TRIM( ADJUSTL( file_time ) ) // '.dat'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OPEN(UNIT = 882, FILE = file_name)

		DO k_ind = 1, N

			WRITE(882,f_d12p6,ADVANCE = 'no')  wno( k_ind )
			WRITE(882,f_d32p17,ADVANCE = 'no')en_spec_V( k_ind )
			WRITE(882,f_d32p17,ADVANCE = 'no')tr_spec_V( k_ind )
			IF ( coupling_status .NE. 0 ) THEN
				WRITE(882,f_d32p17,ADVANCE = 'no')tr_spec_V_self( k_ind )
				WRITE(882,f_d32p17,ADVANCE = 'no')tr_spec_V_intr( k_ind )
			END IF
			IF ( visc_status .EQ. 1 ) THEN
				WRITE(882,f_d32p17,ADVANCE = 'no')visc * laplacian_k( k_ind ) * en_spec_V( k_ind )
			END IF
			IF ( forc_status .EQ. 1 ) THEN
				WRITE(882,f_d32p17,ADVANCE = 'yes')fr_spec( k_ind )
			ELSE
				WRITE(882,f_d32p17,ADVANCE = 'yes')zero
			END IF

		END DO

		CLOSE(882)

  END
! </f>

	SUBROUTINE write_magnetic_energy_spectrum()
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! Write all spectral data in a single file (optional)
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  P  R  I  N   T          O  U  T  P  U  T
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_sp )) &
		            // 'spectral_data_B_t_'// TRIM( ADJUSTL( file_time ) ) // '.dat'

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OPEN(UNIT = 882, FILE = file_name)

		DO k_ind = 1, N

			WRITE(882,f_d12p6,ADVANCE = 'no')  wno( k_ind )
			WRITE(882,f_d32p17,ADVANCE = 'no')en_spec_B( k_ind )
			WRITE(882,f_d32p17,ADVANCE = 'no')tr_spec_B( k_ind )
			WRITE(882,f_d32p17,ADVANCE = 'no')tr_spec_B_self( k_ind )
			IF ( diff_status .EQ. 1 ) THEN
				WRITE(882,f_d32p17,ADVANCE = 'no')tr_spec_B_intr( k_ind )
				WRITE(882,f_d32p17,ADVANCE = 'yes')diff * laplacian_k( k_ind ) * en_spec_B( k_ind )
			ELSE
				WRITE(882,f_d32p17,ADVANCE = 'yes')tr_spec_B_intr( k_ind )
			END IF
			! WRITE(882,f_d32p17,ADVANCE = 'yes')dyn_rate_spec( k_ind )

		END DO

		CLOSE(882)

  END
! </f>

	SUBROUTINE write_kinetic_temporal_data
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! write the data in time, for every timestep
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!   E N E R G Y    V S    T I M E
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF ( t_step .EQ. 0 ) THEN

			file_name = TRIM( ADJUSTL( file_address ) ) // 'energy_V_vs_time.dat'
			!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			OPEN(unit = 4004, file = file_name )
			! File where energy vs time will be written. With additional data

		END IF

		WRITE(4004,f_d8p4,ADVANCE   ='no')  time_now
		WRITE(4004,f_d32p17,ADVANCE ='no')  energy_V
		WRITE(4004,f_d32p17,ADVANCE ='no')  ds_rate_self_V
		WRITE(4004,f_d32p17,ADVANCE ='no')  ds_rate_intr_V
		IF ( visc_status .EQ. 1 ) THEN
			WRITE(4004,f_d32p17,ADVANCE ='no') ds_rate_visc_V
		END IF
		IF ( forc_status .EQ. 1 ) THEN
			WRITE(4004,f_d32p17,ADVANCE ='no') ds_rate_forc
		END IF
		WRITE(4004,f_d32p17,ADVANCE ='no')  ds_rate_net_V
		WRITE(4004,f_d32p17,ADVANCE ='no')  ds_rate_net_V + ds_rate_net_B
		WRITE(4004,f_d32p17,ADVANCE ='no')  energy_tot
		WRITE(4004,f_d32p17,ADVANCE ='yes') skewness

		IF ( t_step .EQ. t_step_total ) THEN
			CLOSE(4004)
		END IF
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	END
! </f>

	SUBROUTINE write_magnetic_temporal_data
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! write the data in time, for every timestep
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!   E N E R G Y    V S    T I M E
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF ( t_step .EQ. 0 ) THEN

			file_name = TRIM( ADJUSTL( file_address ) ) // 'energy_B_vs_time.dat'
			!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			OPEN(unit = 4005, file = file_name )
			! File where energy vs time will be written. With additional data

		END IF

		WRITE(4005,f_d8p4,ADVANCE   ='no')  time_now
		WRITE(4005,f_d32p17,ADVANCE ='no')  energy_B
		WRITE(4005,f_d32p17,ADVANCE ='no')  ds_rate_self_B
		WRITE(4005,f_d32p17,ADVANCE ='no')  ds_rate_intr_B
		IF ( diff_status .EQ. 1 ) THEN
			WRITE(4005,f_d32p17,ADVANCE ='no') ds_rate_diff_B
		END IF
		WRITE(4005,f_d32p17,ADVANCE ='no')  ds_rate_net_B
		! WRITE(4005,f_d32p17,ADVANCE ='no')  dynamo_exp_calc
		WRITE(4005,f_d32p17,ADVANCE ='yes')  dynamo_exp

		IF ( t_step .EQ. t_step_total ) THEN
			CLOSE(4005)
		END IF
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	END
! </f>

	SUBROUTINE print_running_status
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! print the running status of the program, when called during debug
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		IF ( t_step .EQ. 0 ) THEN

			WRITE(*,'(A63)')'-----------------------------------------------------------'
			WRITE(*,'(A63)')'|   TIME     |     ENERGY-V  |    ENERGY-B   |  ENERGY TOT|'
			WRITE(*,'(A63)')'-----------------------------------------------------------'

		END IF

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		WRITE(*,'(A4,F8.4,A4,F12.8,A4,F12.8,A4,F12.8,A4)')'|   ',time_now,' |  '&
		,energy_V,' |  ',energy_B,' |  ',energy_tot,' | '
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		IF ( t_step .EQ. t_step_total ) THEN
			WRITE(*,'(A63)')'-----------------------------------------------------------'
	END IF

	END
! </f>

	SUBROUTINE print_error_timestep()
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! print error if time step chosen is large
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
		WRITE(*,'(A40)')       TRIM( ADJUSTL( 'ERROR: TIME STEP TOO LARGE') )
		WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
		WRITE(*,'(A40,F10.6)') TRIM( ADJUSTL( ' RESET THE TIME STEP (AT MAX) AS :') ),dt_max
		WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
		WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	END
! </f>

	SUBROUTINE print_error_nan()
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL THIS SUBROUTINE TO:
	! print error when NaN is encountered
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
		WRITE(*,'(A40,F8.4)')  TRIM( ADJUSTL( 'ERROR: NAN ENCOUNTERED BEFORE T = ') ), time_now
		WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
		WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
		WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	END
! </f>

	SUBROUTINE write_sim_start()
	! <f
		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		WRITE(*,'(A60)')TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
		WRITE(*,'(A60)')TRIM( ADJUSTL('   E  V  O  L  U  T  I  O  N        S  T  A  R  T  S ') )
		WRITE(*,'(A60)')TRIM( ADJUSTL('IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII') )
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	END
! </f>

	SUBROUTINE write_sim_end()
	! <f
		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		WRITE(*,'(A60)')TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
		WRITE(*,'(A60)')TRIM( ADJUSTL('  E  V  O  L  U  T  I  O  N        E  N  D  S ' ) )
		WRITE(*,'(A60)')TRIM( ADJUSTL('IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII') )
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	END
! </f>

END MODULE system_basicoutput
