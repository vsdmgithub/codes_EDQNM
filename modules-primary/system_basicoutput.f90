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
! MODULE: system_basicoutput
! LAST MODIFIED: 21 JUNE 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! OUTPUT FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_basicoutput
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the outputs from the simulation are done in this module.
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
  ! CHARACTER(LEN =40) ::sub_dir_t

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

    sub_dir_sp  =   'k_data/'
    ! Sub directory name to store spectral data

    ! sub_dir_t   =   't_data/'
    ! Sub directory name to store temporal data

    type_sim    =  'N' // TRIM( ADJUSTL( N_char ) ) // '/'
    ! type_sim    =   'INVISCID_N' // TRIM( ADJUSTL( N_char ) ) // '/'
    ! type of simulation, the data is storing

    CALL get_simulation_name(name_sim)
    ! REF-> <<< system_auxilaries >>>
    ! Creating dated and timed name for the simulation for this particular type

    ! name_sim    =   'test_sim'
    ! Use this to give CUSTOM SIMULATION NAME

    file_address =   TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) //  &
                     TRIM( ADJUSTL( name_sim ) ) // '/'
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

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) )

    ! CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir_t ) ) )

    ! Command to create the main directory and sub directories (name_sim) in the desired path
    ! If exists already, it won't be an error

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
    WRITE(233,"(A3,A20,A2,I8)")      '1.',' No of modes    ',                '= ',N
    WRITE(233,"(A3,A20,A2,ES8.2)")   '2.',' Time step   ',                   '= ',dt
    WRITE(233,"(A3,A20,A2,I8)")      '3.',' Total time steps   ',            '= ',t_step_total
    WRITE(233,"(A3,A20,A2,F8.2)")    '4.',' Total time ',                    '= ',time_total
    WRITE(233,"(A3,A20,A2,F8.6)")    '5.',' Viscosity  ',                    '= ',viscosity
    WRITE(233,"(A3,A20,A2,I8)")      '6.',' No of saves   ',                 '= ',no_of_saves
    WRITE(233,"(A3,A20,A2,F8.3)")    '7.',' Initial energy ',                '= ',init_energy
    WRITE(233,"(A3,A20,A2,F8.2)")    '8.',' Smallest wavenumber',            '= ',min_wno
    WRITE(233,"(A3,A20,A2,F8.2)")    '9.',' Largest wavenumber ',            '= ',max_wno
    WRITE(233,"(A3,A20,A2,I8)")      '10.','Total Triad count ',             '= ',no_of_triads
    WRITE(233,"(A3,A20,A2,F8.2)")    '11.','Localness of triad, cutoff ',    '= ',localness_cutoff_ratio
    WRITE(233,"(A3,A20,A2,I8)")      '12.','Forcing status   ',              '= ',forcing_status

    CLOSE(233)
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    WRITE(*,"(A40)") TRIM(ADJUSTL('--------------------------------------------------------------------'))
    WRITE(*,"(A40)") TRIM(ADJUSTL('------  EDQNM  EQUATION----------------------'))
    WRITE(*,"(A40)") TRIM(ADJUSTL(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'))
    WRITE(*,"(A40)") TRIM(ADJUSTL('-----------PARAMETERS OF SIMULATION------------'))
    WRITE(*,"(A40)") TRIM(ADJUSTL('--------------------------------------------------------------------'))
    WRITE(*,"(A3,A20,A2,I8)")      '1.',' No of modes    ',                '= ',N
    WRITE(*,"(A3,A20,A2,ES8.2)")   '2.',' Time step   ',                   '= ',dt
    WRITE(*,"(A3,A20,A2,I8)")      '3.',' Total time steps   ',            '= ',t_step_total
    WRITE(*,"(A3,A20,A2,F8.2)")    '4.',' Total time ',                    '= ',time_total
    WRITE(*,"(A3,A20,A2,F8.6)")    '5.',' Viscosity  ',                    '= ',viscosity
    WRITE(*,"(A3,A20,A2,I8)")      '6.',' No of saves   ',                 '= ',no_of_saves
    WRITE(*,"(A3,A20,A2,F8.3)")    '7.',' Initial energy ',                '= ',init_energy
    WRITE(*,"(A3,A20,A2,F8.2)")    '8.',' Smallest wavenumber',            '= ',min_wno
    WRITE(*,"(A3,A20,A2,F8.2)")    '9.',' Largest wavenumber ',            '= ',max_wno
    WRITE(*,"(A3,A20,A2,I8)")      '10.','Total Triad count ',             '= ',no_of_triads
    WRITE(*,"(A3,A20,A2,F8.2)")    '11.','Localness of triad, cutoff ',    '= ',localness_cutoff_ratio
    WRITE(*,"(A3,A20,A2,I8)")      '12.','Forcing status   ',              '= ',forcing_status

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

  SUBROUTINE write_temporal_data
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y / E N S T R O P H Y   V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    file_name = TRIM( ADJUSTL( file_address ) ) // 'temporal_data.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    OPEN(unit = 4004, file = file_name )
    ! File where energy vs time will be written. With additional data

    DO t_step = 0, t_step_total

      WRITE(4004,f_d8p4,ADVANCE   ='no')  time      ( t_step )
      WRITE(4004,f_d32p17,ADVANCE ='no')  en_time   ( t_step )
      WRITE(4004,f_d32p17,ADVANCE ='no')  es_time   ( t_step )

      IF ( viscosity .GT. tol_float ) THEN
        WRITE(4004,f_d32p17,ADVANCE ='no')ds_time   ( t_step )
      END IF

      WRITE(4004,f_d32p17,ADVANCE ='yes') sk_time   ( t_step )

    END DO

    CLOSE(4004)
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
      WRITE(*,'(A63)')'|   TIME     |     ENERGY    |   ENSTROPHY   |  DISS  RATE  |'
      WRITE(*,'(A63)')'-----------------------------------------------------------'

    END IF

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(*,'(A4,F8.4,A4,F12.8,A4,F12.8,A4,F12.8,A4)')'|   ',time_now,' |  '&
      ,energy,' |  ',enstrophy,' |  ',dissipation_rate,' | '
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


END MODULE system_basicoutput
