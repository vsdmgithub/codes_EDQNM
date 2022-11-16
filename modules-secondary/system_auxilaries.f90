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
! MODULE: system_auxilaries
! LAST MODIFIED: 15 NOV 2022
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! AUXILARY FUNCTIONS FOR EDQNM CODE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_auxilaries
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This module has subroutines, that are minor functions used everywhere in the EDQNM code.
! These is a completely independent module
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_constants
	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT  NONE

  CONTAINS
! </f>

	SUBROUTINE find_CFL_timestep(dt_max1,dt_max2,delta_t)
! <f
	! CALL this to find a nice rounded of time step based on limits
		IMPLICIT  NONE
		DOUBLE PRECISION,INTENT(IN)::dt_max1,dt_max2
		DOUBLE PRECISION,INTENT(OUT)::delta_t

    delta_t = 10.0D0 ** DBLE( FLOOR( DLOG10( dt_max1 ) ) )

    DO WHILE( delta_t .GE. dt_max2 )
      delta_t = 10.0D0 ** DBLE( FLOOR( DLOG10( hf *delta_t ) ) )
    END DO

    IF( delta_t * 5.0D0 .LT. dt_max2 ) THEN
      delta_t = delta_t * 5.0D0
    END IF

	END
! </f>

	SUBROUTINE step_to_time_convert(step,time,delta_t)
! <f
	! CALL this to convert time step into actual time of simulation
		IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(IN)::step
		DOUBLE PRECISION,INTENT(IN)::delta_t
		DOUBLE PRECISION,INTENT(OUT)::time
		time=DBLE( step ) * delta_t
	END
! </f>

	SUBROUTINE time_to_step_convert(time,step,delta_t)
! <f
	! CALL this to convert time step into actual time of simulation
		IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(OUT)::step
		DOUBLE PRECISION,INTENT(IN)::delta_t
		DOUBLE PRECISION,INTENT(IN)::time
		step  = CEILING( time / delta_t )
	END
! </f>

	SUBROUTINE init_random_seed
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! randomize the seed, to generate random numbers. If not called, everytime it will
	! provide same set of random numbers
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! -----------------------------------------------------------------
    ! THIS PRODUCES NEW RANDOM NUMBERS EVERYTIME CALLED (uniform dist.)
    ! -----------------------------------------------------------------
    INTEGER, DIMENSION(8) :: seed_values
    INTEGER(KIND=4)::seed_size
    ! Declare an assumed shape, dynamic array
    INTEGER, DIMENSION(:),ALLOCATABLE :: seed
    ! gfortran SUBROUTINE to return date and time INFORMATion
    ! from the real time system clock. Works DOwn to milliseconds
    ! and stores the eight return values IN array values.
    CALL DATE_AND_TIME(VALUES=seed_values)
    ! restart the state of the pseuDOranDOm number generator
    ! k = mINimum size of seed (12 on my system)
    seed_size=20
    CALL RANDOM_SEED(size=seed_size)
    ! ALLOCATE memory to seed
    ALLOCATE(seed(seed_size))
    ! assign INFORMATion IN values to seed
    seed(:) = seed_values(:)
    ! seed the ranDOm number generator
    CALL RANDOM_SEED(put=seed)
    ! -----------------------------------------------------------------------

	END
! </f>

	SUBROUTINE get_simulation_name( sim_char )
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! This provides a string  of format
	!     'run_d121120_t104022' for run dated 12/11/2020 timed 10:40:22
	!
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! _________________________
		! LOCAL  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CHARACTER(LEN=8)::year_char,month_char,date_char
		CHARACTER(LEN=8)::hour_char,min_char,sec_char
		INTEGER,DIMENSION(8)::values
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CHARACTER(LEN=*),INTENT(OUT)::sim_char

		CALL DATE_AND_TIME(VALUES=values)
		! Gets the date and time as integer array

		values(1)   =   MOD(values(1),2000)
		! writing only last two digits of year

		WRITE(year_char,'(I2)')			 values(1)
		WRITE(month_char,'(I2.2)')	 values(2)
		WRITE(date_char,'(I2.2)')		 values(3)
		WRITE(hour_char,'(I2.2)')		 values(5)
		WRITE(min_char,'(I2.2)')	   values(6)
		WRITE(sec_char,'(I2.2)')		 values(7)
		! Self-explained

		sim_char    =  'run_D'//TRIM(ADJUSTL(date_char))//TRIM(ADJUSTL(month_char))//&
		TRIM(ADJUSTL(year_char))//'_T'//TRIM(ADJUSTL(hour_char))//&
		TRIM(ADJUSTL(min_char))//TRIM(ADJUSTL(sec_char))

	END
! </f>

  SUBROUTINE find_time_step( t_ref, t_fix )
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
  ! CALL this to fix a timestep from a reference time, with nice decimal places,
  ! e.g 0.00452 ->  0.002, 0.00123 -> 0.001 , 0.008922 -> 0.005
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    INTEGER(KIND=4)::st,dec
    DOUBLE PRECISION,INTENT(IN)::t_ref
    DOUBLE PRECISION,INTENT(OUT)::t_fix

    st  = 0
    dec = 0
    ! Decimal place, starting with zero

    IF ( t_ref .LT. 10.0D0 ) THEN

        DO WHILE( st .NE. 1)
        ! Loop runs till the decimal place is found

            IF ( t_ref > 10.0D0 ** ( - dec ) ) THEN

                t_fix = DBLE( FLOOR( t_ref * ( 10.0D0 ** (dec) ) ) )
                ! First non-zero digit in the 'dec' decimal place

                st    = 1
                ! Found.

            ELSE

                dec   = dec + 1
                ! Exploring the next decimal place

            END IF

        END DO

    END IF

    ! To have cutoff of type 'a * 10(-d)' with a=1,2,5.
    IF ( t_fix .GT. 5 ) THEN

        t_fix = 5.0D0 * ( 10.0D0 ** ( - dec ) )

    ELSE IF ( t_fix .GT. 2 ) THEN

        t_fix = 2.0D0 * ( 10.0D0 ** ( - dec ) )

    ELSE

        t_fix = 1.0D0 * ( 10.0D0 ** ( - dec ) )

    END IF

    END
! </f>

	DOUBLE PRECISION FUNCTION solid_angle( dim0 )
! <f
	! ------------
	! FUNCTION TO: Calculate the solid angle in 'd' dimension
	! -------------
		IMPLICIT NONE
		! _________________________
		! TRANSFER  VARIABLES
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DOUBLE PRECISION,INTENT(IN)::dim0
		DOUBLE PRECISION::dd
		dd = hf * ( dim0 + 1 )
		solid_angle = two * ( ( hf * two_pi ) ** dd ) / DGAMMA( dd )

		RETURN

	END
! </f>
END MODULE system_auxilaries
