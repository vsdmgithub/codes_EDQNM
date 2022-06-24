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
! MODULE: system_basicfunctions
! LAST MODIFIED: 21 June 2022
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SYSTEM PARAMETERS FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>
MODULE system_basicfunctions
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the system parameters, corresponding to the different variant of EDQNM equation we are studying
! is defined here. Note global variables is common amongst these different variants.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_initialcondition
  USE system_basicoutput
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT  NONE

  CONTAINS
! </f>

	SUBROUTINE allocate_edqnm_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       allocate all arrays related to edqnm algorithm
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y     A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( p_ind_min( N, N ), p_ind_max( N, N ) )
    ALLOCATE( kqp_status( N, N, N ), geom_fac( N, N, N ) )
    ALLOCATE( eddy_array( N ) )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( time    (0 : t_step_total) )
    ALLOCATE( en_time (0 : t_step_total) )
    ALLOCATE( es_time (0 : t_step_total) )
    ALLOCATE( ds_time (0 : t_step_total) )
    ALLOCATE( sk_time (0 : t_step_total) )

    DO t_step        = 0 , t_step_total

      time( t_step ) = dt * t_step
      ! Time axis declaration

    END DO

    IF ( ( viscosity .GT. tol_float ) .AND. ( forcing_status .EQ. 1 ) ) THEN

      ALLOCATE( forcer( N ), forcer_template( N ) )
      CALL init_forcing_template
      ! creates forcing arrays

    END IF
    
  END
! </f>

  SUBROUTINE init_edqnm_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to initialize all in-built arrays that need not change during the evolution of system.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION:: wno_p_min,wno_p_max
    DOUBLE PRECISION:: z_f, x_f, y_f

    kqp_status          =  0
    geom_fac            =  zero
    ! Reseting to zero for safety

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !  L I M I T S   O F    I N T E G R A T I O N
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    DO k_ind = 1 , N
    DO q_ind = 1 , N

      wno_p_min   =  DABS( wno( k_ind ) - wno( q_ind ) )
      wno_p_max   =  DABS( wno( k_ind ) + wno( q_ind ) )

      IF ( wno_p_min .LT. wno( 1 ) ) THEN
        p_ind_min( k_ind, q_ind )   =   1
      ELSE IF (find_index_floor( wno_p_min ) .LE. 1) THEN
        p_ind_min( k_ind, q_ind )   =    1
      ELSE
        p_ind_min( k_ind, q_ind )   =   find_index_floor( wno_p_min )
      END IF
      ! FINDING THE MINIMUM OF 'p_ind' FOR EVERY PAIR OF 'k_ind,q_ind'

      IF ( wno_p_max .GT. wno( N-1 ) ) THEN
        p_ind_max( k_ind, q_ind )   =   N
      ELSE IF (find_index_ceiling( wno_p_max ) .GE. N ) THEN
        p_ind_max( k_ind, q_ind )   =   N
      ELSE
        p_ind_max( k_ind, q_ind )   =   find_index_ceiling( wno_p_max )+1
      END IF
      ! FINDING THE MAXIMUM OF 'p_ind' FOR EVERY PAIR OF 'k_ind,q_ind'

    END DO
    END DO

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !  T R I A N G L E     C H E C K   &   C O S I N E S   O F    I T.
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    DO k_ind = 1 , N
    DO q_ind = 1 , N
    DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

      IF ( find_triangle_compatibility(k_ind, q_ind, p_ind) .EQ. 1) THEN

        kqp_status( k_ind, q_ind, p_ind )  =   1
        ! MEANING, THE GIVEN THREE wnoENTUM 'p,k,q' CAN FORM A TRIANGLE.

        CALL find_cosine( k_ind, q_ind, p_ind, z_f )
        CALL find_cosine( k_ind, p_ind, q_ind, y_f )
        CALL find_cosine( q_ind, p_ind, k_ind, x_f )
        ! FINDING COSINES FOR ALL THREE SIDES ONCE IT IS APPROVED TO PARTICIPATE IN THE TRIAD INTERACTION

        geom_fac( k_ind, q_ind, p_ind )  =  (z_f ** thr ) - x_f * y_f
        ! GEOMETRIC FACTOR IN THE E.D.Q.N.M

        no_of_triads = no_of_triads + 1
        ! Counting the triad

      END IF

    END DO
    END DO
    END DO

    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   T  R  I  A  D      D  E  B  U  G
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    DO k_ind = 1, N
    DO q_ind = 1, N
    DO p_ind = 1, N

      IF (p_ind_max(k_ind,q_ind) .GT. N) THEN
          PRINT*,'ERROR IN UPPER BOUNDARY',p_ind_max(k_ind,q_ind)
      END IF

      IF (p_ind_min(k_ind,q_ind) .LT. 1) THEN
          PRINT*,'ERROR IN LOWER BOUNDARY',p_ind_min(k_ind,q_ind),k_ind,q_ind
      END IF

      IF ( kqp_status( k_ind, q_ind, p_ind ) .NE. kqp_status( k_ind, p_ind, q_ind) ) THEN
          PRINT*,'ERROR IN TRIAD at k=',k_ind,' q=',q_ind,' p=',p_ind
      END IF

      IF ( kqp_status( k_ind, q_ind, p_ind ) .NE. kqp_status( q_ind, k_ind, p_ind) ) THEN
          PRINT*,'ERROR IN TRIAD at k=',k_ind,' q=',q_ind,' p=',p_ind
      END IF

    END DO
    END DO
    END DO

  END
! </f>

  SUBROUTINE init_forcing_template
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Initialize a template of forcing similar to the initial condition.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::A0,s_exp

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !  F  O  R  C  I  N  G       T  E  M  P  L  A  T  E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    A0      =   zero
    ! Normalization constant

    s_exp   =   two
    ! Integral scale spectrum exponent

    DO k_ind = 1, N

      A0 = A0 + ( wno( k_ind ) ** s_exp ) * DEXP( - hf * ( wno( k_ind ) / &
                  wno( ind_integral ) ) ** two ) * wno_band( k_ind )

    END DO

    DO k_ind = 1, N

      forcer_template( k_ind )   =   ( wno( k_ind ) ** s_exp )
      forcer_template( k_ind )   =   forcer_template( k_ind ) * DEXP(&
                                     - hf * ( wno( k_ind ) / wno( ind_integral ) ) ** two )

    END DO

    forcer_template     =    forcer_template    /   A0

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

  SUBROUTINE compute_spectral_data
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! TO compute all spectral data and write them
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! S P E C T R A L    D A T A
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    CALL write_spectrum('spectrum',spec) !  ENERGY FILE
    ! REF-> <<< system_basicoutput >>>

    CALL transfer_term

    CALL write_spectrum('transfer',transfer_spec) !  ENERGY TRANSFER FILE
    ! REF-> <<< system_basicoutput >>>

	  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	  !   N  E  T       F  L  U  X
	  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    flux       =   zero
	  DO k_ind = 1, N
      flux( k_ind ) = - SUM( transfer_spec( : k_ind) * wno_band( : k_ind) )
	  END DO

    CALL write_spectrum('flux',flux) !  FLUX  FILE
    ! REF-> <<< system_basicoutput >>>

    IF ( viscosity .GT. tol_float ) THEN

      CALL write_spectrum('dissipation',two * viscosity * laplacian_k * spec) !  DISSIPATION FILE
      ! REF-> <<< system_basicoutput >>>

    END IF

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

  SUBROUTINE compute_temporal_data
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! TO compute all temporal data and write them
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! T E M P O R A L   D A T A
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !  ENERGY,ENSTROPHY, DISSIPATION AND SKEWNESS VS TIME

    energy             = SUM( spec * wno_band )
    en_time(t_step)    = energy

    enstrophy          = SUM( laplacian_k * spec * wno_band )
    es_time(t_step)    = enstrophy

    IF ( viscosity .GT. tol_float ) THEN
      dissipation_rate = two * viscosity * enstrophy
      ds_time(t_step)  = dissipation_rate
    END IF

    skewness           = SUM( ( transfer_spec + forcer ) * laplacian_k * wno_band )
    skewness           = skewness * ( enstrophy ** ( -1.5D0 )) * DSQRT(135.0D0/98.0D0)
    sk_time(t_step)    = skewness

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END
! </f>

  SUBROUTINE transfer_term
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the transfer term for every k
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   T   R   A   N   S   F   E   R       T   E   R   M
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    transfer_spec   =   zero
    ! Reseting the transfer term

    DO k_ind = 1 , N
    DO q_ind = 1 , N
    DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )

    ! The third wnoentum runs through p_min to p_max, determined from the 'k,q'

      IF ( kqp_status( k_ind, q_ind, p_ind )  .EQ. 1 ) THEN
      ! If three wnoentum can form a triangle

        CALL transfer_term_integrand
        ! Getting the integrand term

        transfer_spec( k_ind ) = transfer_spec( k_ind ) + integrand * wno_band( q_ind ) * wno_band( p_ind )
        ! Summation terms over all possible q,p for a given k.

      END IF

    END DO
    END DO
    END DO

  END
! </f>

  SUBROUTINE transfer_term_integrand
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get triad_density
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   I  N  T  E  G  R  A  N  D
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    CALL damping_factor

    integrand  =  (wno( k_ind )**(two)) * spec( q_ind ) - ( wno( q_ind )**(two) ) * spec( k_ind )
    integrand  =  integrand * spec( p_ind ) * geom_fac( k_ind, q_ind, p_ind ) * damping / wno( p_ind )
    ! STITCHING THE INTEGRATION TERMS ONE BY ONE.

  END
! </f>

	SUBROUTINE damping_factor
! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the damping factor for the third wnoent
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   E  D  D  Y            F  R  E  Q  U  E  N  C  Y
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    eddy_array   = laplacian_k * spec * wno_band

    eddy_k       = SUM( eddy_array( : k_ind) ) ** hf
    eddy_q       = SUM( eddy_array( : q_ind) ) ** hf
    eddy_p       = SUM( eddy_array( : p_ind) ) ** hf

    eddy_freq    = eddy_constant * ( eddy_k + eddy_q + eddy_p )

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! D A M P I N G       F A C T O R
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    viscous_freq = viscosity * ( laplacian_k( k_ind ) + &
                    laplacian_k( q_ind ) + laplacian_k( p_ind ) )

    damping      = one - DEXP( -( eddy_freq + viscous_freq ) * time_now )

    damping      = damping / ( eddy_freq + viscous_freq )

  END
! </f>

  SUBROUTINE perform_debug
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check Nan in data
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    DO k_ind = 1, N
      IF ( spec ( k_ind) .NE. spec ( k_ind) ) THEN

        NaN_count =  NaN_count + 1

        PRINT*,"NaN encountered before t=",time_now
        ! IF any NaN is encountered, the loop is exited without any further continuation.

      END IF
    END DO

  END
! </f>

	SUBROUTINE deallocate_edqnm_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       deallocate all edqnm arrays
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y     D  E   A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( p_ind_min, p_ind_max )
    DEALLOCATE( kqp_status, geom_fac )
    DEALLOCATE( forcer, forcer_template )
    DEALLOCATE( eddy_array )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( time, en_time, es_time, ds_time, sk_time )

  END
! </f>

  SUBROUTINE find_cosine( i1, i2, i3 , cosine)
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this calculate the cosine of angle opposite to side with index 'i3'
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! TRANSFER  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(OUT)::cosine
    INTEGER(KIND=4),INTENT(IN)::i1, i2, i3

    cosine  = ( wno( i1 ) ** two ) + ( wno( i2 ) ** two  ) - ( wno( i3 ) ** two )

    cosine  = cosine / ( two * wno( i1 ) * wno( i2 ) )

    ! FINDING COSINE OF ANGLE BETWEEN SIDES 'i1,i2'

  END
! </f>

  INTEGER FUNCTION find_index_floor( wno0 )
! <f
  ! ------------
  ! FUNCTION TO: Calculate the largest index whose wnoentum is smaller than the given wnoentum 'wno0'
  ! -------------
    DOUBLE PRECISION::wno0

    find_index_floor  =   FLOOR ( DLOG( wno0 / wno_base) / log_lambda )

    RETURN

  END
! </f>

  INTEGER FUNCTION find_index_ceiling( wno0 )
! <f
  ! ------------
  ! FUNCTION TO: Calculate the smallest index whose wnoentum is larger than the given wnoentum 'wno0'
  ! -------------
    DOUBLE PRECISION::wno0

    find_index_ceiling  =   CEILING ( DLOG( wno0 / wno_base) / log_lambda )

    RETURN

  END
! </f>

  INTEGER FUNCTION find_triangle_compatibility( i1, i2, i3 )
! <f
  ! ------------
  ! FUNCTION TO: Calculate the compatibility that three wnoentum can form a triangle, 1 is yes, 0 means no.
  ! -------------
    IMPLICIT NONE
    ! _________________________
    ! TRANSFER  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4),INTENT(IN)::i1, i2, i3

    find_triangle_compatibility =  0

    IF ( wno(i1) + wno(i2) .GT. wno(i3) ) THEN
    IF ( wno(i3) + wno(i2) .GT. wno(i1) ) THEN
    IF ( wno(i1) + wno(i3) .GT. wno(i2) ) THEN

      find_triangle_compatibility =  1

    END IF
    END IF
    END IF

  RETURN

  END
! </f>

END MODULE system_basicfunctions
