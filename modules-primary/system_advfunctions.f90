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
! MODULE: system_advfunctions
! LAST MODIFIED: 21 June 2022
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER FOR  EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_advfunctions
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and updates it by a step, using the subroutines
! 1. rk4_algorithm
! 2. time_derivative
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_basicfunctions

  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT NONE
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::flux_pos,flux_neg
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::flux_pos_local,flux_neg_local
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::flux_pos_nonlocal,flux_neg_nonlocal
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  CONTAINS
! </f>

	SUBROUTINE allocate_flux_decomposition_arrays
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       allocate all arrays used in this module
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y     A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( flux_pos( N ), flux_neg( N ) )
    ALLOCATE( flux_pos_local( N ), flux_neg_local( N ) )
    ALLOCATE( flux_pos_nonlocal( N ), flux_neg_nonlocal( N ) )

  END
! </f>

  SUBROUTINE flux_decomposition
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to decompose the net flux at 'k' into positive and negative flux.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   F L U X     D E C O M P O S I T I O N     T   E   R   M
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    flux_pos          = zero
    flux_neg          = zero
    flux_pos_nonlocal = zero
    flux_neg_nonlocal = zero
    flux_pos_local    = zero
    flux_neg_local    = zero
    ! Reseting

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   P  O  S  I  T  I  V  E     F  L  U  X
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    DO k2_ind = 1 , N
    DO k_ind  = k2_ind+1 , N
    DO q_ind  = 1 , k2_ind
    DO p_ind  = 1 , k2_ind

	    IF ( kqp_status( k_ind, q_ind, p_ind )  .EQ. 1 ) THEN
	    ! If three wnoentum can form a triangle

        CALL transfer_term_integrand
        ! Getting the integrand term

        triad_sides =   (/ wno( k_ind ), wno( q_ind ), wno( p_ind ) /)
        ! Array consisting the 3 sides of the triangle

        min_max_ratio   =   MINVAL( triad_sides ) / MAXVAL( triad_sides )
        ! Finding the ratio of minimum to maximum sides of triangle

        IF ( min_max_ratio .LE. localness_cutoff_ratio ) THEN

          flux_pos_nonlocal( k2_ind ) = flux_pos_nonlocal( k2_ind ) &
          + integrand * wno_band( k_ind ) * wno_band( q_ind ) * wno_band( p_ind )
          ! Summation terms over all possible k,q,p for a given k'
					! in the region of nonlocal interactions.

        END IF

        flux_pos( k2_ind ) = flux_pos( k2_ind ) &
        + integrand * wno_band( k_ind ) * wno_band( q_ind ) * wno_band( p_ind )
        ! Summation terms over all possible k,q,p for a given k'.

	    END IF

    END DO
    END DO
    END DO
    END DO

    flux_pos_local = flux_pos - flux_pos_nonlocal
    ! Getting the local flux by subracting

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !  N  E  G  A  T  I  V  E     F  L  U  X
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    DO k2_ind = 1 , N
    DO k_ind  = 1 , k2_ind
    DO q_ind  = k2_ind+1 , N
    DO p_ind  = k2_ind+1 , N

	    IF ( kqp_status( k_ind, q_ind, p_ind )  .EQ. 1 ) THEN
	    ! If three wnoentum can form a triangle

        CALL transfer_term_integrand
        ! Getting the integrand term

        triad_sides =   (/ wno( k_ind ), wno( q_ind ), wno( p_ind ) /)
        ! Array consisting the 3 sides of the triangle

        min_max_ratio   =   MINVAL( triad_sides ) / MAXVAL( triad_sides )
        ! Finding the ratio of minimum to maximum sides of triangle

        IF ( min_max_ratio .LE. localness_cutoff_ratio ) THEN

          flux_neg_nonlocal( k2_ind ) = flux_neg_nonlocal( k2_ind ) &
          + integrand * wno_band( k_ind ) * wno_band( q_ind ) * wno_band( p_ind )
          ! Summation terms over all possible k,q,p for a given k'
					!  in the region of nonlocal interactions.

        END IF

        flux_neg( k2_ind ) = flux_neg( k2_ind ) &
        + integrand * wno_band( k_ind ) * wno_band( q_ind ) * wno_band( p_ind )
        ! Summation terms over all possible k,q,p for a given k'.

      END IF

    END DO
    END DO
    END DO
    END DO

		! flux_neg    =   flux_pos    -   flux
    ! Easy to get, rather than doing the above calculation, if nonlocal local decomposition is not needed.

    flux_neg_local = flux_neg - flux_neg_nonlocal
    ! Getting the local flux by subracting

    CALL write_spectrum('flux_pos'			    ,flux_pos)
    CALL write_spectrum('flux_neg'			    ,flux_neg)
    CALL write_spectrum('flux_pos_local'    ,flux_pos_local)
    CALL write_spectrum('flux_pos_nonlocal'	,flux_pos_nonlocal)
    CALL write_spectrum('flux_neg_local'		,flux_neg_local)
    CALL write_spectrum('flux_neg_nonlocal'	,flux_neg_nonlocal)
    ! REF-> <<< system_basicoutput >>>

  END

  SUBROUTINE deallocate_flux_decomposition_arrays
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
	! To deallocate all the flux decomposition arrays.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y        D  E  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(flux_pos, flux_neg)
    DEALLOCATE(flux_pos_local, flux_neg_local)
    DEALLOCATE(flux_pos_nonlocal, flux_neg_nonlocal)

  END

END MODULE system_advfunctions
