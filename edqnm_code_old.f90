module parameters_mod
    ! -----------------------------------------------------------------
    ! ALL THE PARAMETERS ARE DECLARED HERE
    ! -----------------------------------------------------------------
    implicit none
    integer (kind=4)::N,s_max,s,s_print,i_test_error,no_of_files
    double precision::lb,k0,a,vcs,t_max,t_start,t_end,dt,sc_vcs_time
    double precision,dimension(:),allocatable::k,dk,ku,kl
    character(len=40)::path
    contains
    subroutine initialization_parameters
    ! --------------------------------------------------
    ! DEFINING ALL GLOBAL VARIABLES
    ! ---------------------------------------------------
    implicit none
    integer (kind=4)::i
        ! No of discrete wavenumbers, which will form shells in 3d space with radial width 'dk_i'
        N=61
        allocate(k(N),dk(N),ku(N),kl(N))
        ! Defining the discrete wavenumbers {k_i} which are equispaced in log scale as
        ! EQN:- $k_i=k_0 \lambda ^{i-1}$
        lb=2.0**(1.0/4.0)
        k0=2.0**(-3.0)
        do i=1,N
            k(i)=k0*(lb**(i-1))
        end do
        ! Defining the band for each discrete k_i as k^{u}_i,k^{l}_i  for above and below
        ! EQN:- $\Delta k_i=k_i \ln(\lambda)$
        ! EQN:- $k_i^{l}=\frac{\Delta k_i}{\lambda -1}
        ! EQN:- $k_i^{u}=\frac{\lambda \Delta k_i}{\lambda -1}
        do i=1,N
            dk(i)=k(i)*Log(lb)
            ku(i)=dk(i)*lb/(lb-1)
            kl(i)=dk(i)/(lb-1)
        end do
        ! Defining Viscosity
        vcs=4.4*(10**(-4.0))
        ! Definng constant in eddy viscosity
        a=0.54
        ! Defning time's in integration
        t_start=0.0
        t_end=10.0
        t_max=t_end-t_start
        dt=1.0/(10**(5.0))
        ! Time steps and saving details
        s_max=INT(t_max/dt)
        no_of_files=1000
        s_print=INT(s_max/no_of_files)
        ! Scales in the system
        sc_vcs_time=1.0/(vcs*(k(N)**(2.0))) ! Time scale of smallest eddy
        ! Specifying the folder to save all outputs (this folder has to be created where the .f95 file is)
        path='data/'
        ! Specifying where to look for NaN first
        i_test_error=INT(N/6)
    end subroutine initialization_parameters
    subroutine print_wavenumber
    ! Prints the wavenumbers   along with upper, lower and width of band
    implicit none
        character(len=40)::l1
        integer (kind=4)::i
        l1='wavenumbers'
        open(unit=10,file=trim(adjustl(path))//trim(adjustl(l1))//'.dat',status='unknown')
        do i=1,N
        write(10,101,advance='no')kl(i)
        write(10,101,advance='no')k(i)   
        write(10,101,advance='no')ku(i)
        write(10,101,advance='yes')dk(i)
        end do
        close(10)
    ! --------------
    ! FORMATS
    ! -------------- 
        101 format (f24.10)
    end subroutine print_wavenumber
end module parameters_mod
module function_mod
    use parameters_mod
    implicit none
    ! -------------------------------------------------------------------------------------
    ! ALL THE FUNCTIONS FOR E.D.Q.N.M MODEL ARE DEFINED HERE
    ! -------------------------------------------------------------------------------------
    contains
    subroutine find_maxindex(k_ref,i_ref)
    ! Given a wavenumber, finds the minimum wavenumber lattice which is larger than that
    implicit none
    double precision,intent(in)::k_ref
    integer(kind=4),intent(out)::i_ref
    i_ref=INT(Log(k_ref/k0)/Log(lb))+1
    end subroutine find_maxindex
    subroutine find_minindex(k_ref,i_ref)
    ! Given a wavenumber, finds the maximum wavenumber lattice which is larger than that
    implicit none
    double precision,intent(in)::k_ref
    integer(kind=4),intent(out)::i_ref
    i_ref=INT(Log(k_ref/k0)/Log(lb))
    end subroutine find_minindex
    subroutine integration_limits(i1,i2,i3_min,i3_max)
    ! given k1,k2 the limits of k3 so that three forms a triangle
    implicit none
    integer (kind=4),intent(in)::i1,i2
    double precision::k3_min,k3_max
    integer (kind=4),intent(out)::i3_min,i3_max
    k3_min=ABS(k(i1)-k(i2))
    k3_max=ABS(k(i1)+k(i2))
    if (k3_min<k(1)) then
    i3_min=1
    else
    call find_maxindex(k3_min,i3_min)
    end if
    if (k3_max>k(N)) then
    i3_max=N
    else
    call find_minindex(k3_max,i3_max)
    end if    
    end subroutine integration_limits
    subroutine find_cosine(i1,i2,i3,cosine)
    ! returns the cosine of angle formed by sides k(i1) and k(i2)
    implicit none
    integer(kind=4),intent(in)::i1,i2,i3
    double precision,intent(out)::cosine
    cosine=(k(i1)**(2.0)+k(i2)**(2.0)-k(i3)**(2.0))/(2.0*k(i1)*k(i2))
    end subroutine find_cosine
    subroutine geometric_factor(i1,i2,i3,b)
    implicit none
    ! To calculate the geometric factor 'b' with a triad 'k_1,k_2,k_3'
    integer (kind=4),intent(in)::i1,i2,i3
    double precision,intent(out)::b
    double precision::x,y,z
    call find_cosine(i1,i2,i3,z)
    call find_cosine(i1,i3,i2,y)
    call find_cosine(i2,i3,i1,x)
    b=(k(i2)/k(i1))*(z**(3.0)-x*y) ! conventionally x has to be -x here, but our notation change has another negative sign
    end subroutine geometric_factor
    subroutine damping_factor(i1,i2,i3,en,t,th)
    ! To calculate the eddy damping factor $\theta _{k_1k_2k_3}$
    implicit none
    integer (kind=4),intent(in)::i1,i2,i3
    double precision,dimension(N),intent(in)::en
    double precision,intent(in)::t
    double precision,intent(out)::th
    double precision::m1,m2,m3,m,w
    call eddy_viscosity(i1,en,m1)    
    call eddy_viscosity(i2,en,m2)
    call eddy_viscosity(i3,en,m3)
    m=m1+m2+m3
    w=vcs*(k(i1)**(2.0)+k(i2)**(2.0)+k(i3)**(2.0))
    th=(1.0-exp(-(m+w)*t))/(m+w)
    end subroutine damping_factor
    subroutine eddy_viscosity(i_max,en,m)
    ! To calculate the eddy viscosity for a given wavenumber, with full energy spectrum
    ! EQN:- $\mu_{k_i}=a(\sum_{j=1}^{i}k_i^2 E(k_i)dk_i)^{1/2}
    implicit none
    integer (kind=4),intent(in)::i_max
    double precision,dimension(N),intent(in)::en
    double precision,intent(out)::m
    integer (kind=4)::i
    m=0.0
    do i=1,i_max
        m=m+(k(i)**(2.0))*en(i)*dk(i)
    end do
    m=a*(m**(0.5))
    end subroutine eddy_viscosity
    subroutine energy_factor(i1,i2,i3,en,ef)
    ! Gives the transfer term(partially), which includes only the energies.
    implicit none
    integer (kind=4),intent(in)::i1,i2,i3
    double precision,dimension(N),intent(in)::en
    double precision,intent(out)::ef
    ef=((k(i1)**(2.0))*en(i2)-(k(i2)**(2.0))*en(i1))*en(i3)*k(i1)/(k(i2)*k(i3))
    end subroutine energy_factor
    subroutine triangle_check(i1,i2,i3,ans)
    ! Checks if three momentum forms a triangle
    implicit none
    integer (kind=4),intent(in)::i1,i2,i3
    integer (kind=4),intent(out)::ans
    if (((k(i2)+k(i3))>k(i1)).AND.(ABS(k(i2)-k(i3))<k(i1))) then
        ans=1
    end if
    end subroutine triangle_check
    subroutine transfer_integrand(i1,i2,i3,en,t,tr_int)
    ! Term that has to be summed over two momentum 'k_2,k_3' which forms triangle for given 'k_1'
    implicit none
    integer (kind=4),intent(in)::i1,i2,i3
    double precision,dimension(N),intent(in)::en
    double precision,intent(out)::tr_int
    double precision,intent(in)::t
    double precision::th,b,ef
    call geometric_factor(i1,i2,i3,b)
    call damping_factor(i1,i2,i3,en,t,th)
    call energy_factor(i1,i2,i3,en,ef)
    tr_int=th*b*ef
    end subroutine transfer_integrand
     subroutine transfer_integrated(i1,t,en,tr_int)
!     transfer term 'T' integrated over 'k_2,k_3'
    implicit none
    integer (kind=4),intent(in)::i1
    double precision,intent(in)::t
    double precision,dimension(N),intent(in)::en
    double precision::tr
    double precision,intent(out)::tr_int
    integer (kind=4)::i2,i3,i3_min,i3_max,ct_p,ans
    ans=0
    tr_int=0.0
    ct_p=0
        do i2=1,N
            if (i2/=i1) then
                call integration_limits(i1,i2,i3_min,i3_max)
                do i3=i3_min,i3_max
                    call triangle_check(i1,i2,i3,ans)
                    if (ans==1) then
                        call transfer_integrand(i1,i2,i3,en,t,tr)
                        ct_p=ct_p+1
                        tr_int=tr_int+tr*dk(i2)*dk(i3) 
                   end if
                    ans=0
                end do
            end if
        end do
    end subroutine transfer_integrated
  end module function_mod
module solver_mod
    ! ----------------------------------------------------------------------------------------
    ! SOLVES EDQNM WITH AN INITIAL CONDITION USING RK4
    ! ---------------------------------------------------------------------------------------
    use function_mod
    implicit none
    double precision,dimension(:),allocatable::en_init
    double precision,dimension(:,:),allocatable::en_soln
    ! --------------------------------------------------------------------------------------------
    ! TAKES INITIAL CONDITION FOR 'E(k,0)' AND SOLVES FOR 'E(k,t)' 
    ! --------------------------------------------------------------------------------------------
    contains
    subroutine initial_condition
    implicit none
    integer (kind=4)::i
    double precision::en0,en_thres,en_exp,en_total
    allocate(en_init(N))
    en_total=0.0

!   INITIAL CONDITION CREATED
!---------------------------------------------

!    en_exp=1.0
!    en_thres=2.0**(-60.0)
!    do i=1,N
!           en_init(i)=(k(i)**(en_exp))*exp(-(en_exp*(k(i)**(2.0)))/(2.0*(k(1)**(2.0))))+en_thres
!           en_total=en_total+en_init(i)
!    end do
!    ! gives a spectrum of intial energies. Then we normalize the total energy to 0.5.
!       en0=1.0/(2.0*en_total)
!    ! en0 is the normalization factor we have to multiply to original spectrum.
!    do i=1,N
!           en_init(i)=en0*(k(i)**(en_exp))*exp(-(en_exp*(k(i)**(2.0)))/(2.0*(k(1)**(2.0))))+en_thres
!    end do

    en_thres=2**(-60.0)
    do i=1,N
       en_init(i)=(k(i)**(2.0))*exp(-(k(i)**(2.0))/(1.0))+en_thres
      en_total=en_total+en_init(i)
    end do
    en0=1.0/(2.0*en_total)
    do i=1,N
       en_init(i)=en0*(k(i)**(2.0))*exp(-(k(i)**(2.0))/(1.0))+en_thres
    end do

!   INITIAL CONDITION FROM FILE
!-----------------------------------------------
!    open(unit=10,file='energy99.dat')
!        do i=1,N
!            read(10,'(f32.18)',advance='no')k(i)
!            read(10,'(f32.18)',advance='yes')en_init(i)
!            print*,en_init(i)
!        end do
!    close(10)

   end subroutine initial_condition
    subroutine print_initial_condition
    ! Prints the initial condition as "k(i) en_init(i) entrosphy(i)"
    implicit none
    character(len=40)::l1
    integer (kind=4)::i
    l1='energy_init'
    open(unit=11,file=trim(adjustl(path))//trim(adjustl(l1))//'.dat',status='unknown')
        do i=1,N
            write(11,102,advance='no')k(i) 
            write(11,103,advance='yes')en_init(i)
        end do
    close(11)
    ! --------------
    ! FORMATS
    ! -------------- 
   
        102 format (f24.14)
        103 format (f32.15)
    end subroutine print_initial_condition
     subroutine print_parameters
     implicit none
     double precision::m_N,m_1,eddy_timeN,eddy_time1,en_tot_init
     call eddy_viscosity(N,en_init,m_N)                                                      
     call eddy_viscosity(1,en_init,m_1)
     en_tot_init=SUM(en_init)
     eddy_timeN=1.0/m_N
     eddy_time1=1.0/m_1
    open(unit=23,file=trim(adjustl(path))//'parameters.dat',status='unknown')
        write(23,*)'no of shells= ',N
        write(23,*)'viscosity=',vcs
        write(23,*)'integration time= ',t_max
        write(23,*)'time step= ',dt
        write(23,*)'eddy_time_scale_N= ',eddy_timeN
        write(23,*)'viscous_time_scale_N= ',sc_vcs_time
        write(23,*)'largest_eddy_turn_overtime= ',eddy_time1
         write(23,*)'largest_viscous_time_scale=',(1.0/(vcs*k(1)*k(1)))
        write(23,*)'total initial energy = ',en_tot_init     
    close(23)
    end subroutine
   subroutine rk4_solver
    ! RK4 algorithm
    implicit none
    double precision,dimension(N)::en_temp,den1,den2,den3,den4
    allocate(en_soln(N,s_max))
    en_soln(:,1)=en_init(:)
    do s=1,s_max-1
        en_temp=en_soln(:,s)
        call edqnm_eqn(en_temp,s,den1)
        en_temp=en_soln(:,s)+0.5*den1
        call edqnm_eqn(en_temp,s,den2)
        en_temp=en_soln(:,s)+0.5*den2
        call edqnm_eqn(en_temp,s,den3)
        en_temp=en_soln(:,s)+den3
        call edqnm_eqn(en_temp,s,den4)
        en_soln(:,s+1)=en_soln(:,s)+(den1+2.0*(den2+den3)+den4)/6.0
        ! CHECK for "NaN"
         if(ABS(en_soln(i_test_error,s+1))/=ABS(en_soln(i_test_error,s+1))) then
            print*,"Simulation interupted at time",s
            exit
         end if
         ! Print solution for every 's_print' time steps
        if (MOD(s,s_print)==0) then
            call print_solution(s)    
            print*,INT(s/s_print)
        end if
    end do    
    end subroutine rk4_solver
   subroutine edqnm_eqn(en,s0,den)
    ! To find $dE=dt \times \dot(E(k))$
    implicit none
    double precision,dimension(N),intent(in)::en
    integer (kind=4),intent(in)::s0
    double precision,dimension(N),intent(out)::den
    double precision::t,tr_int2
    integer (kind=4)::i
    t=s0*dt
    do i=1,N                                                                                                                                                            
        call transfer_integrated(i,t,en,tr_int2)
!        tr_int2=0.0
        den(i)=dt*(tr_int2-2.0*vcs*(k(i)**(2.0))*en(i))
    end do
    end subroutine edqnm_eqn
    subroutine print_solution(s)
    ! Print the solution for s=124 time step as file name 'data/energy124.dat' with format 'k(i) en(124,i) k^2 en(124,i)'
    implicit none
    integer(kind=4),intent(in)::s
    character(len=40)::l1,c
    integer (kind=4)::i
        l1='energy'
        write (c,21),INT(s/s_print)
        open(unit=13,file=trim(adjustl(path))//trim(adjustl(l1))//trim(adjustl(c))//'.dat',status='unknown')
            do i=1,N
                write(13,104,advance='no')k(i) 
                write(13,105,advance='yes')en_soln(i,s)
            end do
        close(13)
    ! --------------
    ! FORMATS
    ! --------------
        21 format(I5)
        104 format (f24.14)
        105 format (f32.15)
    end subroutine print_solution
end module solver_mod
program edqnm_code
use solver_mod
implicit none
character(len=4)::c1
call initialization_parameters
call print_wavenumber
call initial_condition
call print_initial_condition
call print_parameters
print*,"Enter 'Y' to simulate"
read*,c1
if (c1=='y') then
    call rk4_solver
    print*,"Simulation Completed"
else
    print*,"Simulation Terminated"
end if
end program edqnm_code
