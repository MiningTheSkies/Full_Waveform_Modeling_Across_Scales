  program diffusion

  implicit none

  include "constants.h"

! number of timesteps
  integer, parameter :: NSTEP = 30000
 
! time step in seconds
  double precision, parameter :: DT = 100000000. ! s
 
! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.
 
! model parameters  (SI)
  double precision, parameter :: LENGTH = 3.0d+03 ! m
  double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
  double precision, parameter :: THERMALCONDUCTIVITY = 10.0d-01 ! cal/m/s/K
  double precision, parameter :: HEATCAPACITY = 0.3d+03 ! cal/kg/K

  integer ispec,i,j,iglob,itime

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLL) :: xigll
 
! weights
  double precision, dimension(NGLL) :: wgll
 
! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLL,NGLL) :: hprime

! anchors
  double precision, dimension(NSPEC) :: x1,x2

! global grid points
  double precision, dimension(NGLOB) :: x

! material properties
  double precision, dimension(NGLL,NSPEC) :: rho,heat_capacity,thermal_conductivity

! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

! local mass matrix
  double precision mass_local

! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

! temperature and temperature time derivative
  double precision, dimension(NGLOB) :: temperature,dtemperature_dt

! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

! time marching
  double precision deltat,deltatover2
  double precision dh,diffusivity,time_step

! end fluxes
  double precision flux_1,flux_NGLOB

! end temperatures
  double precision temperature_1,temperature_NGLOB

! derivatives
  double precision dtdx,flux,templ,temp(NGLL)

! movie
  character(len=50) moviefile

! declaring some variables
  integer iglob1, iglob2, iglobj
  double precision, dimension(NGLOB) :: rhs_global
  double precision rhs_local
  double precision :: stiffness_local

!++++++++++++++++++++++++++++++++++++++++++++++++++

  call define_derivative_matrix(xigll,wgll,hprime)

! evenly spaced anchors between 0 and 1
  do ispec = 1,NSPEC
    x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
    x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
  enddo

! set up the mesh properties
  do ispec = 1,NSPEC
    do i = 1,NGLL
      rho(i,ispec) = DENSITY
      thermal_conductivity(i,ispec) = THERMALCONDUCTIVITY
      heat_capacity(i,ispec) = HEATCAPACITY
      dxidx(i,ispec) = 2. / (x2(ispec)-x1(ispec))
      jacobian(i,ispec) = (x2(ispec)-x1(ispec)) / 2.
    enddo
  enddo

! set up local to global numbering
  iglob = 1
  do ispec = 1,NSPEC
    do i = 1,NGLL
      if(i > 1) iglob = iglob+1
      ibool(i,ispec) = iglob
    enddo
  enddo

! get the global grid points
  do ispec = 1,NSPEC
    do i = 1,NGLL
      iglob = ibool(i,ispec)
      x(iglob) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
    enddo
  enddo


! ! declaring some variables
!   integer iglob1, iglob2, iglobj
!   double precision, dimension(NGLOB) :: rhs_global
!   double precision stiffness_local


! calculate the global mass matrix 'mass_global'
! put your codes here

  ! set up iglob, then local mass, then global mass
  ! 2 loops, first to NSPEC, second to NGLL
  ! local mass is w * rho * cp * J
  ! directly add the local mass to the global mass at that iglob element
  mass_global(:) = 0.
  do ispec = 1, NSPEC
    do i = 1, NGLL
      iglob = ibool(i,ispec)
      mass_local = wgll(i) * rho(i, ispec) * heat_capacity(i, ispec) * jacobian(i, ispec)
      mass_global = mass_global(iglob) + mass_local
    enddo
  enddo



! estimate the time step 'time_step'
! put your codes here

  ! define deltat as the DT that is pre-defined, then deltatover2 as DT /2
  deltat = DT
  deltatover2 = DT / 2

  print *,'time step estimate: ',time_step,' seconds'



! set up the boundary conditions ... I need to fix this and add more using slide 7. ***DON'T FORGET TO DO THIS!***
! put your codes here

  ! initial temperatures
  temperature_1 = 10
  temperature_NGLOB = 0



! initialize
  temperature(:) = 0.
  dtemperature_dt(:) = 0.

  do itime = 1,NSTEP



! update temperature
! put your codes here

    ! set left bound temp
    temperature(1) = temperature_1
    ! set right bound temp
    temperature(NGLOB) = temperature_NGLOB
    ! set temps in between - this included in the loop later too
    temperature(:) = temperature(:) + deltatover2 * dtemperature_dt(:)


    ! I'm not really sure where I'm supposed to define the rhs section since it doesn't correspond to any of the "put your codes here" sections, so I'll just leave it here

    ! define rhs stiffness matrix

    ! initialize rhs_global
    rhs_global(:) = 0.

    ! loop over NPEC, then NGLL
    do ispec = 1,NSPEC
      do i = 1,NGLL

        ! initialize rhs_local
        rhs_local = 0.
        iglob = ibool(i,ispec)

        ! loop over NGLL for setting stiffness & rhs
        do j = 1,NGLL
          iglobj = ibool(j,ispec)
          ! stiffness
          stiffness_local = sum( wgll(:) * thermal_conductivity(:,ispec) * hprime(i,:) * hprime(j,:) * &
            dxidx(:,ispec) * dxidx(:,ispec) * jacobian(:, ispec))
          rhs_local = rhs_local - stiffness_local * temperature(iglobj)
        enddo
        ! set left & right iglobs as 1 and the numbter of GLL points
        iglob1 = ibool(1,ispec)
        iglob2 = ibool(NGLL,ispec)
        ! special cases for first and last points... first point has rhs_local - [other] while last point has rhs_local + [other]
        if(ispec == 1 .and. i == 1) then
          rhs_local = rhs_local - thermal_conductivity(i,ispec) * sum(temperature(iglob1:iglob2) * hprime(:,i) * dxidx(i, ispec))
        else if (ispec == NSPEC .and. i == NGLL) then
          rhs_local = rhs_local + thermal_conductivity(i,ispec) * sum(temperature(iglob1:iglob2) * hprime(:,i) * dxidx(i, ispec))
        endif
        ! add rhs_local to the rhs_global at particular iglob element
        rhs_global(iglob) = rhs_global(iglob) + rhs_local
      enddo
    enddo
    ! update temperature rate with rhs_global divided by the mass_gobal
    dtemperature_dt(:) = rhs_global(:) / mass_global(:)
    ! update the old temp with the new temperature rate modified by half of DT
    temperature(:) = temperature(:) + deltatover2 * dtemperature_dt(:)






! write out snapshots
    if(mod(itime-1,1000) == 0) then
      write(moviefile,10) itime
10    format('snapshot',i5.5)
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(temperature(iglob))
      enddo
      close(10)
    endif

  enddo ! end time loop

  end program diffusion
