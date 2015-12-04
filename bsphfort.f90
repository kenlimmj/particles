! ===================================== !
! Specifies what double precision is    !
! ===================================== !
module precision
  implicit none
  integer, parameter :: WP = kind(1.0d0)
end module

! ===================================== !
! Timing and number of particles        !
! ===================================== !
module state
  use precision
  implicit none
  real(WP) :: time, dt, tfin
  integer :: n ! Number of particles
  integer :: iter, itermax, frame
  real(WP) :: mass
end module

! ===================================== !
! Constant values needed for SPH        !
! ===================================== !
module constants
  use precision
  implicit none
  real(WP), parameter :: GRAVITY = -9.81_WP
  real(WP), parameter :: PI = 3.14159265_WP
  real(WP), parameter :: BOX_SIZE = 1.0_WP
  real(WP), parameter :: PARTICLE_DENSITY = 1000.0_WP
  real(WP), parameter :: h = 0.05_WP    ! KERNEL_SIZE
  real(WP), parameter :: BULK_MODULUS = 1000.0_WP
  real(WP), parameter :: E = 0.75  ! Resistution Coefficient

  real(WP), parameter :: ARTIFICIAL_PRESSURE_STRENGTH = 1.0e-5_WP
  real(WP), parameter :: ARTIFICIAL_PRESSURE_RADIUS = 0.3_WP * h
  integer , parameter :: ARTIFICIAL_PRESSURE_POWER = 4

  real(WP), parameter :: VISCOSITY = 1.0e-1_WP

  real(WP) :: PRESSURE_RADIUS_FACTOR

end module constants

! ===================================== !
! Arrays related to Position and        !
! Velocity of particles                 !
! ===================================== !
module posvel
  use precision
  implicit none
  real(WP),dimension(:),allocatable ::  px,  py,  pz  ! Current position
  real(WP),dimension(:),allocatable ::  vx,  vy,  vz  ! Velocity
  real(WP),dimension(:),allocatable ::  vlx,  vly,  vlz  ! Velocity half time step behind
  real(WP),dimension(:),allocatable ::  rho           ! Density
end module

! ===================================== !
! Arrays related to quickly identifying !
! nearest neighbors                     !
! ===================================== !
module neighbors
  use precision
  use constants
  implicit none
  integer,dimension(:),allocatable :: nc
  integer,dimension(:,:),allocatable :: nbs
end module neighbors

! ===================================== !
! Arrays related to forces              !
! ===================================== !
module forces
  use precision
  implicit none
  real(WP),dimension(:),allocatable :: fx, fy, fz  ! Total force
  logical, parameter :: USE_SURFACE_TENSION = .false.
end module forces

module kernels
  implicit none
  contains
    ! ================================================= !
    ! Pressure Kernel                                   !
    ! ================================================= !
    real(WP) function pressurekernel(i,j)
      use constants
      use posvel
      use state
      implicit none
      integer, intent(in) :: i,j
      real(WP) :: dx, dy, dz, r2, c, h2

      dx = px(i) - px(j)
      dy = py(i) - py(j)
      dz = pz(i) - pz(j)
      r2 = dx*dx + dy*dy + dz*dz
      h2 = h*h
      c = 45.0_WP*BULK_MODULUS*0.5/PI*mass

      if(r2 .ge. h2) then
        pressurekernel = 0.0_WP
      else
        pressurekernel = c / h**5 &
                        *(1.0_WP - sqrt(r2/h2))**2 &
                        /sqrt(r2/h2) &
                        *(rho(i)+rho(j)-2.0_WP*PARTICLE_DENSITY) &
                        /(rho(i)*rho(j))
      end if

      return
    end function pressurekernel

    ! ============================================ !
    ! 6th degree polynomial kernel for viscosity   !
    ! ============================================ !
    real(WP) function visckernel(i,j)
      use constants
      use posvel
      use state
      implicit none
      integer, intent(in) :: i,j
      real(WP) :: dx, dy, dz, r2, c, h2

      dx = px(i) - px(j)
      dy = py(i) - py(j)
      dz = pz(i) - pz(j)
      r2 = dx*dx + dy*dy + dz*dz
      h2 = h*h
      c = -45.0_WP*VISCOSITY/PI*mass

      if(r2 .ge. h2) then
        visckernel = 0.0_WP
      else
        visckernel = c / h**5*(1.0_WP - sqrt(r2/h2))/(rho(j))
      end if

      return
    end function visckernel

    ! ===================================== !
    !  Density kernel
    ! ===================================== !
    real(WP) function densitykernel(i,j)
      use constants
      use posvel
      use state
      implicit none
      integer, intent(in) :: i,j
      real(WP) :: dx, dy, dz, r2, c, h2

      dx = px(i) - px(j)
      dy = py(i) - py(j)
      dz = pz(i) - pz(j)
      r2 = dx*dx + dy*dy + dz*dz
      h2 = h*h
      c = 315.0_WP/64.0_WP/PI*mass

      if(r2 .ge. h2) then
        densitykernel = 0.0_WP
      else
        densitykernel = c/h**9 *(h2-r2)**3
      end if

      return
    end function densitykernel

end module kernels


! ===================================== !
! Compute the pressure radius factor    !
! ===================================== !
subroutine computePressureRadiusFactor
  use constants
  implicit none
  real(WP) :: r, c

  r = ARTIFICIAL_PRESSURE_RADIUS
  c = 315.0_WP / (64.0_WP * PI)

  PRESSURE_RADIUS_FACTOR = h**9 / &
    (c * (h*h - r*r)**3)

  return
end subroutine computePressureRadiusFactor

! ===================================== !
! Run the simulation                    !
! ===================================== !
program sph_run
  use state
  use posvel
  implicit none

  time = 0.0_WP
  dt = 0.0_WP ! NEED SMART WAY TO CALCULATE DT FOR STABILITY?

  call init

  iter = 0

  do while (time .lt. tfin)
    iter = iter + 1
    if(iter .gt. itermax) then
      print*, "Exceeded maximum iterations. Exiting..."
      stop
    end if

    call step
    time = time + dt

    if(mod(iter,20) .eq. 0) then
      print*,"Iteration: ",iter," Time: ",time
      print*,"Max U: ", maxval(vx(:))," Max V: ",maxval(vy(:))," Max W: ",maxval(vz(:))
      print*,"  "
    end if

    if(mod(iter,10) .eq. 0) then
      call particle_write
    end if

  end do

  call deinit

end program sph_run

! ================================================ !
! Read in command line arguments and initialize    !
! ================================================ !
subroutine init
  use state
  use constants
  use posvel
  use neighbors
  use forces
  implicit none
  character(30) :: tfin_string, dt_string, itermax_string, init_file
  integer :: i, unit
  real, dimension(:), allocatable :: data

  ! Read in command line arguments for tfin, dt, max_iter, init_file
  call get_command_argument(4,init_file)
  if(len_trim(init_file) .eq. 0) then
    print*, "Incorrect number of command line arguments"
    print*, "Correct usage is sphfort (tfin) (dt) (max_iter) (init_file)"
    print*, "Exiting..."
    stop
  end if
  call get_command_argument(1,tfin_string)
  read (tfin_string, *) tfin

  call get_command_argument(2,dt_string)
  read (dt_string, *) dt

  call get_command_argument(3,itermax_string)
  read (itermax_string, *) itermax

  ! Find length of initial condition file (number of particles)
  unit = 20
  open(unit,file=trim(init_file),status="old",action="read")
  read(unit,*) n    ! First line will be number of particles

  ! Allocate everything
  allocate(px(n)); px = 0.0_WP
  allocate(py(n)); py = 0.0_WP
  allocate(pz(n)); pz = 0.0_WP
  allocate(vx(n)); vx = 0.0_WP
  allocate(vy(n)); vy = 0.0_WP
  allocate(vz(n)); vz = 0.0_WP
  allocate(vlx(n)); vlx = 0.0_WP
  allocate(vly(n)); vly = 0.0_WP
  allocate(vlz(n)); vlz = 0.0_WP
  allocate(rho(n)); rho = 0.0_WP
  allocate(nc(n)); nc = 0
  allocate(nbs(n,1)); nbs = 0
  allocate(fx(n)); fx = 0.0_WP
  allocate(fy(n)); fy = 0.0_WP
  allocate(fz(n)); fz = 0.0_WP

  ! Populate location and velocity from initial condition file
  allocate(data(6))
  do i = 1, n
    read(unit,*) data(:)
    px(i) = data(1)
    py(i) = data(2)
    pz(i) = data(3)
    vx(i) = data(4)
    vy(i) = data(5)
    vz(i) = data(6)
  end do
  close(unit)
  deallocate(data)

  call computePressureRadiusFactor

  call mass_change

  ! Print out initial particle positions
  frame = 0
  call particle_write

  return
end subroutine init

subroutine mass_change
  use posvel
  use state
  use constants
  implicit none
  integer :: i
  real(WP) :: rho0, rho2s, rhos

  mass = 1
  rho0 = PARTICLE_DENSITY
  rho2s = 0.0_WP
  rhos = 0.0_WP

  call compute_density

  do i = 1, n
    rho2s = rho2s + rho(i)*rho(i)
    rhos = rhos + rho(i)
  end do
  mass = mass*rho0*rhos/rho2s

  return
end subroutine mass_change

! ===================================== !
! Deallocate all arrays in the modules    !
! ===================================== !
subroutine deinit
  use state
  use constants
  use posvel
  use neighbors
  use forces
  implicit none

  ! Deallocate everything
  deallocate(px)
  deallocate(py)
  deallocate(pz)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(vlx)
  deallocate(vly)
  deallocate(vlz)
  deallocate(nc)
  deallocate(nbs)
  deallocate(fx)
  deallocate(fy)
  deallocate(fz)
  deallocate(rho)

  return
end subroutine deinit


subroutine particle_write
  use posvel
  use state
  implicit none
  integer :: i, fid
  character(64) :: buffer, start,filename

!  write(buffer,"(ES12.3)") time
  start = "particles_"
  write(buffer,"(a, i5.5)") trim(start), frame
  filename = trim(adjustl(buffer))
  open(fid,file=filename,status="REPLACE",action="WRITE")
  write(fid,*) n
  do i = 1, n
    write(fid,FMT="(6(ES22.15,1x))") px(i), py(i), pz(i), vx(i), vy(i), vz(i)
  end do

  frame = frame + 1

  return
end subroutine

! ===================================== !
! One step of SPH                       !
! ===================================== !
subroutine step
  implicit none

  call neighbor_find
  call compute_density
  call compute_pressure
  call compute_visc
  call ext_forces
  call posvel_update

end subroutine step


subroutine neighbor_find
  use state
  use neighbors
  use posvel
  implicit none
  integer :: i, j, q, p
  integer, dimension(:,:,:),allocatable :: part_count
  integer, dimension(:,:,:,:), allocatable :: binpart
  integer :: nbinx, nbiny, nbinz, max_part_guess
  integer :: binx, biny, binz
  integer :: cbinx, cbiny, cbinz
  integer :: stx, sty, stz
  real(WP) :: y_height, distx, disty, distz, tdist

  deallocate(nbs)
  y_height = maxval(py(:))
  nbinx = ceiling(BOX_SIZE/h) - 1
  nbiny = ceiling(y_height/h) - 1
!  nbinz = ceiling(BOX_SIZE/h) - 1
nbinz = 1


  max_part_guess = 10000*n/(nbinx*nbiny*nbinz)
  allocate(part_count(nbinx,nbiny,nbinz))
  allocate(binpart(max_part_guess,nbinx,nbiny,nbinz)); binpart = 0

  part_count = 0

  do i = 1, n
    binx = NINT((px(i))/(BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((py(i))/(y_height+1.0e-15_WP)*(real(nbiny,WP)-1.0_WP)) + 1
    binz = NINT((pz(i))/(BOX_SIZE)*(real(nbinz,WP)-1.0_WP)) + 1
    part_count(binx, biny, binz) = part_count(binx, biny, binz) + 1
    binpart(part_count(binx,biny,binz),binx,biny,binz) = i
  end do

  max_part_guess = 27*maxval(part_count(:,:,:))
  allocate(nbs(n,max_part_guess)); nbs = 0
  nc = 0
  do i = 1, n
    binx = NINT((px(i))/(BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((py(i))/(y_height+1.0e-15_WP)*(real(nbiny,WP)-1.0_WP)) + 1
    binz = NINT((pz(i))/(BOX_SIZE)*(real(nbinz,WP)-1.0_WP)) + 1
    q = 0
    do stx = -1, 1
      cbinx = binx + stx
      if(cbinx.lt.1 .or. cbinx.gt.nbinx) cycle
      do sty = -1, 1
        cbiny = biny + sty
        if(cbiny.lt.1 .or. cbiny.gt.nbiny) cycle
        do stz = -1, 1
          cbinz = binz + stz
          if(cbinz.lt.1 .or. cbinz.gt.nbinz) cycle
          do j = 1, max_part_guess
            if (binpart(j,cbinx,cbiny,cbinz) .eq. 0) exit
            p = binpart(j,cbinx,cbiny,cbinz)
            distx = px(i) - px(p)
            disty = py(i) - py(p)
            distz = pz(i) - pz(p)
            tdist = sqrt(distx*distx + disty*disty + distz*distz)
            if(tdist .le. h .and. p.ne.i) then
              q = q + 1
              nbs(i,q) = binpart(j,cbinx,cbiny,cbinz)
            end if
          end do
        end do
      end do
    end do
    nc(i) = q
  end do

  deallocate(binpart)
  deallocate(part_count)

  return
end subroutine neighbor_find

subroutine compute_density
  use state
  use constants
  use neighbors
  use kernels
  use posvel
  implicit none
  integer :: i, j, l, nb
  real(WP) :: rhosum

  rho = 315.0_WP/64.0_WP/PI/h**3*mass

  do i = 1, n
    rhosum = 0.0_WP

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      rhosum = rhosum + densitykernel(i,nb)

    end do
  rho(i) = rho(i) + rhosum

  end do

  return
end subroutine compute_density

subroutine compute_pressure
  use constants
  use neighbors
  use posvel
  use forces
  use kernels
  use state
  implicit none
  integer :: i, j, l, nb
  real(WP) :: sx, sy, sz
  real(WP) :: dx, dy, dz
  real(WP) :: k, scorr

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    ! Compute pressure
    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
!      scorr = 0.0_WP
!!      ! Surface Tension
!!      if(USE_SURFACE_TENSION) then
!!      scorr = -ARTIFICIAL_PRESSURE_STRENGTH &
!!              * (poly6kernel(i,nb) * PRESSURE_RADIUS_FACTOR) &
!!              ** ARTIFICIAL_PRESSURE_POWER
!!      end if

      dx = px(i) - px(nb)
      dy = py(i) - py(nb)
      dz = pz(i) - pz(nb)

      k = pressurekernel(i,nb)
      sx = sx + k*dx
      sy = sy + k*dy
      sz = sz + k*dz

    end do

    fx(i) = sx
    fy(i) = sy
    fz(i) = sz

  end do

  return
end subroutine compute_pressure

subroutine compute_visc
  use constants
  use neighbors
  use posvel
  use forces
  use kernels
  use state
  implicit none
  integer :: i, j, l, nb
  real(WP) :: sx, sy, sz
  real(WP) :: dx, dy, dz
  real(WP) :: k, scorr

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    ! Compute pressure
    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)

      dx = vx(i) - vx(nb)
      dy = vy(i) - vy(nb)
      dz = vz(i) - vz(nb)

      k = visckernel(i,nb)
      sx = sx + k*dx
      sy = sy + k*dy
      sz = sz + k*dz

    end do

    fx(i) = fx(i) + sx
    fy(i) = fy(i) + sy
    fz(i) = fz(i) + sz

  end do

  return
end subroutine compute_visc

subroutine ext_forces
  use forces
  use constants
  implicit none

  ! Only gravity for now
  fy = fy + GRAVITY

  return
end subroutine ext_forces

! Leapfrog time integration
! From Bindel 2014 SPH code
subroutine posvel_update
  use posvel
  use state
  use forces
  implicit none

  if(iter.eq.1) then
    vlx = vx
    vly = vy
    vlz = vz

    vlx = vlx + 0.5_WP*dt * fx
    vly = vly + 0.5_WP*dt * fy
    vlz = vlz + 0.5_WP*dt * fz

    vx = vx + dt*fx
    vy = vy + dt*fy
    vz = vz + dt*fz

    px = px + dt*vlx
    py = py + dt*vly
    pz = pz + dt*vlz

  else

    vlx = vlx + dt * fx
    vly = vly + dt * fy
    vlz = vlz + dt * fz

    vx = vlx
    vy = vly
    vz = vlz

    vx = vx + 0.5_WP*dt*fx
    vy = vy + 0.5_WP*dt*fy
    vz = vz + 0.5_WP*dt*fz

    px = px + dt*vlx
    py = py + dt*vly
    pz = pz + dt*vlz
  end if

  call enforce_BCs

  return
end subroutine posvel_update

subroutine enforce_BCs
  use constants
  use posvel
  use state
  implicit none
  integer :: i

  do i = 1, n

    if(px(i) .lt. 0.0_WP) then
      call damp_reflect(i,0.0_WP,px(i),vx(i),vlx(i))
    end if

    if(px(i) .gt. BOX_SIZE) then
      call damp_reflect(i,BOX_SIZE,px(i),vx(i),vlx(i))
    end if

    if(py(i) .lt. 0.0_WP) then
      call damp_reflect(i,0.0_WP,py(i),vy(i),vly(i))
    end if

    if(py(i) .gt. BOX_SIZE) then
      call damp_reflect(i,BOX_SIZE,py(i),vy(i),vly(i))
    end if

    if(pz(i) .lt. 0.0_WP) then
      call damp_reflect(i,0.0_WP,pz(i),vz(i),vlz(i))
    end if

    if(pz(i) .gt. BOX_SIZE) then
      call damp_reflect(i,BOX_SIZE,pz(i),vz(i),vlz(i))
    end if

  end do

  return
end subroutine enforce_BCs

! Taking from 2014 Bindel SPH code
subroutine damp_reflect(i,wall,pos,vel, vel_lag)
  use constants
  use posvel
  implicit none
  integer, intent(in) :: i
  real(WP), intent(in) :: wall
  real(WP), intent(inout) :: pos, vel, vel_lag
  real(WP) :: tbounce

  if(vel .eq. 0.0_WP) return

  ! Scale back distance to get when collision occurred
  tbounce = (pos-wall)/vel
  px(i) = px(i) -(1.0_WP-E)*tbounce*vx(i)
  py(i) = py(i) -(1.0_WP-E)*tbounce*vy(i)
  pz(i) = pz(i) -(1.0_WP-E)*tbounce*vz(i)

  ! Reflect
  pos = 2.0_WP * wall-pos
  vel = -vel
  vel_lag = -vel_lag

  ! Damp all velocities
  vx(i) = vx(i) * E
  vy(i) = vy(i) * E
  vz(i) = vz(i) * E
  vlx(i) = vlx(i) * E
  vly(i) = vly(i) * E
  vlz(i) = vlz(i) * E

  return
end subroutine damp_reflect




