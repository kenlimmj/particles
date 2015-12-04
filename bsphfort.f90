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
  real(WP), parameter :: GRAVITY = -9.8_WP
  real(WP), parameter :: PI = 3.14159265_WP
  real(WP), parameter :: BOX_SIZE =1.0_WP
  real(WP), parameter :: PARTICLE_DENSITY = 1000.0_WP
  real(WP), parameter :: h = 0.05_WP      ! KERNEL_SIZE
  real(WP), parameter :: h2 = h*h         ! KERNEL_SIZE**2
  real(WP), parameter :: h3 = h*h2         ! KERNEL_SIZE**2
  real(WP), parameter :: h5 = h3*h2     ! KERNEL_SIZE**5
  real(WP), parameter :: h9 = h5*h2*h2     ! KERNEL_SIZE**9
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
  integer,dimension(:,:,:),allocatable :: part_count
  integer,dimension(:,:,:,:),allocatable :: binpart

  ! Expected max particles in single bin
  integer, parameter :: nbinx = floor(BOX_SIZE/h) - 1
  integer, parameter :: nbiny = floor(BOX_SIZE/h) - 1
  integer, parameter :: nbinz = floor(BOX_SIZE/h) - 1
  integer :: max_part_guess
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
      real(WP) :: dx, dy, dz, r2, c

      dx = px(i) - px(j)
      dy = py(i) - py(j)
      dz = pz(i) - pz(j)
      r2 = dx*dx + dy*dy + dz*dz

      if(r2 .gt. h2) then
        pressurekernel = 0.0_WP
      else
        c = 45.0_WP*mass/PI/h5/rho(j)/rho(i)*(1.0_WP-sqrt(r2/h2))
        pressurekernel = c * BULK_MODULUS*0.5_WP &
                        *(rho(i)+rho(j)-2.0_WP*PARTICLE_DENSITY) &
                        *(1.0_WP - sqrt(r2/h2))/sqrt(r2/h2)
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
      real(WP) :: dx, dy, dz, r2, c

      dx = px(i) - px(j)
      dy = py(i) - py(j)
      dz = pz(i) - pz(j)
      r2 = dx*dx + dy*dy + dz*dz

      if(r2 .gt. h2) then
        visckernel = 0.0_WP
      else
        c = 45.0_WP*mass/PI/h5/rho(j)/rho(i)*(1.0_WP-sqrt(r2/h2))
        visckernel = -c * VISCOSITY
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
      real(WP) :: dx, dy, dz, r2, c

      dx = px(i) - px(j)
      dy = py(i) - py(j)
      dz = pz(i) - pz(j)
      r2 = dx*dx + dy*dy + dz*dz

      if(r2 .gt. h2) then
        densitykernel = 0.0_WP
      else
        c = (315.0_WP/64.0_WP/PI)*mass/h9
        densitykernel = c*(h2-r2)*(h2-r2)*(h2-r2)
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

  PRESSURE_RADIUS_FACTOR = h9 / &
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
  dt = 0.0_WP

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

    if(mod(iter,100) .eq. 0) then
      print*,"Iteration: ",iter," Time: ",time
      print*,"Max U: ", maxval(vx(:))," Max V: ",maxval(vy(:))," Max W: ",maxval(vz(:))
      print*,"  "
    end if

    if(mod(iter,100) .eq. 0) then
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
  open(unit,file=trim(init_file),action="read")
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
  allocate(fx(n)); fx = 0.0_WP
  allocate(fy(n)); fy = 0.0_WP
  allocate(fz(n)); fz = 0.0_WP


  ! Guess for allocation that at most particles
  ! will be concentrated 100 times above even
  max_part_guess = 100*n/(nbinx*nbiny*nbinz)
  allocate(part_count(nbinx,nbiny,nbinz))
  allocate(binpart(max_part_guess,nbinx,nbiny,nbinz)); binpart = 0
  allocate(nbs(n,max_part_guess)); nbs = 0

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
  call neighbor_find
  call compute_density
  rho0 = PARTICLE_DENSITY
  rho2s = 0.0_WP
  rhos = 0.0_WP

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
  deallocate(part_count)
  deallocate(binpart)
  deallocate(nbs)

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
  fid = 0
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
  use posvel
  use state
  implicit none
  real(WP) :: t1, t2, t3, t4, t5, t6, t7

  call cpu_time(t1)
  call neighbor_find
  call cpu_time(t2)
  call compute_density
  call cpu_time(t3)
  call compute_pressure
  call cpu_time(t4)
  call compute_visc
  call cpu_time(t5)
  call ext_forces
  call cpu_time(t6)
  call posvel_update
  call cpu_time(t7)

!  print*, "neighbor time = ", t2-t1
!  print*, "density time = ", t3-t2
!  print*, "pressure time = ", t4-t3
!  print*, "visc time = ", t5-t4
!  print*, "ext forces time = ", t6-t5
!  print*, "posvel time = ", t7-t6
!  print*, "total time = ", t7-t1
!  print*, " "
!  print*, " "


end subroutine step


subroutine neighbor_find
  use state
  use neighbors
  use posvel
  implicit none
  integer :: i, j, q, p
  integer :: binx, biny, binz
  integer :: cbinx, cbiny, cbinz
  integer :: stx, sty, stz
  real(WP) :: y_height, distx, disty, distz, tdist

  part_count = 0
  binpart = 0

  do i = 1, n
    binx = NINT((px(i))/(BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((py(i))/(BOX_SIZE)*(real(nbiny,WP)-1.0_WP)) + 1
    binz = NINT((pz(i))/(BOX_SIZE)*(real(nbinz,WP)-1.0_WP)) + 1
    part_count(binx, biny, binz) = part_count(binx, biny, binz) + 1
    binpart(part_count(binx,biny,binz),binx,biny,binz) = i
  end do

  nbs = 0
  nc = 0
  do i = 1, n
    binx = NINT((px(i))/(BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((py(i))/(BOX_SIZE)*(real(nbiny,WP)-1.0_WP)) + 1
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
            if(tdist .lt. h .and. p.ne.i) then
              q = q + 1
              nbs(i,q) = binpart(j,cbinx,cbiny,cbinz)
            end if
          end do
        end do
      end do
    end do
    nc(i) = q
  end do

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
  real(WP) :: rhosum, c, r2
  real(WP) :: dx, dy, dz, kpd

  rho = (315.0_WP/64.0_WP/PI)*mass/h3
  c = 315.0_WP/64.0_WP/PI*mass/h9
  do i = 1, n
    rhosum = 0.0_WP

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)

      dx = px(i) - px(nb)
      dy = py(i) - py(nb)
      dz = pz(i) - pz(nb)
      r2 = dx*dx + dy*dy + dz*dz
      kpd = h2 - r2
      rhosum = rhosum + c * kpd * kpd * kpd

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
  real(WP) :: c, r2, m, q


  c = 45.0_WP*mass/PI/h5*BULK_MODULUS*0.5_WP

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
      r2 = dx*dx + dy*dy + dz*dz
      m = sqrt(r2/h2)
      q = 1.0_WP - m


      k = c / rho(nb) * q &
          * (rho(i)+rho(nb)-2.0_WP*PARTICLE_DENSITY) &
          * q / m


      sx = sx + k*dx
      sy = sy + k*dy
      sz = sz + k*dz

    end do

    fx(i) = sx/rho(i)
    fy(i) = sy/rho(i)
    fz(i) = sz/rho(i)

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
  real(WP) :: dvx, dvy, dvz
  real(WP) :: dx, dy, dz
  real(WP) :: k, c, r2

  c = -45.0_WP*mass/PI/h5*VISCOSITY

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    ! Compute viscous
    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)

      dvx = vx(i) - vx(nb)
      dvy = vy(i) - vy(nb)
      dvz = vz(i) - vz(nb)
      dx = px(i) - px(nb)
      dy = py(i) - py(nb)
      dz = pz(i) - pz(nb)
      r2 = dx*dx + dy*dy + dz*dz

      k = c / rho(nb) * (1.0_WP-sqrt(r2/h2))

      sx = sx + k*dvx
      sy = sy + k*dvy
      sz = sz + k*dvz

    end do

    fx(i) = fx(i) + sx/rho(i)
    fy(i) = fy(i) + sy/rho(i)
    fz(i) = fz(i) + sz/rho(i)

  end do

  return
end subroutine compute_visc

subroutine ext_forces
  use forces
  use constants
  use posvel
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

  if(abs(vel) .le. 1.0e-7_WP) return

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




