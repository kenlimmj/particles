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
  integer :: iter, itermax
end module

! ===================================== !
! Constant values needed for SPH        !
! ===================================== !
module constants
  use precision
  implicit none
  integer,  parameter :: ITERATIONS_PER_STEP = 10
  real(WP), parameter :: GRAVITY = 1.0_WP
  real(WP), parameter :: PI = 3.14159265_WP
  real(WP), parameter :: BOX_SIZE = 1.0_WP
  real(WP), parameter :: PARTICLE_DENSITY = 6378.0_WP * 100.0_WP
  real(WP), parameter :: PARTICLE_MASS = 1.0_WP
  real(WP), parameter :: FORCE_EPSILON = 600.0_WP
  real(WP), parameter :: KERNEL_SIZE = 0.1_WP

  real(WP), parameter :: ARTIFICIAL_PRESSURE_STRENGTH = 1.0e-5_WP
  real(WP), parameter :: ARTIFICIAL_PRESSURE_RADIUS = 0.3_WP * KERNEL_SIZE
  integer , parameter :: ARTIFICIAL_PRESSURE_POWER = 4

  real(WP), parameter :: ARTIFICIAL_VISCOSITY = 1.0e-2_WP
  real(WP), parameter :: VORTICITY_COEFFICIENT = 1.0e-4_WP
  logical,  parameter :: USE_VORTICITY_CONFINEMENT = .true.

  real(WP) :: PRESSURE_RADIUS_FACTOR

  real(WP),dimension(:), allocatable :: lm ! Lambda

end module constants

! ===================================== !
! Arrays related to Position and        !
! Velocity of particles                 !
! ===================================== !
module posvel
  use precision
  implicit none
  real(WP),dimension(:),allocatable ::  ox,  oy,  oz  ! Original position
  real(WP),dimension(:),allocatable ::  px,  py,  pz  ! Current position
  real(WP),dimension(:),allocatable :: cpx, cpy, cpz  ! Candidate position
  real(WP),dimension(:),allocatable ::  vx,  vy,  vz  ! Velocity
  real(WP),dimension(:),allocatable :: cvx, cvy, cvz ! Candidate velocity
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
  integer, parameter :: num_of_bins = ceiling(2*BOX_SIZE/KERNEL_SIZE) ! Bins in each direction
end module neighbors

! ===================================== !
! Arrays related to forces              !
! ===================================== !
module forces
  use precision
  implicit none
  real(WP),dimension(:),allocatable :: fx, fy, fz  ! Total force
  real(WP),dimension(:),allocatable :: dpx, dpy, dpz  ! Pressure
  real(WP),dimension(:),allocatable :: vox, voy, voz  ! Vorticity
end module forces

module kernels
  implicit none
  contains
    ! ===================================== !
    ! 6th degree polynomial kernel          !
    ! ===================================== !
    real(WP) function poly6kernel(i,j)
      use constants
      use posvel
      implicit none
      integer, intent(in) :: i,j
      real(WP) :: dx, dy, dz, r, c

      dx = cpx(i) - cpx(j)
      dy = cpy(i) - cpy(j)
      dz = cpz(i) - cpz(j)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      c = 315.0_WP / (64.0_WP * PI)

      if(r .gt. KERNEL_SIZE) then
        poly6kernel = 0.0_WP
      else
        poly6kernel = c * (KERNEL_SIZE*KERNEL_SIZE - r*r)**3 &
                      / KERNEL_SIZE**9
      end if

      return
    end function poly6kernel

    ! ===================================== !
    !   Viscosity kernel                      !
    ! ===================================== !
    real(WP) function viscositykernel(i,j)
      use constants
      use posvel
      implicit none
      integer, intent(in) :: i,j
      real(WP) :: dx, dy, dz, r, c

      dx = cpx(i) - cpx(j)
      dy = cpy(i) - cpy(j)
      dz = cpz(i) - cpz(j)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      c = 7.5_WP / PI

      if(r .gt. KERNEL_SIZE) then
        viscositykernel = 0.0_WP
      else
        viscositykernel = c * (KERNEL_SIZE*KERNEL_SIZE - r*r)**3 &
                          / KERNEL_SIZE**9
      end if

      return
    end function viscositykernel
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

  PRESSURE_RADIUS_FACTOR = 1.0_WP / &
    (c * (KERNEL_SIZE*KERNEL_SIZE - r*r)**3 &
    / KERNEL_SIZE**9)

  return
end subroutine computePressureRadiusFactor



! ===================================== !
! Spiky gradient kernel                 !
! ===================================== !
subroutine spikykernel(i,j, gv)
  use constants
  use posvel
  implicit none
  integer, intent(in) :: i,j
  real(WP), dimension(3), intent(out) :: gv
  real(WP) :: dx, dy, dz, r, c
  real(WP) :: f

  dx = cpx(i) - cpx(j)
  dy = cpy(i) - cpy(j)
  dz = cpz(i) - cpz(j)
  r = sqrt(dx*dx + dy*dy + dz*dz)
  c = 15.0_WP / PI
  if(r .lt. 1.0e-4_WP) r = 1.0e-4_WP

  if( r .gt. KERNEL_SIZE) then
    gv(1) = 0.0_WP
    gv(2) = 0.0_WP
    gv(3) = 0.0_WP
  else
    f = c / KERNEL_SIZE**6 * (KERNEL_SIZE - r)**2 / r
    gv(1) = f * dx
    gv(2) = f * dy
    gv(3) = f * dz
  end if

  return
end subroutine spikykernel

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

    if(mod(itermax,10) .eq. 0) then
      print*,"Iteration: ",iter," Time: ",time
      print*,"Max U: ", maxval(vx(:))," Max V: ",maxval(vy(:))," Max W: ",maxval(vz(:))
      print*,"  "
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
  allocate(lm(n))
  allocate(ox(n))
  allocate(oy(n))
  allocate(oz(n))
  allocate(px(n))
  allocate(py(n))
  allocate(pz(n))
  allocate(cpx(n))
  allocate(cpy(n))
  allocate(cpz(n))
  allocate(vx(n))
  allocate(vy(n))
  allocate(vz(n))
  allocate(cvx(n))
  allocate(cvy(n))
  allocate(cvz(n))
  allocate(nc(n))
  allocate(nbs(n,1))
  allocate(fx(n))
  allocate(fy(n))
  allocate(fz(n))
  allocate(dpx(n))
  allocate(dpy(n))
  allocate(dpz(n))
  allocate(vox(n))
  allocate(voy(n))
  allocate(voz(n))

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


  return
end subroutine init

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
  deallocate(lm)
  deallocate(ox)
  deallocate(oy)
  deallocate(oz)
  deallocate(px)
  deallocate(py)
  deallocate(pz)
  deallocate(cpx)
  deallocate(cpy)
  deallocate(cpz)
  deallocate(vx)
  deallocate(vy)
  deallocate(vz)
  deallocate(cvx)
  deallocate(cvy)
  deallocate(cvz)
  deallocate(nc)
  deallocate(nbs)
  deallocate(fx)
  deallocate(fy)
  deallocate(fz)
  deallocate(dpx)
  deallocate(dpy)
  deallocate(dpz)
  deallocate(vox)
  deallocate(voy)
  deallocate(voz)

  return
end subroutine deinit


! ===================================== !
! One step of SPH                       !
! ===================================== !
subroutine step
  use constants, only : ITERATIONS_PER_STEP, USE_VORTICITY_CONFINEMENT
  implicit none
  integer :: subiter

  ! Apply the forces
  call force_apply

  ! Compute candidate velocities and positions
  call posvel_candidate

  ! Find neighbors and compute lambda
  call neighbor_find

  do subiter = 1, ITERATIONS_PER_STEP
    ! Compute lamnda
    call compute_lambda

    ! Jiggle the particles
    call particle_jiggle
  end do

  ! Confine vorticity, update vel, and calculate force
  if(USE_VORTICITY_CONFINEMENT) then
    call compute_omega
    call vorticity_force
  end if

  ! Update particle positions due to vorticity
  call viscosity_update

  ! Update particle position
  call pos_update

end subroutine step



subroutine force_apply
  use forces
  use constants
  implicit none

  ! Apply forces
  fx(:) = vox(:)
  fy(:) = GRAVITY
  fy(:) = fy(:) + voy(:)
  fz(:) = voz(:)

  return
end subroutine force_apply

subroutine posvel_candidate
  use forces
  use posvel
  use state
  implicit none

  ! Update candidate velocity
  cvx(:) = vx(:) + fx(:) * dt
  cvy(:) = vy(:) + fy(:) * dt
  cvz(:) = vz(:) + fz(:) * dt

  ! TODO VELOCITY DAMPING (IF USED) GOES HERE
  cpx(:) = px(:) + cvx(:) * dt
  cpy(:) = py(:) + cvy(:) * dt
  cpz(:) = pz(:) + cvz(:) * dt

  return
end subroutine posvel_candidate

subroutine neighbor_find
  use state
  use neighbors
  use posvel
  implicit none
  integer :: i, j, q
  integer, dimension(num_of_bins, num_of_bins, num_of_bins) :: part_count
  integer, dimension(:,:,:,:), allocatable :: binloc
  integer :: binx, biny, binz, max_part_guess
  integer :: cbinx, cbiny, cbinz
  integer :: stx, sty, stz

  deallocate(nbs)
  max_part_guess = 1000*n/num_of_bins**3
  allocate(binloc(max_part_guess,num_of_bins,num_of_bins, num_of_bins))

  part_count = 0.0
  do i = 1, n
    binx = cpx(i)/((2*BOX_SIZE)/num_of_bins) + num_of_bins/2 + 1
    biny = cpy(i)/((2*BOX_SIZE)/num_of_bins) + num_of_bins/2 + 1
    binz = cpz(i)/((2*BOX_SIZE)/num_of_bins) + num_of_bins/2 + 1
    part_count(binx, biny, binz) = part_count(binx, biny, binz) + 1
    binloc(part_count(binx,biny,binz),binx,biny,binz) = i
  end do

  max_part_guess = 27*maxval(part_count(:,:,:))
  allocate(nbs(n,max_part_guess))
  do i = 1, n
    binx = cpx(i)/((2*BOX_SIZE)/num_of_bins) + num_of_bins/2 + 1
    biny = cpy(i)/((2*BOX_SIZE)/num_of_bins) + num_of_bins/2 + 1
    binz = cpz(i)/((2*BOX_SIZE)/num_of_bins) + num_of_bins/2 + 1
    q = 0
    do stx = -1, 1
      do sty = -1, 1
        do stz = -1, 1
          do j = 1, max_part_guess
            cbinx = binx + stx
            cbiny = biny + sty
            cbinz = binz + stz
            if (binloc(j,cbinx,cbiny,cbinz) .eq. 0) exit
            q = q + 1
            nbs(i,q) = binloc(j,cbinx,cbiny,cbinz)
          end do
        end do
      end do
    end do
    nc(i) = q
  end do

  deallocate(binloc)

  return
end subroutine neighbor_find

subroutine compute_lambda
  use state
  use constants
  use neighbors
  use kernels
  implicit none
  integer :: i, j, l, nb
  real(WP) :: rho, numer, denom
  real(WP) :: sx, sy, sz
  real(WP), dimension(3) :: gv

  do i = 1, n
    rho = 0.0_WP
    numer = 0.0_WP
    denom = 0.0_WP

    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      rho = rho + PARTICLE_MASS * poly6kernel(i,nb)

      call spikykernel(i,nb, gv)
      sx = sx + gv(1)
      sy = sy + gv(2)
      sz = sz + gv(3)

      ! I think this is the same as what is in the C version
      denom = denom + (gv(1)*gv(1)+ gv(2)*gv(2) + gv(3)*gv(3)) &
                      /(PARTICLE_DENSITY*PARTICLE_DENSITY)
    end do
    numer = rho / PARTICLE_DENSITY - 1.0

    denom = denom + (sx*sx + sy*sy + sz*sz) &
                    /(PARTICLE_DENSITY*PARTICLE_DENSITY)

    lm(i) = -numer / (denom + FORCE_EPSILON)

  end do

  return
end subroutine compute_lambda

subroutine particle_jiggle
  use constants
  use neighbors
  use posvel
  use forces
  use kernels
  use state
  implicit none
  integer :: i, j, l, nb
  real(WP) :: sx, sy, sz
  real(WP), dimension(3) :: gv
  real(WP) :: k, scorr

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    ! Compute pressure
    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      scorr = 0.0_WP
      ! Surface Tension
      ! TODO ENABLE/DISABLE SURFACE TENSION
      scorr = -ARTIFICIAL_PRESSURE_STRENGTH &
              * (poly6kernel(i,nb) * PRESSURE_RADIUS_FACTOR) &
              ** ARTIFICIAL_PRESSURE_POWER

      call spikykernel(i,nb, gv)
      k = lm(i) + lm(nb) + scorr
      sx = sx + gv(1) * k
      sy = sy + gv(2) * k
      sz = sz + gv(3) * k
    end do

    dpx(i) = sx / PARTICLE_DENSITY
    dpy(i) = sy / PARTICLE_DENSITY
    dpz(i) = sz / PARTICLE_DENSITY
  end do

    ! Update candidate position
    cpx(:) = cpx(:) + dpx(:)
    cpy(:) = cpy(:) + dpy(:)
    cpz(:) = cpz(:) + dpz(:)

    ! Update candidate velocity
    cvx(:) = (cpx(:) - px(:)) / dt
    cvy(:) = (cpy(:) - py(:)) / dt
    cvz(:) = (cpz(:) - pz(:)) / dt

  ! Resolve collisions
  do i = 1, n

    ! Floor
    if(cpy(i) .lt. 0.0_WP) then
      cpy(i) = 0.0_WP
      cvy(i) = -cvy(i)
    end if

    ! Left Wall
    if(cpx(i) .lt. -BOX_SIZE) then
      cpx(i) = -BOX_SIZE
      cvx(i) = -cvx(i)
    end if

    ! Right Wall
    if(cpx(i) .gt. BOX_SIZE) then
      cpx(i) = BOX_SIZE
      cvx(i) = -cvx(i)
    end if

    ! Back Wall
    if(cpz(i) .lt. -BOX_SIZE) then
      cpz(i) = -BOX_SIZE
      cvz(i) = -cvz(i)
    end if

    ! Front Wall
    if(cpz(i) .gt. BOX_SIZE) then
      cpz(i) = BOX_SIZE
      cvz(i) = -cvz(i)
    end if

  end do

  return
end subroutine particle_jiggle

subroutine confine_vorticity
  use forces
  implicit none

  ! STILL NEED TO INCLUDE THIS

  return
end subroutine

subroutine compute_omega
  use state
  use posvel
  use forces
  use neighbors
  implicit none
  integer :: i, j, l, nb
  real(WP) :: sx, sy, sz
  real(WP) :: wx, wy, wz
  real(WP), dimension(3) :: gv

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      wx = cvx(nb) - cvx(i)
      wy = cvy(nb) - cvy(i)
      wz = cvz(nb) - cvz(i)

      call spikykernel(i,nb, gv)
      sx = sx + wy * gv(3) - wz * gv(2)
      sy = sy + wz * gv(1) - wx * gv(3)
      sz = sz + wx * gv(2) - wy * gv(1)
    end do

    ! Reuse dp array
    dpx(i) = sx
    dpy(i) = sy
    dpz(i) = sz
  end do

  return
end subroutine compute_omega

subroutine vorticity_force
  use state
  use constants
  use forces
  use neighbors
  implicit none
  integer :: i, j, l, nb
  real(WP) :: sx, sy, sz
  real(WP) :: dpmag
  real(WP), dimension(3) :: gv

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      call spikykernel(i,nb, gv)
      dpmag = sqrt(dpx(nb)*dpx(nb) + dpy(nb)*dpy(nb) + dpz(nb)*dpz(nb))
      gv(1) = gv(1) * dpmag / PARTICLE_DENSITY
      gv(2) = gv(2) * dpmag / PARTICLE_DENSITY
      gv(3) = gv(3) * dpmag / PARTICLE_DENSITY

      sx = sx + gv(1)
      sy = sy + gv(2)
      sz = sz + gv(3)
    end do

    dpmag = sqrt(sx*sx + sy*sy + sz*sz)
    if( dpmag .lt. 1.0e-5_WP) dpmag = 1.0e-5_WP

    sx = sx / dpmag
    sy = sy / dpmag
    sz = sz / dpmag

    vox(i) = (sy*dpz(i) - sz*dpy(i)) * VORTICITY_COEFFICIENT
    voy(i) = (sz*dpx(i) - sx*dpz(i)) * VORTICITY_COEFFICIENT
    voz(i) = (sx*dpy(i) - sy*dpx(i)) * VORTICITY_COEFFICIENT

  end do

  return
end subroutine vorticity_force

subroutine viscosity_update
  use state
  use constants
  use posvel
  use neighbors
  use kernels
  implicit none
  integer :: i, j, l, nb
  real(WP) :: sx, sy, sz
  real(WP) :: dx, dy, dz, k

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      dx = cvx(nb) - cvx(i)
      dy = cvy(nb) - cvy(i)
      dz = cvz(nb) - cvz(i)

      k = viscositykernel(i,nb)
      sx = sx + k * dx
      sy = sy + k * dy
      sz = sz + k * dz
    end do

    cvx(i) = cvx(i) + ARTIFICIAL_VISCOSITY * sx
    cvy(i) = cvy(i) + ARTIFICIAL_VISCOSITY * sy
    cvz(i) = cvz(i) + ARTIFICIAL_VISCOSITY * sz

  end do

  return
end subroutine viscosity_update

subroutine pos_update
  use posvel
  implicit none

  ! Velocities
  vx = cvx
  vy = cvy
  vz = cvz

  ! Positions
  px = cpx
  py = cpy
  pz = cpz

  return
end subroutine pos_update




