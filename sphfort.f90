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
  integer,  parameter :: ITERATIONS_PER_STEP = 20
  real(WP), parameter :: GRAVITY = 1.0_WP
  real(WP), parameter :: PI = 3.14159265_WP
  real(WP), parameter :: BOX_SIZE = 0.5_WP
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
  logical, parameter :: USE_INITIAL_VELOCITY_DAMPING = .false.
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
  logical, parameter :: USE_SURFACE_TENSION = .true.
end module forces

module kernels
  implicit none
  contains
    ! ================================================= !
    ! 6th degree polynomial kernel for density          !
    ! ================================================= !
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

      if(r .ge. KERNEL_SIZE) then
        poly6kernel = 0.0_WP
      else
        poly6kernel = c * (KERNEL_SIZE*KERNEL_SIZE - r*r)**3 &
                      / KERNEL_SIZE**9
      end if

      return
    end function poly6kernel

    ! ============================================ !
    ! 6th degree polynomial kernel for viscosity   !
    ! ============================================ !
    real(WP) function poly6visckernel(i,j)
      use constants
      use posvel
      implicit none
      integer, intent(in) :: i,j
      real(WP) :: dx, dy, dz, r, c

      dx = vx(i) - vx(j)
      dy = vy(i) - vy(j)
      dz = vz(i) - vz(j)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      c = 315.0_WP / (64.0_WP * PI)

      if(r .ge. KERNEL_SIZE) then
        poly6visckernel = 0.0_WP
      else
        poly6visckernel = c * (KERNEL_SIZE*KERNEL_SIZE - r*r)**3 &
                      / KERNEL_SIZE**9
      end if

      return
    end function poly6visckernel

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


      if(r .ge. KERNEL_SIZE) then
        viscositykernel = 0.0_WP
      else
        viscositykernel = c / KERNEL_SIZE**3 &
                            * (-r**3/(2.0_WP*KERNEL_SIZE**3) &
                              + r*r/(KERNEL_SIZE*KERNEL_SIZE) &
                              + KERNEL_SIZE/(2.0_WP*r) - 1.0_WP)


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

  PRESSURE_RADIUS_FACTOR = KERNEL_SIZE**9 / &
    (c * (KERNEL_SIZE*KERNEL_SIZE - r*r)**3)

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
  c = 15.0/PI
  if(r .lt. 1.0e-4_WP) r = 1.0e-4_WP

  if( r .ge. KERNEL_SIZE) then
    gv(1) = 0.0_WP
    gv(2) = 0.0_WP
    gv(3) = 0.0_WP
  else
    f = c * (KERNEL_SIZE-r)**3/KERNEL_SIZE**6
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

    if(mod(iter,1) .eq. 0) then
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
  open(unit,file=trim(init_file),status="old",action="read")
  read(unit,*) n    ! First line will be number of particles

  ! Allocate everything
  allocate(lm(n)); lm = 0.0_WP
  allocate(px(n)); px = 0.0_WP
  allocate(py(n)); py = 0.0_WP
  allocate(pz(n)); pz = 0.0_WP
  allocate(cpx(n)); cpx = 0.0_WP
  allocate(cpy(n)); cpy = 0.0_WP
  allocate(cpz(n)); cpz = 0.0_WP
  allocate(vx(n)); vx = 0.0_WP
  allocate(vy(n)); vy = 0.0_WP
  allocate(vz(n)); vz = 0.0_WP
  allocate(cvx(n)); cvx = 0.0_WP
  allocate(cvy(n)); cvy = 0.0_WP
  allocate(cvz(n)); cvz = 0.0_WP
  allocate(nc(n)); nc = 0
  allocate(nbs(n,1)); nbs = 0
  allocate(fx(n)); fx = 0.0_WP
  allocate(fy(n)); fy = 0.0_WP
  allocate(fz(n)); fz = 0.0_WP
  allocate(dpx(n)); dpx = 0.0_WP
  allocate(dpy(n)); dpy = 0.0_WP
  allocate(dpz(n)); dpz = 0.0_WP
  allocate(vox(n)); vox = 0.0_WP
  allocate(voy(n)); voy = 0.0_WP
  allocate(voz(n)); voz = 0.0_WP

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

  ! Print out initial particle positions
  call particle_write

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


subroutine particle_write
  use posvel
  use state
  implicit none
  integer :: i, fid
  character(64) :: buffer, filename

  write(buffer,"(ES12.3)") time
  filename = "particles_"//trim(adjustl(buffer))
  open(fid,file=filename,status="REPLACE",action="WRITE")
  write(fid,*) n
  do i = 1, n
    write(fid,FMT="(6(ES22.15,1x))") px(i), py(i), pz(i), vx(i), vy(i), vz(i)
  end do

  return
end subroutine

! ===================================== !
! One step of SPH                       !
! ===================================== !
subroutine step
  use constants, only : ITERATIONS_PER_STEP, USE_VORTICITY_CONFINEMENT
  use posvel
  implicit none
  integer :: subiter

  ! Apply the forces
  call force_apply
!  print*,'after force', maxval(cvx), cpx(maxloc(cvx)), cpy(maxloc(cvx)),cpz(maxloc(cvx))
  ! Compute candidate velocities and positions
  call posvel_candidate
!  print*,'after posvel', maxval(cvx), cpx(maxloc(cvx)), cpy(maxloc(cvx)),cpz(maxloc(cvx))
  ! Find neighbors and compute lambda
  call neighbor_find
!  print*,'after neighbor', maxval(cvx), cpx(maxloc(cvx)), cpy(maxloc(cvx)),cpz(maxloc(cvx))
  do subiter = 1, ITERATIONS_PER_STEP
    ! Compute lamnda
    call compute_lambda
!  print*,'after lambda', maxval(cvx), maxval(cvy), maxval(cvz)
    ! Jiggle the particles
    call particle_jiggle
!  print*,'after jiggle', maxval(cvx), maxval(cvy), maxval(cvz)
  end do
  call update_candidate_vel

  ! Confine vorticity, update vel, and calculate force
  if(USE_VORTICITY_CONFINEMENT) then
!  print*,'vorticity stuff', maxval(cvx), cpx(maxloc(cvx)), cpy(maxloc(cvx)),cpz(maxloc(cvx))
    call compute_omega
    call vorticity_force
  end if
!  print*,'before visc', maxval(cvx), maxval(cvy), maxval(cvz)
  ! Update particle positions due to vorticity
  call viscosity_update
!  print*,'update with visc', maxval(cvx), maxval(cvy), maxval(cvz)

  ! Update particle position
  call pos_update
!  print*,'update pos', maxval(cvx), cpx(maxloc(cvx)), cpy(maxloc(cvx)),cpz(maxloc(cvx))

end subroutine step



subroutine force_apply
  use forces
  use constants
  implicit none

  ! Apply forces
  fx(:) = vox(:)
  fy(:) = -GRAVITY
  fy(:) = fy(:) + voy(:)
  fz(:) = voz(:)

  return
end subroutine force_apply

subroutine posvel_candidate
  use forces
  use posvel
  use state
  implicit none
  integer :: i
  real(WP) :: mag

  ! Update candidate velocity
  cvx(:) = vx(:) + fx(:) * dt
  cvy(:) = vy(:) + fy(:) * dt
  cvz(:) = vz(:) + fz(:) * dt

  if(USE_INITIAL_VELOCITY_DAMPING .and. time<1.0_WP) then
    do i = 1,n
      mag = sqrt(cvx(i)*cvx(i) + cvy(i)*cvy(i) + cvz(i)*cvz(i))
      if(mag .gt. 2.0_WP) then
        cvx(i) = cvx(i)/mag
        cvy(i) = cvy(i)/mag
        cvz(i) = cvz(i)/mag
      end if
    end do
  end if

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
  integer :: i, j, q, p
  integer, dimension(:,:,:),allocatable :: part_count
  integer, dimension(:,:,:,:), allocatable :: binpart
  integer :: nbinx, nbiny, nbinz, max_part_guess
  integer :: binx, biny, binz
  integer :: cbinx, cbiny, cbinz
  integer :: stx, sty, stz
  real(WP) :: y_height, distx, disty, distz, tdist

  deallocate(nbs)
  y_height = maxval(cpy(:))
  nbinx = ceiling(BOX_SIZE*2.0_WP/KERNEL_SIZE) - 1
!  nbiny = ceiling(y_height/KERNEL_SIZE) - 1
!  nbinz = ceiling(BOX_SIZE*2.0_WP/KERNEL_SIZE) - 1
nbiny = 1
nbinz = 1
!  max_part_guess = 10000*n/(nbinx*nbiny*nbinz)
  max_part_guess = n
  allocate(part_count(nbinx,nbiny,nbinz))
  allocate(binpart(max_part_guess,nbinx,nbiny,nbinz)); binpart = 0

  part_count = 0
!  print*, 'num of bins', nbinx, nbiny, nbinz
  do i = 1, n
    binx = NINT((cpx(i)-(-BOX_SIZE))/(2.0_WP*BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((cpy(i)-0.0_WP)/(y_height+1.0e-15_WP)*(real(nbiny,WP)-1.0_WP)) + 1
    binz = NINT((cpz(i)-(-BOX_SIZE))/(2.0_WP*BOX_SIZE)*(real(nbinz,WP)-1.0_WP)) + 1
    part_count(binx, biny, binz) = part_count(binx, biny, binz) + 1
    binpart(part_count(binx,biny,binz),binx,biny,binz) = i
  end do

  max_part_guess = 27*maxval(part_count(:,:,:))
  allocate(nbs(n,max_part_guess)); nbs = 0
  nc = 0
  do i = 1, n
    binx = NINT((cpx(i)-(-BOX_SIZE))/(2.0_WP*BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((cpy(i)-0.0_WP)/(y_height+1.0e-15_WP)*(real(nbiny,WP)-1.0_WP)) + 1
    binz = NINT((cpz(i)-(-BOX_SIZE))/(2.0_WP*BOX_SIZE)*(real(nbinz,WP)-1.0_WP)) + 1
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
            distx = cpx(i) - cpx(p)
            disty = cpy(i) - cpy(p)
            distz = cpz(i) - cpz(p)
            tdist = sqrt(distx*distx + disty*disty + distz*distz)
            if(tdist .le. KERNEL_SIZE .and. p.ne.i) then
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
    numer = rho / PARTICLE_DENSITY - 1.0_WP

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
      if(USE_SURFACE_TENSION) then
      scorr = -ARTIFICIAL_PRESSURE_STRENGTH &
              * (poly6kernel(i,nb) * PRESSURE_RADIUS_FACTOR) &
              ** ARTIFICIAL_PRESSURE_POWER
      end if

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
!  print*,"before update", maxval(cvx), cpx(maxloc(cvx)), cpy(maxloc(cvx)),cpz(maxloc(cvx))
    ! Update candidate position
    cpx(:) = cpx(:) + dpx(:)
    cpy(:) = cpy(:) + dpy(:)
    cpz(:) = cpz(:) + dpz(:)

!  print*,"before collision", maxval(cvx), cpx(maxloc(cvx)), cpy(maxloc(cvx)),cpz(maxloc(cvx))

  ! Resolve collisions
  do i = 1, n

    ! Floor
    if(cpy(i) .lt. 0.0_WP) then
      cpy(i) = 0.0_WP
    end if

    ! Left Wall
    if(cpx(i) .lt. -BOX_SIZE) then
      cpx(i) = -BOX_SIZE
    end if

    ! Right Wall
    if(cpx(i) .gt. BOX_SIZE) then
      cpx(i) = BOX_SIZE
    end if

    ! Back Wall
    if(cpz(i) .lt. -BOX_SIZE) then
      cpz(i) = -BOX_SIZE
    end if

    ! Front Wall
    if(cpz(i) .gt. BOX_SIZE) then
      cpz(i) = BOX_SIZE
    end if

  end do

  return
end subroutine particle_jiggle

subroutine update_candidate_vel
  use posvel
  use state
  use constants
  implicit none
  integer :: i

  cvx = (cpx - px)/ dt
  cvy = (cpy - py)/ dt
  cvz = (cpz - pz)/ dt

  ! Resolve collisions
  do i = 1, n

    ! Floor
    if(cpy(i) .lt. 0.0_WP) then
      cvy(i) = -cvy(i)
    end if

    ! Left Wall
    if(cpx(i) .lt. -BOX_SIZE) then
      cvx(i) = -cvx(i)
    end if

    ! Right Wall
    if(cpx(i) .gt. BOX_SIZE) then
      cvx(i) = -cvx(i)
    end if

    ! Back Wall
    if(cpz(i) .lt. -BOX_SIZE) then
      cvz(i) = -cvz(i)
    end if

    ! Front Wall
    if(cpz(i) .gt. BOX_SIZE) then
      cvz(i) = -cvz(i)
    end if
  end do

  return
end subroutine update_candidate_vel

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
  real(WP) :: dvx, dvy, dvz, k
  real(WP) :: r

  vx = cvx
  vy = cvy
  vz = cvz

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

!    print*,i, vx(i), cpx(i), cpy(i), cpz(i)

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      dvx = vx(nb) - vx(i)
      dvy = vy(nb) - vy(i)
      dvz = vz(nb) - vz(i)

      k = poly6visckernel(i,nb)
      r = (cpx(i)-cpx(nb))**2+(cpy(i)-cpy(nb))**2+(cpz(i)-cpz(nb))**2
      r = sqrt(r)
!      print*,r-KERNEL_SIZE,k, i, nb, k*dvx, dvx, vx(nb), vx(i), px(nb),px(i)
      sx = sx + k * dvx
      sy = sy + k * dvy
      sz = sz + k * dvz
    end do

!    print*,'damn', i, sx, sy, sz
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




