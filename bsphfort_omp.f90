! ===================================== !
! Specifies what double precision is    !
! ===================================== !
module precision

  implicit none
  integer, parameter :: WP = 8
end module

! ===================================== !
! Timing and number of particles        !
! ===================================== !
module state
  use precision
  implicit none
  real(WP) :: time, dt, tfin
  integer :: n ! Number of particles
  integer :: procs
  integer :: iter, frame
  real(WP) :: frame_freq
  real(WP) :: mass
  INTEGER :: my_id
end module

! ===================================== !
! Constant values needed for SPH        !
! ===================================== !
module constants
  use precision
  implicit none
  real(WP), parameter :: GRAVITY = -9.8_WP
  real(WP), parameter :: PI = 3.14159265_WP
  real(WP) :: BOX_SIZE
  real(WP), parameter :: PARTICLE_DENSITY = 1000.0_WP
  real(WP) :: h, h2, h3, h4, h5, h9       ! KERNEL_SIZE
  real(WP), parameter :: BULK_MODULUS = 1000.0_WP
  real(WP), parameter :: E = 0.75  ! Resistution Coefficient
  real(WP), parameter :: VISCOSITY = 1.0e-2_WP
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
  !DIR$ attributes align:64 :: px, py, pz
  !DIR$ attributes align:64 :: vx, vy, vz
  !DIR$ attributes align:64 :: vlx, vly, vlz
  !DIR$ attributes align:64 :: rho
  !DIR$ ASSUME_ALIGNED px: 64
  !DIR$ ASSUME_ALIGNED py: 64
  !DIR$ ASSUME_ALIGNED pz: 64
  !DIR$ ASSUME_ALIGNED vx: 64
  !DIR$ ASSUME_ALIGNED vy: 64
  !DIR$ ASSUME_ALIGNED vz: 64
  !DIR$ ASSUME_ALIGNED vlx: 64
  !DIR$ ASSUME_ALIGNED vly: 64
  !DIR$ ASSUME_ALIGNED vlz: 64
  !DIR$ ASSUME_ALIGNED rho: 64

end module

! ===================================== !
! Arrays related to quickly identifying !
! nearest neighbors                     !
! ===================================== !
module neighbors
  use precision
  use constants
  implicit none
  integer,dimension(:,:,:,:),allocatable :: binpart
  integer,dimension(:,:,:),allocatable :: part_count
  integer,dimension(:,:),allocatable :: nbs
  integer,dimension(:),allocatable :: nc

  ! Expected max particles in single bin
  integer :: nbinx, nbiny, nbinz
  integer :: max_part_guess
  !DIR$ attributes align:64 :: nc, nbs, part_count, binpart
  !DIR$ ASSUME_ALIGNED nc: 64
  !DIR$ ASSUME_ALIGNED nbs: 64
  !DIR$ ASSUME_ALIGNED part_count: 64
  !DIR$ ASSUME_ALIGNED binpart: 64
end module neighbors

! ===================================== !
! Arrays related to forces              !
! ===================================== !
module forces
  use precision
  implicit none
  real(WP),dimension(:),allocatable :: fx, fy, fz  ! Total force
  !DIR$ attributes align:64 :: fx, fy, fz
  !DIR$ ASSUME_ALIGNED fx: 64
  !DIR$ ASSUME_ALIGNED fy: 64
  !DIR$ ASSUME_ALIGNED fz: 64
end module forces

! ===================================== !
! Timers for diagnostics                !
! ===================================== !
module timers
  use precision
  implicit none
  real(WP) :: neigh_t, dens_t, force_t, ef_t, ps_t, tot_t
end module timers

! ===================================== !
! Run the simulation                    !
! ===================================== !
program sph_run
  use state
  use posvel
  use timers
  use omp_lib
  implicit none
  real(WP) :: init_start, init_stop
  real(WP) :: write_start, write_stop, write_sum
  real(WP) :: deinit_start, deinit_stop
  real(WP) :: total_time
  real(WP) :: twrite


  init_start = omp_get_wtime()
  time = 0.0_WP
  dt = 0.0_WP

  call init
  twrite = frame_freq
  init_stop = omp_get_wtime()

  !$OMP PARALLEL NUM_THREADS(procs) default(shared)
  iter = 0
  time = 0.0_WP
  write_sum = 0.0_WP
  do while (time .lt. tfin)
    !$OMP MASTER
    iter = iter + 1
    !$OMP END MASTER

    ! Simulate the fluid flow
    call step
    !$OMP MASTER
    time = time + dt
    !$OMP END MASTER

    !$OMP MASTER
      if(time .ge. twrite) then
        write_start = omp_get_wtime()
        twrite = twrite + frame_freq
        print*,"Iteration: ",iter," Time: ",time
        print*,"Max U: ", maxval(vx(:))," Max V: ",maxval(vy(:))," Max W: ",maxval(vz(:))
        print*,"  "
        call particle_write
        write_stop = omp_get_wtime()
        write_sum = write_sum + (write_stop-write_start)
      end if
    !$OMP END MASTER
    !$OMP BARRIER
  end do
  !$OMP END PARALLEL

  deinit_start = omp_get_wtime()
  call deinit
  deinit_stop = omp_get_wtime()



  total_time = deinit_stop-init_start
  write(*,"(A)")                 "Timing Information:                Total Time              Time Per Step"
  write(*,"(A)")                 "---------------------             ------------            ---------------"
  write(*,"(A,ES25.11)")         "Initialization:        ",init_stop-init_start
  write(*,"(A,ES25.11,ES25.11)") "Neighbor Calculation:  ",neigh_t,neigh_t/iter
  write(*,"(A,ES25.11,ES25.11)") "Density Calculation:   ",dens_t,dens_t/iter
  write(*,"(A,ES25.11,ES25.11)") "Force Calculation:     ",force_t,force_t/iter
  write(*,"(A,ES25.11,ES25.11)") "Ext. Force Calculation:",ef_t,ef_t/iter
  write(*,"(A,ES25.11,ES25.11)") "Posvel Calculation:    ",ps_t,ps_t/iter
  write(*,"(A,ES25.11,ES25.11)") "File Writing:          ",write_sum,write_sum/iter
  write(*,"(A,ES25.11)")         "Deinitialization Time: ",deinit_stop-deinit_start
  write(*,"(A,ES25.11,ES25.11)") "Total Time:            ",total_time,total_time/iter

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
  use timers
  implicit none
  character(30) :: init_file, proc_string
  integer :: i, unit
  real, dimension(:), allocatable :: data

  ! Read in command line argument init_file
  call get_command_argument(1,init_file)
  call get_command_argument(2,proc_string)
  if(len_trim(proc_string) .eq. 0) then
    print*, "Incorrect number of command line arguments"
    print*, "Correct usage is ./bsphfort (init_file) (procs)"
    print*, "Exiting..."
    stop
  end if

  read (proc_string, *) procs

  ! Find length of initial condition file (number of particles)
  unit = 20
  open(unit,file=trim(init_file),action="read")
  read(unit,*) n            ! First line will be number of particles
  read(unit,*) h            ! First line will be number of particles
  h2=h*h; h4=h2*h2; h3=h2*h; h5=h3*h2; h9=h5*h3*h
  read(unit,*) tfin         ! Second line is finish time of simulation
  read(unit,*) dt           ! Third line is time step
  read(unit,*) frame_freq   ! Fourth line is viz frame frequency
  read(unit,*) BOX_SIZE   ! Fifth line is BOX_SIZE

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


  ! Neighbor finding size creations, each bin is ~h
    nbinx = floor(BOX_SIZE/h) - 1
    nbiny = floor(BOX_SIZE/h) - 1
    nbinz = floor(BOX_SIZE/h) - 1

 ! Need to find good way to estimate this so array not needlessly large
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

  ! Print out initial particle positions
  frame = 0
  call particle_write

  ! Initialize timers
  neigh_t = 0.0_WP
  dens_t = 0.0_WP
  force_t = 0.0_WP
  ef_t = 0.0_WP
  ps_t = 0.0_WP
  tot_t = 0.0_WP

  call mass_change

  return
end subroutine init

subroutine mass_change
  use posvel
  use state
  use constants
  use neighbors
  use forces
  implicit none
  real(WP) :: rho0, rho2s, rhos
  integer :: i, j, q, p, nb, l
  integer :: binx, biny, binz, tbins, bin
  integer :: cbinx, cbiny, cbinz
  integer :: stx, sty, stz
  real(WP) :: distx, disty, distz, tdist2
  real(WP) :: rhosum, c, r2
  real(WP) :: dx, dy, dz, kpd

  mass = 1

  do i = 1, n
    binx = NINT((px(i))/(BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((py(i))/(BOX_SIZE)*(real(nbiny,WP)-1.0_WP)) + 1
    binz = NINT((pz(i))/(BOX_SIZE)*(real(nbinz,WP)-1.0_WP)) + 1
    part_count(binx, biny, binz) = part_count(binx, biny, binz) + 1
    binpart(part_count(binx,biny,binz),binx,biny,binz) = i
  end do

  tbins = nbinx*nbiny*nbinz
  do bin = 0, tbins-1
    binz = bin / nbinx / nbiny + 1
    biny = mod(bin / nbinx, nbinz) + 1
    binx = mod(bin,nbinx) + 1
    do i = 1, part_count(binx,biny,binz)
      p = binpart(i,binx,biny,binz)
      q = 0
      do stz = -1, 1
        cbinz = binz + stz
        if(cbinz.lt.1 .or. cbinz.gt.nbinz) cycle
        do sty = -1, 1
          cbiny = biny + sty
          if(cbiny.lt.1 .or. cbiny.gt.nbiny) cycle
          do stx = -1, 1
            cbinx = binx + stx
            if(cbinx.lt.1 .or. cbinx.gt.nbinx) cycle
            do j = 1, part_count(cbinx,cbiny, cbinz)
              nb = binpart(j,cbinx,cbiny,cbinz)
              distx = px(p) - px(nb)
              disty = py(p) - py(nb)
              distz = pz(p) - pz(nb)
              tdist2 = distx*distx + disty*disty + distz*distz
              if(tdist2 .lt. h2 .and. nb.ne.p) then
                q = q + 1
                nbs(p,q) = nb
              end if
            end do
          end do
        end do
      end do
      nc(p) = q
    end do
  end do

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
  deallocate(fx)
  deallocate(fy)
  deallocate(fz)
  deallocate(rho)
  deallocate(nc)
  deallocate(nbs)
  deallocate(part_count)
  deallocate(binpart)

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
  use constants
  use state
  use timers
  use omp_lib
  implicit none
  real(WP) :: t1, t2, t3, t4, t5, t6, t7


  t1 = omp_get_wtime()
  call neighbor_find
  t2 = omp_get_wtime()
  call compute_density
  t3 = omp_get_wtime()
  call compute_forces
  t5 = omp_get_wtime()
  call ext_forces
  t6 = omp_get_wtime()
  call posvel_update
  t7 = omp_get_wtime()

  neigh_t = neigh_t  + (t2-t1)/procs
  dens_t  = dens_t   + (t3-t2)/procs
  force_t = force_t  + (t5-t3)/procs
  ef_t    = ef_t     + (t6-t5)/procs
  ps_t    = ps_t     + (t7-t6)/procs

end subroutine step


subroutine neighbor_find
  use state
  use neighbors
  use posvel
  use omp_lib
  implicit none
  integer :: i, j, q, p, nb, l
  integer :: binx, biny, binz, tbins, bin
  integer :: cbinx, cbiny, cbinz
  integer :: mnbinx, mnbiny, mnbinz
  integer :: stx, sty, stz
  integer :: step
  real(WP) :: distx, disty, distz, tdist2

  !$OMP WORKSHARE
  part_count = 0
  binpart = 0
  nbs = 0
  nc = 0
  !$OMP END WORKSHARE

  !$OMP SINGLE
  do i = 1, n
    binx = NINT((px(i))/(BOX_SIZE)*(real(nbinx,WP)-1.0_WP)) + 1
    biny = NINT((py(i))/(BOX_SIZE)*(real(nbiny,WP)-1.0_WP)) + 1
    binz = NINT((pz(i))/(BOX_SIZE)*(real(nbinz,WP)-1.0_WP)) + 1
    part_count(binx, biny, binz) = part_count(binx, biny, binz) + 1
    binpart(part_count(binx,biny,binz),binx,biny,binz) = i
  end do
  !$OMP END SINGLE

  tbins = nbinx*nbiny*nbinz
  step = tbins/procs/100
  !$OMP DO private(binz,biny,binx,p,nb,cbinx,cbiny,cbinz,distx,disty,distz,tdist2) &
  !$OMP SCHEDULE(DYNAMIC,step)
  do bin = 0, tbins-1
    binz = bin / nbinx / nbiny + 1
    biny = mod(bin / nbinx, nbinz) + 1
    binx = mod(bin,nbinx) + 1
    do i = 1, part_count(binx,biny,binz)
      p = binpart(i,binx,biny,binz)
      q = 0
      do stz = -1, 1
        cbinz = binz + stz
        if(cbinz.lt.1 .or. cbinz.gt.nbinz) cycle
        do sty = -1, 1
          cbiny = biny + sty
          if(cbiny.lt.1 .or. cbiny.gt.nbiny) cycle
          do stx = -1, 1
            cbinx = binx + stx
            if(cbinx.lt.1 .or. cbinx.gt.nbinx) cycle
            do j = 1, part_count(cbinx,cbiny, cbinz)
              nb = binpart(j,cbinx,cbiny,cbinz)
              distx = px(p) - px(nb)
              disty = py(p) - py(nb)
              distz = pz(p) - pz(nb)
              tdist2 = distx*distx + disty*disty + distz*distz
              if(tdist2 .lt. h2 .and. nb.ne.p) then
                q = q + 1
                nbs(p,q) = nb
              end if
            end do
          end do
        end do
      end do
      nc(p) = q
    end do
  end do
  !$OMP END DO
  return
end subroutine neighbor_find

subroutine compute_density
  use state
  use constants
  use neighbors
  use posvel
  implicit none
  integer :: i, j, l, nb
  real(WP) :: rhosum, c, r2
  real(WP) :: dx, dy, dz, kpd

  !$OMP WORKSHARE
  rho = (315.0_WP/64.0_WP/PI)*mass/h3
  !$OMP END WORKSHARE
  c = 315.0_WP/64.0_WP/PI*mass/h9

  !$OMP DO private(rhosum,nb,dx,dy,dz,r2,kpd)
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
  !$OMP END DO

  return
end subroutine compute_density

subroutine compute_forces
  use constants
  use neighbors
  use posvel
  use forces
  use state
  implicit none
  integer :: i, j, l, nb
  real(WP) :: sx, sy, sz
  real(WP) :: dx, dy, dz
  real(WP) :: dvx, dvy, dvz
  real(WP) :: kpres, kvisc
  real(WP) :: cpres, cvisc, r2, m, q

  cpres = 45.0_WP*mass/PI/h5*BULK_MODULUS*0.5_WP
  cvisc = -45.0_WP*mass/PI/h5*VISCOSITY

  !$OMP DO private(sx,sy,sz,nb,dx,dy,dz,dvx,dvy,dvz,r2,m,q,kpres,kvisc)
  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    ! Compute pressure
    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      dx = px(i) - px(nb)
      dy = py(i) - py(nb)
      dz = pz(i) - pz(nb)
      dvx = vx(i) - vx(nb)
      dvy = vy(i) - vy(nb)
      dvz = vz(i) - vz(nb)

      r2 = dx*dx + dy*dy + dz*dz
      m = sqrt(r2/h2)
      q = 1.0_WP - m

      kpres = cpres / rho(nb) * q &
          * (rho(i)+rho(nb)-2.0_WP*PARTICLE_DENSITY) &
          * q / m

      kvisc = cvisc / rho(nb) * q

      sx = sx + kpres*dx+kvisc*dvx
      sy = sy + kpres*dy+kvisc*dvy
      sz = sz + kpres*dz+kvisc*dvz

    end do
    fx(i) = sx
    fy(i) = sy
    fz(i) = sz

  end do
  !$OMP END DO

  return
end subroutine compute_forces

subroutine ext_forces
  use forces
  use constants
  use posvel
  implicit none

  !$OMP WORKSHARE
  fx = fx/rho
  fy = fy/rho + GRAVITY
  fz = fz/rho
  !$OMP END WORKSHARE

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

    !$OMP WORKSHARE
    vlx = vx
    vly = vy
    vlz = vz
    !$OMP END WORKSHARE

    !$OMP WORKSHARE
    vlx = vlx + 0.5_WP*dt * fx
    vly = vly + 0.5_WP*dt * fy
    vlz = vlz + 0.5_WP*dt * fz
    !$OMP END WORKSHARE

    !$OMP WORKSHARE
    px = px + dt*vlx
    py = py + dt*vly
    pz = pz + dt*vlz

    vx = vx + dt*fx
    vy = vy + dt*fy
    vz = vz + dt*fz
    !$OMP END WORKSHARE

  else

    !$OMP WORKSHARE
    vlx = vlx + dt * fx
    vly = vly + dt * fy
    vlz = vlz + dt * fz

    px = px + dt*vlx
    py = py + dt*vly
    pz = pz + dt*vlz
    !$OMP END WORKSHARE

    !$OMP WORKSHARE
    vx = vlx
    vy = vly
    vz = vlz
    !$OMP END WORKSHARE

    !$OMP WORKSHARE
    vx = vx + 0.5_WP*dt*fx
    vy = vy + 0.5_WP*dt*fy
    vz = vz + 0.5_WP*dt*fz
    !$OMP END WORKSHARE
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

  !$OMP DO
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
  !$OMP END DO

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




