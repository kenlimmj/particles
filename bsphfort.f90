! ===================================== !
! Specifies what double precision is    !
! ===================================== !
module precision

  implicit none
  integer, parameter :: WP = 8
end module

module parallel
  use precision
  use mpi
  integer :: irank,ierr,nproc,nprocx,nprocy,nprocz
  real(WP),dimension(:),allocatable :: sendbuf
  real(WP),dimension(:,:),allocatable :: recvbuf
  integer,dimension(:),allocatable :: sendbuf2
  integer,dimension(:),allocatable :: recvbuf2

end module
! ===================================== !
! Timing and number of particles        !
! ===================================== !
module state
    use parallel
  use precision
  implicit none
  real(WP) :: time, dt, tfin
  integer :: n ! Number of particles
  integer :: n_ ! Number of particles
  integer :: iter, frame
  real(WP) :: frame_freq
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
  real(WP),dimension(:),allocatable,target ::  px,  py,  pz  ! Current position
  real(WP),dimension(:),allocatable,target ::  vx,  vy,  vz  ! Velocity
  real(WP),dimension(:),allocatable,target ::  vlx,  vly,  vlz  ! Velocity half time step behind
  real(WP),dimension(:),allocatable ::  rho           ! Density
  integer ,dimension(:),allocatable,target ::  ghost  
  integer                                  ::  vars=9
  type ptr
    real(WP), pointer :: p(:)  
  end type ptr
  type(ptr),dimension(9)                ::  dat
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

  real(WP) :: neigh_t, dens_t, force_t, ef_t, ps_t,comm_t, tot_t

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

  call MPI_Init(ierr)
  call MPI_COMM_RANK( MPI_COMM_WORLD, irank, ierr )

  init_start = omp_get_wtime()
  time = 0.0_WP
  dt = 0.0_WP

  call init
  twrite = frame_freq
  init_stop = omp_get_wtime()

  iter = 0
  do while (time .lt. tfin)
    iter = iter + 1

    ! Simulate the fluid flow
    call step
    time = time + dt

    if(time .ge. twrite) then
      write_start = omp_get_wtime()
      twrite = twrite + frame_freq
      !print*,"Iteration: ",iter," Time: ",time
      !print*,"Max U: ", maxval(vx(:))," Max V: ",maxval(vy(:))," Max W: ",maxval(vz(:))
      !print*,"  "
      call particle_write
      write_stop = omp_get_wtime()
      write_sum = write_sum + (write_stop-write_start)
    end if

  end do

  deinit_start = omp_get_wtime()
  call deinit
  deinit_stop = omp_get_wtime()


  total_time = deinit_stop-init_start

  if (irank.eq.0)then
    write(*,"(A)")                 "Timing Information:                Total Time              Time Per Step"
    write(*,"(A)")                 "---------------------             ------------            ---------------"
    write(*,"(A,ES25.11)")         "Initialization:        ",init_stop-init_start
    write(*,"(A,ES25.11,ES25.11)") "Neighbor Calculation:  ",neigh_t,neigh_t/iter
    write(*,"(A,ES25.11,ES25.11)") "Density Calculation:   ",dens_t,dens_t/iter
    write(*,"(A,ES25.11,ES25.11)") "Force Calculation:     ",force_t,force_t/iter
    write(*,"(A,ES25.11,ES25.11)") "Ext. Force Calculation:",ef_t,ef_t/iter
    write(*,"(A,ES25.11,ES25.11)") "Posvel Calculation:    ",ps_t,ps_t/iter
    write(*,"(A,ES25.11,ES25.11)") "Communcation:          ",comm_t,comm_t/iter
    write(*,"(A,ES25.11)")         "Deinitialization Time: ",deinit_stop-deinit_start
    write(*,"(A,ES25.11,ES25.11)") "Total Time:            ",total_time,total_time/iter
  endif
  call MPI_FINALIZE(ierr)

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
  use parallel
  implicit none
  character(30) :: init_file
  integer :: i, unit
  real, dimension(:), allocatable :: data
  character(64) :: buffer

  ! Read in command line argument init_file
  call get_command_argument(1,init_file)
  if(len_trim(init_file) .eq. 0) then
    !print*, "Incorrect number of command line arguments"
    !print*, "Correct usage is ./bsphfort (init_file)"
    !print*, "Exiting..."
    stop
  end if


  ! Find length of initial condition file (number of particles)
  unit = 20+irank
  write(buffer,"(a, i5.5)") trim(init_file),irank
  open(unit,file=trim(adjustl(buffer)),action="read")
  read(unit,*) n_            ! First line will be number of particles
  read(unit,*) n            ! First line will be number of particles
  read(unit,*) h            ! First line will be number of particles
  h2=h*h; h3=h2*h; h5=h3*h2; h9=h5*h3*h
  read(unit,*) tfin         ! Second line is finish time of simulation
  read(unit,*) dt           ! Third line is time step
  read(unit,*) frame_freq   ! Fourth line is viz frame frequency
  read(unit,*) nprocx
  read(unit,*) nprocy
  read(unit,*) nprocz
  read(unit,*) BOX_SIZE   ! Fifth line is BOX_SIZE
  nproc=nprocx*nprocy*nprocz

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
  allocate(ghost(n)); ghost = 0

  ! Neighbor finding size creations, each bin is ~2h
  nbinx = floor(BOX_SIZE/h) - 1
  nbiny = floor(BOX_SIZE/h) - 1
  nbinz = floor(BOX_SIZE/h) - 1

 ! Need to find good way to estimate this so array not needlessly large
  max_part_guess = 500*n_/(nbinx*nbiny*nbinz)
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

  dat(1)%p  => px
  dat(2)%p  => py
  dat(3)%p  => pz
  dat(4)%p  => vx
  dat(5)%p  => vy
  dat(6)%p  => vz
  dat(7)%p  => vx
  dat(8)%p  => vy
  dat(9)%p  => vz

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
  use parallel
  implicit none
  integer :: i
  real(WP) :: rho0, rho2s, rhos, rhos_,rho2s_

  mass = 1
  call communicate
  call neighbor_find
  call compute_density
  rho0 = PARTICLE_DENSITY
  rho2s = 0.0_WP
  rhos = 0.0_WP

  do i = 1, n
    if (ghost(i).eq.0) rho2s = rho2s + rho(i)*rho(i)
    if (ghost(i).eq.0) rhos = rhos + rho(i)
  end do

  call MPI_ALLREDUCE(rhos,rhos_,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(rho2s,rho2s_,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  mass = mass*rho0*rhos_/rho2s_
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
  deallocate(ghost)

  return
end subroutine deinit


subroutine particle_write
  use posvel
  use state
  implicit none
  integer :: i, fid
  character(64) :: buffer, start,filename,catcommand

  write(buffer,"(a, i5.5)") 'buf',irank
  filename = trim(adjustl(buffer))
  fid = irank


  open(fid,file=filename,status="REPLACE",action="WRITE")
  do i = 1, n
    if (ghost(i).eq.0) write(fid,FMT="(6(ES22.15,1x))") px(i), py(i), pz(i), vx(i), vy(i), vz(i)
  end do
  close(fid)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (irank .eq. 0)then
    start = "particles_"
    write(buffer,"(a, i5.5)") trim(start), frame
    filename = trim(adjustl(buffer))
    
    
    open(nproc,file=filename,status="REPLACE",action="WRITE")
    write(nproc,*) n_
    close(nproc)
    
    
    do i = 0,nproc-1
      write(buffer,"(a,a, i5.5,a,a,i5.5)") 'cat ','buf',i,' >> ',trim(start),frame
      catcommand = trim(adjustl(buffer))
      call system(catcommand)
    enddo
  endif

  frame = frame + 1

  return
end subroutine

subroutine communicate
  use parallel
  use posvel
  use neighbors
  use state
  use forces
  implicit none
  integer :: i,j,k,nn,proc,proc_orig,proc_recv,n_new,maxsend,maxsend_,disp,v, total_count,sub_count
  integer,dimension(0:nproc-1) :: sendcount,recvcount,flag
  real(WP) :: epsx,epsy,epsz
  type(ptr),dimension(9) :: dat_new
  real(WP),dimension(:),allocatable,target :: px_new,py_new,pz_new,vx_new,vy_new,vz_new,vlx_new,vly_new,vlz_new
  integer,dimension(:),allocatable,target :: ghost_new
  integer,dimension(0:7) :: balance
  integer :: missing=0,fac=1

  !print*,irank, 'begmin    ',minval(px),minval(py),minval(pz),minval(vx),minval(vy),minval(vz)
  !print*,irank, 'begmax    ',maxval(px),maxval(py),maxval(pz),maxval(vx),maxval(vy),maxval(vz)

  balance=0
  sendcount=0
  n_new=0
  do nn = 1,n
    proc_orig = nprocx*nprocy*floor((pz(nn)/(BOX_SIZE/nprocz)))    &
              + nprocx*       floor((py(nn)/(BOX_SIZE/nprocy)))    &
              +               floor((px(nn)/(BOX_SIZE/nprocx)))    
    flag=0
    do i = -1,1
      do j = -1,1
        do k = -1,1
          epsx=fac*h*i
          epsy=fac*h*j
          epsz=fac*h*k

          
          proc = nprocx*nprocy*floor(((pz(nn)+epsx)/(BOX_SIZE/nprocz)))    &
               + nprocx*       floor(((py(nn)+epsy)/(BOX_SIZE/nprocy)))    &
               +               floor(((px(nn)+epsz)/(BOX_SIZE/nprocx)))    
          if ((proc>=0).and.(proc<nproc)) then
            if ((ghost(nn).eq.0).and.(flag(proc).eq.0))then
              balance(proc) = balance(proc)+1
              flag(proc)=1
              sendcount(proc)=sendcount(proc)+1
            endif
          endif
        enddo
      enddo
    enddo
!    if(sum(flag)>0 .and. ghost(nn).eq.0) print*,flag(-2) 
  enddo
  !print*,'anac', irank,n,n_,'bal',balance




  n_new=sendcount(irank)

  maxsend=maxval(sendcount)


  call MPI_ALLREDUCE(maxsend,maxsend_,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_ALLTOALL(sendcount(0),1,MPI_INTEGER,recvcount(0),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

  !if(irank.eq.0) print*,irank,'recv', recvcount

  !print*,irank,'send    ', sendcount

  allocate(sendbuf(nproc*maxsend_))
  allocate(recvbuf(nproc*maxsend_,1:vars))
  allocate(sendbuf2(nproc*maxsend_))
  allocate(recvbuf2(nproc*maxsend_))

  allocate(px_new   (sum(recvcount))); px_new   =0.0_WP
  allocate(py_new   (sum(recvcount))); py_new   =0.0_WP
  allocate(pz_new   (sum(recvcount))); pz_new   =0.0_WP
  allocate(vx_new   (sum(recvcount))); vx_new   =0.0_WP
  allocate(vy_new   (sum(recvcount))); vy_new   =0.0_WP
  allocate(vz_new   (sum(recvcount))); vz_new   =0.0_WP
  allocate(vlx_new   (sum(recvcount))); vx_new   =0.0_WP
  allocate(vly_new   (sum(recvcount))); vy_new   =0.0_WP
  allocate(vlz_new   (sum(recvcount))); vz_new   =0.0_WP
  allocate(ghost_new(sum(recvcount))); ghost_new=0
  
  dat_new(1)%p  => px_new
  dat_new(2)%p  => py_new
  dat_new(3)%p  => pz_new
  dat_new(4)%p  => vx_new
  dat_new(5)%p  => vy_new
  dat_new(6)%p  => vz_new
  dat_new(7)%p  => vlx_new
  dat_new(8)%p  => vly_new
  dat_new(9)%p  => vlz_new

  do v = 1,vars
    sendcount=0
    do nn = 1,n
      proc_orig = nprocx*nprocy*floor((pz(nn)/(BOX_SIZE/nprocz)))    &
                + nprocx*       floor((py(nn)/(BOX_SIZE/nprocy)))    &
                +               floor((px(nn)/(BOX_SIZE/nprocx)))    
      flag=0
      do i = -1,1
        do j = -1,1
          do k = -1,1
            epsx=fac*h*i
            epsy=fac*h*j
            epsz=fac*h*k
           

            proc = nprocx*nprocy*floor(((pz(nn)+epsx)/(BOX_SIZE/nprocz)))    &
                 + nprocx*       floor(((py(nn)+epsy)/(BOX_SIZE/nprocy)))    &
                 +               floor(((px(nn)+epsz)/(BOX_SIZE/nprocx)))    
            if ((proc>=0).and.(proc<nproc)) then
              if ((ghost(nn).eq.0).and.(flag(proc).eq.0))then
                flag(proc)=1
                sendcount(proc)=sendcount(proc)+1
                if(v.eq.1 .and. proc.eq.proc_orig)then
                    sendbuf2(proc*maxsend_+sendcount(proc)) = 0
                elseif(v.eq.1 .and. proc.ne.proc_orig)then
                    sendbuf2(proc*maxsend_+sendcount(proc)) = 1
                endif
                sendbuf(proc*maxsend_+sendcount(proc)) = dat(v)%p(nn)
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    if (size(recvbuf)>0)then
      call MPI_ALLTOALL(sendbuf,maxsend_,MPI_REAL8,recvbuf(1,v),maxsend_,MPI_REAL8,MPI_COMM_WORLD,ierr)
      call MPI_ALLTOALL(sendbuf2,maxsend_,MPI_INTEGER,recvbuf2,maxsend_,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  enddo

  

  n_new=0
  do i=0,nproc-1
      do j=1,recvcount(i)
        disp=(i)*maxsend_+j
        n_new=n_new+1
        
        do v = 1,vars
          dat_new(v)%p(n_new)=recvbuf(disp,v)
        enddo
        ghost_new(n_new)=recvbuf2(disp)
      enddo
  enddo


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
  deallocate(ghost)

  allocate(px   (size(px_new))); px   =px_new
  allocate(py   (size(px_new))); py   =py_new
  allocate(pz   (size(px_new))); pz   =pz_new
  allocate(vx   (size(px_new))); vx   =vx_new
  allocate(vy   (size(px_new))); vy   =vy_new
  allocate(vz   (size(px_new))); vz   =vz_new
  allocate(vlx  (size(px_new))); vlx  =vlx_new
  allocate(vly  (size(px_new))); vly  =vly_new
  allocate(vlz  (size(px_new))); vlz  =vlz_new
  allocate(rho  (size(px_new))); rho  =0.0_WP
  allocate(nc   (size(px_new))); nc   =0
  allocate(fx   (size(px_new))); fx   =0.0_WP
  allocate(fy   (size(px_new))); fy   =0.0_WP
  allocate(fz   (size(px_new))); fz   =0.0_WP
  allocate(ghost(size(px_new))); ghost=ghost_new

  allocate(nbs  (size(px_new),max_part_guess)); nbs = 0

  dat(1)%p  => px
  dat(2)%p  => py
  dat(3)%p  => pz
  dat(4)%p  => vx
  dat(5)%p  => vy
  dat(6)%p  => vz
  dat(7)%p  => vlx
  dat(8)%p  => vly
  dat(9)%p  => vlz

  deallocate(px_new)
  deallocate(py_new)
  deallocate(pz_new)
  deallocate(vx_new)
  deallocate(vy_new)
  deallocate(vz_new)
  deallocate(ghost_new)
  deallocate(sendbuf)
  deallocate(recvbuf)
  deallocate(sendbuf2)
  deallocate(recvbuf2)


  n=size(px)
  
!  sub_count=0
!  do i=1,n
!    if (ghost(i).eq.0) sub_count=sub_count+1
!  enddo
!  call MPI_ALLREDUCE(sub_count,total_count,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!  if (irank.eq.0) print*,sub_count,total_count,'  total',n_
!!  vx=5.0_WP
!  print*,sub_count,size(px), irank, 'min    ',minval(px),minval(py),minval(pz),minval(vx),minval(vy),minval(vz)
!  print*,sub_count,size(px), irank, 'max    ',maxval(px),maxval(py),maxval(pz),maxval(vx),maxval(vy),maxval(vz)
  return
end subroutine communicate


! ===================================== !
! One step of SPH                       !
! ===================================== !
subroutine step
  use constants
  use state
  use timers
  use omp_lib
  implicit none
  real(WP) :: t1, t2, t3, t4, t5, t6, t7, t8

  t1 = omp_get_wtime()
  if (n>0) call neighbor_find
  t2 = omp_get_wtime()
  call compute_density
  t3 = omp_get_wtime()
  call compute_forces
  t5 = omp_get_wtime()
  call ext_forces
  t6 = omp_get_wtime()
  call posvel_update
  t7 = omp_get_wtime()
  call communicate
  t8 = omp_get_wtime()

  neigh_t = neigh_t  + t2-t1
  dens_t  = dens_t   + t3-t2
  force_t = force_t  + t5-t3
  ef_t    = ef_t     + t6-t5
  ps_t    = ps_t     + t7-t6
  comm_t  = comm_t   + t8-t7

end subroutine step


subroutine neighbor_find
  use state
  use neighbors
  use posvel
  use omp_lib
  implicit none
  integer :: i, j, q, p, nb
  integer :: binx, biny, binz
  integer :: cbinx, cbiny, cbinz
  integer :: stx, sty, stz, housemates
  real(WP) :: distx, disty, distz, tdist2

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

  ! Left wall
  binx = 1
    do binz = 1, nbinz
      do biny = 1, nbiny
        do i = 1, part_count(binx,biny,binz)
          p = binpart(i,binx,biny,binz)
          binpart(i,binx,biny,binz) = -1
          do stz = -1, 1
            cbinz = binz + stz
            if(cbinz.lt.1 .or. cbinz.gt.nbinz) cycle
            do sty = -1, 1
              cbiny = biny + sty
              if(cbiny.lt.1 .or. cbiny.gt.nbiny) cycle
              do stx = 0, 1
                cbinx = binx + stx
                do j = 1, part_count(cbinx,cbiny, cbinz)
                  nb = binpart(j,cbinx,cbiny,cbinz)
                  if(nb .eq. -1) cycle
                  distx = px(p) - px(nb)
                  disty = py(p) - py(nb)
                  distz = pz(p) - pz(nb)
                  tdist2 = distx*distx + disty*disty + distz*distz
                  if(tdist2 .lt. h2) then
                    nc(p) = nc(p) + 1
                    nbs(p,nc(p)) = nb
                    nc(nb) = nc(nb) + 1
                    nbs(nb,nc(nb)) = p
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  ! Right wall
  binx = nbinx
    do binz = 1, nbinz
      do biny = 1, nbiny
        do i = 1, part_count(binx,biny,binz)
          p = binpart(i,binx,biny,binz)
          binpart(i,binx,biny,binz) = -1
          do stz = -1, 1
            cbinz = binz + stz
            if(cbinz.lt.1 .or. cbinz.gt.nbinz) cycle
            do sty = -1, 1
              cbiny = biny + sty
              if(cbiny.lt.1 .or. cbiny.gt.nbiny) cycle
              do stx = -1, 0
                cbinx = binx + stx
                do j = 1, part_count(cbinx,cbiny, cbinz)
                  nb = binpart(j,cbinx,cbiny,cbinz)
                  if(nb .eq. -1) cycle
                  distx = px(p) - px(nb)
                  disty = py(p) - py(nb)
                  distz = pz(p) - pz(nb)
                  tdist2 = distx*distx + disty*disty + distz*distz
                  if(tdist2 .lt. h2) then
                    nc(p) = nc(p) + 1
                    nbs(p,nc(p)) = nb
                    nc(nb) = nc(nb) + 1
                    nbs(nb,nc(nb)) = p
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  ! Bottom wall
  biny = 1
    do binz = 1, nbinz
      do binx = 2, nbinx-1
        do i = 1, part_count(binx,biny,binz)
          p = binpart(i,binx,biny,binz)
          binpart(i,binx,biny,binz) = -1
          do stz = -1, 1
            cbinz = binz + stz
            if(cbinz.lt.1 .or. cbinz.gt.nbinz) cycle
            do sty = 0, 1
              cbiny = biny + sty
              do stx = -1, 1
                cbinx = binx + stx
                do j = 1, part_count(cbinx,cbiny, cbinz)
                  nb = binpart(j,cbinx,cbiny,cbinz)
                  if(nb .eq. -1) cycle
                  distx = px(p) - px(nb)
                  disty = py(p) - py(nb)
                  distz = pz(p) - pz(nb)
                  tdist2 = distx*distx + disty*disty + distz*distz
                  if(tdist2 .lt. h2) then
                    nc(p) = nc(p) + 1
                    nbs(p,nc(p)) = nb
                    nc(nb) = nc(nb) + 1
                    nbs(nb,nc(nb)) = p
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
  end do

  ! Top wall
  biny = nbiny
    do binz = 1, nbinz
      do binx = 2, nbinx-1
        do i = 1, part_count(binx,biny,binz)
          p = binpart(i,binx,biny,binz)
          binpart(i,binx,biny,binz) = -1
          do stz = -1, 1
            cbinz = binz + stz
            if(cbinz.lt.1 .or. cbinz.gt.nbinz) cycle
            do sty = -1, 0
              cbiny = biny + sty
              do stx = -1, 1
                cbinx = binx + stx
                do j = 1, part_count(cbinx,cbiny, cbinz)
                  nb = binpart(j,cbinx,cbiny,cbinz)
                  if(nb .eq. -1) cycle
                  distx = px(p) - px(nb)
                  disty = py(p) - py(nb)
                  distz = pz(p) - pz(nb)
                  tdist2 = distx*distx + disty*disty + distz*distz
                  if(tdist2 .lt. h2) then
                    nc(p) = nc(p) + 1
                    nbs(p,nc(p)) = nb
                    nc(nb) = nc(nb) + 1
                    nbs(nb,nc(nb)) = p
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
  end do

  ! Forward wall
  binz = 1
    do biny = 2, nbiny-1
      do binx = 2, nbinx-1
        do i = 1, part_count(binx,biny,binz)
          p = binpart(i,binx,biny,binz)
          binpart(i,binx,biny,binz) = -1
          do stz = 0, 1
            cbinz = binz + stz
            do sty = -1, 1
              cbiny = biny + sty
              do stx = -1, 1
                cbinx = binx + stx
                do j = 1, part_count(cbinx,cbiny, cbinz)
                  nb = binpart(j,cbinx,cbiny,cbinz)
                  if(nb .eq. -1) cycle
                  distx = px(p) - px(nb)
                  disty = py(p) - py(nb)
                  distz = pz(p) - pz(nb)
                  tdist2 = distx*distx + disty*disty + distz*distz
                  if(tdist2 .lt. h2) then
                    nc(p) = nc(p) + 1
                    nbs(p,nc(p)) = nb
                    nc(nb) = nc(nb) + 1
                    nbs(nb,nc(nb)) = p
                  end if
                end do
              end do
            end do
          end do
        end do
    end do
  end do

  ! Back Wall
  binz = nbinz
    do biny = 2, nbiny-1
      do binx = 2, nbinx-1
        do i = 1, part_count(binx,biny,binz)
          p = binpart(i,binx,biny,binz)
          binpart(i,binx,biny,binz) = -1
          do stz = -1, 0
            cbinz = binz + stz
            do sty = -1, 1
              cbiny = biny + sty
              do stx = -1, 1
                cbinx = binx + stx
                do j = 1, part_count(cbinx,cbiny, cbinz)
                  nb = binpart(j,cbinx,cbiny,cbinz)
                  if(nb .eq. -1) cycle
                  distx = px(p) - px(nb)
                  disty = py(p) - py(nb)
                  distz = pz(p) - pz(nb)
                  tdist2 = distx*distx + disty*disty + distz*distz
                  if(tdist2 .lt. h2) then
                    nc(p) = nc(p) + 1
                    nbs(p,nc(p)) = nb
                    nc(nb) = nc(nb) + 1
                    nbs(nb,nc(nb)) = p
                  end if
                end do
              end do
            end do
          end do
        end do
    end do
  end do

  ! Rest
  do binz = 2, nbinz-1
    do biny = 2, nbiny-1
      do binx = 2, nbinx-1
        do i = 1, part_count(binx,biny,binz)
          p = binpart(i,binx,biny,binz)
          binpart(i,binx,biny,binz) = -1
          do stz = -1, 1
            cbinz = binz + stz
            do sty = -1, 1
              cbiny = biny + sty
              do stx = -1, 1
                cbinx = binx + stx
                do j = 1, part_count(cbinx,cbiny, cbinz)
                  nb = binpart(j,cbinx,cbiny,cbinz)
                  if(nb .eq. -1) cycle
                  distx = px(p) - px(nb)
                  disty = py(p) - py(nb)
                  distz = pz(p) - pz(nb)
                  tdist2 = distx*distx + disty*disty + distz*distz
                  if(tdist2 .lt. h2) then
                    nc(p) = nc(p) + 1
                    nbs(p,nc(p)) = nb
                    nc(nb) = nc(nb) + 1
                    nbs(nb,nc(nb)) = p
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

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

  rho = (315.0_WP/64.0_WP/PI)*mass/h3
  c = 315.0_WP/64.0_WP/PI*mass/h9
  do i = 1, n
    rhosum = 0.0_WP

    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      if(nb .le. i) cycle

      dx = px(i) - px(nb)
      dy = py(i) - py(nb)
      dz = pz(i) - pz(nb)
      r2 = dx*dx + dy*dy + dz*dz
      kpd = h2 - r2
      rhosum = rhosum + c * kpd * kpd * kpd
      rho(nb) = rho(nb) + c * kpd * kpd * kpd

    end do

  rho(i) = rho(i) + rhosum
  end do
  !print*,  'minr   ',minval(rho)
  !print*,  'maxr   ',maxval(rho)

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


  fx = 0.0_WP
  fy = 0.0_WP
  fz = 0.0_WP

  cpres = 45.0_WP*mass/PI/h5*BULK_MODULUS*0.5_WP
  cvisc = -45.0_WP*mass/PI/h5*VISCOSITY

  do i = 1, n
    sx = 0.0_WP
    sy = 0.0_WP
    sz = 0.0_WP

    ! Compute pressure
    l = nc(i)
    do j = 1, l
      nb = nbs(i,j)
      if(nb.le.i) cycle
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
      fx(nb) = fx(nb) - kpres*dx-kvisc*dvx
      sy = sy + kpres*dy+kvisc*dvy
      fy(nb) = fy(nb) - kpres*dy-kvisc*dvy
      sz = sz + kpres*dz+kvisc*dvz
      fz(nb) = fz(nb) - kpres*dz-kvisc*dvz

    end do
    fx(i) = fx(i)+sx
    fy(i) = fy(i)+sy
    fz(i) = fz(i)+sz


  end do

  return
end subroutine compute_forces

subroutine ext_forces
  use forces
  use constants
  use posvel
  use parallel
  implicit none

  fx = fx/rho
  fy = fy/rho
  fz = fz/rho

  ! Only gravity for now
  fy = fy + GRAVITY

!  if (irank.eq.0) print*,  'min    ',minval(fx),minval(fy),minval(fz)
!  if (irank.eq.0) print*,  'max    ',maxval(fx),maxval(fy),maxval(fz)
  return
end subroutine ext_forces

! Leapfrog time integration
! From Bindel 2014 SPH code
subroutine posvel_update
  use posvel
  use state
  use forces
  use constants
  implicit none
  integer :: i
  real(WP) :: epss=0.000001

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

  do i=1,n
    if (px(i)<0.0_WP) px(i)=0.0_WP+epss
    if (py(i)<0.0_WP) py(i)=0.0_WP+epss
    if (pz(i)<0.0_WP) pz(i)=0.0_WP+epss
    if (px(i)>BOX_SIZE) px(i)=BOX_SIZE-epss
    if (py(i)>BOX_SIZE) py(i)=BOX_SIZE-epss
    if (pz(i)>BOX_SIZE) pz(i)=BOX_SIZE-epss

  enddo
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




