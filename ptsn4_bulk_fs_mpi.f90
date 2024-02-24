PROGRAM hamiltonian
  IMPLICIT NONE
  include 'mpif.h'

  INTEGER :: i, j, k, l, it, base, core, job ,mpi_err, mpi_size, num_procs, mpi_rank
  Real*8, Parameter :: pi = ACos(-1.D0)
  Integer, Dimension(:), Allocatable :: rx, ry, rz, orbx, orby
  Real*8, Dimension(:), Allocatable :: k1, k2,k3, hopr, hopi
  Complex*16, Dimension(:,:,:), Allocatable :: hamil
  Complex*16, Parameter :: iota = (0.D0, 1.D0)
  Integer, Parameter :: norb = 84 ! number of orbitals
  Integer, Parameter :: num = 882000!7490880 ! number of lines in hr file
  Integer, Parameter :: kkm = 110
  real*8, Parameter :: alpha = 1.D0!0.9789
  real*8, Parameter :: beta = 0.5!0.553
  real*8, Parameter :: gama = 1.D0!0.9842

  integer, parameter :: ns = (kkm+1)**2
  integer, allocatable :: pardis(:,:)
  real*8, dimension(:,:), allocatable :: res, res_red

  Integer,Parameter :: nmax=norb
  Integer,Parameter :: Nm=nmax
  Integer,Parameter :: LDA=Nm
  Integer :: INFO
  Integer,Parameter ::LWORK=2*Nm-1
  Real*8 :: RWORK(3*Nm-2)
  Complex*16 :: Bk(ns, LDA,Nm)
  Complex*16 :: WORK(LWORK)
  real*8 :: Wk(ns, Nm)
  External :: ZHEEV

  ! Open the data files
  Open(12, File="harald_hr.dat", Status="Unknown") 
  Open(13, File="kmesh_file.dat", Status="Unknown")
  !Open(13, File="k-data.dat", Status="Unknown")
  Open(15, File="fermi_surface_0kx_mpi.dat", Status="Unknown")
  
  
  Allocate(rx(num))
  Allocate(ry(num))
  Allocate(rz(num))
  Allocate(k1(ns))
  Allocate(k2(ns))
  Allocate(k3(ns))  
  Allocate(orbx(num))
  Allocate(orby(num))
  Allocate(hopr(num))
  Allocate(hopi(num))
  Allocate(hamil(ns, norb, norb))
  allocate(res(ns,Nm))
  allocate(res_red(ns,Nm))
  
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
print*, '1. reading data file'
 !------------------- Read the data from the files------------------------!
   Do i = 1, num
     Read(12, *) rx(i), ry(i), rz(i), orbx(i), orby(i), hopr(i), hopi(i)
   End Do
   
   Do i = 1, ns
     Read(13, *) k1(i), k2(i)
   End Do

   ! Close the files
   Close(12)
   Close(13)

  !---------------------------starting mpi-----------------------------!
  CALL MPI_Init(mpi_err)
  CALL MPI_Comm_size(MPI_COMM_WORLD,mpi_size,mpi_err)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

  mpi_rank=mpi_rank+1

  ALLOCATE(pardis(mpi_size,ns))
  base=ns/mpi_size

  core=1
  job=1
  do i=1,mpi_size
    do k=1,ns
        pardis(i,k)=0
    enddo
  enddo



  do i=1,ns
    pardis(core,job)=i
    job=job+1

    if (job.gt.base) then
      core=core+1
      job=1
    endif
  enddo
  if (mpi_rank == 1) then
    write(*,*) "1. reading data files ..."
  endif

!----------------- Fourier transformation of hr file (parallel do)------------------------! 

  if (mpi_rank == 1) then
    write(*,*) "3. starting Fourier transformation ..."
  endif
  hamil = (0.d0, 0.d0)
  it=1
  do while ((pardis(mpi_rank,it).gt.0).and.(it.le.ns))   
    j = pardis(mpi_rank,it)

    Do i = 1, num
      hamil(j, orbx(i), orby(i)) = hamil(j, orbx(i), orby(i)) + &
      (hopr(i) + iota*hopi(i)) * exp(iota*((alpha*(0.142)*pi - beta*k1(j))* rx(i) + &
      (alpha*(0.142)*pi + beta*k1(j)) * ry(i) +  gama*k2(j)*rz(i))) 
    End Do
    Bk(j,:,:) = hamil(j,:,:)
    !-------------diagonalising the hamiltonian---------------!
    Call ZHEEV('V','U',Nm,Bk(j,:,:),LDA,Wk(j,:),WORK,LWORK,RWORK,INFO)

    If( INFO.Gt.0 ) Then
    Write(*,*)'The algorithm failed to compute eigenvalues.'
    Stop
    End If

    res(j,:) = Wk(j,:)
    it=it+1
  end do
  
  !-------------------------------- Ending MPI---------------------------------!
  call mpi_allreduce(res, res_red, norb*ns,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, mpi_err)
  CALL MPI_Barrier(MPI_COMM_WORLD,mpi_err)

  !write(*,*) "forking data from different cores to a sing file ..."


  if (mpi_rank == 1) then
    write(*,*) "4. full Fourier transformation is done!"
    do j = 1, ns
      ! write(15,*)  -1/pi*aimag(res_red(i)) 
      write(15,*) k1(j),k2(j),res_red(j,:)
    enddo 
    write(*,*) "5. data generated successfully!"
  end if   
  
  !deallocate(pardis)
  call MPI_Finalize(mpi_err)
  close(15)
    
  print*, 'data generated successfully!'
   
   


END PROGRAM hamiltonian
