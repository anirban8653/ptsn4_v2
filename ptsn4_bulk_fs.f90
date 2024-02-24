PROGRAM hamiltonian
  IMPLICIT NONE

  INTEGER :: i, j
  Real*8, Parameter :: pi = ACos(-1.D0)
  Integer, Dimension(:), Allocatable :: rx, ry, rz, orbx, orby
  Real*8, Dimension(:), Allocatable :: k1, k2,k3, hopr, hopi
  Complex*16, Dimension(:,:,:), Allocatable :: hamil
  Complex*16, Parameter :: iota = (0.D0, 1.D0)
  Integer, Parameter :: norb = 84 ! number of orbitals
  Integer, Parameter :: num = 882000!7490880 ! number of lines in hr file
  Integer, Parameter :: kkm = 50
  integer, parameter :: ns = (kkm+1)**2


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
  Open(15, File="fermi_surface_0kx.dat", Status="Unknown")
  
  
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
  
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
print*, '1. reading data file'
! Read the data from the files
   Do i = 1, num
     Read(12, *) rx(i), ry(i), rz(i), orbx(i), orby(i), hopr(i), hopi(i)
   End Do
   
   Do i = 1, ns
     Read(13, *) k1(i), k2(i)
   End Do

   ! Close the files
   Close(12)
   Close(13)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  hamil = (0.D0, 0.D0)
  print*, '3. starting fourier transformation'
  Do j = 1,  ns

    Do i = 1, num
      hamil(j, orbx(i), orby(i)) = hamil(j, orbx(i), orby(i)) + &
      (hopr(i) + iota*hopi(i)) * exp(iota*((0.9789*0.284*pi-0.553*k1(j))* rx(i) + &
      (0.9789*0.284*pi+0.553*k1(j)) * ry(i) +  0.9842*k2(j)*rz(i))) 
    End Do
    Bk(j,:,:) = hamil(j,:,:)
    !-------------diagonalising the hamiltonian---------------!
    Call ZHEEV('V','U',Nm,Bk(j,:,:),LDA,Wk(j,:),WORK,LWORK,RWORK,INFO)

    If( INFO.Gt.0 ) Then
    Write(*,*)'The algorithm failed to compute eigenvalues.'
    Stop
    End If
    
    write(15,*) k1(j),k2(j),Wk(j,:)
    write(*,*) j
  enddo
   
  Close(15) 
   
  print*, 'data generated successfully!'
   
   


END PROGRAM hamiltonian
