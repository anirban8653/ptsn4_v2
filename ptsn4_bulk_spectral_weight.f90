PROGRAM hamiltonian
   IMPLICIT NONE

   INTEGER :: i, j, k, l
   Complex*16 :: trace
   Real*8, Parameter :: pi = ACos(-1.D0)
   Real*8, Parameter :: afb = 0.005
   Real*8, Parameter :: x = 0.D0
   Real*8, Parameter :: vimp = 100.D0
   Integer, Dimension(:), Allocatable :: rx, ry, rz, orbx, orby
   Real*8, Dimension(:), Allocatable :: kx, ky, kz, hopr, hopi!, kxp, kzp
   Complex*16, DIMENSION(:,:), ALLOCATABLE :: identity
   Complex*16, Dimension(:,:,:), Allocatable :: hamil,Gk
   Complex*16, Parameter :: iota = (0.D0, 1.D0)
   Integer, Parameter :: norb = 84 ! number of orbitals
   Integer, Parameter :: num = 882000!7490880 ! number of lines in hr file
   Integer, Parameter :: kkm = 2
   

   ! Open the data files
   Open(12, File="harald_hr.dat", Status="Unknown") 
   Open(13, File="kmesh_file.dat", Status="Unknown")
   !Open(14, File="kmesh_file_xz.dat", Status="Unknown")
   Open(15, File="A_ptsn4_bulk.dat", Status="Unknown")
   
   Allocate(identity(norb,norb))
   Allocate(rx(num))
   Allocate(ry(num))
   Allocate(rz(num))
   Allocate(kx((kkm+1)**3))
   Allocate(ky((kkm+1)**3))
   Allocate(kz((kkm+1)**3))
   !Allocate(kxp((kkm+1)**2))
   !Allocate(kzp((kkm+1)**2))
   Allocate(orbx(num))
   Allocate(orby(num))
   Allocate(hopr(num))
   Allocate(hopi(num))
   Allocate(hamil((kkm+1)**3, norb, norb))
   Allocate(Gk((kkm+1)**3, norb, norb))
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
print*, '1. reading data file'
! Read the data from the files
   Do i = 1, num
     Read(12, *) rx(i), ry(i), rz(i), orbx(i), orby(i), hopr(i), hopi(i)
   End Do
   
   Do i = 1, (kkm+1)**3
     Read(13, *) kx(i), kz(i), ky(i)
   End Do

   ! Close the files
   Close(12)
   Close(13)
   Close(14)
   
 print*, '2. identity matrix'
  !Initialize the identity matrix
   Do i = 1, norb
     Do j = 1, norb
       If (i == j) Then
         identity(i,j) = 1.0
       Else
         identity(i,j) = 0.0
       End If
     End Do
   End Do
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   hamil = (0.D0, 0.D0)
   print*, '3. starting fourier transformation'
   Do j = 1,  (kkm+1)**2
  
     Do i = 1, num
       hamil(j, orbx(i), orby(i)) = hamil(j, orbx(i), orby(i)) + &
       (hopr(i) + iota*hopi(i)) * exp(iota*(kx(((kkm)/2 + 1)+ (j - 1)*(kkm+1)) * rx(i) + &
       ky(((kkm)/2 + 1) + (j - 1)*(kkm+1)) * ry(i) + kz(((kkm)/2 + 1)  + (j - 1)*(kkm+1)) * rz(i))) 
     End Do
  
   Gk(j,:,:) = inv((x + iota * afb) * identity(:,:) - hamil(j,:,:))
 
  !calculating the trace of greens function
  trace = (0.D0,0.D0)
  do i = 1, norb
  trace = trace + Gk(j,i,i)
  end do 
  
  !writing the output file
   write(15,*)  -1/pi*aimag(trace) 
   
  write(*,*) j
   enddo
   
  Close(15) 
   
  print*, 'data generated successfully!'
   
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Contains

Function inv(A) result(Ainv)
  Complex*16, dimension(:,:), intent(in) :: A
  Complex*16, dimension(size(A,1),size(A,2)) :: Ainv

  Complex*16, dimension(size(A,1)) :: work  ! work array for LAPACK
  Integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  Integer :: n, info

  ! External procedures defined in LAPACK
  External ZGETRF
  External ZGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  Call ZGETRF(n, n, Ainv, n, ipiv, info)

  If (info /= 0) then
     Stop 'Matrix is numerically singular!'
  End If

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  Call ZGETRI(n, Ainv, n, ipiv, work, n, info)

  If (info /= 0) then
     Stop 'Matrix inversion failed!'
  End If
  Return
End Function inv

END PROGRAM hamiltonian
