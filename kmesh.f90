PROGRAM GenerateDataFile
   INTEGER :: num_divisions
   REAL :: pi, delta, x,y
   INTEGER :: i,j,k
   CHARACTER(20) :: file_name
  
   !read the num_division value
   Write(*,*) "put the num division value"
   READ(*,*) num_divisions

   ! Calculate the value of pi
   pi = ACOS(-1.0)

   ! Calculate the step size (delta) for each division
   delta = (2.0 * pi) / REAL(num_divisions)

   ! Open the file for writing
   WRITE(file_name, '(A)') 'kmesh_file.dat'
   OPEN(UNIT=10, FILE=file_name, STATUS='REPLACE')

   ! Generate and write the data to the file
   DO i = 1, num_divisions + 1
   DO j = 1, num_divisions + 1
   ! DO k = 1, num_divisions + 1
      x = -pi + (REAL(i - 1) * delta)
      y = -pi + (REAL(j - 1) * delta)
      ! z = -pi + (REAL(k - 1) * delta)
     WRITE(10, *) x, y
   ! ENDDO
   ENDDO
   ENDDO

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'Data file generated successfully!'

END PROGRAM GenerateDataFile


