PROGRAM GenerateDataFile
   INTEGER, PARAMETER :: num_divisions = 99
   REAL :: pi, delta, delta1, delta2, x, y, zeta
   INTEGER :: i
   CHARACTER(20) :: file_name
  

   ! Calculate the value of pi
   pi = ACOS(-1.0)

   ! Calculate the step size (delta) for each division
   zeta = 0.3297 * 2 * pi
   delta = zeta / REAL(num_divisions)

   ! Open the file for writing
   file_name = 'kpath_data.dat'
   OPEN(UNIT=10, FILE=file_name, STATUS='unknown')

   ! Generate and write the first set of data to the file
   DO i = 1, num_divisions + 1
      x = 0 + REAL(i - 1) * delta
      y = 0 - REAL(i - 1) * delta
      z = 0
      WRITE(10, *) x, y, z
   END DO

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'First set of data generated successfully!'

!--------------------------------------------------------------------!
   ! Calculate the step size (delta) for the second set of data
   delta1 = (0.67 * 2 * pi - zeta) / REAL(num_divisions)
   delta2 = (2*zeta) / REAL(num_divisions)

   ! Reopen the file for appending
   OPEN(UNIT=10, FILE=file_name, STATUS='old', POSITION='append')

   ! Generate and append the second set of data to the file
   DO i = 1, num_divisions + 1
      x = zeta + REAL(i - 1) * delta1
      y = -zeta + REAL(i - 1) * delta2
      z = 0
      WRITE(10, *) x, y, z
   END DO

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'Second set of data appended successfully!'

!-------------------------------------------------------------------!
   ! Calculate the step size (delta) for the Third set of data
   delta1 = (0.67 * 2 * pi - 0.5 * 2 * pi) / REAL(num_divisions)
   delta2 = (0.5 * 2* pi - zeta) / REAL(num_divisions)

   ! Reopen the file for appending
   OPEN(UNIT=10, FILE=file_name, STATUS='old', POSITION='append')

   ! Generate and append the second set of data to the file
   DO i = 1, num_divisions + 1
      x = (2*pi - zeta) - REAL(i - 1) * delta1
      y = zeta + REAL(i - 1) * delta2
      z = 0
      WRITE(10, *) x, y, z
   END DO

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'Third set of data appended successfully!'


!-------------------------------------------------------------------!
   ! Calculate the step size (delta) for the Fourth set of data
   delta1 = (0.5 * 2 * pi) / REAL(num_divisions)
   delta2 = (0.5 * 2* pi ) / REAL(num_divisions)

   ! Reopen the file for appending
   OPEN(UNIT=10, FILE=file_name, STATUS='old', POSITION='append')

   ! Generate and append the second set of data to the file
   DO i = 1, num_divisions + 1
      x = 0.5 * 2* pi  - REAL(i - 1) * delta1
      y = 0.5 * 2* pi  - REAL(i - 1) * delta2
      z = 0  + REAL(i - 1) * delta2
      WRITE(10, *) x, y, z
   END DO

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'Fourth set of data appended successfully!'



!-------------------------------------------------------------------!
   ! Calculate the step size (delta) for the Fifth set of data
   delta1 = (0.5 * 2 * pi) / REAL(num_divisions)
   delta2 = (0.5 * 2* pi ) / REAL(num_divisions)

   ! Reopen the file for appending
   OPEN(UNIT=10, FILE=file_name, STATUS='old', POSITION='append')

   ! Generate and append the second set of data to the file
   DO i = 1, num_divisions + 1
      x = 0.0  + REAL(i - 1) * delta1
      y = 0.0  + REAL(i - 1) * delta2
      z = 0.5 * 2* pi
      WRITE(10, *) x, y, z
   END DO

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'Fifth set of data appended successfully!'
   
   

!-------------------------------------------------------------------!
   ! Calculate the step size (delta) for the Sixth set of data
   delta1 = ((2 * pi - zeta) - 0.5 * 2 * pi) / REAL(num_divisions)
   delta2 = (0.5 * 2* pi - zeta ) / REAL(num_divisions)

   ! Reopen the file for appending
   OPEN(UNIT=10, FILE=file_name, STATUS='old', POSITION='append')

   ! Generate and append the second set of data to the file
   DO i = 1, num_divisions + 1
      x = 0.5 * 2* pi  + REAL(i - 1) * delta1
      y = 0.5 * 2* pi  - REAL(i - 1) * delta2
      z = 0.5 * 2* pi
      WRITE(10, *) x, y, z
   END DO

   ! Close the file
   CLOSE(10)

   WRITE(*, *) 'Sixth set of data appended successfully!'   


!-------------------------------------------------------------------!
   ! Calculate the step size (delta) for the Seventh set of data
delta1 = ((2 * pi - zeta) - zeta) / REAL(num_divisions)
delta2 = (2 * zeta ) / REAL(num_divisions)

! Reopen the file for appending
OPEN(UNIT=10, FILE=file_name, STATUS='old', POSITION='append')

! Generate and append the second set of data to the file
DO i = 1, num_divisions + 1
   x = ( 2* pi - zeta)  - REAL(i - 1) * delta1
   y = zeta  - REAL(i - 1) * delta2
   z = 0.5 * 2* pi
   WRITE(10, *) x, y, z
END DO

! Close the file
CLOSE(10)

WRITE(*, *) 'Seventh set of data appended successfully!'   

END PROGRAM GenerateDataFile
