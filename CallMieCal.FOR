      PROGRAM CallMieCal
C
      REAL I1,I2,I3,I4,I6
      REAL I5 
      REAL O1
C  
      READ(5,'(F20.0)') I1
      READ(5,'(F20.0)') I2
      READ(5,'(F20.0)') I3
      READ(5,'(F20.0)') I4
      READ(5,'(F20.0)') I5
      READ(5,'(F20.0)') I6 
      

C
C *** Call the Mie Calibration (MC) Subroutine ***
C
	CALL MIECAL (I1,I2,I3,I4,I5,I6,O1)
      WRITE(*,*) O1  
      STOP 'CallMieCal - completed.'
C
      END
