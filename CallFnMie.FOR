      PROGRAM CallFnMie
C
      REAL I1
      REAL O1,THETADEG,PI
C  
  
C
C *** Call the Mie Calibration (MC) Subroutine ***

	READ(5,'(F20.0)') THETA
	PI=ACOS(-1.0)		
	O1=FNMIE(THETA)
      WRITE(*,*) O1  
      STOP 'CallFnMie - completed.'
C
      END
