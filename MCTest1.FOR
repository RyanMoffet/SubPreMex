C *** The "Gateway Routine" ***
C     this subroutine gets the data to and from matlab
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (integer) Replace integer by integer*8 on the DEC Alpha and the
C     SGI 64-bit platforms
C
      integer plhs(*), prhs(*)
      integer mxCreateFull,mxGetPr
      integer RMP, RRP, RIP, WLP, DPP, RP
C-----------------------------------------------------------------------
C
      integer nlhs, nrhs
      integer mxGetM, mxGetN, mxIsNumeric
      integer m, n, size
	real*8 REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,QEXT
C     Check for proper number of arguments. 
      if(nrhs .ne. 5) then
         call mexErrMsgTxt('Six inputs required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgTxt('One output required.')
      endif

C     Get the size of the input array.
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
      size = m*n

C     Check to insure the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input must be a number.')
      endif

C     Create matrix for the return argument.
      plhs(1) = mxCreateFull(m,n,0)
      RMP = mxGetPr(prhs(1))
      RRP = mxGetPr(prhs(2))
      RIP = mxGetPr(prhs(3))
      WLP = mxGetPr(prhs(4))
	DPP = mxGetPr(prhs(5))
      RP = mxGetPr(plhs(1))
      call mxCopyPtrToReal8(RMP,REFMED,size)
	call mxCopyPtrToReal8(RRP,REFREAL,size)
	call mxCopyPtrToReal8(RIP,REFIMAG,size)
	call mxCopyPtrToReal8(WLP,WAVELEN,size)
	call mxCopyPtrToReal8(DPP,DIAMETER,size)



C     Call the computational subroutine.
      CALL TSTBH (REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,QEXT)
      
C     Load the data into y_pr, which is the output to MATLAB
      call mxCopyReal8ToPtr(QEXT,RP,size)     

      return
      end
C *** Computational Subroutine ***
 
C **********************************************************************
      SUBROUTINE TSTBH (REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,QEXT)
C
C Returns scattered intensity per dTHETA at zenith THETA
C over azimuth PHI range specified by LGEOM with ILLUMination
C
C  ILLUM   Beams            |   LGEOM  Axis        PHI Range
C    1    single            |     1     on         -PI to +PI
C    2    dual,incoherent   |     2     off  PHI0-DPHI to PHI0+DPHI
C    3    dual,coherent     |
C
      REAL*8 REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,X
      REAL*8 THETA
      REAL*8 QSCA,QEXT,QBCK
      COMPLEX REFREL,S1,S2,SI1,SI2
C
C
      X=3.14159265*DIAMETER*REFMED/WAVELEN
	THETA=90
      CALL BHMIE (X,REFREL,THETA,S1,S2,SI1,SI2,QSCA,QEXT,QEXT)
C
      RETURN
      END
C **********************************************************************
      SUBROUTINE BHMIE (X,REFREL,THETA,S1,S2,SI1,SI2,QSCA,QEXT,QBCK)
C
C Subroutine BHMIE calculates scattering irradiance
C functions S1 and S2 for a given size parameter X,
C relative complex refractive index REFREL and scattering
C angle THETA.
C
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
      REAL*8 AMU,THETA,PI0,PI1,PI,TAU
      REAL*8 X,X1,QSCA,QSCA1,QEXT,QEXT1,QBCK,QBCK1
      COMPLEX REFREL,REFREL1,Y,D1,XI1,XI,S1,S2,SI1,SI2,S1B
      COMPLEX D(3000),AN(3000),BN(3000)
      EQUIVALENCE (D,BN)
      DATA X1,REFREL1/0.0,(0.0,0.0)/
C
      IF(X.NE.X1 .OR. REFREL.NE.REFREL1) THEN
C
C Series terminated after NSTOP terms
C       
        XSTOP=X+4.*X**(1./3.)+2.0
        NSTOP=NINT(XSTOP)
        IF (NSTOP.GT.3000) THEN
          WRITE (6,*) 'ERROR IN SUBROUTINE BHMIE'
          WRITE (6,*) 'ARRAY D() IS NOT LARGE ENOUGH'
          WRITE (6,*) 'TRY DIMENSION OF GREATER THAN',NSTOP
          STOP
        ENDIF
C
C Logarithmic derivative D(J) calculated by DOWNWARD
C recurrence beginning with initial value 0.0+I*0.0
C at J = NMX
C
        Y=X*REFREL
        NMX=NINT(MAX(XSTOP,CABS(Y)))+15
        D1=CMPLX(0.0,0.0)
        DO 100 N=NMX,2,-1
          RN=REAL(N)
          D1=(RN/Y)-(1./(D1+RN/Y))
  100     D(MIN(NSTOP,N-1))=D1
C
C Riccati-Bessel functions with real argument X
C calculated by upward recurrence
C
        DX=DBLE(X)
        PSI1=COS(DX)
        PSI=SIN(DX)
        CHI1=-SIN(X)
        CHI=COS(X)
        APSI=REAL(PSI)
        XI=CMPLX(APSI,-CHI)
        QSCA=0.0
        QEXT=0.0
        S1B=CMPLX(0.0,0.0)
        SGNN=1.
C
        DO 200 N=1,NSTOP
C
          DN=DBLE(N)
          RN=REAL(N)
          PSI0=PSI1
          PSI1=PSI
          PSI=(2.D0*DN-1.D0)*PSI1/DX-PSI0
          CHI0=CHI1
          CHI1=CHI
          CHI=(2.*RN-1.)*CHI1/X-CHI0
          APSI1=APSI
          APSI=REAL(PSI)
          XI1=XI
          XI=CMPLX(APSI,-CHI)
          AN(N)=((D(N)/REFREL+RN/X)*APSI-APSI1)
     1         /((D(N)/REFREL+RN/X)*XI-XI1)
          BN(N)=((D(N)*REFREL+RN/X)*APSI-APSI1)
     1         /((D(N)*REFREL+RN/X)*XI-XI1)
C
          SGNN=-SGNN
          QSCA=QSCA+(2.*RN+1.)*(CABS(AN(N))**2+CABS(BN(N))**2)
          QEXT=QEXT+(2.*RN+1.)*REAL(AN(N)+BN(N))
          S1B=S1B+SGNN*(2.*RN+1.)*(BN(N)-AN(N))
C
  200   CONTINUE
C
        QSCA1=(2./X**2)*QSCA
        QEXT1=(2./X**2)*QEXT
        QBCK1=(1./X**2)*CABS(S1B)**2
        X1=X
        REFREL1=REFREL
C
      ENDIF
C
      QSCA=QSCA1
      QEXT=QEXT1
      QBCK=QBCK1
C
C Angular dependent factors PI and TAU
C calculated by upward recurrence
C
      AMU=COS(THETA)
      PI1=0.0
      PI=1.0
      S1=CMPLX(0.0,0.0)
      S2=CMPLX(0.0,0.0)
      SI1=CMPLX(0.0,0.0)
      SI2=CMPLX(0.0,0.0)
      SGNN=1.
C
      DO 300 N=1,NSTOP
C
        RN=REAL(N)
        FN=(2.*RN+1.)/(RN*(RN+1.))
        SGNN=-SGNN
        TAU=RN*AMU*PI-(RN+1.)*PI1
        S1=S1+FN*(BN(N)*TAU+AN(N)*PI)
        S2=S2+FN*(AN(N)*TAU+BN(N)*PI)
        SI1=SI1+SGNN*FN*(BN(N)*TAU-AN(N)*PI)
        SI2=SI2+SGNN*FN*(AN(N)*TAU-BN(N)*PI)
C
        PI0=PI1
        PI1=PI
        PI=((2.*RN+1.)/RN)*AMU*PI1-(RN+1.)*PI0/RN
C
  300 CONTINUE
C
      RETURN
      END
