C *** The "Gateway Routine" ***
C     this subroutine gets the data to and from matlab
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on the DEC Alpha and the
C     SGI 64-bit platforms
C
c      mwpointer plhs, prhs
c      mwpointer mxGetPr, mxCreateFull
c      mwpointer REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,LGEOM,R
C-----------------------------------------------------------------------
C

      integer nlhs, nrhs
	integer plhs(1), prhs
      integer mxGetM, mxGetN, mxIsNumeric
	integer mxCreateDoubleMatrix
      integer m, n, size
	integer LGEOM,LG
	real REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER
	real R,RM,RR,RI,WL,DP
C     Check for proper number of arguments. 
      if(nrhs .ne. 6) then
         call mexErrMsgTxt('Six inputs required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgTxt('One output required.')
      endif

C     Get the size of the input array.
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      size = m*n

C     Check to insure the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input must be a number.')
      endif

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(m,n,0)
      RMP = mxGetPr(prhs(1))
      RRP = mxGetPr(prhs(2))
      RIP = mxGetPr(prhs(3))
      WLP = mxGetPr(prhs(4))
	DPP = mxGetPr(prhs(5))
      LGP = mxGetPr(prhs(6))
      RP= mxGetPr(plhs(1))
      call mxCopyPtrToReal8(RMP,RM,size)
	call mxCopyPtrToReal8(RRP,RR,size)
	call mxCopyPtrToReal8(RIP,RI,size)
	call mxCopyPtrToReal8(WLP,WL,size)
	call mxCopyPtrToReal8(DPP,DP,size)
C	call mxCopyPtrToReal8(LGP,LG,1)


C     Call the computational subroutine.%VAL(RM),%VAL(RR),%VAL(RI),%VAL(WL),%VAL(DP),%VAL(LG),%VAL(R))
      CALL MC(RM,RR,RI,WL,DP,LG,R) 
      

C     Load the data into y_pr, which is the output to MATLAB
      call mxCopyReal8ToPtr(R,RP,size)     

      return
      end
C *** Computational Subroutine ***
      SUBROUTINE MC(REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,LGEOM,R)
C
C Program SCATTER is designed to calculate scattering cross section
C (cm2) vs particle diameter (um) for a sphere illuminated by a laser
C beam for various scattering geometries.
C Formulae by  W.W. Szymanski (1986), init. program by S. Palm (1986).
C Final corrections by A. Majerowicz (1986).
C Corrections to off-axis calculations and modification of integration
C tolerance conditions by M. Stolzenburg (1990).
C
      PARAMETER ( ndim = 15 )
C
      REAL PI
      REAL BETA1,BETA2,THETA0,PHI0,THETALO,THETAHI
      REAL WAVELEN,DIAMETER,REFMED,REFREAL,REFIMAG,R
      REAL TOLR,PRECCAL,PRECSTO
      REAL AUX(ndim),X,RV,FSYM,FACTOR,RREF(101)
      INTEGER ILLUM,LGEOM,LSTEP,NUMSTEP,NPREC,LDBG
      COMPLEX REFREL
      CHARACTER*80 FILNAM
      CHARACTER*20 INPUT
      CHARACTER*1 FLAG(101)
      LOGICAL LOPEX
C
      COMMON/SCATPI/PI
      COMMON/MIE/X,REFREL,ILLUM,BETA1,BETA2,THETA0,PHI0
C
      EXTERNAL FNMIE
C
      PI=ACOS(-1.0)
C
      
C Calculates scattering cross section of particles over specified
C collection geometry (LGEOM) with specified illumination (ILLUM)
C
C  ILLUM   Beams            |   LGEOM  Axis        PHI Range
C    1    single            |     1     on         -PI to +PI
C    2    dual,incoherent   |     2     off  PHI0-DPHI to PHI0+DPHI
C    3    dual,coherent     |
C
C *** This if statement is for the off/on axis geometries ***
      IF (LGEOM.EQ.1) THEN
	  BETA1DEG=7.2
	  BETA2DEG=172.8
        BETA1=BETA1DEG*PI/180.
	  BETA2=BETA2DEG*PI/180.
        THETALO=BETA1
        THETAHI=BETA2
      ELSE
        BETA1DEG=0.
	  BETA2DEG=73.
        BETA1=BETA1DEG*PI/180.
	  BETA2=BETA2DEG*PI/180. 
        THETA0DEG=90
        THETA0=THETA0DEG*PI/180.
        PHI0DEG=90.
        PHI0=PHI0DEG*PI/180.
        THETALO=THETA0-BETA2
        THETAHI=THETA0+BETA2
      ENDIF

      ILLUM=1

C ****   Calculate Machine Precision ****
C
      NPREC=0
   80 NPREC=NPREC+1
      IF (1.+1./2.**NPREC.GT.1.) GOTO 80
      PRECCAL=1./2.**(NPREC-1)
      NPREC=0
   85 NPREC=NPREC+1
      PRECSTO=1.+1./2.**NPREC
      IF (PRECSTO.GT.1.) GOTO 85
      PRECSTO=1./2.**(NPREC-1)
      TOLR=1.0*10.**-5.
C *** LDBG is for showing the output of QATR ***
      LDBG=0
C
C **** USE SYMMETRY ABOUT THETA=90 TO MINIMIZE INTEGRATION RANGE ****
C
      FSYM=1.
      IF (ABS(THETAHI+THETALO-PI).LT.1.E-5) THEN
        THETAHI=PI/2.
        IF (ILLUM.EQ.1) THEN
          ILLUM=2
        ELSE
          FSYM=2.
        ENDIF
      ENDIF
C
C **** CALCULATE PARAMETERS ****
C
      FACTOR=FSYM*(0.5E-4*WAVELEN/PI)**2
      REFREL=CMPLX(REFREAL,REFIMAG)/REFMED
      X=PI*DIAMETER*REFMED/WAVELEN
      VOL=PI*DIAMETER**3/6.
C
C ********************
C Start of Integration
C ********************
C
      PRM=DIAMETER
      REFREL=CMPLX(REFREAL,REFIMAG)/REFMED
      CALL QATR (THETALO,THETAHI,0.,TOLR,ndim,LDBG,FNMIE,R,IER,AUX)

C
C      IF (IER.GT.ndim) WRITE(6,605)
      R=R*FACTOR
      RV=R/VOL*1.E8
C      WRITE(6,580) PRM,R,RV,IER
C
C Finish up.
C
      RETURN
	END
C **********************************************************************
      SUBROUTINE QATR (XL,XU,EPS,ERR,NDIM,LDBG,FCT,Y,IER,AUX)
C
C Subroutine to integrate a function FCT(x) using Romberg quadrature
C XL,XU = lower, upper limits of integral
C EPS = absolute error bound for the numerical integration
C ERR = relative error bound for the numerical integration
C NDIM = maximum number of bisections of integration domain
C LDBG = flag for debugging output
C FCT = user defined function or integrand to be integrated
C Y = returned value of integral
C IER = error flag or number of bisections done
C   1<=IER<=NDIM: number of bisections done, convergence criterion met
C       IER>NDIM: convergence criterion not met, NDIM bisections done
C
C DELT1 = estimated error for previous iteration
C DELT2 = estimated error for current iteration
C ICNV = no. of consecutive iterations ABS(DELT2) has decreased
C JCNV = no. of consecutive iterations ABS(DELT2) has been
C        within at least one of the error bounds EPS or ERR
C Convergence is accepted when (JCNV>=NCNV) or (ICNV>=MCNV and JCNV>0).
C When integrating over highly oscillatory Mie scattering patterns false
C apparent convergence sometimes occurs when only two consecutive
C iterations are considered.  More complex convergence criteria applied
C here virtually eliminates the possibility of acceptance of a false
C convergence.
C
      DIMENSION AUX(1)
C
      I=1
      AUX(1) = (FCT(XL)+FCT(XU)) / 2.
      DELT2 = 0.
      H = XU-XL
      IF (H.EQ.0.0) GOTO 8
C
      E = EPS/ABS(H)
      MCNV = 5
      ICNV = 0
      NCNV = 3
      JCNV = 0
      JJ = 1
      DO 7 I=2,NDIM
        RJJ = REAL(JJ)
        SM = 0.
        DO 3 J=1,JJ
          X = XL + ((REAL(J)-0.5)/RJJ) * H
    3     SM = SM + FCT(X)
        AUX(I) = (AUX(I-1) + SM/RJJ) / 2.
        Y = AUX(1)
        Q = 1.
        DO 4 II=I-1,1,-1
          Q = 4. * Q
    4     AUX(II) = AUX(II+1) + (AUX(II+1)-AUX(II)) / ( Q-1. )
        DELT1 = DELT2
        DELT2 = Y-AUX(1)
        IF (LDBG.EQ.2) WRITE(6,500) I,DELT2,ABS(DELT2/AUX(1))
        IF (ABS(DELT2).GT.ABS(DELT1)) ICNV = -1
        ICNV = ICNV+1
        IF (ABS(DELT2).GT.E .AND. ABS(DELT2/AUX(1)).GT.ERR) JCNV = -1
        JCNV = JCNV+1
        IF ((ICNV.GE.MCNV.AND.JCNV.GT.0) .OR. (JCNV.GE.NCNV)) GOTO 8
    7   JJ= 2 * JJ
      I=NDIM+1
    8 IER=I
      Y = H*AUX(1)
      IF (LDBG.EQ.1) WRITE(6,500) I,DELT2,ABS(DELT2/AUX(1))
C
      RETURN
  500 FORMAT (' Qatr Iter = ',I3,5X,'Err = ',1PE11.4,
     +    5X,'Rel Err = ',E10.4)
      END
C **********************************************************************
      REAL FUNCTION FNMIE (THETA)
C
C Returns scattered intensity per dTHETA at zenith THETA
C over azimuth PHI range specified by LGEOM with ILLUMination
C
C  ILLUM   Beams            |   LGEOM  Axis        PHI Range
C    1    single            |     1     on         -PI to +PI
C    2    dual,incoherent   |     2     off  PHI0-DPHI to PHI0+DPHI
C    3    dual,coherent     |
C
      INTEGER LGEOM,ILLUM
      REAL PI,THETA,BETA1,BETA2,THETA0,PHI0,DPHI,DPHI1
      REAL X,I1,I2,FD
      REAL QSCA,QEXT,QBCK
      COMPLEX REFREL,S1,S2,SI1,SI2
C
      COMMON/SCATPI/PI
      COMMON/MIE/X,REFREL,ILLUM,LGEOM,BETA1,BETA2,THETA0,PHI0
C
      CALL BHMIE (X,REFREL,THETA,S1,S2,SI1,SI2,QSCA,QEXT,QBCK)
C
      IF(ILLUM.EQ.1) THEN
        I1=CABS(S1)**2
        I2=CABS(S2)**2
      ELSEIF(ILLUM.EQ.2) THEN
        I1=CABS(S1)**2+CABS(SI1)**2
        I2=CABS(S2)**2+CABS(SI2)**2
      ELSE
        I1=CABS(S1+SI1)**2
        I2=CABS(S2+SI2)**2
      ENDIF
C
      IF(LGEOM.EQ.1) THEN
        FNMIE=(I1+I2)*PI*SIN(THETA)
      ELSE
        DPHI=(COS(BETA2)-COS(THETA)*COS(THETA0))
     1   /(SIN(THETA)*SIN(THETA0))
        DPHI=ACOS(SIGN(MIN(ABS(DPHI),1.),DPHI))
        FD=0.25*(SIN(2.*(PHI0+DPHI))-SIN(2.*(PHI0-DPHI)))
        IF(ABS(THETA-THETA0).LT.BETA1) THEN
          DPHI1=(COS(BETA1)-COS(THETA)*COS(THETA0))
     1     /(SIN(THETA)*SIN(THETA0))
          DPHI1=ACOS(SIGN(MIN(ABS(DPHI1),1.),DPHI1))
          DPHI=DPHI-DPHI1
          FD=FD-0.25*(SIN(2.*(PHI0+DPHI1))-SIN(2.*(PHI0-DPHI1)))
        ENDIF
        FNMIE=(I1*(DPHI-FD)+I2*(DPHI+FD))*SIN(THETA)
      ENDIF
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
      REAL AMU,THETA,PI0,PI1,PI,TAU
      REAL X,X1,QSCA,QSCA1,QEXT,QEXT1,QBCK,QBCK1
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
