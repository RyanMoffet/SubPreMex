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
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
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
      REAL*8 THETA
      REAL*8 X,QSCA,QSCA1,QEXT,QEXT1,QBCK,QBCK1
      COMPLEX REFREL,S1,S2,SI1,SI2
      qsca=1
	qbck=1
	qext=1
	s1=1
	s2=1
	si1=1
	si2=1
	return
	end
