C *** The "Gateway Routine" ***
C     this subroutine gets the data to and from matlab
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------
C     (integer) Replace integer by integer*8 on the DEC Alpha and the
C     SGI 64-bit platforms
C
      integer plhs(*), prhs(*)
      integer mxCreateFull,mxGetPr
C	integer LGP
      integer RMP, RRP, RIP, WLP, DPP, LGP, RP
C-----------------------------------------------------------------------
C
      integer nlhs, nrhs
      integer mxGetM, mxGetN, mxIsNumeric
      integer m, n, size
	integer LGEOM
	real*8 REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,R,LG
C	real RP,RMP,RRP,RIP,WLP,DPP
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
      plhs(1) = mxCreateFull(m,n,0)
      RMP = mxGetPr(prhs(1))
      RRP = mxGetPr(prhs(2))
      RIP = mxGetPr(prhs(3))
      WLP = mxGetPr(prhs(4))
	DPP = mxGetPr(prhs(5))
      LGP = mxGetPr(prhs(6))
      RP= mxGetPr(plhs(1))
      call mxCopyPtrToReal8(RMP,REFMED,size)
	call mxCopyPtrToReal8(RRP,REFREAL,size)
	call mxCopyPtrToReal8(RIP,REFIMAG,size)
	call mxCopyPtrToReal8(WLP,WAVELEN,size)
	call mxCopyPtrToReal8(DPP,DIAMETER,size)
	call mxCopyPtrToReal8(LGP,LG,size)


C     Call the computational subroutine.%VAL(RM),%VAL(RR),%VAL(RI),%VAL(WL),%VAL(DP),%VAL(LG),%VAL(R))
C      CALL MC(RM,RR,RI,WL,DP,LG,R)
c      CALL MC(%VAL(REFMED),%VAL(REFREAL),%VAL(REFIMAG),%VAL(WAVELEN), 
c     +        %VAL(DIAMETER),%VAL(LG),%VAL(R))
      CALL add(REFMED,REFREAL,REFIMAG,WAVELEN, 
     +        DIAMETER,LG,R)
      
C      RP=plhs(1)
C     Load the data into y_pr, which is the output to MATLAB
      call mxCopyReal8ToPtr(R,RP,size)     

      return
      end
C *** Computational Subroutine ***
C *** Computational Subroutine ***
      SUBROUTINE add(REFMED,REFREAL,REFIMAG,WAVELEN,DIAMETER,LG,R)

      REAL*8 PI,REFMED,REFREAL,REFIMAG
      REAL*8 WAVELEN,DIAMETER,R,LG
      

	LGR=REAL(LG)

	R=REFMED+REFREAL+REFIMAG+REAL(WAVELEN)+DBLE(DIAMETER)+LG

	END
