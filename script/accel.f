C--------------------------------------------------------------------
C
C Subroutine ACCEL --> Calls subroutines that calulate potential
C                      energy and interatomic forces. Converts
C                      forces to accelerations. 
C                                               
C  by    CS Carmer & B Weiner                                               
C        LATEST MODIFICATION      11/ 7/92
C  Modification by H.B. Zhang 2013.11.12 SCFCRT issue
C--------------------------------------------------------------------
 
      SUBROUTINE ACCEL(N,M,X,Y,Z,AX,AY,AZ,SX,SY,SZ,ZPOT,IPOT,
     + LTEST,ITYP,NKEEP,NIMAGE,NRECT,LBOUND,LXROT,LMOMEN,
     + ISNAP,IMOPCALL,NBD)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER QND
      DOUBLE PRECISION M
      LOGICAL LTEST,LAST,LBOUND,LXROT,LMOMEN,LGHOST,LCOMPOT
      LOGICAL NOPOT,DEBUGM, LSCFCRT
      PARAMETER (ND = 2000)
      PARAMETER (QND = 500)
      DIMENSION X(ND),Y(ND),Z(ND)
      DIMENSION AX(ND),AY(ND),AZ(ND),M(ND)
      DIMENSION FFX(ND),FFY(ND),FFZ(ND),FMOP(3*QND),FZIN(3*ND)
      DIMENSION ITYP(ND),IQTYP(ND)           
      DIMENSION XQC(ND),YQC(ND),ZQC(ND)                             
      DIMENSION CMOP(3,QND)
      COMMON /RESTART/ LAST,ICF,IRUN,IDENRC
C      CHARACTER  MOLNAM*20,CLINE*60,CLINE1*60,CSTAT*10
      COMMON /FFILES/ MOLNAM,IPOS
      COMMON /DYNINF/ NEIG,NBAS
      COMMON /DYNFIX/ JFIX(ND),NOTFIX
      COMMON /STEPS/  DT,NTIME,NSNAP
      COMMON /COMBINED/ NQC,ICOL(ND)
      COMMON /LNKS/ LCOMPOT,LGHOST,IGHOST(ND),ILINK(ND),IPTR(ND)
      COMMON /NTMS/ EPS,NAT,ITMAX
      COMMON /DYNZIN1/ NAA,IGIT,IXS,IOPTYP,MXGEOM,ICHG,
     +                 NCHARG,IBASIS,EK,ELED(10)
      COMMON /DYNZIN2/ LRESET,LCONEX
C ********************************************************************
C In order to change force from Hartrees/angstrom to (atomic mass units)
C *(angstroms)/(femto second)**2 one must multiply by 0.262549
C To change form Bohrs to angstroms one must divide by 0.52917725
C To change form Kcal to Hartrees one must divide by  627.5
C---------------------------------------------------------------------
      PARAMETER (HCONV = 0.2625499961D0,BOHR=0.52917725D0,
     + ACALH=6.2750499D2)
C---------------------------------------------------------------------
C  Electron volt to Hartree conversion factor (multiply)
C---------------------------------------------------------------------
      PARAMETER (EVHAR = 3.674903D-02)
C ********************************************************************
      DATA IACCEL /-1/ , IDEN /1/, DEBUGM/.FALSE./
      
      NOPOT = .TRUE.
      ICOMPOT = 0
      ZPOT=0
 
      DO 10 I=1,N
        AX(I)=0.0D0
        AY(I)=0.0D0
        AZ(I)=0.0D0
   10 CONTINUE   
   
C**********************MIRROR QANTUM ATOMS ONTO EMPIRICAL ATOMS******   
      IF(IMOPCALL == 100) 
     + CALL QUANTUM2EMPIRICAL(N,M,X,Y,Z,SX,SY,SZ,ITYP,ICOL,IPOT
     + ,LBOUND,NBD)
 
C*********************************************************************
C LENNARD-JONES POTENTIAL   ===>  IPOT = 2
C*********************************************************************
 
      IF (IPOT.EQ.2) THEN
          NOPOT = .FALSE.
          ZPOT = 0.D0
          DO 30 I=1,(N-1)
          DO 20 J=(I+1),N
          DELX=X(I)-X(J)
          DELY=Y(I)-Y(J)
          DELZ=Z(I)-Z(J)
C ********************************************************************
C
C DETERMINE IN WHICH CELL THE PERTINENT ATOMS RESIDE
C
C ********************************************************************
          IF (LBOUND) CALL PERIODIC(DELX,DELY,DELZ,SX,SY,SZ)
 
          R=SQRT(DELX*DELX+DELY*DELY+DELZ*DELZ)
          
	    IF (R.LT.10) THEN
	      IF (R.LT.0.1) R=0.1
		  CALL POTENTAL(R,FORCE,POT)
            FORCE = FORCE/M(I)
		ELSE
		  FORCE =0
		  POT=0
	    ENDIF
          AX(I)=AX(I)+FORCE*DELX
          AY(I)=AY(I)+FORCE*DELY
          AZ(I)=AZ(I)+FORCE*DELZ
 
          AX(J)=AX(J)-FORCE*DELX
          AY(J)=AY(J)-FORCE*DELY
          AZ(J)=AZ(J)-FORCE*DELZ
 
          ZPOT=ZPOT+POT
 
   20   CONTINUE
   30   CONTINUE
      ENDIF
C*********************************************************************
C ZINDO POTENTIAL   ===>   IPOT = 1
C*********************************************************************
 

C*********************************************************************
C MOPAC POTENTIAL   ===>   IPOT = 3
C*********************************************************************
 
  140 IF (IPOT.EQ.3) THEN
         NOPOT = .FALSE.
         IF (ICOMPOT.NE.0) GOTO 150
         NKEEP = N
 
         IF (NIMAGE.NE.0) THEN
           CALL IMAGE (X,Y,Z,NEXTRA,SX,SY,SZ,N,ICOL,ITYP,IPRINT,NIMAGE,
     +                 NRECT,IPOT)
           NKEEP = NEXTRA
         ENDIF
 
         DO II = 1,NKEEP
            CMOP(1,II) = X(II)
            CMOP(2,II) = Y(II)
            CMOP(3,II) = Z(II)
         ENDDO
 
  150    IF (ICOMPOT.NE.0) NKEEP = NQC
         DO I = 1,NKEEP
           FMOP(I) = 0.0D0
           FMOP(I+NKEEP) = 0.0D0
           FMOP(I+2*NKEEP) = 0.0D0
         ENDDO

         CALL getmopacversion(NVERSION)	 
	   WRITE(61,*) 'MOPAC VERSION IS',NVERSION
	   WRITE(61,*)
	   WRITE(61,*) 'NQC=',NQC
C**********************************************************	     
         WRITE(*,*) '*********RUNNING MOPAC******'
         IF( (NVERSION.EQ.7).OR.(DEBUGM.EQV..TRUE.) )
     +     CALL callmopac(NQC,XQC,YQC,ZQC,IQTYP,FMOP,ZPOT)
           LSCFCRT = .False.
           CALL testSCFCRAT(LSCFCRT)
           if( LSCFCRT .eq. .True. ) then
               write(*,*) 'SCFCRT will be changed to 10D-2'
               CALL callmopac_SCF2(NQC,XQC,YQC,ZQC,IQTYP,FMOP,ZPOT)
               write(*,*) 'SCFCRT will be changed back from 10D-2'
           end if
         WRITE(*,*) '*********END MOPAC******'
C**********************************************************	     
         IF(DEBUGM.EQV..TRUE.) THEN  
	     WRITE(61,*) 'AFTER'	     
	     WRITE(61,*) 'ZPOT 1 =',ZPOT
           WRITE(61,*) 'nkeep=',NKEEP
	     do I = 1,NKEEP
	       wriTE(61,*) 'X=',X(I),Y(I),Z(I)
	     enddo
	     write(61,*) '' 
	     do I = 1,(3*NKEEP)
	       write(61,*) 'FMOP=',FMOP(I)
	     enddo
	     wriTE(61,*) '' 
	   ENDIF

         ZPOT = ZPOT*EVHAR
 
C*************FOR COMBINED POTENTIAL**********************************
         IF (ICOMPOT.NE.0) THEN
           IF(ICOMPOT.EQ.7) IPOT=7
           IF(ICOMPOT.EQ.10) IPOT=10
           IF(ICOMPOT.EQ.11) IPOT=11
           IF(ICOMPOT.EQ.15) IPOT=15
           GOTO 160
         ENDIF
C*********************************************************************
 
         FSUM = 0.D0
         DO  I=1,NKEEP
C************************* ADDED 11/21/91 ****************************
           FMOP(I) = -FMOP(I) * HCONV/ACALH
           FMOP(I+NKEEP) = -FMOP(I+NKEEP) * HCONV/ACALH
           FMOP(I+2*NKEEP) = -FMOP(I+2*NKEEP) * HCONV/ACALH
C**********************************************************************
           FSUMI = FMOP(I)**2 + FMOP(I+NKEEP)**2 + FMOP(I+2*NKEEP)**2
           FSUM = FSUM + FSUMI
           FSUMI = DSQRT(FSUMI)
         ENDDO
 
         FSUM = DSQRT(FSUM)

         NLIM = N
         IF (LTEST) NLIM = NKEEP
         DO  I=1,NLIM
           IF (JFIX(I).EQ.0) THEN
              AX(I) = FMOP(I)/M(I)
              AY(I) = FMOP(I+NKEEP)/M(I)
              AZ(I) = FMOP(I+2*NKEEP)/M(I)
           ELSE
              AX(I) = 0.0D0
              AY(I) = 0.0D0
              AZ(I) = 0.0D0
           ENDIF
         ENDDO
      ENDIF
 
C*********************************************************************
C STILLINGER-WEBER TYPE             IPOT = 4 OR IPOT = 5
C POTENTIAL FUNCTION
C                     
C MPB DIAMOND POTENTIAL             IPOT = 8
C J.Appl.Phys. 70,543 (1991)
C
C EMPIRICAL DIAMOND POTENTIAL BASED ON AM1   IPOT = 9
C*********************************************************************
 
      IF (IPOT.EQ.4.OR.IPOT.EQ.5.OR.IPOT.EQ.8.OR.IPOT.EQ.9.
     +    OR.IPOT.EQ.12) THEN
        NOPOT = .FALSE.
        NN = N
        IF (NIMAGE.NE.0) THEN
          CALL IMAGE (X,Y,Z,NEXTRA,SX,SY,SZ,N,ICOL,ITYP,
     +               IPRINT,NIMAGE,NRECT,IPOT)
          NN = NEXTRA
        ENDIF
C*********************************************************************
C SUBROUTIN  WEBER RETURNS THE INTERATOMIC FORCES AS DETERMINED
C BY THE STILLINGER-WEBER POTENTIAL FOR SILICON          (IF IPOT = 4
C BY THE HOUFENG DAI'S PARAMETERIZATION OF ZINDO SILICON (IF IPOT = 5
C BY MPB POTENTIAL FOR DIAMOND                           (IF IPOT = 8
C*********************************************************************
      IF(IPOT.EQ.4.OR.IPOT.EQ.5.OR.IPOT.EQ.8) THEN
      CALL WEBER(NN,M,X,Y,Z,FFX,FFY,FFZ,ZPOT,IPOT,ITYP,LBOUND,SX,SY,SZ)
C*********************************************************************
C SUB AM1C RETURNS THE INTERATOMIC FORCES AS DETERMINED
C BY EMPIRICAL DIAMOND POTENTIAL BASED ON AM1  7/92      (IF IPOT = 9
C*********************************************************************
      ELSEIF (IPOT.EQ.9) THEN
        CALL AM1C(NN,M,X,Y,Z,FFX,FFY,FFZ,ZPOT,IPOT,LBOUND,SX,SY,SZ)
      ELSEIF (IPOT.EQ.12) THEN
        CALL EMPNI(N,M,X,Y,Z,FFX,FFY,FFZ,ZPOT,IPOT,ITYP,LBOUND,SX,SY,SZ)
      ENDIF
 
        DO I=1,N 
          IF (JFIX(I).EQ.1) FFX(I) = 0.0D0
          IF (JFIX(I).EQ.1) FFY(I) = 0.0D0
          IF (JFIX(I).EQ.1) FFZ(I) = 0.0D0
          AX(I)=AX(I)+FFX(I)/M(I)
          AY(I)=AY(I)+FFY(I)/M(I)
          AZ(I)=AZ(I)+FFZ(I)/M(I)
        ENDDO
      ENDIF
C*********************************************************************
C COMBINED POTENTIAL
C ZINDO/HUOFENGDAI                 IPOT = 6
C*********************************************************************
 

C*********************************************************************
C COMBINED POTENTIALS
C HFD/MOPAC AM1                   IPOT = 7
C CHEP1/MOPAC AM1                 IPOT = 10
C EMPNI/MOPCA                     IPOT = 11
C*********************************************************************
 
  160 IF( (IPOT.EQ.7).OR.(IPOT.EQ.10).OR.(IPOT.EQ.11).OR.
     +    (IPOT.EQ.15) ) THEN
        NOPOT = .FALSE.  
        IF (ICOMPOT.NE.0) GOTO 170
C*********************************************************************
C FIRST FIND ALL FORCES-->ACCEL DUE TO EMPIRICAL POTENTIAL
C*********************************************************************
      IF(IPOT.EQ.7) THEN
       CALL WEBER(N,M,X,Y,Z,FFX,FFY,FFZ,ZPOT,IPOT,ITYP,LBOUND,SX,SY,SZ)
      ELSEIF(IPOT.EQ.11) THEN
       CALL EMPNI(N,M,X,Y,Z,FFX,FFY,FFZ,ZPOT,IPOT,ITYP,LBOUND,SX,SY,SZ)
      ELSEIF(IPOT.EQ.15) THEN
       CALL EMPNI(N,M,X,Y,Z,FFX,FFY,FFZ,ZPOT,IPOT,ITYP,LBOUND,SX,SY,SZ)
      ENDIF
C*********************************************************************
C IDENTIFY QUANTUM CHEMICAL ATOMS ---> LOOK AT ICOL()
C*********************************************************************
        NQC=0
        DO I=1,N
          IF (ABS(ICOL(I)).EQ.3) THEN
            NQC       = NQC +1
            IPTR(NQC) = I
            XQC(NQC)  = X(I)
            YQC(NQC)  = Y(I)
            ZQC(NQC)  = Z(I)
            IQTYP(NQC)= ITYP(I)
          ENDIF
        ENDDO
        IF(IACCEL.EQ.-1) WRITE(61,*) 'NQC FROM ACCEL =',NQC
        IF (NQC.EQ.0) GOTO 180
C*********************************************************************
C IDENTIFY EMPIRICAL ATOMS THAT ARE NEAREST NEIGHBOR TO THE QC CLUSTER
C CONVERT THESE TO SATURATING QC ATOMS (GHOST ATOMS)
C*********************************************************************
        IF(LGHOST) THEN
          WRITE(61,*) 'CALCULATING LINK LIST'
          CALL GHOST (N,NQC,XQC,YQC,ZQC,IPTR,ILINK,IQTYP,ICOL,IPOT,ITYP)
        ELSE
          CALL LNK(X,Y,Z,NQC,XQC,YQC,ZQC,IQTYP,IPOT,LBOUND,SX,SY,SZ)
        ENDIF
C*********************************************************************
C     WRITE CURRENT QC COORDS TO UNIT 64 (CHEM3D FORMAT)
C*********************************************************************
        IF (MOD(IACCEL,NTIME).EQ.0.OR.IACCEL.LE.0) THEN
          CALL PRT3D(NQC,XQC,YQC,ZQC,IQTYP,DT,ISNAP,NTIME)
          IF (LGHOST) STOP
        ENDIF

C*********************************************************************
C SEND COORDINATES OF QC ATOMS TO THE QUANTUM CHEMICAL PRGRAM
C RETURN FORCES IN FMOP(I),FMOP(I+NQC),FMOP(I+2*NQC)
C*********************************************************************
        DO I=1,NQC
          CMOP(1,I) = XQC(I)
          CMOP(2,I) = YQC(I)
          CMOP(3,I) = ZQC(I)
        ENDDO

        IF (IPOT.EQ.7) ICOMPOT = 7
        IF (IPOT.EQ.10) ICOMPOT = 10
        IF (IPOT.EQ.11) ICOMPOT = 11
        IF (IPOT.EQ.15) ICOMPOT = 15
        IPOT = 3
        GOTO 140
C*********************************************************************
C CONVERT QC FORCES TO UNITS OF AMU*ANGSTROM/FS**2             
C EMPIRICAL FORCES ALREADY CONVERTED
C*********************************************************************
  170   QCVERT =  HCONV/ACALH
        DO I=1,NQC
          II = IPTR(I)
          FMOP(I)       = -FMOP(I)*QCVERT
          FMOP(I+NQC)   = -FMOP(I+NQC)*QCVERT
          FMOP(I+2*NQC) = -FMOP(I+2*NQC)*QCVERT
        ENDDO

C*********************************************************************
C SUM THE FORCES DUE TO QM AND EP POTENTIALS                     
C*********************************************************************
        DO I = 1,NQC-NIMAGE
          II = IPTR(I)
          FFX(II) =  FFX(II) + FMOP(I)
          FFY(II) =  FFY(II) + FMOP(I+NQC)
          FFZ(II) =  FFZ(II) + FMOP(I+2*NQC)
        ENDDO
 
  180   IDUMMY = 0
        FXSUM=0.D0
        FYSUM=0.D0
        FZSUM=0.D0
        DO I=1,N
          IF (JFIX(I).EQ.1) FFX(I) = 0.0D0
          IF (JFIX(I).EQ.1) FFY(I) = 0.0D0
          IF (JFIX(I).EQ.1) FFZ(I) = 0.0D0
          AX(I)=AX(I)+FFX(I)/M(I)
          AY(I)=AY(I)+FFY(I)/M(I)
          AZ(I)=AZ(I)+FFZ(I)/M(I)
        ENDDO
      ENDIF
 
      IF (NOPOT) THEN
         WRITE(61,*) 'STOPPING: IPOT = ',IPOT
         STOP
      ENDIF

      IACCEL = IACCEL + 1
      RETURN
      END
C--------------------------------------------------------------------
C
C Subroutine POTENTAL --> determines interatomic forces based on
C                         a Lennard-Jones potential
C
C--------------------------------------------------------------------
      SUBROUTINE POTENTAL(R,FORCE,POT)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M
      DOUBLE PRECISION MORSE_EXP,MORSE_TERM       
C********************CONVERSION FACTORS*******************************

C************CONVERT ENERGY FROM EV TO UNITS OF HARTREE***************
      EPSILON1=3.67502D-2
C*********************************************************************  
C***********************CONVERT 1 EV TO KCAL/MOL**********************
      ENERGY_CONST=23.06D0
C*********************************************************************  
C**CONVERT FORCES IN KCAL/MOL/A TO AMU*ANGSTROM / FEMPTOSECOND**2
C**********************************************************************
      FORCE_CONST=0.0004184D0
C*********************************************************************  
C**CONVERT FORCES IN EV/A TO UNITS OF AMU*ANGSTROM / FEMPTOSECOND**2
C**********************************************************************  
      FCONV=ENERGY_CONST*FORCE_CONST
C********************************************************************   
C      IF(R.GT.2.1) THEN
C	   ALPHA=1.2
C	   Ro=2.2
C	   D=2.08   
C            RI=R-Ro	                
C            MORSE_EXP = DEXP( -ALPHA*RI )
C            MORSE_TERM =  1 - MORSE_EXP   
C  
C 
C*********************************************************************
C Calculate the two-body potential due to I and J
C*********************************************************************
C          
C            POT = (D*MORSE_TERM*MORSE_TERM - 1.942)*EPSILON1
C            FORCE=(-2*D*MORSE_TERM*MORSE_EXP*ALPHA)*FCONV
C	    write(*,*) 'MORSE','R=',R,'FORCE=',FORCE,'POT=',POT
C      ELSE
      EPSILON=0.8
      DELTA=1.961
       
      RI=DELTA/R
      RI3=RI*RI*RI
      RI6=RI3*RI3
      RI12=RI6*RI6
 
      G=24.D0*RI*RI6*(2.D0*RI6-1.D0)
      FORCE=G*RI*EPSILON*FCONV
      POT=4.D0*(RI12-RI6)*EPSILON*EPSILON1
C      	    write(*,*) 'LJ','R=',R,'FORCE=',FORCE,'POT=',POT
C      ENDIF
C      
C      if(abs(FORCE) > 1e-1)write(*,*) 'FORCE=',FORCE,POT,'POT'
      RETURN
      END
      

      SUBROUTINE QUANTUM2EMPIRICAL(NCUR,M,X,Y,Z,SXCUR,SYCUR,SZCUR
     +	,ITYP,ICOL,IPOT,LBOUND,N)	 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M, MO
      PARAMETER (ND = 2000)
      COMMON /DYNFIX/ JFIX(ND),NOTFIX	
      DIMENSION X(ND),Y(ND),Z(ND)
      DIMENSION XO(ND),YO(ND),ZO(ND),ITYPO(ND),MO(ND)
      DIMENSION ITYP(ND),ICOL(ND),M(ND)
	LOGICAL LBOUND
C*******COMPUTE NUMBER OF ORIGINAL - "REAL" ATOMS & SYSTEM BODERS
	N1=NCUR-8*N
  	SX=SXCUR/3
	SY=SYCUR/3
	SZ=0;

	WRITE(61,*) 'N,NCUR,SX,SY,SZ,SXCUR,SYCUR,SZCUR,LBOUND'
	WRITE(61,*)  N,NCUR,SX,SY,SZ,SXCUR,SYCUR,SZCUR,LBOUND	
	WRITE(61,*) '---------------------------------'
	
	NO=0
      DO I=1,N1
C*************************RENORM CORDINATES***************************	
	  IF (ITYP(I).NE.1) THEN
	     NO=NO+1
	     XO(NO)=X(I)
	     YO(NO)=Y(I)
	     ZO(NO)=Z(I)
	     ITYPO(NO)=ITYP(I)
	     MO(NO)=M(I)	
	  ENDIF
      ENDDO 
	
	IF (NO.NE.N) THEN
	  WRITE(61,*) 'ERROR IN QUANTUM2EMPIRICAL'
	  WRITE(61,*) 'INITIAL BOUNDARY NUMBER IS',N,'WHILE COMPUTED number IS',NO	   	
	  STOP
	ENDIF	
     
C*********ADD -SX -SY CLUSTER*****************************************	 
	DO  I=(N1+1),(N1+N)
        M(I) = MO(I-N1)
        XNEW = XO(I-N1)-SX
	  YNEW = YO(I-N1)-SY
	  ZNEW = ZO(I-N1)		
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW	
	  ITYP(I) = ITYPO(I-N1)
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C*********************************************************************

C*********ADD -SX CLUSTER*********************************************	 
	DO  I=(N1+N+1),(N1+2*N)
        M(I) = MO(I-(N1+N))
        XNEW = XO(I-(N1+N))-SX
	  YNEW = YO(I-(N1+N))
	  ZNEW = ZO(I-(N1+N))	
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW
	  ITYP(I) = ITYPO(I-(N1+N))
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C**********************************************************************		

C*********ADD -SX +SY CLUSTER******************************************	 
	DO  I=(N1+2*N+1),(N1+3*N)
        M(I) = MO(I-(N1+2*N))
        XNEW = XO(I-(N1+2*N))-SX
	  YNEW = YO(I-(N1+2*N))+SY
	  ZNEW = ZO(I-(N1+2*N))
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW
	  ITYP(I) = ITYPO(I-(N1+2*N))
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C**********************************************************************		

C*********ADD -SY CLUSTER**********************************************	 
	DO  I=(N1+3*N+1),(N1+4*N)
        M(I) = MO(I-(N1+3*N))
        XNEW = XO(I-(N1+3*N))
	  YNEW = YO(I-(N1+3*N))-SY
	  ZNEW = ZO(I-(N1+3*N))
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW	
	  ITYP(I) = ITYPO(I-(N1+3*N))
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C**********************************************************************	

C*********ADD +SY CLUSTER**********************************************	 
	DO  I=(N1+4*N+1),(N1+5*N)
        M(I) = MO(I-(N1+4*N))
        XNEW = XO(I-(N1+4*N))
	  YNEW = YO(I-(N1+4*N))+SY
	  ZNEW = ZO(I-(N1+4*N))
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW
	  ITYP(I) = ITYPO(I-(N1+4*N))
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C**********************************************************************	

C*********ADD +SX +SY CLUSTER******************************************
	DO  I=(N1+5*N+1),(N1+6*N)
        M(I) = MO(I-(N1+5*N))
	  XNEW = XO(I-(N1+5*N))+SX
	  YNEW = YO(I-(N1+5*N))+SY
	  ZNEW = ZO(I-(N1+5*N))
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW
	  ITYP(I) = ITYPO(I-(N1+5*N))
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C**********************************************************************	

C*********ADD +SX -SY CLUSTER******************************************
	DO  I=(N1+6*N+1),(N1+7*N)
        M(I) = MO(I-(N1+6*N))
	  XNEW = XO(I-(N1+6*N))+SX
	  YNEW = YO(I-(N1+6*N))-SY
	  ZNEW = ZO(I-(N1+6*N))
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW
	  ITYP(I) = ITYPO(I-(N1+6*N))
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C**********************************************************************

C*********ADD +SX CLUSTER**********************************************
	DO  I=(N1+7*N+1),(N1+8*N)
        M(I) = MO(I-(N1+7*N))
	  XNEW = XO(I-(N1+7*N))+SX
	  YNEW = YO(I-(N1+7*N))
	  ZNEW = ZO(I-(N1+7*N))
        IF(LBOUND) CALL BOUNDARY(XNEW,YNEW,ZNEW,SXCUR,SYCUR,SZCUR)
        X(I) = XNEW
	  Y(I) = YNEW
	  Z(I) = ZNEW	
	  ITYP(I) = ITYPO(I-(N1+7*N))
C*********ALL ATOMS ARE FIXED AND ASIGNED TO THE MORSE (15) POTENTIAL*
	  ICOL(I) = 15
	  JFIX(I)	= 1
      ENDDO
C**********************************************************************
      N1=N1+8*N	
	SX=SX*3
	SY=SY*3	
C**********************************************************************	 
      RETURN
      END

C=======================================================================
C H.B. Zhang 
C to see whether SCFCRT is satisfied or not
C=======================================================================
      subroutine testSCFCRAT(LSCFCRT)
      implicit none
      
      logical :: LSCFCRT
      character(44) :: aline, SCFFail
      
      LSCFCRT = .False.
      SCFFail = ' FOR SOME REASON THE SCF CALCULATION FAILED.'
      
      open(unit=4444,file='mopac.out')
      do while (.True.)
         read(4444,'(A)',end=1234) aline
         if( aline .eq. SCFFail ) then
             LSCFCRT = .True.
             goto 1234
         end if
      end do
            
1234  close(4444)      end subroutine testSCFCRAT
