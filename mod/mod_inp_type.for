*******************************************************************************
*                                 MODCON                                      *
*******************************************************************************

C NEXT MINOR CHANGE: ITERATE USING GUSTAFFSON WHEN GETTING PRESS FROM DENS & T
C                   WITHOUT COLMAS !!!!!!!!!!
C NEXT MINOR CHANGE: IMPROVE THE CONDITION FOR DETERMINING P(TAU_0), I.E.THE
C                  TOPMOST POINT OF THE ATMOSPHERE (FROM BOEHM-VITENSE CHAP. 9)
C nEXT MINOR CHENGE: DO NOT DETERMINE RMU IF RHO IS NOT KNOWN, BUT PRESSURE IS,
C                    SINCE RHO WILL GET DETERMINED BY GUSTAFFSON ANYWAY.
C NEXT CHANGE: FOR ATMOSPHERES TABULATED TOO SPARESLY, INTRODUCE AN 
C              INTERPOLATION JUST AFTER READING THE ATMOSPHERE, I.E. SHORTLY
C              AFTER MODINP IS CALLED.
C NEXT MAJOR CHANGE: USE MANY FREQUENCIES TO GET THE ROSSELAND MEAN OPACITY
C                    THIS SHOULD NOT BE TOO DIFFICULT, SINCE PUTTING J=0 INTO
C                    ABSKO GIVES ROSSELAND OPACITY !!!!!!!!!!!!!!!!!!!!!!!
C                    ONE WILL JUST HAVE TO SEE TO IT THAT AT THE END, THE 
C                    KAPPA AT THE PROPER WAVELENGTH IS CALCULATED AGAIN.

C This code reads an input model of any number of height steps and, if wanted, 
C interpolates between them to create a model with equidistant steps in 
C log(tau). If any important atmospheric quantities are missing in the input
C atmosphere, then the code calculates them at all the height points.
C At least the two variables TEMP and PRESS have to be
C specified at input. Depending on the variable, the interpolation is linear
C in TLOG (for Z, TEMP, BFIELD, VTURB, VEL), or exponential (for PEL, PRESS,
C KAPPA, DENS). Depending on what the input variables are, different parts 
C of the code may have to be changed. THE INPUT ROUTINE SHOULD BE CHECKED EVERY 
C TIME. 

C IBOT: The number of height (depth) points of the output atmosphere

	PARAMETER (NDIM=1024)

	REAL      Z(NDIM), VEL(NDIM), PRESS(NDIM), BFIELD(NDIM)
	REAL      DENS(NDIM), TEMP(NDIM), KAPPA(NDIM), NUMDENS(NDIM)
	REAL      TAU(NDIM), TLOG(NDIM), VTURB(NDIM), PEL(NDIM)
	REAL      GAMMA(NDIM), CHII(NDIM), RMU(NDIM), COLMAS(NDIM)
	REAL      DENSLN(NDIM), RHO(NDIM), PEL2(NDIM), KAPPA2(NDIM)
	REAL      TREF(10), RMUREF(10), NEL(NDIM), DUMMY(NDIM)
	REAL      VELHOR(NDIM), VT(NDIM), discard
        real      tlogtp, gr_fac
	CHARACTER MODNAM*55, BOUNDS*1, OUTNAM*25, FOROUT*1, WAVANS*1
	CHARACTER TYPINP*1, DIREC*25, gravnew*1
	LOGICAL   ZGOOD

	COMMON/IBOT/IBOT
	COMMON/OUTNAM/OUTNAM
	COMMON/MODNAM/MODNAM
	COMMON/TBOUND/ TUPPER, TLOWER
	COMMON/NDD/NDD
	COMMON/CONST/GASCON, GRAVAC
	COMMON/NVEL/NVEL, NGAMMA, NCHII

C Commons used by INITAB.
	COMMON/NEWT/NEWT
	COMMON/COUTR/NTO,NTPO(10)

        save vturb, z

	NDD=NDIM

C Set constants

	PI=3.141592653589793238462643
C Gas constant in cgs units: erg/mol/K
	GASCON=8.314510E+07                       
C Solar accelaration due to gravity in  cgs units: cm/sec**2
	GRAVAC=27397                             
C Boltzmann's constant (erg/K)'
        BOLTZK=1.380658E-16 

C Set all input parameters, read filenames and upper and lower bounds 
C from input file. forout normally 'A'. bounds normally 'T'
C 
c        TYPINP='2'
c	FOROUT='A'
	NRAY=1
c        BOUNDS='T'
        WAVANS='Y'
	WAVLEN = 5000

C Read the name of the input model, the number of height points and the
C value of the microturbulence. If the microturbulence is input zero, then
C the code reads the file VTURB.DAT for the microturbulence.

        OPEN (UNIT=31, FILE='input_para_type.dat', STATUS='old')
c        READ (31,1020) MODNAM
c	READ (31,1020) OUTNAM
	READ (31,1020) TYPINP
	READ (31,1020) FOROUT
	READ (31,1020) BOUNDS
	READ (31,*) TUPPER, TLOWER
	IF (TUPPER.GT.TLOWER) THEN
	   THELP=TUPPER
	   TUPPER=TLOWER
	   TLOWER=THELP
	END IF
        READ (31,*) taulog
        READ (31,*) discard
        READ (31,*) GR_FAC
	READ (31,*) NOUT
	READ (31,*) NRAY
        CLOSE (UNIT=31)


c added by C.Norris October '14
	MODNAM = 'inputfile'
c'F3_vhydro_012000_0'
c'testfile2'
	OUTNAM = 'outputfile'
c	TUPPER = -4.5
c	TLOWER = 2.5
c	IF (TUPPER.GT.TLOWER) THEN
c	   THELP=TUPPER
c	   TUPPER=TLOWER
c	   TLOWER=THELP
c	END IF
c
c	taulog = -9
c	VT = 4.1
c	GR_FAC = 1.00
c       1.00 for G2
c       0.73 for F3
c	NOUT = 256
c	NRAY = 512

        GRAVAC = GRAVAC*GR_FAC
        print *, 'Reading 2-d model with ', nray, ' rays '
        print *, 'Interpolation uses ', tupper, tlower, nout

C Beginning of loop over rays, i.e. for 2-D data, where the atmosphere along
C each ray is determined one after the other (and the data file is read one
C ray at a time).

	DO IRAY=1,NRAY

C Read the data from input model atmosphere. Depending on what the input
C variables are this part of the code may have to be changed.
c	   print *, 'at line 135: rho = ', rho(1),rho(40)
	CALL MODINP (TLOG, Z, TEMP, PEL, KAPPA, PRESS, VTURB, 
     1               DENS, BFIELD, VEL, BOUNDS, GAMMA, CHII, 
     2               NEL, NUMDENS, TAU, COLMAS, VELHOR, TYPINP, 
     3               IRAY, NRAY)

        IF (BOUNDS.EQ.'N') THEN
	  NOUT = IBOT
	END IF
C Find out which variables have to be determined, could streamline the following
c if it is only used for muram input and kurucz output as this will slow things c down considerably.. 

	NVTURB=0
	NTLOG=0
	NPEL=0
	NKAPPA=0
	NDENS=0
	NPRESS=0
	NZ=0
	NTAU=0
	NNEL=0
	NCOLMAS=0
	NNUMDENS=0
	NVEL=0
	NGAMMA=0
	NCHII=0
	DO I=1, IBOT
	  IF (TLOG(I).NE.0.0) NTLOG=NTLOG+1
	  IF (PEL(I).NE.0.0) NPEL=NPEL+1
	  IF (KAPPA(I).NE.0.0) NKAPPA=NKAPPA+1
	  IF (DENS(I).NE.0.0) NDENS=NDENS+1
	  IF (NUMDENS(I).NE.0.0) NNUMDENS=NNUMDENS+1
	  IF (PRESS(I).NE.0.0) NPRESS=NPRESS+1
	  IF (VTURB(I).NE.0.0) NVTURB=NVTURB+1
	  IF (COLMAS(I).NE.0.0) NCOLMAS=NCOLMAS+1
	  IF (Z(I).NE.0.0) NZ=NZ+1
	  IF (TAU(I).NE.0.0) NTAU=NTAU+1
	  IF (NEL(I).NE.0.0) NNEL=NNEL+1
	  IF (VEL(I).NE.0.0) NVEL=NVEL+1
	  IF (GAMMA(I).NE.0.0) NGAMMA=NGAMMA+1
	  IF (CHII(I).NE.0.0) NCHII=NCHII+1
	END DO
	NGRAV=0
c	print *, 'at line 180:'
c	print *, 'NVTURB,NTLOG,NPEL,NKAPPA =', NVTURB,NTLOG,NPEL,NKAPPA
c	print *, 'NDENS,NPRESS,NZ,NTAU =', NDENS,NPRESS,NZ,NTAU
c	print *, 'NNEL,NCOLMAS,NNUMDENS =', NNEL,NCOLMAS,NNUMDENS
c	print *, 'NVEL,NGAMMA,NCHII =', NVEL,NGAMMA,NCHII

C Create the electron pressure if the electron number density is known.
C or vice versa.. (for Kurucz output, NEL is needed, after all)

	IF ((NPEL.EQ.0).AND.(NNEL.NE.0)) THEN
	  DO I=1, IBOT
	    PEL(I)=NEL(I)*BOLTZK*TEMP(I)
	    NPEL=NPEL+1
	  END DO
	  PRINT*,'ELECTRON PRESSURE CREATED FROM ELECTRON DENSITY'
	ELSE IF ((NPEL.NE.0).AND.(NNEL.EQ.0)) THEN
	  DO I=1, IBOT
	    NEL(I)=PEL(I)/(BOLTZK*TEMP(I))
	    NNEL=NNEL+1
	  END DO
	   PRINT*,'ELECTRON DENSITY CREATED FROM ELECTRON PRESSURE'
	END IF

C Create the gas pressure if the total number density is known.

	IF ((NPRESS.EQ.0).AND.(NNUMDENS.NE.0)) THEN
	  DO I=1, IBOT
	    PRESS(I)=NUMDENS(I)*BOLTZK*TEMP(I)
	    NPRESS=NPRESS+1
	  END DO
	  PRINT*,'GAS PRESSURE CREATED FROM TOTAL NUMBER DENSITY'
	END IF

C Initialize the abscont package if necessary

	IF ((NPEL.EQ.0).OR.(NKAPPA.EQ.0).OR.(NDENS.EQ.0)
     1     .OR.(WAVANS.EQ.'Y').OR.(NPRESS.EQ.0)) 
     2     THEN
	  LREAD=9
	  LWRITE=10
	  LSCR1=98
	  LSCR2=99
	  IOUTS=0
	  IJ=1
60	  IF (IJ.EQ.1) THEN
	    DIREC='RAD:'
	    NDIREC=13
	    IJ=IJ+1
c	    print *, ij
c	    print *, ' Trying ', DIREC(1:NDIREC)//'modconinp.dat'
	  ELSE IF (IJ.EQ.2) THEN
	    DIREC=' '
	    NDIREC = 0
c	    NDIREC = 1
	    IJ=IJ+1
c	    print *, ij
c            print *, ' Trying ', 'modconinp.dat'
	  ELSE IF (IJ.EQ.3) THEN
	    DIREC='INV:'
	    NDIREC=1
	    IJ=IJ+1
	    print *, ij
            print *, ' Trying ', DIREC(1:NDIREC)//'modconinp.dat'
	  ELSE
	    PRINT*,'modconinp.dat WAS NOT FOUND ANYWHERE'
	    STOP
	  END IF
	  OPEN (UNIT=LREAD,FILE=DIREC(1:NDIREC)//'modconinp.dat', 
     1         STATUS='OLD',ERR=60)
C          IF (IOUTS.EQ.1) OPEN (UNIT=LWRITE, FILE='modconout.dat',        !VMS
C     1         STATUS='NEW',SHARED)                                       !VMS
          IF (IOUTS.EQ.1) OPEN (UNIT=LWRITE, FILE='modconout.dat',       !UNIX
     1         STATUS='NEW')                                      !UNIX
	  OPEN (UNIT=LSCR1,STATUS='SCRATCH', FORM='UNFORMATTED')
	  OPEN (UNIT=LSCR2,STATUS='SCRATCH', FORM='UNFORMATTED')
	  CALL INITAB(LREAD,LWRITE,LSCR1,LSCR2,IOUTS)
	  NEWT=2
	  NTO=0
	  NTPO(1)=1
	END IF

C Create the microturbulence velocity if it is not known at input.
	IF ((NVTURB.EQ.0).AND.((FOROUT.EQ.'S').OR.(FOROUT.EQ.'A'))) THEN
	  DO I=1, IBOT
	    VTURB(I)=VT(I)*1.0E+5         ! [cm/sec]
	    NVTURB=NVTURB+1
	  END DO
	END IF
	
c	print *, 'at line 265', z(1), temp(1), press(1), dens(1)
C If gas law is to be used to derive density or pressure, then determine the 
C proper RMU values by interpolating in the table of RMU(TEMP). Note that this
C table is only valid for the sun! It has been taken from the HSRA.
C Although superseded for many purposes, the following is still useful (e.g. to
C obtain the starting point of an iteration).

C NOTE THAT IT IS MORE EXACT TO USE DENSITY DETERMINED FROM ABSKO IF PRESSURE
C IS KNOWN !!!!!!!!!!!!!!!!!!!

	IF ((NPRESS.EQ.0).OR.(NDENS.EQ.0)) THEN
	  TREF(1)=4000.
	  RMUREF(1)=1.3
	  TREF(2)=5000.
	  RMUREF(2)=1.298
	  TREF(3)=6000.
	  RMUREF(3)=1.296
	  TREF(4)=7000.
	  RMUREF(4)=1.293
	  TREF(5)=8000.
	  RMUREF(5)=1.285
	  TREF(6)=9000.
	  RMUREF(6)=1.27
	  TREF(7)=10000.
	  RMUREF(7)=1.24
	  TREF(8)=11000.
	  RMUREF(8)=1.2
	  DO I=1, IBOT
	    DO J=1, 8
	      IF (TEMP(I).GE.TREF(8)) THEN
		RMU(I)=RMUREF(7)+(RMUREF(8)-RMUREF(7))
     1                *(TEMP(I)-TREF(7))/(TREF(8)-TREF(7))
	      ELSE IF ((J.EQ.1).AND.(TEMP(I).LE.TREF(1))) THEN
	        RMU(I)=RMUREF(1)
	      ELSE IF ((TEMP(I).LT.TREF(J)).AND.(TEMP(I).GE.TREF(J-1)))
     1                 THEN
		RMU(I)=RMUREF(J-1)+(RMUREF(J)-RMUREF(J-1))
     1                *(TEMP(I)-TREF(J-1))/(TREF(J)-TREF(J-1))
	      END IF
	    END DO
	  END DO
	END IF

C Determine TLOG if it is not known, but TAU is known.

	ITERBIG=0
300	IF ((NTAU.NE.0).AND.(NTLOG.EQ.0)) THEN
	  DO I=1, IBOT
	    IF (TAU(I).EQ.0.0) TAU(I)=1.02E-10
	    TLOG(I)=LOG10(TAU(I))
	    NTLOG=NTLOG+1
	  END DO
	  PRINT*,'LOG(TAU) CREATED FROM TAU'
	END IF

C Determine the height scale from the column mass and the gas density

	IF ((NZ.EQ.0).AND.(NCOLMAS.NE.0).AND.(NDENS.NE.0)) THEN
	  Z0=-999.
	  Z(1)=0.0
	  DO I=2, IBOT
	    Z(I)=Z(I-1)-(COLMAS(I)-COLMAS(I-1))*2./(DENS(I)+DENS(I-1))
	  END DO
	  IF (NTLOG.NE.0) THEN
	    I=1
	    DO WHILE (TLOG(I).LT.0.0)
	      I=I+1
	    END DO
	    IF (I.EQ.1) THEN
	      PRINT*,'PROBLEMS: TLOG(I) IS NEVER NEGATIVE'
	      DO J=1, IBOT
	        PRINT*,'TLOG,Z',J,TLOG(J),Z(J)
	      END DO
	    END IF
	    DELZ=Z(I)+(Z(I)-Z(I-1))*TLOG(I)/(TLOG(I-1)-TLOG(I))
	  ELSE 
	    DELZ=0.0
	  END IF
	  DO I=1, IBOT
	    Z(I)=Z(I)-DELZ
	    NZ=NZ+1
	  END DO
	  PRINT*,'HEIGHT SCALE CREATED FROM COLMAS & DENSITY'
	END IF

C Convert depth scale into height scale

	IF (Z(2).GT.Z(1)) THEN
	  DO I=1,IBOT
	    Z(I)=-Z(I)
	  END DO
	END IF

C Determine P(tau), PEL(tau), KAPPA(tau) and if need DENS(tau) if only T(tau) 
C is given, using the method proposed by Boehm Vitense.

	IF ((NPRESS.EQ.0).AND.(NDENS.EQ.0).AND.(NZ.EQ.0).AND.
     1      (NKAPPA.EQ.0).AND.(NTLOG.NE.0)) THEN
	  IF (NGRAV.EQ.0) THEN
	    PRINT*,'IS NON-SOLAR GRAVITATIONAL ACCELARATION REQUIRED ',
     1             '(Y/N, DEF=N)'
	    READ 1020,GRAVNEW
	    IF (GRAVNEW.EQ.'Y') THEN
	      PRINT*,'ENTER NEW LOG(G) VALUE (SOLAR VALUE=4.47)'
	      READ *, GRAVAC
	      GRAVAC=10.**GRAVAC
	    END IF
	  END IF

          CALL PTAUCR (TEMP, TLOG, PRESS, PEL, DENS, KAPPA, NPRESS, 
     1                 NPEL, NDENS, NKAPPA, NDD)
	  PRINT*,'PRESSURE, DENSITY,PEL, KAPPA CREATED FROM TEMP & TLOG'
	END IF
c	print *, 'at line 378', z(1), temp(1), press(1), dens(1)
C THE FOLLOWING ROUTES ARE CONSIDERABLY WORSE THAN THE ONE USED ABOVE.
CC For the case when only T(tau) is input, determine z(tau), T(z) and P(z)
CC iteratively within subr PZCREA
C
C	IF ((NPRESS.EQ.0).AND.(NDENS.EQ.0).AND.(NZ.EQ.0).AND.
C     1      (NKAPPA.EQ.0).AND.(NTLOG.NE.0)) 
C	  IF (NGRAV.EQ.0) THEN
C	    PRINT*,'NON-SOLAR GRAVITATIONAL ACCELARATION REQUIRED ',
C     1             '(Y/N, DEF=N)'
C	    READ 9000,GRAVNEW
C	    IF (GRAVNEW.EQ.'Y') THEN
C	      PRINT*,'ENTER NEW LOG(G) VALUE (SOLAR VALUE=4.47)'
C	      READ *, GRAVAC
C             GRAVAC=10.**GRAVAC
C	    END IF
C	  END IF
C     2      CALL PZCREA (TEMP, TLOG, RMU, Z, PRESS, PEL, KAPPA, DENS, 
C     3                   NZ, NPRESS, NPEL, NKAPPA, NDENS, NDD)
C
CC Alternative route. Work via column mass density (COLMAS). This leads to a 
CC simpler equation for the pressure (cf. Mihalas, 1978, Stellar Atmospheres, 
CC p. 170
C
C	IF ((NPRESS.EQ.0).AND.(NDENS.EQ.0).AND.(NZ.EQ.0).AND.
C     1      (NKAPPA.EQ.0).AND.(NTLOG.NE.0)) 
C	  IF (NGRAV.EQ.0) THEN
C	    PRINT*,'NON-SOLAR GRAVITATIONAL ACCELARATION REQUIRED ',
C     1             '(Y/N, DEF=N)'
C	    READ 9000,GRAVNEW
C	    IF (GRAVNEW.EQ.'Y') THEN
C	      PRINT*,'ENTER NEW LOG(G) VALUE (SOLAR VALUE=4.47)'
C	      READ *, GRAVAC
C             GRAVAC=10.**GRAVAC
C	    END IF
C	  END IF
C     2      CALL CPCREA (TEMP, TLOG, PRESS, PEL, KAPPA, 
C     3                   NPRESS, NPEL, NKAPPA, NDD)

C Create pressure if both pressure and density are unknown, but T and z are 
C known.
C THIS PART WILL HAVE TO BE CHANGED TO TAKE INTO ACCOUNT THE FACT THAT THE
C DENS IS BETTER DETERMINED FROM PEABS AND RMU IS THEN BETTER OBTAINED FROM 
C THE NEW DENS AND PRESSURE, ITERATING UNTIL CONVERGENCE.

	IF ((NPRESS.EQ.0).AND.(NDENS.EQ.0).AND.(NZ.NE.0)) THEN
	  IF ((NKAPPA.EQ.0).AND.(NTAU.EQ.0).AND.(NCOLMAS.EQ.0)) THEN
	    PRINT*,'ENTER DESIRED PRESSURE (CGS) AT Z=0. IF NO ',
     1             'GRID POINT LIES AT Z=0, ENTER P AT' 
	    PRINT*,'FIRST GRID POINT JUST BELOW Z=0.  IF INPUT ',
     1             'P(Z=0)=0, THEN P(Z=0)=1.31E5 IS USED'
	    PRINT*,'IF Z=0 IS NOT AMONG THE LISTED POINTS, ENTER P ',
     1             'AT LOWEST LISTED POINT'
	    READ*, PRESS0
	    IF (PRESS0.EQ.0.0) PRESS0=1.31E5
	    ZGOOD=.TRUE.
	    CALL HYDST (TEMP, PRESS, PRESS0, Z, RMU, ZGOOD, ICRIT, 
     1                  NDD)
	    NPRESS=IBOT
	    PRINT*,'PRESSURE CREATED FROM Z AND TEMPERATURE'
	  ELSE IF (((NKAPPA.NE.0).OR.(NTAU.NE.0)).AND.(NCOLMAS.EQ.0)) THEN
	    PRINT*,'SINCE KAPPA(Z) AND/OR TAU(Z) ARE KNOWN TOGETHER ',
     1             'WITH T(Z), THERE IS NO'
	    PRINT*,'FREEDOM IN THE CHOICE OF P(Z_0). AT PRESENT THE ',
     1             'CODE CANNOT HANDLE THIS CASE'
	    STOP
	  END IF
	END IF

cC Determine the column mass density for known z and density (only required 
cC for some purposes)
c
c	print *, 'NCOLMAS, NZ, NDENS ', NCOLMAS, NZ, NDENS
c	IF ((NCOLMAS.EQ.0).AND.(NZ.NE.0).AND.(NDENS.NE.0)) THEN
c	  COLMAS(1)=2.E-4
c	  COLMAS(1)=8.E-4
c	  NCOLMAS=1
c	  DO I=2, IBOT
c	    COLMAS(I)=COLMAS(I-1)+(DENS(I-1)+DENS(I))/2.*(Z(I-1)-Z(I))
c	    NCOLMAS=NCOLMAS+1
c	  END DO
c	END IF
cC THE COLMAS AT THE UPPER END OF THE ATMOSPHERE (AS DETERMINED ABOVE) CAN BE 
cC IMPROVED, E.G. BY EXTRAPOLATING FROM BELOW !!!!!!!!!!!!!!!!!!!!!

C Call the ABSCONT routine package to determine the density from Temp and 
C Pressure. This means that
C PEL and the absorption coefficient are also determined. If not already known
C they are stored in the relevant variables.
c	print *,'at line 474 rho = ',rho(1),rho(40)
	IF ((NDENS.EQ.0).AND.(NPRESS.NE.0)) THEN
	  NPEL2=NPEL
	  NKAPPA2=0
	  IF (NPEL.NE.0) THEN
	    DO I=1, IBOT
	      PEL2(I)=PEL(I)
	      KAPPA2(I)=0.0
	    END DO
	  END IF

C	  IF ((WAVANS.EQ.'N').AND.(WAVLEN.EQ.0.0)) WAVLEN=5000.0
	  IF (WAVLEN.EQ.0.0) WAVLEN=5000.0
	  CALL PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL2, KAPPA2, RHO, 
     1                NPEL2, NKAPPA2, 1, IBOT, NDD)
	  DO I=1, IBOT
	    DENS(I)=RHO(I)
	    NDENS=NDENS+1
	  END DO
	  PRINT*,'GAS DENSITY CREATED FROM PRESSURE AND TEMPERATURE ',
     1           'USING ABSKO'
	  IF (NPEL.EQ.0) THEN
	    DO I=1, IBOT
	      PEL(I)=PEL2(I)
	    END DO
	    NPEL=IBOT
	    PRINT*,'ELECTRON PRESSURE CREATED FROM TEMP. & GAS PRESSURE'
	  END IF
	  IF (NKAPPA.EQ.0) THEN
	    DO I=1, IBOT
	      KAPPA(I)=KAPPA2(I)
	    END DO
	    NKAPPA=IBOT
	    PRINT*,'ABS COEFF CREATED FROM ELECTRON PRESSURE & TEMP.'
	  END IF
	END IF
C Determine the column mass density for known z and density (only required
C for some purposes)  !!! THIS IS WRONG,I think!!! thermal pressure????

        IF ((NCOLMAS.EQ.0)) then
          if ((NZ.NE.0).AND.(NDENS.NE.0)) THEN
            if (npress .ne. 0) then
              colmas(1) = press(1) / gravac
            else 
              COLMAS(1)=2.E-4
            end if
            NCOLMAS=1
            DO I=2, IBOT
              COLMAS(I)=COLMAS(I-1)+(DENS(I-1)+DENS(I))/2.*(Z(I-1)-Z(I))
              NCOLMAS=NCOLMAS+1
            END DO
c	    print *, 'COLMAS CALCULATED FROM GAS PRESSURE & Z'
	   end if
        END IF
c	print *, 'At line 526:'
c	print *, 'z,temp,press,dens = ',z(1), temp(1), press(1), dens(1)
c	print *, '1: colmas, rho, dummy = ', colmas(1), rho(1), dummy(1)
c	print *, '40:colmas, rho, dummy =', colmas(40), rho(40), dummy(40)
c	print *, '1: pel2, pel = ',  pel2(1), pel(1)
c	print *, '40:pel2, pel = ',  pel2(40), pel(40)

C THE COLMAS AT THE UPPER END OF THE ATMOSPHERE (AS DETERMINED ABOVE) CAN BE
C IMPROVED, E.G. BY EXTRAPOLATING FROM BELOW !!!!!!!!!!!!!!!!!!!!!

C Determine gas density from the column mass (COLMAS). This is slightly 
C complicated, since it requires a differentiation, so that accuracy is lost. 
C To regain this DENS is smoothed.
	IF ((NDENS.EQ.0).AND.(NPRESS.EQ.0).AND.(NCOLMAS.NE.0).AND.
     1           (NZ.NE.0)) THEN
	  DO I=2, IBOT-1
	    DENS(I)=ABS((COLMAS(I-1)-COLMAS(I+1))/(Z(I-1)-Z(I+1)))
	  END DO
C In the following is an alternative route in which the density DENS is simply
C smoothed directly using running means, but this gives large errors 
C	  DENS(2)=(DENS(2)+DENS(3))/2.
C	  DENS(IBOT-1)=(DENS(IBOT-2)+DENS(IBOT-1))/2.
C	  DO I=3, IBOT-2
C	    DENS(I)=(DENS(I-1)+DENS(I)+DENS(I+1))/3.
C	  END DO
C	  DENS(1)=DENS(2)+(DENS(2)-DENS(3))/(Z(2)-Z(3))*(Z(1)-Z(2))
C	  DENS(IBOT)=DENS(IBOT-1)+(DENS(IBOT-1)-DENS(IBOT-2))
C     1             /(Z(IBOT-1)-Z(IBOT-2))*(Z(IBOT)-Z(IBOT-1))
C So we try smoothing the logarithmic density. 
	  DO I=2, IBOT-1
	    DENSLN(I)=LOG(DENS(I))
	  END DO
	  DENSLN(2)=(DENSLN(2)+DENSLN(3))/2.
	  DENSLN(IBOT-1)=(DENSLN(IBOT-2)+DENSLN(IBOT-1))/2.
	  DO I=3, IBOT-2
	    DENSLN(I)=(DENSLN(I-1)+DENSLN(I)+DENSLN(I+1))/3.
	  END DO
	  DENSLN(1)=DENSLN(2)+(DENSLN(2)-DENSLN(3))/(Z(2)-Z(3))
     1             *(Z(1)-Z(2))
	  DENSLN(IBOT)=DENSLN(IBOT-1)+(DENSLN(IBOT-1)-DENSLN(IBOT-2))
     1             /(Z(IBOT-1)-Z(IBOT-2))*(Z(IBOT)-Z(IBOT-1))
	  DO I=1, IBOT
	    DENS(I)=EXP(DENSLN(I))
	    NDENS=NDENS+1
	  END DO
	  PRINT*,'DENSITY CREATED FROM COLUMN MASS AND Z'
	ELSE IF ((NDENS.EQ.0).AND.(NPRESS.EQ.0).AND.((NCOLMAS.EQ.0).OR.
     1           (NZ.EQ.0))) THEN
	  IF (NCOLMAS.NE.0) THEN
C Determine the pressure by integrating (trivially) the 
C Eq. of hydrostat. Equil.
C To determine the const. of Integr. the pressure at some height must be
C input by the user!
	    DO I=1, IBOT
	      PRESS(I)=GRAVAC*COLMAS(I)
	    END DO
	    PRINT*,'ENTER DEPTH INDEX FOR WHICH YOU WANT TO ENTER GAS ',
     1             'PRESSURE'
	    READ*,INDEX
	    PRINT*,'ENTER PRESSURE AT INDEX:',INDEX
	    READ*,PRESS0
	    IF (PRESS0.EQ.0.) PRESS0=1.31E5
	    PRESDIF=PRESS0-PRESS(INDEX)
	    DO I=1, IBOT
	      PRESS(I)=PRESS(I)+PRESDIF
	      IF (PRESS(I).LE.0.0) THEN
	        PRINT*,'PRESSURE AT LEVEL',I,' IS ZERO OR NEGATIVE:',
     1                  PRESS(I)
	        STOP
	      END IF
	    END DO
	    NPRESS=IBOT
	    PRINT*,'P CREATED FROM COL MASS!'
	  ELSE
	    PRINT*,'DENSITY AND PRESSURE ARE NOT KNOWN AND NO ',
     1             'COLUMN MASS IS KNOWN'
	    PRINT*,'INSUFFICIENT INFORMATION TO CALCULATE THE UNKNOWN ',
     1             'ATMOSPHERIC VARIABLES'
	    STOP
	  END IF
	END IF

C Create the gas pressure if not known at input (density and temperature req.)

	IF ((NPRESS.EQ.0).AND.(NDENS.NE.0)) THEN
	  IF (NCOLMAS.NE.0) THEN
C Determine the pressure by integrating (trivially) the Eq. of hydrostat. Equil.
C To determine the const. of Integr., get rho from Gustaffson and compare with
C (correct) density read from file. Iterate until convergence.
	    DO I=1, IBOT
	      PRESS(I)=GRAVAC*COLMAS(I)
	      IF (DENS(I).LE.5.*DENS(1)) II=I
	    END DO
	    IF (II.LT.3) II=3
	    PRINT*,'II=',II,' COLMAS MUST HAVE CONVERGED TO A ',
     1             'REASONABLE VALUE BY THIS MANY Z STEPS'
	    PRESS1=DENS(II)*TEMP(II)*GASCON/RMU(II)
	    CONST1=PRESS1-PRESS(II)
	    CALL PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL2, KAPPA2, 
     1                  RHO, 0, 0, II, II, NDD)
C	PRINT*,'II,TLOG,DENS,RHO',II,TLOG(II),DENS(II),RHO(II)
400	    IF (ABS((RHO(II)-DENS(II))/DENS(II)).GT.0.001) THEN
	      CONST2=DENS(II)/RHO(II)*(GRAVAC*COLMAS(II)+CONST1)
     1              -GRAVAC*COLMAS(II)
	      PRESS2=PRESS(II)+CONST2
	      DUMMY(II)=PRESS2
	      CALL PEABS (TEMP, DUMMY, WAVANS, WAVLEN, PEL2, KAPPA2, 
     1                    RHO, 0, 0, II, II, NDD)
	      CONST1=CONST2
C	PRINT*,'DENS,RHO,CONST2', DENS(II),RHO(II),CONST2
	      GOTO 400
	    END IF
	    DO I=1, IBOT
	      PRESS(I)=GRAVAC*COLMAS(I)+CONST2
	    END DO
	    NPRESS=IBOT
	    PRINT*,'GAS PRESSURE CREATED FROM COLMAS, DENSITY & TEMP'
	    IF ((NPEL.EQ.0).AND.(NNEL.EQ.0)) THEN
	      DO I=1, IBOT
	        PEL(I)=PEL2(I)
	      END DO
	      NPEL=IBOT
	      PRINT*,'EL. PRESS. CREATED FROM PRESS & TEMP'
	    END IF
	  ELSE
C Determine pressure if colmas is not known, but dens and temp are. 
C HERE ALSO ITERATE TO GET A BETTER VALUE OF PRESS: 1ST DETERMINE NEW 
C DENS USING THE NEW PRESSURE (FROM PEABS) AND CHANGE THE PRESSURE UNTIL THE 
C DERIVED DENS EQUALS THE ORIGINAL INPUT VALUE !!!!!!!!!!!!!!!!!!!!!
C SEE HOW THINGS HAVE BEEN SOLVED JUST ABOVE !!!!!!!!!!!!!!!!!!!!!!

	    DO I=1,IBOT
	      PRESS(I)=DENS(I)*TEMP(I)*GASCON/RMU(I)
	      NPRESS=NPRESS+1
	    END DO
	    PRINT*,'GAS PRESSURE CREATED FROM DENSITY AND TEMPERATURE'
	  END IF
	END IF
c	print *, 'At line 654:'
c	print *, 'z,temp,press,dens = ',z(1), temp(1), press(1), dens(1)
c	print *, '1: colmas, rho, dummy = ', colmas(1), rho(1), dummy(1)
c	print *, '40:colmas, rho, dummy =', colmas(40), rho(40), dummy(40)
c	print *, '1: pel2, pel = ',  pel2(1), pel(1)
c	print *, '40:pel2, pel = ',  pel2(40), pel(40)
	
C Alternative technique for determining the gas pressure from density by 
C integrating the eq. of hydrostatic equilibrium for the case of known density. 
C
C	IF ((PRESS(1).EQ.0.0).AND.(DENS(1).NE.0.0)) THEN
C	  IF (RMU(IBOT).EQ.0.0) RMU(IBOT)=1.28  !!!too simple ????
C	  PRESS(IBOT)=DENS(IBOT)*TEMP(IBOT)*GASCON/RMU(IBOT)
C	  DO I=IBOT-1,1, -1
C	    PRESS(I)=PRESS(I+1)-GRAVAC*(Z(I+1)-Z(I))
C     1              *(DENS(I)+DENS(I+1))/2.
C	  END DO
C	END IF

C Create the opacity if both Z and TLOG are known.

	IF ((NKAPPA.EQ.0).AND.(NZ.NE.0).AND.(NTLOG.NE.0)) THEN
	  KAPPA(1)=(10.**TLOG(2)-10.**TLOG(1))/ABS(Z(2)-Z(1))/DENS(1)
	  NKAPPA=NKAPPA+1
	  DO I=2, IBOT-1
	    KAPPA(I)=(10.**TLOG(I+1)-10.**TLOG(I-1))/ABS(Z(I+1)-Z(I-1))
     1              /DENS(I)
	    NKAPPA=NKAPPA+1
	  END DO
	  KAPPA(IBOT)=ABS(10.**TLOG(IBOT)-10.**TLOG(IBOT-1))
     1               /ABS(Z(IBOT)-Z(IBOT-1))/DENS(IBOT)
	  NKAPPA=NKAPPA+1
	  PRINT *,'KAPPA CREATED FROM Z AND TLOG'
	END IF

C Call the ABSCONT routine package to calculate the electron pressure and
C the absorption coefficient if not already known. Loop over temperature and 
C gas pressure. 
C Now in principle ABSCONT can also give the density if not already known.
C However, it is already determined from ABSCONT further up and this second
C line of defense has so far not been needed, so that it is commented out
C again and would require further programming to implement here.
c	print *, 'at line 706:'
c	print *, 'rho = ', rho(1), rho(40)
c	print *, '1: pel2, pel = ',  pel2(1), pel(1)
c	print *, '40:pel2, pel = ',  pel2(40), pel(40)
	IF ((NPEL.EQ.0).OR.(NKAPPA.EQ.0).OR.(WAVANS.EQ.'Y').OR.
     1    (NDENS.EQ.0)) THEN
C	  NKAPPA2=NKAPPA
C	  IF (NDENS.EQ.0) NKAPPA2=0
C	  IF ((WAVANS.EQ.'N').AND.(NKAPPA2.EQ.0)) WAVLEN=5000.0
	IF ((WAVANS.EQ.'N').AND.(NKAPPA.EQ.0)) WAVLEN=5000.0
	NKAPPA=0
	  CALL PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL, KAPPA, RHO, 
     1                NPEL, NKAPPA, 1, IBOT, NDD)
c	  print *, 'PEABS run at line 716', rho(1), rho(40)
	  IF (NPEL.EQ.0) THEN
	    NPEL=IBOT
c	    PRINT*,'ELECT PRESSURE CREATED FROM TEMPERATURE & GAS PRESS'
	  END IF
	  IF (NKAPPA.EQ.0) THEN
C	    DO I=1, IBOT
C	      KAPPA(I)=KAPPA2(I)
C	    END DO
	    NKAPPA=IBOT
c	    PRINT*,'ABS COEFF CREATED FROM ELECT PRESSURE & TEMPERATURE'
	  END IF
          IF (NNEL.EQ.0) THEN
            DO I=1, IBOT
              NEL(I)=PEL(I)/(BOLTZK*TEMP(I))
              NNEL=NNEL+1
            END DO
c            PRINT*,'ELECTRON DENSITY CREATED FROM ELECTRON PRESSURE' 
          END IF
c	  print *, 'ibot = ', ibot
C	  IF (NDENS.EQ.0) THEN
C	    DO I=1,IBOT
C	      DENS(I)=RHO(I)
C	    END DO
C	    NDENS=IBOT
C	    PRINT*,'DENSITY CREATED FROM PRESS & TEMPERAT (USING ABSKO)'
C	  END IF
	END IF
c	print *, 'at line 745', z(1), temp(1), press(1), dens(1)
c	print *, 'At line 745:'
c	print *, '1:z,temp,press,dens = ',z(1), temp(1), press(1), dens(1)
c	print *, '40:z,temp,press,dens=',z(40),temp(40),press(40),dens(40)
c	print *, '1: colmas, rho, dummy = ', colmas(1), rho(1), dummy(1)
c	print *, '40:colmas, rho, dummy =', colmas(40), rho(40), dummy(40)
c	print *, 'pel(1-5)', pel(1),pel(2),pel(3),pel(4),pel(5)
c	print *, 'pel(6-10)', pel(6),pel(7),pel(8),pel(9),pel(10)
c	print *, '1: pel2, pel = ',  pel2(1), pel(1)
c	print *, '10:pel2, pel = ',  pel2(10), pel(10)
c	print *, '35:pel2, pel = ',  pel2(35), pel(35)
c	print *, '38:pel2, pel = ',  pel2(38), pel(38)
c	print *, '39:pel2, pel = ',  pel2(39), pel(39)
c	print *, '40:pel2, pel = ',  pel2(40), pel(40)
c, rho(1), dummy(1), pel2(1), pel(1)
C Close the files created by the ABSCONT package

	CLOSE (UNIT=LREAD)
	IF (IOUTS.EQ.1) CLOSE (UNIT=LWRITE)
	CLOSE (UNIT=LSCR1)
	CLOSE (UNIT=LSCR2)

C Determine Z (height scale from known optical depth, absorp. coeff. and density
C Also shift the Z scale such that Z=0 corresponds to TLOG=0.

	IF ((NZ.EQ.0).AND.(NTLOG.NE.0).AND.(NDENS.NE.0).AND.
     1     (NKAPPA.NE.0)) THEN
	  Z(1)=0.0
	  TAU(1)=10**TLOG(1)
	  DO I=2, IBOT
	    TAU(I)=10**TLOG(I)
	    Z(I)=Z(I-1)-2./(KAPPA(I)*DENS(I)+KAPPA(I-1)*DENS(I-1))
     1          *ABS(TAU(I)-TAU(I-1))
	    NZ=NZ+1
	  END DO
	  I=1
	  DO WHILE (TLOG(I).LT.0.0)
	    I=I+1
	  END DO
	  DELZ=Z(I)+(Z(I)-Z(I-1))*TLOG(I)/(TLOG(I-1)-TLOG(I))
	  DO I=1, IBOT
	    Z(I)=Z(I)-DELZ
	  END DO
	  PRINT*,'Z SCALE CREATED FROM TLOG, DENSITY AND ABSORPT COEFF'
	END IF

C Integrate KAPPA to get the optical depth TAU.

!	IF ((((NTLOG.EQ.0).OR.(WAVANS.EQ.'Y')).AND.(NKAPPA.NE.0).AND.
!     1        (NDENS.NE.0).AND.(NZ.NE.0)).AND.(FOROUT.NE.'M')) THEN
	IF ((((NTLOG.EQ.0).OR.(WAVANS.EQ.'Y')).AND.(NKAPPA.NE.0).AND.
     1        (NDENS.NE.0).AND.(NZ.NE.0))) THEN
	  IF ((NZ.NE.0).AND.(NTLOG.EQ.0)) THEN
            if (iray .eq. 1) then
              TLOG(1) = taulog
	      IF ((TLOG(1).LT.-10.).OR.(TLOG(1).GE.1.)) THEN
	        PRINT*,'ILLEGAL VALUE OF TLOG(1)'
	        STOP
	      END IF
              tlogtp = tlog(1)
            else
              tlog(1) = tlogtp
            end if 
	  END IF
	  TAU(1)=10.**TLOG(1)
	  DO I=2,IBOT
	    TAU(I)=TAU(I-1)+(KAPPA(I)*DENS(I)+KAPPA(I-1)*DENS(I-1))/2.
     1            *ABS(Z(I)-Z(I-1))
	    TLOG(I)=LOG10(TAU(I))
	    NTLOG=NTLOG+1
	  END DO
c	  PRINT*,'TLOG CREATED FROM Z, DENSITY & ABSORPTION COEFF'
	ELSE IF ((NTLOG.EQ.0).AND.((NKAPPA.EQ.0).OR.(NDENS.EQ.0).OR.
     1          (NZ.EQ.0))) THEN
	  PRINT*,'SORRY, CANNOT CALCULATE NEW TLOG SCALE, SINCE NO ',
     1           'HEIGHT SCALE IS AVAILABLE!!!!!!'
	END IF

C If Z value was determined from the column mass density at a time when TLOG was
C not yet known, then shift the Z-scale now so that Z=0 corresponds to TLOG=0.

	IF ((NTLOG.NE.0).AND.(NZ.NE.0).AND.(Z0.EQ.-999.)) THEN
	  I=1
	  DO WHILE (TLOG(I).LT.0.0)
	    I=I+1
	  END DO
	  IF ((I.EQ.1).OR.(I.GT.IBOT)) THEN
	    PRINT*,'PROBLEMS: TLOG(I) IS NEVER OR ALWAYS NEGATIVE'
	    DO J=1, IBOT
	      PRINT*,'TLOG,Z',J,TLOG(J),Z(J)
	    END DO
	  END IF
	  DELZ=Z(I)+(Z(I)-Z(I-1))*TLOG(I)/(TLOG(I-1)-TLOG(I))
	  DO I=1, IBOT
	    Z(I)=Z(I)-DELZ
	  END DO
	END IF

C Check if all the variables required for the output have been calculated. If 
C this is not the case then go back to label 300 to run through the code again
C and determine the remaining the variables.

	IF (((FOROUT.EQ.'S').AND.(NVTURB.EQ.0)).OR.((FOROUT.NE.'M').AND.
     1      (NTLOG.EQ.0)).OR.(NPEL.EQ.0).OR.(NKAPPA.EQ.0).OR.
     2      (NDENS.EQ.0).OR.(NPRESS.EQ.0).OR.(NZ.EQ.0)) THEN
	  PRINT*,('*', I=1,78)
	  PRINT*,'NOT ALL VARIABLES DETERMINED IN THIS SINGLE RUN. ',
     1           'GOING THROUGH CODE AGAIN'
	  PRINT*,('*', i=1, 78)
	  ITERBIG=ITERBIG+1
	  IF (ITERBIG.LT.4) THEN
	    GOTO 300
	  ELSE
	    PRINT*,'EVEN AFTER THREE FULL ITERATIONS NOT ALL ',
     1             'QUANTITIES ARE KNOWN'
	    PRINT*,'NTLOG,NPEL,NKAPPA,NDENS,NPRESS,NZ,NVTURB'
            PRINT*, NTLOG,NPEL,NKAPPA,NDENS,NPRESS,NZ,NVTURB
	    STOP
	  END IF
	END IF

C Interpolate the results linearly in order to create an output model with
C the number of height points input above (NOUT)
c	print *,'at line 834', z(1), temp(1), press(1), dens(1)
	IF (BOUNDS.NE.'N') THEN
          IF (BOUNDS.EQ.'T') THEN
c	    IF (FOROUT.EQ.'A') THEN
c	      atmlow=alog10(colmas(ibot))
c	    ELSE
c	      atmlow=tlog(ibot)
c	    END IF
c           IF (atmlow.LT.TLOWER) THEN
c	    print *, 'ibot = ', ibot, TLOG(Ibot)
            IF (TLOG(IBOT).LT.TLOWER) THEN

!	     	  print *, atmlow, tlower, tupper
	     	  print *, tlog(ibot), tlower, tupper

              PRINT *,' YOUR LOWER TLOG BOUNDARY IS TOO LARGE, THE ',
     1              ' MODEL DOES NOT GO SO DEEP DOWN'
              PRINT 1010,' PLEASE ENTER NEW TLOG BOUNDS : '
              READ *, TUPPER, TLOWER
              IF (TUPPER.GT.TLOWER) THEN
                THELP=TUPPER
                TUPPER=TLOWER
                TLOWER=THELP
              END IF
            END IF
          END IF
c	  print *, 'at line 856', z(1), temp(1), press(1), dens(1)
	  IF (FOROUT.EQ.'A') THEN
c need to convert... 
	    CALL INTPOL(TLOG, Z, TEMP, nel, KAPPA, PRESS, VTURB,
     1               colmas, BFIELD, VEL, BOUNDS, GAMMA, CHII, NOUT)
	  ELSE
	    CALL INTPOL(TLOG, Z, TEMP, PEL, KAPPA, PRESS, VTURB, 
     1                 DENS, BFIELD, VEL, BOUNDS, GAMMA, CHII, NOUT)
	  END IF
	END IF

C Write the final model.
c	print *, 'at line 867',  z(1), temp(1), press(1), dens(1)
	IF (FOROUT.EQ.'M') THEN
	  CALL MODOUTM (Z, TEMP, PRESS, DENS, TLOG, PEL, KAPPA, NOUT,
     1                  IRAY, NRAY)
	ELSE IF (FOROUT.EQ.'K') THEN
	  CALL STEFFOUT (TLOG, Z, TEMP, PRESS, DENS, PEL, KAPPA, 
     1                   VEL, VELHOR, WLOUT, NOUT, IRAY, NRAY)
	ELSE IF (FOROUT.EQ.'I') THEN
	  CALL INTOUT (TLOG, TEMP, PRESS,PEL, KAPPA, GRAVAC,
     1			NOUT, IRAY, NRAY)
	ELSE IF (FOROUT.EQ.'A') THEN
	  CALL KUROUT (COLMAS, TEMP, PRESS, NEL, KAPPA, VTURB, 
     1	               NOUT, IRAY, NRAY)
	ELSE IF (FOROUT.EQ.'T') THEN
	  WLOUT=WAVLEN
	  CALL TLOGOUT (TLOG, Z, NOUT, WLOUT, IRAY, NRAY)
	ELSE
	  WLOUT=WAVLEN
	  CALL MODOUT (TLOG, Z, TEMP, PEL, KAPPA, PRESS, VTURB, DENS,
     1                BFIELD, VEL, GAMMA, CHII, NOUT, WLOUT, IRAY, NRAY)
	END IF
	
	END DO
C The above is the end of loop over rays (only important for 2-D data)

1010	FORMAT (A,$)	
1020	FORMAT (A)	
1030	FORMAT (A,A,$)	
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE PTAUCR (TEMP, TLOG, PRESS, PEL, DENS, KAPPA, NPRESS, 
     1                 NPEL, NDENS, NKAPPA, NDD)

	PARAMETER (NDIM=1024, NNN=5)

	REAL      PRESS(NDIM), TEMP(NDIM), KAPPA(NDIM), TAU(NDIM)
	REAL      TLOG(NDIM), PEL(NDIM), DENS(NDIM), RHO(NDIM)
	REAL      POLD(NDIM), TLBEG(NNN), PBEG(NNN)
	CHARACTER WAVANS*1
	
	COMMON/IBOT/IBOT
	COMMON/CONST/GASCON, GRAVAC

	DATA TLBEG/-7.00, -5.91, -5.05, -4.06, -3.00/
	DATA PBEG /6.32E0,1.06E2,3.18E2,1.01E3,3.46E3/

	IF (NDD.NE.NDIM) THEN
	  PRINT*,'DIMENSIONS NOT COMPATIBLE: SUBR. PTCREA'
	  STOP
	END IF

C Preliminaries

	DO I=1, IBOT
          TAU(I)=10**TLOG(I)
	END DO

C Interpolate in table to get a very first estimate of the pressure at the 
C top of the atmosphere.

	DO J=1, NNN
	  IF (TLOG(1).GE.TLBEG(NNN)) THEN
	    PRESS(1)=PBEG(NNN)
	  ELSE IF (TLOG(1).LE.TLBEG(1)) THEN
	    PRESS(1)=PBEG(1)
	  ELSE IF (((TLOG(1).LT.TLBEG(J)).AND.(TLOG(1).GE.TLBEG(J-1)))
     1            .AND.(J.NE.1)) THEN
	    PRESS(1)=PBEG(J-1)+(PBEG(J)-PBEG(J-1))*(TLOG(1)-TLBEG(J-1))
     1              /(TLBEG(J)-TLBEG(J-1))
	  END IF
	END DO

C Improve the value of the pressure at the first grid point (top of 
C atmosphere) by iteration.

	WAVANS='Y'
	WAVLEN=5000.0
	ICON=0
	DO WHILE (ABS(PRESS(1)-POLD(1)).GT.0.01*PRESS(1))
	  POLD(1)=PRESS(1)
	  CALL PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL, KAPPA, RHO, 
     1                NPEL, NKAPPA, 1, 1, NDD)
	  PRESS(1)=GRAVAC*TAU(1)/KAPPA(1)
	  ICON=ICON+1
	END DO
	PRINT*,'STEP I=1 CONVERGED AFTER', ICON,' ITERATIONS'
	DENS(1)=RHO(1)

C Extrapolate pressure to next lower height in atmosphere. Then iterate to 
C improve this value.

	DO I=2, IBOT
C The next line has been changed: from a * to a +, which should be much better!
C	  PRESS(I)=PRESS(I-1)*GRAVAC/KAPPA(I-1)*(TAU(I)-TAU(I-1))
	  PRESS(I)=PRESS(I-1)+GRAVAC/KAPPA(I-1)*(TAU(I)-TAU(I-1))
	  ICON=0
	  DO WHILE (ABS(PRESS(I)-POLD(I)).GT.0.01*PRESS(I))
	    POLD(I)=PRESS(I)
	    CALL PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL, KAPPA, RHO, 
     1                  NPEL, NKAPPA, I, I, NDD)
	    PRESS(I)=PRESS(I-1)+GRAVAC*2./(KAPPA(I)+KAPPA(I-1))
     1              *(TAU(I)-TAU(I-1))
	    ICON=ICON+1
	  END DO
	  PRINT*,'STEP I=', I,' CONVERGED AFTER',ICON,'ITERATIONS'
	  DENS(I)=RHO(I)
	END DO
	NKAPPA=IBOT
	NPEL=IBOT
	NPRESS=IBOT
	NDENS=IBOT
	RETURN
	END
	

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE CPCREA (TEMP, TLOG, PRESS, PEL, KAPPA, 
     1                   NPRESS, NPEL, NKAPPA, NDD)

	PARAMETER (NDIM=1024, NNN=11)

	REAL      PRESS(NDIM), TEMP(NDIM), KAPPA(NDIM), RHO(NDIM)
	REAL      TAU(NDIM), TLOG(NDIM), PEL(NDIM), COLMAS(NDIM)
	REAL      POLD(NDIM), CMOLD(NDIM), TLREF(NNN), CMREF(NNN)
	CHARACTER WAVANS*1
	
	COMMON/IBOT/IBOT
	COMMON/CONST/GASCON, GRAVAC

	DATA TLREF/-7., -5.91, -5.05, -4.06, -3., -2., -1., 0., 1., 
     1              2., 3./
	DATA CMREF/1.E-3, 4.733E-3, 1.2585E-2, 3.7898E-2, 0.1270,
     1             0.4635, 1.6618, 4.7575, 6.8458, 9.5562, 18.2678/

	IF (NDD.NE.NDIM) THEN
	  PRINT*,'DIMENSIONS NOT COMPATIBLE: SUBR. CPCREA'
	  STOP
	END IF

C Preliminaries

	DO I=1, IBOT
          TAU(I)=10**TLOG(I)
	END DO
	I=1
	DO WHILE (TLOG(I).LT.0.0)
	  I=I+1
	END DO
	I0=I

C Initial guess for COLMAS(tau) taken from the HSRA (in variables TLREF and 
C CMREF)

	DO I=1, IBOT
	  DO J=1, NNN		
	    IF (TLOG(I).GE.TLREF(NNN)) THEN
	      COLMAS(I)=CMREF(NNN-1)+(CMREF(NNN)-CMREF(NNN-1))
     1                 *(TLOG(I)-TLREF(NNN-1))/(TLREF(NNN)-TLREF(NNN-1))
	    ELSE IF ((J.EQ.1).AND.(TLOG(I).LE.TLREF(1))) THEN
	      COLMAS(I)=CMREF(1)+(CMREF(2)-CMREF(1))
     1                 *(TLOG(I)-TLREF(1))/(TLREF(2)-TLREF(1))
	    ELSE IF ((TLOG(I).LT.TLREF(J)).AND.(TLOG(I).GE.TLREF(J-1)))
     1               THEN
	      COLMAS(I)=CMREF(J-1)+(CMREF(J)-CMREF(J-1))
     1                 *(TLOG(I)-TLREF(J-1))/(TLREF(J)-TLREF(J-1))
	    END IF
	  END DO
	END DO

	PRINT*,'ENTER DESIRED PRESSURE (CGS) AT TAU=1. IF P0=0, THEN ',
     1         '1.3E5 IS USED'
	READ*, PRESS0
	IF (PRESS0.EQ.0.0) PRESS0=1.31E5
	PRATIO=PRESS0/1.31E5
	DO I=1, IBOT
	  COLMAS(I)=COLMAS(I)*PRATIO
	END DO

	NSTEP=50
	DO N=1, NSTEP
	  DO I=1, IBOT
	    CMOLD(I)=COLMAS(I)
	    POLD(I)=PRESS(I)
	  END DO

C Integrate hydrostatic equilibrium to obtain p(colmas)

	  DO I=1, IBOT
	    PRESS(I)=GRAVAC*COLMAS(I)
	  END DO

C Calculate electron pressure and absorption coefficient

	  WAVANS='Y'
	  WAVLEN=5000.0
	  CALL PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL, KAPPA, RHO, 
     1                NPEL, NKAPPA, 1, IBOT, NDD)

C Determine new COLMAS-scale. Keep Colmas(1) fixed. This does not matter, since
C PRESS is scaled anyway to conform to its value at tau=0.

	  DO I=2, IBOT
	    COLMAS(I)=COLMAS(I-1)-2./(KAPPA(I)+KAPPA(I-1))
     1               *ABS(TAU(I)-TAU(I-1))
	  END DO

C Check for convergence

	  IF ((ABS(COLMAS(1)-CMOLD(1)).LT.ABS(0.01*COLMAS(1))).AND.
     1        (ABS(COLMAS(IBOT)-CMOLD(IBOT)).LT.ABS(0.01*COLMAS(IBOT)))
     2        .AND.(ABS(PRESS(1)-POLD(1)).LT.0.01*PRESS(1)).AND.
     3        (ABS(PRESS(IBOT)-POLD(IBOT)).LT.0.01*PRESS(IBOT)))
     4        GOTO 100
	END DO

	PRINT*,'FAILED TO CONVERGE AFTER', NSTEP,' STEPS'
	STOP
100	PRINT*,'CONVERGED AFTER',N,' STEPS'

	NPEL=IBOT
	NPRESS=IBOT
	NKAPPA=IBOT	

	RETURN
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE PZCREA (TEMP, TLOG, RMU, Z, PRESS, PEL, KAPPA, DENS, 
     1                     NZ, NPRESS, NPEL, NKAPPA, NDENS, NDD)

	PARAMETER (NDIM=1024, NNN=11)

	REAL      Z(NDIM), PRESS(NDIM), TEMP(NDIM), KAPPA(NDIM)
	REAL      TAU(NDIM), TLOG(NDIM), PEL(NDIM), RMU(NDIM)
	REAL      DENS(NDIM), TLOGTRY(NDIM), TAUCOR(NDIM)
	REAL      TAUOLD(NDIM), TAUTRY(NDIM), ZOLD(NDIM), RHO(NDIM)
	REAL      TLREF(NNN), ZREF(NNN), ZINIT(NDIM)
	LOGICAL   ZGOOD
	CHARACTER WAVANS*1
	
	COMMON/IBOT/IBOT
	COMMON/CONST/GASCON, GRAVAC

	DATA TLREF/-7., -5.91, -5.05, -4.06, -3., -2., -1., 0., 1., 
     1              2., 3./
	DATA ZREF/1.02E8, 7.20E7, 6.54E7, 5.43E7, 4.20E7, 2.83E7, 
     1            1.38E7, 0.00E0,-6.38E6,-1.41E7,-3.12E7/

	IF (NDD.NE.NDIM) THEN
	  PRINT*,'DIMENSIONS NOT COMPATIBLE: SUBR. PZCREA'
	  STOP
	END IF

C Preliminaries

	DO I=1, IBOT
          TAU(I)=10**TLOG(I)
	  TAUTRY(I)=TAU(I)
	END DO
	I=1
	DO WHILE (TLOG(I).LT.0.0)
	  I=I+1
	END DO
	I0=I

C Initial guess for z(tau) determined by interpolating between selected values 
C taken from the HSRA (in variables TLREF and ZREF)

	DO I=1, IBOT
	  DO J=1, NNN		
	    IF (TLOG(I).GE.TLREF(NNN)) THEN
	      ZINIT(I)=ZREF(NNN-1)+(ZREF(NNN)-ZREF(NNN-1))
     1            *(TLOG(I)-TLREF(NNN-1))/(TLREF(NNN)-TLREF(NNN-1))
	    ELSE IF ((J.EQ.1).AND.(TLOG(I).LE.TLREF(1))) THEN
	      ZINIT(I)=ZREF(1)+(ZREF(2)-ZREF(1))*(TLOG(I)-TLREF(1))
     1            /(TLREF(2)-TLREF(1))
	    ELSE IF ((TLOG(I).LT.TLREF(J)).AND.(TLOG(I).GE.TLREF(J-1)))
     1              THEN
	      ZINIT(I)=ZREF(J-1)+(ZREF(J)-ZREF(J-1))
     1            *(TLOG(I)-TLREF(J-1))/(TLREF(J)-TLREF(J-1))
	    END IF
	  END DO
	  Z(I)=ZINIT(I)
	END DO

	PRINT*,'ENTER PRESSURE (CGS) AT Z=0 (=TAU=1). IF P0=0, THEN ',
     1         ' 1.3E5 IS USED'
	READ*, PRESS0
	IF (PRESS0.EQ.0.0) PRESS0=1.31E5

C	INEW=2
	DO N=1, 20*IBOT

C Integrate hydrostatic equilibrium to obtain p(z)

	  IP=0
	  ZGOOD=.TRUE.
	  CALL HYDST (TEMP, PRESS, PRESS0, Z, RMU, ZGOOD, ICRIT, NDD)

	  IF (.NOT.ZGOOD) THEN
	    PRINT*,'PRESSURE AT TOP OF ATMOSPHERE IS TOO SMALL FOR ',
     1             'ABSKO. LOWER UPPER'
	    PRINT*,'LIMIT ON ATMOSPHERE, I.E. CHOOSE LARGER LOG(TAU)',
     1             'FOR FIRST POINT'
	    PRINT*,'OR ANOTHER PRESSURE'
	    IBAD=IBOT
	    GOTO 50
	  ELSE IF (ZGOOD) THEN
	    PRINT*,'ITERATION STEP NUMBER', N
	  END IF	  

C Calculate electron pressure and absorption coefficient. For log(tau)>0,
C these have to be calculated right from the top downwards, since changes in 
C Z change the pressure, thus changing PEL and KAPPA....... (the whole mess is 
C due to the fact that although TLOG is integrated from the top down, the 
C pressure is integrated from Z=0 (i.e. TLOG=0) outwards)

	  WAVANS='Y'
	  WAVLEN=5000.0
	  CALL PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL, KAPPA, RHO, 
     1                NPEL, NKAPPA, 1, IBOT, NDD)

C Calculate the density

	  DO I=1,IBOT
	    DENS(I)=RMU(I)*PRESS(I)/TEMP(I)/GASCON
	  END DO

C Store the Z-scale for use in determining the new one.

	  DO I=1, IBOT
	    ZOLD(I)=Z(I)
	    TAUOLD(I)=TAUTRY(I)
	  END DO

C Calculate the new optical depth scale.

	  TAUTRY(1)=TAU(1)
	  TLOGTRY(1)=LOG10(TAUTRY(1))
	  DO I=2, IBOT
	    TAUTRY(I)=TAUTRY(I-1)+ABS(Z(I)-Z(I-1))
     1               *(KAPPA(I)*DENS(I)+KAPPA(I-1)*DENS(I-1))/2.
	    TLOGTRY(I)=LOG10(TAUTRY(I))	    
	  END DO

C Calculate the new Z scale, by correcting the old one in view of differences
C between the current tau scale (TAUTRY) and the true one (TAU)

	  Z(1)=ZOLD(1)
	  DO I=2, IBOT
	    TAUCOR(I)=(TLOG(I)-TLOG(I-1))/(TLOGTRY(I)-TLOGTRY(I-1))
	    Z(I)=Z(I-1)+(ZOLD(I)-ZOLD(I-1))*TAUCOR(I)
C	    Z(I)=Z(I-1)+(ZOLD(I)-ZOLD(I-1))*(TLOG(I)-TLOG(I-1))
C     1          /(TLOGTRY(I)-TLOGTRY(I-1))
	  END DO

	  ZSHIFT=Z(I0)+(Z(I0)-Z(I0-1))*TLOG(I0)/(TLOG(I0-1)-TLOG(I0))
	  DO I=1, IBOT
	    Z(I)=Z(I)-ZSHIFT
	  END DO

C Check for convergence of the TAUTRY values to the input TAU values.

	  IRUN=1
	  I=1
	  PRINT *,'I,TAUTRY,TAU, DT',I,TAUTRY(I),TAU(I), 
     1             (TAUTRY(I)-TAU(I))/TAU(I)
	  DO I=2, IBOT
	    IF (ABS(TAUTRY(I)-TAU(I)).LE.0.01*ABS(TAU(I)))
C	    IF (ABS(TAUTRY(I)-TAU(I)).LE.0.03*ABS(TAU(I)))
     1         IRUN=IRUN+1
	    PRINT *,'I,TAUTRY,TAU',I,TAUTRY(I),TAU(I), 
     1              (TAUTRY(I)-TAU(I))/TAU(I)
	  END DO
	  IF (IRUN.EQ.IBOT) GOTO 100

50	  DO I=1, IBOT
	    IF (((TAUTRY(I).LT.TAUOLD(I)).AND.(TAUOLD(I).LT.TAU(I))).OR.
     1         ((TAUTRY(I).GT.TAUOLD(I)).AND.(TAUOLD(I).GT.TAU(I))).OR.
     2         ((ABS(TAUTRY(I)-TAUOLD(I)).LT.0.02*TAUTRY(I)).AND.
     3         (ABS(TAUTRY(I)-TAU(I)).GT.0.05*TAU(I)))) IBAD=IBAD+1
	  END DO
	  IF (IBAD.GT.IBOT/2) THEN
	    PRINT*,'CHANGING P0=PRESSURE(Z_0). OLD P0=',PRESS0
C	    PRESS0=PRESS0*TAU(I0)/TAUTRY(I0)  !PRODUCED TOO BIG JUMPS: OVERSHOT
	    PFACT=SQRT(TAU(I0)/TAUTRY(I0))
	    PRESS0=PRESS0*PFACT
	    PRINT*,'NEW P0=',PRESS0
	    IF ((PFACT.GT.1.03).OR.(PFACT.LT.0.97).OR.(.NOT.ZGOOD)) THEN
	      DO I=1, IBOT
	        Z(I)=ZINIT(I)
	      END DO
	    END IF
	  END IF
	  IBAD=0

	END DO

	PRINT*,'FAILED TO CONVERGE AFTER', N-1,' STEPS'
	STOP
100	PRINT*,'CONVERGED AFTER',N,' STEPS'

	NZ=IBOT
	NPEL=IBOT
	NPRESS=IBOT
	NKAPPA=IBOT	
	NDENS=IBOT

	RETURN
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE HYDST (TEMP, PRESS, PRESS0, Z, RMU, ZGOOD, ICRIT, 
     1                    NDD)

C Upon input of temperature, height scale, mean molecular weight and pressure at
C a given height the routine returns the pressure at every height in the 
C atmosphere (it integrates the equations of hydrostatic equilibrium)

	PARAMETER (NDIM=1024)

	LOGICAL ZGOOD
	REAL Z(NDIM), TEMP(NDIM), PRESS(NDIM), RMU(NDIM)
	REAL F(NDIM)

	COMMON/IBOT/IBOT
	COMMON/CONST/GASCON, GRAVAC

	IF (NDD.NE.NDIM) THEN
	  PRINT*,'DIMENSIONS NOT COMPATIBLE: SUBR. HYDST'
	  STOP
	END IF

	ZGOOD=.TRUE.

	I=1
	DO WHILE ((I.LT.IBOT).AND.(Z(I).GT.0.0))
	  I=I+1
	END DO
	I0=I
	PRESS(I0)=PRESS0
	DO I=1,IBOT
          F(I)=GRAVAC*RMU(I)/GASCON/TEMP(I)
	END DO

	SUMINT=0.0
	DO I=I0,2,-1
c	  SUMINT=SUMINT+0.5*(F(I-1)+F(I))*(Z(I-1)-Z(I))
c	  PRESS(I-1)=PRESS(I0)*EXP(-SUMINT)
	  PRESS(I-1)=PRESS(I)*EXP(-0.5*(F(I-1)+F(I))*(Z(I-1)-Z(I)))
C	  IF (PRESS(I-1).LT.0.01) THEN       
	  IF (PRESS(I-1).LT.1.E-10) THEN       
	    PRINT *,'TOO LOW PRESSURE, P=',PRESS(I-1),'  HEIGHT=',I-1
	    ZGOOD=.FALSE.
	    ICRIT=I-1
	    RETURN
	  END IF
	END DO

	SUMINT=0.0
	DO I=I0, IBOT-1
	  SUMINT=SUMINT+0.5*(F(I+1)+F(I))*(Z(I+1)-Z(I))
	  PRESS(I+1)=PRESS(I0)*EXP(-SUMINT)
	  IF (PRESS(I+1).GT.5.0E8) THEN
	    PRINT *,'TOO HIGH PRESSURE, P=',PRESS(I+1),'  HEIGHT=',I+1
	    ZGOOD=.FALSE.
	    ICRIT=I+1
	    RETURN
	  END IF
	END DO

	RETURN
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE MODINP (TLOG, Z, TEMP, PEL, KAPPA, PRESS, VTURB, 
     1                     DENS, BFIELD, VEL, BOUNDS, GAMMA, CHII, 
     2                     NEL, NUMDENS, TAU, COLMAS, VELHOR, TYPINP, 
     3                     IRAY, NRAY)

C Reads the input model. Depending on which variables are known a goto is 
C carried out that to the position with the relevant input. The code assumes 
C that the input atmosphere is written in the order of increasing log(tau).
C
C Atmospheric quantities:
C       COLMAS = column mass
C       COLMASLG = logarithm of column mass
C       DENS   = total mass density [g/cm^3]
C       KAPPA  = absorption coefficient [1/g]
C       NEL    = electron number density [1/cm^3]
C       NUMDENS = total number density
C	NH     = H number density
C       PEL    = electron pressure [cgs]
C       PELG   = logarithm of electron pressure [cgs]
C       PGPTOT = ratio of gas to total pressure [cgs]
C       PRESS  = gas pressure [cgs]
C       PRESSLG = logarithm of gas pressure [cgs]
C       PTOT   = total pressure = gas + turbulent pressure [cgs]
C       TAU    = optical depth
C       TEMP   = temperature [K]
C       TLOG   = logarithm of optical depth
C       VTURB  = microturbulent velocity [cm/s]
C       Z      = height [cm]
C       ZKM    = height [km]
C       BFIELD = magnetic field strength [G]
C       CHII   = azimuthal angle of the magnetic field [degrees]
C       GAMMA  = inclination angle of the magnetic field to the line of sight

	PARAMETER (NDIM=1024)

	REAL      TLOG(NDIM), Z(NDIM), TEMP(NDIM), PEL(NDIM)
	REAL      BFIELD(NDIM), KAPPA(NDIM), PRESS(NDIM), VTURB(NDIM)
	REAL      DENS(NDIM), VEL(NDIM), GAMMA(NDIM), CHII(NDIM)
	REAL      COLMAS(NDIM), COLMASLG(NDIM), PELG(NDIM), ZKM(NDIM)
	REAL      PRESSLG(NDIM), TAU(NDIM), RMU(NDIM), PGPTOT(NDIM)
	REAL      PTOT(NDIM), NEL(NDIM), DENSLN(NDIM), VELHOR(NDIM)
	REAL      NUMDENS(NDIM), NH(NDIM), NP(NDIM), FOO(NDIM)
	CHARACTER MODNAM*55, BOUNDS*1, TYPINP*1, junktext*50
	real	  grav, Teff
        integer   n, ndepth

	COMMON/IBOT/IBOT
	COMMON/MODNAM/MODNAM
	COMMON/TBOUND/ TUPPER, TLOWER
	COMMON/NDD/NDD
	COMMON/DD/DD(NDIM)
	COMMON/CONST/GASCON, GRAVAC
	common/Teff/Teff

        save itot
        save ndepth

	IF (NDD.NE.NDIM) STOP 'SUBR. MODINP: ARRAY DIMENSIONS ARE WRONG'

C Initialize some variables
	III=0
	DO JJJ=1, NDIM
          Z(JJJ)=0.0
	  PELG(JJJ)=0.0
	  PRESSLG(JJJ)=0.0
	  PTOT(JJJ)=0.0
	  PGPTOT(JJJ)=0.0
	  tlog(jjj) = 0.0
	  pel(jjj) = 0.0
	  kappa(jjj) = 0.0
	  dens(jjj) = 0.0
	  colmas(jjj) = 0.0
	  colmaslg(jjj) = 0.0
	  presslg(jjj) = 0.0
	  press(jjj) = 0.0
	  tau(jjj) = 0.0
	  nel(jjj) = 0.0
          zkm(jjj) = 0.0
	END DO


	IF (IRAY.EQ.1) OPEN(UNIT=11, FILE=MODNAM,STATUS='OLD')

	IF (TYPINP.EQ.'K') THEN
	  ITOT=55
        else if (TYPINP .eq. '2') then
          if (iray .eq. 1) then
             read (11,*) n
             if (n .ne. nray) then
                 print *, 'file has ', n, 'rays vs ', nray, 'input'
             end if
             read (11,*) ndepth
          end if 
c          read (11,*) jray
          itot = ndepth
	else 
	  READ (11,*) ITOT
	END IF

c        print *, iray, itot

	IF ((BOUNDS.EQ.'H').AND.(ITOT.LT.IBOT)) THEN
	  STOP 'YOUR FILE DOES NOT HAVE SO MANY HEIGHT STEPS'
	END IF

	IF (TYPINP.EQ.'S') THEN
C goto 100 for input of standard model (correct format for STOPRO)
	  GOTO 100
	ELSE IF (TYPINP.EQ.'T') THEN
C goto 90 for input of model with only tlog, temp known (but model is written in
C standard input format, i.e. that z lies between tlog & temp, but is only 
C considered to be a dummy).
	  GOTO 90
	ELSE IF (TYPINP.EQ.'K') THEN
C goto 80 for input of Kurucz's models
	  GOTO 80
	ELSE IF (TYPINP.EQ.'R') THEN
C goto 21 for input of new Fontenla models (ApJ 2006)
	  GOTO 21
	ELSE IF (TYPINP.EQ.'M') THEN
C goto 4 for Maltby et al. Reference queit sun atmosphere
	  GOTO 4
	ELSE IF (TYPINP.EQ.'U') THEN
C goto 5 for Maltby et al. Sunspot Umbra M atmosphere
	  GOTO 5
	ELSE
C goto 19 for input of old Fontenla models (FAL)
C goto 14 for input of MURAM MODELS
C goto 1 for Holweger-Mueller photosphere model
C goto 2 for z, temp
C goto 3 for Anderson's radiative equilibrium solar atmosphere
C goto 4 for Maltby et al. Reference queit sun atmosphere
C goto 5 for Maltby et al. Sunspot Umbra M atmosphere
C goto 6 for Basri+Marcy K2 atmosphere
C goto 7 for Tom Ayres' models
C goto 8 for Lemaire et al. (1983, A+A ) models
C goto 9 for models written with this code in the MHD output format
C goto 10 for MALTX.
C goto 11 for models with TLOG, Z and TEMP known
C goto 12 for Obridko-Staude umbral models
C goto 13 for Ding-Fang penumbral  model
C goto 14 for models with Z, PRESS, DENS and TEMP known (SPRVAL of Oski).
C goto 15 for Jo's format of the Nordlund-Stein granulation model.
C goto 16 for Mathias Steffen's granulation models (Z,T,P,-V, from bottom to top
C goto 17 for input of first part of Multi atmosphere format (without H pops).
C goto 18 for input of z, T, P and B etc. (e.g. for canopies with P-balance)
C goto 19 for Fontenla models (depth instead of height)
C goto 21 for Fontenla models (ApJ 2006)
C goto 22 for Sim & Jordan (2005) model
c goto 80 for Kurucz's models
C goto 88 for tests (can be changed at will)
          GOTO 14
c          GOTO 22
C goto 99 for standard STOPRO format, but with PEL & KAPPA to be recalculated
	END IF

C Input of Holweger-Mueller (HOLMUL) model atmosphere
1	PRINT*,'READING THE HOLWEGER-MUELLER MODEL'
	PRINT*,'**********************************'
	DO I=1, ITOT
	  READ (11,*) TLOG(I),TEMP(I), PRESS(I),PEL(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of models for which only Z and T are known, i.e. the whole hydrostatics
C must be calculated (practically the same as MODEL4, except that no external
C atmosphere has to be assumed). P_g must be given at the lower boundary.
2	PRINT*,'READING A MODEL FOR WHICH ONLY Z & T ARE KNOWN'
	PRINT*,'**********************************************'
	DO I=1, ITOT
	  READ (11,*) Z(I),TEMP(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of Anderson's radiative equilibrium model with CO
3	PRINT*,'READING ANDERSONS RE MODEL'
	PRINT*,'***************************'
	DO I=1, ITOT
	  READ (11,*) COLMASLG(I), TLOG(I), ZKM(I), TEMP(I), PELG(I), 
     1                PRESSLG(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of Maltby et al. Reference atmosphere
4	PRINT*,'READING THE MALTBY ET AL. QUIET SUN MODEL'
	PRINT*,'*****************************************'
	DO I=1, ITOT
	  READ (11,*) ZKM(I), COLMAS(I), TAU(I), TEMP(I), NEL(I), 
     1                PTOT(I),PGPTOT(I), DENS(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of Maltby et al. Sunspot M or E atmosphere
5	PRINT*,'READING A MALTBY ET AL. UMBRAL MODEL'
	PRINT*,'************************************'
	DO I=1, ITOT
	  READ (11,*) ZKM(I), TAU(I), TEMP(I), PTOT(I), PGPTOT(I), 
     1                DENS(I), NEL(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of Basri+Marcy K2 atmosphere
6	PRINT*,'READING THE BASRI & MARCY K2 MODEL'
	PRINT*,'*********************************'
	DO I=1, ITOT
	  READ (11,*) TEMP(I), PEL(I), PRESS(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	OPEN (UNIT=12, FILE='K2LOGM', STATUS='OLD')
	  READ (12,*)
	  DO I=1, ITOT, 3
	    READ (12,*) COLMASLG(I)
	  END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=12)
	OPEN (UNIT=13, FILE='K2TAU5000', STATUS='OLD')
	  READ (13,*)
	  DO I=1, ITOT, 3
	    READ (13,*) TAU(I)
	  END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=13)
C Interpolation	
	DO I=1, ITOT, 3
	  TLOG(I)=LOG10(TAU(I))
	END DO
	DO I=1, ITOT
	  PRESSLG(I)=LOG(PRESS(I))
	END DO
	DO I=1, ITOT
	  IF (MOD((I-1),3).EQ.1) THEN
	    COLMASLG(I)=COLMASLG(I-1)+(COLMASLG(I+2)-COLMASLG(I-1))
     1                 *(PRESSLG(I)-PRESSLG(I-1))
     2                 /(PRESSLG(I+2)-PRESSLG(I-1))
	    TLOG(I)=TLOG(I-1)+(TLOG(I+2)-TLOG(I-1))
     1                 *(PRESSLG(I)-PRESSLG(I-1))
     2                 /(PRESSLG(I+2)-PRESSLG(I-1))
	  ELSE IF (MOD((I-1),3).EQ.2) THEN
	    COLMASLG(I)=COLMASLG(I-2)+(COLMASLG(I+1)-COLMASLG(I-2))
     1                 *(PRESSLG(I)-PRESSLG(I-2))
     2                 /(PRESSLG(I+1)-PRESSLG(I-2))
	    TLOG(I)=TLOG(I-2)+(TLOG(I+1)-TLOG(I-2))
     1                 *(PRESSLG(I)-PRESSLG(I-2))
     2                 /(PRESSLG(I+1)-PRESSLG(I-2))
	  END IF
	  COLMAS(I)=10**COLMASLG(I)
	END DO
	    
C Density determination (only for Basri and Marcy)
	DO I=1, ITOT
	  RMU(I)=1.24
	  DENS(I)=RMU(I)*PRESS(I)/TEMP(I)/GASCON
	END DO
	IF (III.EQ.0) GOTO 140

C Tom Ayres' models
7       PRINT*,'READING A TOM AYRES MODEL'
	PRINT*,'*************************'
	DO I=1, ITOT
	  READ(11,*) COLMAS(I), TEMP(I), DENS(I), TAU(I), NEL(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Lemaire et al. (1983, A+A ) models
8       PRINT*,'READING THE LEMAIRE ET AL. MODEL'
	PRINT*,'********************************'
	DO I=1, ITOT
	  READ(11,*) ZKM(I), COLMAS(I), TEMP(I), PRESS(I), DENS(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of models written with this code in the MHD output format
9	PRINT*,'READING A MODEL WRITTEN IN THE MHD OUTPUT FORMAT'
	PRINT*,'************************************************'
	DO I=1, ITOT
	  READ(11,*) TLOG(I), Z(I), TEMP(I), PRESS(I), DENSDUM, PEL(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140
C Input of MALTX.
10      DO I=1, ITOT
	  READ(11,*) TLOG(I), Z(I), TEMP(I), PRESS(I), DENS(I), PEL(I),
     1               KAPPA(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of models with TLOG, Z and TEMP known
11	PRINT*,'READING A MODEL WITH TLOG, Z, TEMP ONLY'
	PRINT*,'***************************************'
	DO I=1, ITOT
	  READ(11,*) TLOG(I), Z(I), TEMP(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of Obridko-Staude sunspot models: Z, TAU, TEMP, PEL, PRESS, DENS  known
12	PRINT*,'READING A OBRIDKO-STAUDE MODEL'
	PRINT*,'******************************'
	DO I=1, ITOT
	  READ(11,*) Z(I), TAU(I), TEMP(I), PEL(I), PRESS(I), DENS(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of DING-FANG penumbral model: COLMAS, Z, TEMP, NEL, TAU_5000 known
13	PRINT*,'READING THE DING-FANG MODEL'
	PRINT*,'***************************'
	DO I=1, ITOT
	  READ(11,*) COLMAS(I), ZKM(I), TEMP(I), NEL(I), TAU(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Models with Z, PRESS, DENS and TEMP known (SPRVAL of Oski/muram)
14	discard = 1
c	PRINT*,'READING A MODEL WITH Z,P,DENS AND T ONLY'
c	PRINT*,'****************************************'
	DO I=1, ITOT
c        READ(11,*)  TEMP(i), PRESS(I), DENS(i)
c         READ(11,*)  Z(I), TEMP(I), PRESS(I), DENS(I)
          read(11,*) ztemp,temp(i),press(i),dens(i)
          z(i) = 1e+05*ztemp
c          z(i) = 6.25e+04*ztemp
c	  print *, ztemp, temp(i), press(i), dens(i)
	  
	END DO
c	print *,  ztemp
c	print *, 'vals = ', z(1),temp(1),press(1),dens(1)
c	print *, z(500),temp(500),press(500),dens(500)
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Jo's format of the Nordlund-Stein granulation model.
15      PRINT*,'READING A NORDLUND-STEIN MODEL IN JO BRULS FORMAT'
	PRINT*,'**************************************************'
	DO I=1, ITOT
	  READ(11,*) Z(I), COLMAS(I), TEMP(I), NEL(I), VEL(I), VTURB(I)
C	  Z(I)=6.5E7-2.5E6*(I-1)
C	  COLMAS(I)=10**COLMASLG(I)
	  VEL(I)=VEL(I)*1.E5
	  VTURB(I)=VTURB(I)*1.E5
	END DO
	PRINT*,'NORDLUND-STEIN MODEL IN JO BRULS'' FORMAT READ'
	IF (Z(ITOT).NE.-5.E7) THEN
	  PRINT*,'PROBLEMS WITH READING NORDLUND STEIN MODEL'
	  PRINT*,'ITOT=', ITOT,'Z(ITOT)=',Z(ITOT),' INSTEAD OF -500KM'
	END IF
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Models with Z, PRESS, DENS and TEMP known (Single ray of Steffen's model)
C VELHOR is the horizontal velocity, VEL is the vertical velocity
16	PRINT*,'READING A SINGLE RAY OF STEFFENS MODEL'
	PRINT*,'***************************************'
	ITOT= 61
	READ(11,*)
	DO I=ITOT, 1, -1
	  READ(11,*) Z(I), TEMP(I), PRESS(I), VEL(I), VELHOR(I)
	  VEL(I)=-VEL(I)
	  VTURB(I)=1.E5
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C First part of Multi input atmosphere, Han's format (changed somewhat!!!!)
17	PRINT*,'READING A MODEL IN HANS FORMAT'
	PRINT*,'*******************************'
	READ(11,*)
	DO I=1, ITOT
	  READ(11,*) COLMAS(I), TEMP(I), PELDUMMY, VEL(I), VTURB(I)
	  COLMAS(I)=10**COLMAS(I)
C	  VEL(I)=-VEL(I)
C	  VTURB(I)=1.E5
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Model with Z, T, P, B etc. at input (e.g. for canopies with pressure balance)
18	PRINT*,'READING A MODEL WITH ONLY T,P,Z)'
	PRINT*,'********************************'
	READ(11,*)
	DO I=1, ITOT
	  READ(11,*) TLOGDUMMY, Z(I), TEMP(I), PRESS(I), PELDUMMY, 
     1               KAPPADUMMY, DENSDUMMY, BFIELD(I), VTURB(I), VEL(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Fontenla models A, C, F and P (depth instead of height)    
 19	PRINT*,'READING A FONTENLA ET AL. MODEL'
	PRINT*,'*******************************'

c	XXX SHOULD REALLY include density determination here.. 
c	      dens=numdens*wghtmole+nel*9.108e-28;

c	      if(i > 0) abund[i]=pow((double)10,abund[i])*ascale;
c      	      wghtmole+=abund[i]*amass[i];
c    wghtmole*=1.660e-24;
	DO I = 1, ITOT
	   READ(11,*) IDUMMY, ZKMNEG, TEMP(I), NEL(I), RPROTON,
     1                RH1, RH2, RHE1, RHE2, RHEII, RHEIII
	   ZKM(I) = -ZKMNEG
	   NUMDENS(I) = NEL(I)+RPROTON+RH1+RH2+RHE1+RHE2+RHEII+RHEIII
	ENDDO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C New RISE Fontenla models @@@@@@BEFORE LAURENCE MESSING
 20     PRINT*,'READING A FONTENLA AVRETT (RISE) MODEL'
        PRINT*,'**************************************'
        DO I = 1, ITOT
           READ(11,*) IDUMMY, ZKMNEG, TEMP(I), NEL(I), NH(I),
     1                VTURB(I), VP, BH_MIN, VEL(I)
           ZKM(I) = -ZKMNEG
           NUMDENS(I) = NEL(I)+NH(I)
	   VTURB(I) = 1.e+05 * VTURB(I)
	   VEL(I) = 1.e+05 * VEL(I)
        ENDDO

        IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
        IF (III.EQ.0) GOTO 140

C New RISE Fontenla models 
c 20     PRINT*,'READING A FONTENLA AVRETT (RISE) MODEL'
c        PRINT*,'**************************************'
c        DO I = 1, ITOT
c           READ(11,*) IDUMMY, ZCM, TEMP(I), NEL(I), NP(I), 
c     1			FOO(I), NH(I), VTURB(I), VP, BH_MIN, VEL(I)
c           ZKM(I) = ZCM /10000.0
c           NUMDENS(I) = NP(I)+NH(I)
c	   VTURB(I) = 1.e+05 * VTURB(I)
c	   VEL(I) = 1.e+05 * VEL(I)
c        ENDDO
c
c        IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
c        IF (III.EQ.0) GOTO 140
c
C Kurucz's models

C New Fontenla models C, E, F, H, P, R, S (ApJ 2006)
 21     PRINT*,'READING A FONTENLA AVRETT  MODEL (ApJ 2006)'
        PRINT*,'**************************************'
        DO I = 1, ITOT
          READ(11,*) PRESS(I), TEMP(I), ZKM(I), NEL(I), NP(I), 
     1              FOO(I), VTURB(I)
C           ZKM(I) = -ZKMNEG
           NUMDENS(I) = NEL(I)+NP(I)
	   VTURB(I) = 1.e+05 * VTURB(I)
C	   VEL(I) = 1.e+05 * VEL(I)
        ENDDO

        IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
        IF (III.EQ.0) GOTO 140

C Stuart Sim and Carole Jordan's Eps Eri models based on Thatcher
 22     print *, 'READING ONE OF THE SIM-JORDAN (2005) MODELS'
        print *, '**************************************'
        do i = 1, itot
           read(11,*) colmaslg(i),temp(i),nel(i),nh(i),vturb(i)
           colmas(i) = 10**colmaslg(i)
           temp(i) = 10**temp(i)
           nel(i) = 10**nel(i)
           nh(i) = 10**nh(i)
           vturb(i) = 1.e+05*vturb(i)
       end do
       IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
       IF (III.EQ.0) GOTO 140

80	PRINT*,'READING A KURUCZ RE MODEL'
	PRINT*,'*************************'
	read(11,'(a)') junktext
	print *, junktext
	print *, ' Enter Teff and log gravity '
	read(*,*) teff, grav
	gravac = 10**(grav)
        DO IP=1, 21
          READ (11,*)
        END DO
        read(11,'(a)') junktext
        print *, junktext
        print *, ' Enter number of lines for the model '
        read(*,*) itot
	itot = itot-1

	READ (11,*)
C Skip over first line of Kurucz's model, since it seems to be ill.
        DO I=1, ITOT
	  READ(11,*) COLMAS(I), TEMP(I), PRESS(I), NEL(I), 
     1               abross, ACCRAD,  VTURB(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Just to test things. This input can be changed to please.
88	PRINT*,'READING A TEST INPUT (MAY BE CHANGED AT WILL)'
	PRINT*,'*********************************************'
	DO I=1, ITOT
	  READ(11,*) zneg, PRESS(I), TEMP(I), PEL(I), dummy, 
     1               DENS(I), KAPPA(I)
CCCCC , PRESS(I), PELGAGA, 
CCCCC     1               KAPPA(I), DENS(I), BFIELD(I), VTURB(I), VEL(I)
          z(i)=-zneg
	END DO
	PRINT*,'**********************************'
	PRINT*,'READING IN TEST INPUT !!!!!!!!!!!!'
	PRINT*,'**********************************'
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of T(tau) model, for which the rest of the quantities are to be 
C determined.
90	PRINT*,'READING A MODEL FOR WHICH ONLY T(TAU) IS KNOWN'
	PRINT*,'**********************************************'
	DO I=1, ITOT
	  READ(11,*) TLOG(I), ZDUMMY, TEMP(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Atmospheres in standard STOPRO format, but for which the el. pressure and the
C opacity have to be recalculated.

99	PRINT*,'READING STANDARD INPUT FORMAT (PEL & KAP DUMMIES)'
	PRINT*,'*************************************************'
	DO I=1, ITOT
	  READ(11,*) TLOG(I), Z(I), TEMP(I), PRESS(I), PELDUMMY, 
     1               KAPPADUMMY, DENS(I), BFIELD(I), VTURB(I), VEL(I)
	END DO
	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (III.EQ.0) GOTO 140

C Input of models created with MODEL3.FOR or MODEL4.FOR after Sept. 1986 !!!!!!
C (only for interpolation purposes)
100	PRINT*,'READING A MODEL IN THE OLD OUTPUT FORMAT'
	PRINT*,'****************************************'
	DO I=1, ITOT
	  READ(11,*, ERR=110,IOSTAT=IERR) TLOG(I), Z(I), TEMP(I), 
     1        PRESS(I), PEL(I), KAPPA(I), DENS(I), BFIELD(I), VTURB(I),
     2        VEL(I), GAMMA(I), CHII(I)

	END DO
110	IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	IF (IERR.NE.0) THEN
	  OPEN (UNIT=11, FILE=MODNAM,STATUS='OLD')
	  READ (11,*) ITOT
	  IF ((BOUNDS.EQ.'H').AND.(ITOT.LT.IBOT)) THEN
	    STOP 'YOUR FILE DOES NOT HAVE SO MANY HEIGHT STEPS'
	  END IF
	  DO I=1,ITOT
	  READ(11,*, ERR=120,IOSTAT=IER2) TLOG(I), Z(I), TEMP(I), 
     1        PRESS(I), PEL(I), KAPPA(I), DENS(I), BFIELD(I), VTURB(I),
     2        VEL(I), GAMMA(I)
	    CHII(I)=0.0
	  END DO
120	  IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	END IF
	IF (IER2.NE.0) THEN
	  OPEN (UNIT=11, FILE=MODNAM,STATUS='OLD')
	  READ (11,*) ITOT
	  IF ((BOUNDS.EQ.'H').AND.(ITOT.LT.IBOT)) THEN
	    STOP 'YOUR FILE DOES NOT HAVE SO MANY HEIGHT STEPS'
	  END IF
	  DO I=1,ITOT
	  READ(11,*, ERR=130,IOSTAT=IERR) TLOG(I), Z(I), TEMP(I), 
     1        PRESS(I), PEL(I), KAPPA(I), DENS(I), BFIELD(I), VTURB(I),
     2        VEL(I)
	    GAMMA(I)=0.0
	    CHII(I)=0.0
	  END DO
130        IF (IRAY.EQ.NRAY) CLOSE  (UNIT=11)
	END IF

140	IF (BOUNDS.EQ.'T' .and. TYPINP.EQ.'S') THEN
	  DO I=1, ITOT
	    IF (TLOWER.LE.TLOG(I)) IBOT=I
	  END DO
	  IF (TUPPER.LT.TLOG(1)) 
     1       STOP 'YOUR UPPER TLOG BOUNDARY IS TOO HIGH'
	  IF (TLOWER.GT.TLOG(ITOT)) 
     1       STOP 'YOUR LOWER TLOG BOUNDARY IS TOO DEEP'
	ELSE IF (BOUNDS.EQ.'T' .and. TYPINP.EQ.'K') THEN
	   print *, ' Cannot check boundaries yet '
	   ibot = itot
	ELSE IF (BOUNDS.NE.'H') THEN
	  IBOT=ITOT
	END IF

C Change non-standard formats for variables into the formats standard for 
C STOKES input, etc.

	IF ((ZKM(1).NE.0.0).AND.(Z(1).EQ.0.0)) THEN
	  DO I=1, ITOT
	    Z(I)=ZKM(I)*1.0E+05
	  END DO
	END IF

	IF ((COLMASLG(1).NE.0.0).AND.(COLMAS(1).EQ.0.0)) THEN
	   DO I=1,ITOT
	      COLMAS(I)=10**COLMASLG(I)
	   END DO
	END IF

	IF ((PELG(1).NE.0.0).AND.(PEL(1).EQ.0.0)) THEN
	  DO I=1, ITOT
	    PEL(I)=10**PELG(I)
	  END DO
	END IF

	IF ((PRESSLG(1).NE.0.0).AND.(PRESS(1).EQ.0.0)) THEN
	  DO I=1, ITOT
	    PRESS(I)=10**PRESSLG(I)
	  END DO
	END IF

	IF ((PTOT(1).NE.0.0).AND.(PGPTOT(1).NE.0.0).AND.
     1      (PRESS(1).EQ.0.0)) THEN
	  DO I=1, ITOT
	    PRESS(I)=PGPTOT(I)*PTOT(I)
	  END DO
	END IF

	RETURN
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE MODOUTM (Z, TEMP, PRESS, DENS, TLOG, PEL, KAPPA, 
     1                      NOUT, IRAY, NRAY)

	PARAMETER (NDIM=1024)

	REAL Z(NDIM), TEMP(NDIM), PRESS(NDIM), DENS(NDIM), TLOG(NDIM)
	REAL PEL(NDIM), KAPPA(NDIM)

	CHARACTER OUTNAM*25

	COMMON/IBOT/IBOT
	COMMON/OUTNAM/OUTNAM
	COMMON/NDD/NDD
	COMMON/DD/DD(NDIM)

	IF (NDD.NE.NDIM)STOP 'SUBR. MODOUTM: ARRAY DIMENSIONS ARE WRONG'

C	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='NEW')          !VMS
	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='UNKNOWN')      !UNIX
	
	IF ((NRAY.GT.1).AND.(IRAY.EQ.1)) 
     1    WRITE (21,*) NRAY, 'NO. OF RAYS = No. of atmos'
	WRITE (21,*) NOUT, '=NO. OF DEPTH PTS', IRAY, '=RAY NUMBER'
	DO I=1, NOUT
	  WRITE (21,9000) TLOG(I), Z(I), TEMP(I), PRESS(I), 
     1                    DENS(I), PEL(I), KAPPA(I)
	END DO
	WRITE (21,*)
	WRITE (21,*)' LOG(TAU)    Z        TEMP      P_GAS      DENS',
     1              '      P_ELEC       KAPPA'

	IF (IRAY.EQ.NRAY) CLOSE (UNIT=21)

C Increased accuracy for pressure has been included below (for accurate B)
9000	FORMAT (2X,F6.2,1X,1PE11.4,1X,0PF7.1,1X,1PE11.4,3(1X,1PE10.3))
C9000	FORMAT (2X,F6.2,1X,1PE11.4,1X,0PF7.1,3(1X,1PE10.3))
	RETURN
	END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE MODOUT (TLOG, Z, TEMP, PEL, KAPPA, PRESS, VTURB, 
     1                     DENS, BFIELD, VEL, GAMMA, CHII, NOUT, WLOUT,
     2                     IRAY, NRAY)

C   The F.T. model is written onto a file. The model has NOUT height points.
C   The variables are explained in the main program.
C   Variables:	TLOG:log. to base 10 of continuum optical depth.
C		Z:	Depth of F.T. atmosphere at TLOG
C		T:	Temp. of F.T. atmosphere
C		PEL:Electr. pressure of F.T. atmosphere
C		KAPPA:Absorption coefficient of F.T. atmosphere

	PARAMETER (NDIM=1024)

	REAL      TLOG(NDIM), Z(NDIM), TEMP(NDIM), PEL(NDIM)
	REAL      BFIELD(NDIM), KAPPA(NDIM), PRESS(NDIM), VTURB(NDIM)
	REAL      DENS(NDIM), VEL(NDIM), GAMMA(NDIM), CHII(NDIM)
	CHARACTER OUTNAM*25

	COMMON/IBOT/IBOT
	COMMON/OUTNAM/OUTNAM
	COMMON/NDD/NDD
	COMMON/NVEL/NVEL, NGAMMA, NCHII

	IF (NDD.NE.NDIM) STOP 'SUBR. MODOUT: ARRAY DIMENSIONS ARE WRONG'

C	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='NEW')          !VMS
	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='UNKNOWN')      !UNIX
	WEIGHT=1.

	IF ((NRAY.GT.1).AND.(IRAY.EQ.1)) 
     1    WRITE (21,*) NRAY, 'NO. OF RAYS = No. of atmos'
	WRITE (21,9010) NOUT, WEIGHT,WLOUT,'=WAVLEN FOR TAU & KAPPA', 
     1                  IRAY, '=RAY NUMBER'
        if ((ngamma.ne.0).or.(nchii.ne.0)) then
          do i=1,nout
            write(21,9000) TLOG(I), Z(I), TEMP(I), PRESS(I),
     1                     PEL(I), KAPPA(I), DENS(I), BFIELD(I),
     2                     VTURB(I), VEL(I), GAMMA(I), CHII(I)
          end do
        elseif (nvel.eq.0) then
          do i=1,nout
	    write(21,9100) TLOG(I), Z(I), TEMP(I), PRESS(I), 
     1                     PEL(I), KAPPA(I), DENS(I), BFIELD(I), 
     2                     VTURB(I), VEL(I)
          end do
        else
          do i=1,nout
            write(21,9200) TLOG(I), Z(I), TEMP(I), PRESS(I), 
     1                     PEL(I), KAPPA(I), DENS(I), BFIELD(I), 
     2                     VTURB(I), VEL(I)
          end do
        end if
	WRITE (21,*)
	IF ((NGAMMA.NE.0).OR.(NCHII.NE.0)) THEN
	  WRITE (21,*) ' LOG(TAU)    Z          TEMP    PRESS  ',
     1                 '  EL PRESS   KAPPA     DENS    ',
     2                 'B   MICROT  V  GAM    CHI'
          WRITE (21,*) '            [CM]         [K]    [CGS]  ',
     1                 '   [CGS]    [GM^-1]  [G/CM^3]  ',
     2                 '[G]     [CM/S]    [DEGR]'
	ELSE
	  WRITE (21,*) ' LOG(TAU)    Z          TEMP    PRESS  ',
     1                 '  EL PRESS   KAPPA     DENS    ',
     2                 'B   MICROT  V'
          WRITE (21,*) '            [CM]         [K]    [CGS]  ',
     1                 '   [CGS]    [GM^-1]  [G/CM^3]  ',
     2                 '[G]     [CM/S]'
	END IF


	IF (IRAY.EQ.NRAY) CLOSE (UNIT=21)

9000	FORMAT (1X,F7.3,1X,1PE12.5,1X,0PF8.1,1X,1PE9.3,1X,1PE9.3,1X,
     1          1PE9.3,1X,1PE9.3,1X,0PF3.1,1X,1PE7.1,1X,1PE9.2,1X,
     2          0PF6.1,1X,0PF6.1)
9100	FORMAT (1X,F7.3,1X,1PE12.5,1X,0PF8.1,1X,1PE9.3,1X,1PE9.3,1X,
     1          1PE9.3,1X,1PE9.3,1X,0PF3.1,1X,1PE7.1,1X,0PF2.0)
9200	FORMAT (1X,F7.3,1X,1PE12.5,1X,0PF8.1,1X,1PE9.3,1X,1PE9.3,1X,
     1          1PE9.3,1X,1PE9.3,1X,0PF3.1,1X,1PE7.1,1X,1PE9.2)
9010    FORMAT (I6,7X,F2.0, 7X,F7.0,A, I5,A)
	RETURN
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE TLOGOUT (TLOG, Z, NOUT, WLOUT, IRAY, NRAY)

C   The F.T. model is written onto a file. The model has NOUT height points.
C   The variables are explained in the main program.
C   Variables:	TLOG:log. to base 10 of continuum optical depth.
C		Z:	Depth of F.T. atmosphere at TLOG
C		T:	Temp. of F.T. atmosphere
C		PEL:Electr. pressure of F.T. atmosphere
C		KAPPA:Absorption coefficient of F.T. atmosphere

	PARAMETER (NDIM=1024)

	REAL      TLOG(NDIM), Z(NDIM)
	CHARACTER OUTNAM*25

	COMMON/IBOT/IBOT
	COMMON/OUTNAM/OUTNAM
	COMMON/NDD/NDD

	IF (NDD.NE.NDIM) STOP 'SUBR. MODOUT: ARRAY DIMENSIONS ARE WRONG'

C	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='NEW')          !VMS
	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='UNKNOWN')      !UNIX
	WEIGHT=1.

	IF ((NRAY.GT.1).AND.(IRAY.EQ.1)) THEN
c       1    WRITE (21,*) NRAY, 'NO. OF RAYS = No. of atmos'
	   WRITE (21,*) NRAY, 'NO. OF RAYS = No. of atmos'
	   WRITE (21,9010) NOUT, WEIGHT,WLOUT,'=WAVLEN FOR TAU & KAPPA', 
     1                  IRAY, '=RAY NUMBER'
	   WRITE (21,*)
	   WRITE (21,*) ' LOG(TAU)    Z '
	END IF
	do i=1,nout
	   write(21,9000) TLOG(I), Z(I)
        end do

c	WRITE (21,*)
c	WRITE (21,*) ' LOG(TAU)    Z '

	IF (IRAY.EQ.NRAY) CLOSE (UNIT=21)

9000	FORMAT (1X,F7.3,1X,1PE12.5,1X,0PF8.1,1X,1PE9.3,1X,1PE9.3,1X,
     1          1PE9.3,1X,1PE9.3,1X,0PF3.1,1X,1PE7.1,1X,1PE9.2,1X,
     2          0PF6.1,1X,0PF6.1)
9100	FORMAT (1X,F7.3,1X,1PE12.5,1X,0PF8.1,1X,1PE9.3,1X,1PE9.3,1X,
     1          1PE9.3,1X,1PE9.3,1X,0PF3.1,1X,1PE7.1,1X,0PF2.0)
9200	FORMAT (1X,F7.3,1X,1PE12.5,1X,0PF8.1,1X,1PE9.3,1X,1PE9.3,1X,
     1          1PE9.3,1X,1PE9.3,1X,0PF3.1,1X,1PE7.1,1X,1PE9.2)
9010    FORMAT (I6,7X,F2.0, 7X,F7.0,A, I5,A)
	RETURN
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	SUBROUTINE STEFFOUT (TLOG, Z, TEMP, PRESS, DENS, PEL, 
     1                     KAPPA, VEL, VELHOR, WLOUT, NOUT, IRAY, NRAY)

	PARAMETER (NDIM=1024)

	REAL      TLOG(NDIM), Z(NDIM), TEMP(NDIM), PEL(NDIM)
	REAL      BFIELD(NDIM), KAPPA(NDIM), PRESS(NDIM)
	REAL      DENS(NDIM), VEL(NDIM), VELHOR(NDIM)
	CHARACTER OUTNAM*25

	COMMON/OUTNAM/OUTNAM
	COMMON/NDD/NDD

	IF (NDD.NE.NDIM) STOP 'SUBR. MODOUT: ARRAY DIMENSIONS ARE WRONG'

C	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='NEW')          !VMS
	IF (IRAY.EQ.1) OPEN (UNIT=21, FILE=OUTNAM, STATUS='UNKNOWN')      !UNIX

	GAMMAGAS=4./3.
	WRITE(21,*) IRAY, NOUT, NRAY,
     1            ' RAY NO., NO. OF Z PTS, NO. OF X PTS'
	DO I=1, NOUT
	  VTOT=SQRT(VEL(I)**2+VELHOR(I)**2)
	  C_S=SQRT(GAMMAGAS*PRESS(I)/DENS(I))
	  WRITE (21,9100) Z(I), VTOT, C_S, VTOT/C_S
	END DO
9100    FORMAT(1X,4(1PE10.3))
	IF (IRAY.EQ.NRAY) CLOSE(21)
	RETURN
	END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE INTOUT (TLOG, TEMP, PRESS, PEL, KAPPA, GRAVAC, 
     1			   NOUT, IRAY, NRAY)

C   The F.T. model is written onto a file. The model has NOUT height points.
C   The variables are explained in the main program.
C   Variables:  TLOG:log. to base 10 of continuum optical depth.
C               T:      Temp. of F.T. atmosphere
C               PEL:Electr. pressure of F.T. atmosphere
C               KAPPA:Absorption coefficient of F.T. atmosphere
c
c   Little new output routine for format needed for intense/spectrum
c					(ycu, May 97)
c
c   Input and common block stuff
c
	parameter (ndim=1024, nmod=25)

	real	  tlog(ndim), Temp(ndim), Pel(ndim), Press(ndim)
	real	  kappa(ndim)
	real	  gravac
	real	  Teff
	integer   nout, iray, nray, ndd
	character outnam*25

	common/outnam/outnam
	common/ndd/ndd
	common/Teff/Teff
c
c    Local and output vars (fix he abundance for now.. I am sure that this
c    is actually read somewhere in the program! !!!!!! BODGE BODGE !!!!!!!
c
	integer	  i
	real	  loggrav, heab

	heab = 0.098
	loggrav = log10(gravac)

        if (ndd.ne.ndim) STOP 'SUBR. INTOUT: ARRAY DIMENSIONS ARE WRONG'

        if (iray.eq.1) open (unit=21, file=outnam, status='unknown')   

c    Should really set up the optical depth grid suitable for intense/spectrum. 
c    Originally, it had 25 depth points starting only at log tau = -4. I think, 
c    however, that the spacing is not crucial. Will therefore try and take the 
c    last 25 depth points of models calculated here to avoid further interpolation.

c	write (21, 9300) Teff, loggrav, heab
	write (21, 9300) Teff, gravac, heab

	do i = nout-nmod+1, nout
	   write (21, 9310)
     1		 10**tlog(i), Temp(i), Press(i), Pel(i), kappa(i)
	end do

	if (iray .eq. nray) close (unit=21)

 9300	format(1p3e12.5)
 9310	format(1p5e12.5)

	return
	end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE KUROUT (COLMAS, TEMP, PRESS, NEL, KAPPA, VTURB, 
     1			   NOUT, IRAY, NRAY)

C   The F.T. model is written onto a file. The model has NOUT height points.
C   The variables are explained in the main program.
c
c   Little new output routine for format needed for Kurucz output
c
c   Input and common block stuff
c
	INTEGER   NDIM
        PARAMETER (NDIM=1024)

        REAL      COLMAS(NDIM), TEMP(NDIM), NEL(NDIM), PRESS(NDIM)
        REAL     VTURB(NDIM), KAPPA(NDIM), X
	INTEGER   IRAY, NRAY, NDD, IM
        CHARACTER OUTNAM*25

        COMMON/OUTNAM/OUTNAM
        COMMON/NDD/NDD

        IF (NDD.NE.NDIM) STOP 'SUBR. INTOUT: ARRAY DIMENSIONS ARE WRONG'
c	PRINT *, "print statement output"
        IF (IRAY.EQ.1) then
           OPEN (UNIT=21, FILE=OUTNAM, STATUS='UNKNOWN')
           WRITE (21,*) nray
           WRITE (21,*) nout
        end if
        X = 0.5
c
	DO IM = 1, NOUT
           WRITE (21, 9300) COLMAS(IM), TEMP(IM), PRESS(IM), NEL(IM), 
     1                      KAPPA(IM), X, abs(VTURB(IM))
        END DO

        IF (IRAY .EQ. NRAY) CLOSE (UNIT=21)

 9300   FORMAT ( (1pe15.8, 0pf9.1, 1p5e10.3))

        RETURN
        END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	SUBROUTINE INTPOL (TLOG, Z, TEMP, PEL, KAPPA, PRESS, VTURB, 
     1                     DENS, BFIELD, VEL, BOUNDS, GAMMA, CHII, NOUT)

C Depending on the variable, the interpolation is linear  in TLOG (for Z, TEMP,
C BFIELD, VTURB, VEL), or logarithmic (for PEL, PRESS,  KAPPA, DENS). 

C IBOT:    the number of points per step
C ITOT:    Total number of time steps

	PARAMETER (NDIM=1024)

	REAL      TLOG(NDIM), Z(NDIM), TEMP(NDIM), PEL(NDIM)
	REAL      BFIELD(NDIM), KAPPA(NDIM), PRESS(NDIM), VTURB(NDIM)
	REAL      DENS(NDIM), VEL(NDIM), ATLOG(NDIM), AZ(NDIM)
	REAL      ATEMP(NDIM), AVEL(NDIM), ABFIELD(NDIM), AVTURB(NDIM)
	REAL      ALOGKAPPA(NDIM), ALOGPEL(NDIM), ALOGPRESS(NDIM)
	REAL      ALOGDENS(NDIM), GAMMA(NDIM), AGAMMA(NDIM)
	REAL      CHII(NDIM), ACHII(NDIM)
	CHARACTER BOUNDS*1
	CHARACTER TLOGFILE*25
	
	COMMON/IBOT/IBOT
	COMMON/TBOUND/ TUPPER, TLOWER
	COMMON/NDD/NDD

	IF (NDD.NE.NDIM) STOP 'SUBR. INTPOL: ARRAY DIMENSIONS ARE WRONG'

	DO I=1, IBOT
	  ATLOG(I)=TLOG(I)
	  AZ(I)=Z(I)
	  ATEMP(I)=TEMP(I)
  	  ABFIELD(I)=BFIELD(I)
	  AVTURB(I)=VTURB(I)
	  AVEL(I)=VEL(I)
	  AGAMMA(I)=GAMMA(I)
	  ACHII(I)=CHII(I)
	  IF ((PEL(I).NE.0.0).AND.(KAPPA(I).NE.0.0)) THEN
	    ALOGPEL(I)=LOG(PEL(I))
	    ALOGPRESS(I)=LOG(PRESS(I))
	    ALOGKAPPA(I)=LOG(KAPPA(I))
	    ALOGDENS(I)=LOG(DENS(I))
	  ELSE IF (I.LE.IBOT) THEN
	    IF (PEL(I).EQ.0.0) THEN
	      PRINT*,'ELECTRON PRESSURE IS ZERO. I=', I
	    ELSE
	      PRINT*,'ABSRPTION COEFFICIENT IS ZERO. I=', I
	    END IF
	    STOP 
	  END IF
	END DO

	IF (BOUNDS.EQ.'F') THEN
	   TLOGFILE = "tlogfile.dat"
	ELSE IF (BOUNDS.EQ.'T') THEN
	  DELTL=ABS(TLOWER-TUPPER)/REAL(NOUT-1)
	  TLOG(1)=TUPPER
	ELSE
	  DELTL=(ATLOG(IBOT)-ATLOG(1))/REAL(NOUT-1)
	END IF

	IF (BOUNDS.EQ.'F') THEN
	   OPEN(UNIT=41, FILE=TLOGFILE,STATUS='OLD')
	   DO I=1, NOUT
	      READ(41,*) TLOG(I)
	   END DO
	ELSE
	   DO I=2, NOUT
	      TLOG(I)=TLOG(1)+REAL(I-1)*DELTL
	   END DO
	END IF
c	PRINT *, TLOG(1),' ',TLOG(NOUT)
	
	CLOSE(41)
c	PRINT *,'OLD TLOG(IBOT):',ATLOG(IBOT),
c     1        '  NEW TLOG(NOUT):',TLOG(NOUT)

	IF (BOUNDS.EQ.'T') THEN
	  IF (TUPPER.EQ.ATLOG(1)) THEN
	    IBEG=2
	  ELSE
	    IBEG=1
	  END IF
	  IF (TLOWER.EQ.ATLOG(IBOT)) THEN
	    IEND=NOUT-1
	  ELSE
	    IEND=NOUT
	 END IF
	ELSE IF (BOUNDS.EQ.'F') THEN
	   IBEG = 1
	   IEND = NOUT
	ELSE
	  IBEG=2
	  IEND=NOUT-1
	END IF

	J=1
	DO I=IBEG, IEND
	  DO WHILE (TLOG(I).GT.ATLOG(J))
	    J=J+1
	  END DO
	  X=(TLOG(I)-ATLOG(J-1))/(ATLOG(J)-ATLOG(J-1))
	  Z(I)=AZ(J-1)*(1.-X)+AZ(J)*X
	  TEMP(I)=ATEMP(J-1)*(1.-X)+ATEMP(J)*X
	  BFIELD(I)=ABFIELD(J-1)*(1.-X)+ABFIELD(J)*X
	  VTURB(I)=AVTURB(J-1)*(1.-X)+AVTURB(J)*X
	  VEL(I)=AVEL(J-1)*(1.-X)+AVEL(J)*X
	  GAMMA(I)=AGAMMA(J-1)*(1.-X)+AGAMMA(J)*X
	  CHII(I)=ACHII(J-1)*(1.-X)+ACHII(J)*X
	  PEL(I)=EXP(ALOGPEL(J-1)*(1.-X)+ALOGPEL(J)*X)
	  PRESS(I)=EXP(ALOGPRESS(J-1)*(1.-X)+ALOGPRESS(J)*X)
	  KAPPA(I)=EXP(ALOGKAPPA(J-1)*(1.-X)+ALOGKAPPA(J)*X)
	  DENS(I)=EXP(ALOGDENS(J-1)*(1.-X)+ALOGDENS(J)*X)
	END DO

	IF (IEND.EQ.NOUT-1) THEN
	  TLOG(NOUT)=ATLOG(IBOT)
	  Z(NOUT)=AZ(IBOT)
	  TEMP(NOUT)=ATEMP(IBOT)
	  BFIELD(NOUT)=ABFIELD(IBOT)
	  VTURB(NOUT)=AVTURB(IBOT)
	  VEL(NOUT)=AVEL(IBOT)
	  GAMMA(NOUT)=AGAMMA(IBOT)
	  CHII(NOUT)=ACHII(IBOT)
	  PEL(NOUT)=EXP(ALOGPEL(IBOT))
	  PRESS(NOUT)=EXP(ALOGPRESS(IBOT))
	  KAPPA(NOUT)=EXP(ALOGKAPPA(IBOT))
	  DENS(NOUT)=EXP(ALOGDENS(IBOT))
	END IF

	RETURN
	END	

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      FUNCTION LENSTR(STRING)
 
C Function LENSTR determines the length of a string
 
      CHARACTER  STRING*(*)
      INTEGER    I

      I=1
10    IF (STRING(I:I).EQ.' ') GOTO 20
      I=I+1
      IF (I.GT.LEN(STRING)) GOTO 20
      GOTO 10
 
20    LENSTR=I-1
 
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	SUBROUTINE PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL, KAPPA,       !MODCON
     1                    RHO, NPEL, NKAPPA, IBEG, IEND, NDD)           !MODCON
C	SUBROUTINE PEABS (TEMP, PRESS, WAVANS, WAVLEN, PEL, KAPPA,       !INVERT
C     1                    RHO, NPEL, NKAPPA, IBEG, IEND)                !INVERT

C Calls the PEMAKE and ABSKO subr. to determine electron pressure and 
C continuum absorption coefficient. Requires file MODCONINP to function 
C properly. Since passing on the wavelength is such a mess, a test has been 
C introduced using TESTWL. TESTWL is filled with WLOUT, and it is tested in 
C ABSKO whether TESTWL corresponds to the wavelength for the ISETA and JJ 
C passed on to ABSKO.

C If this subr. is called from MODCON, then comment out all lines ending with
C !INVERT, if called by INVERT, comment out lines with !MODCON

C Written by S.K. Solanki, version of may 1993

	PARAMETER (NDIM=1024)                                            !MODCON
C        INCLUDE 'invpar.inc'                                           !INVERT

	REAL      PRESS(NDIM), TEMP(NDIM), KAPPA(NDIM), PEL(NDIM)
	REAL      RHO(NDIM)
	CHARACTER WAVANS*1, CONANS*1

	COMMON/IBOT/IBOT
	COMMON/NEWT/NEWT
C	COMMON/COUTR/NTO, NTPO(10)
	COMMON/TESTWL/TESTWL

	IF (NDIM.NE.NDD) THEN                                           !MODCON
	  PRINT*,'INCONSISTENCY IN ARRAY DIMENSIONS: SUBR. PEABS'       !MODCON
	  STOP                                                          !MODCON
	END IF                                                          !MODCON

	IF ((IBEG.LT.1).OR.(IBEG.GT.IEND).OR.(IEND.GT.IBOT)) THEN
	  PRINT*,'PROBLEMS WITH IBEG & IEND IN SUBR. PEABS. IBEG,IEND',
     1            IBEG, IEND
	  STOP
	END IF

	  IF (WAVANS.EQ.'N') THEN
C ISETA: number of wavelength set (in the internal representation of ABSKO)
	    ISETA=2                 
C JJ: number of wavelength in the wavelength set (internal to ABSKO)
	    JJ=5                   
	    WLOUT=5000.
	  ELSE
	    IF ((WAVLEN.GT.4150.).AND.(WAVLEN.LE.7100.)) ISETA=2
	    IF ((WAVLEN.GT.7100.).AND.(WAVLEN.LE.13250.)) ISETA=3
	    IF ((WAVLEN.GT.15000.).AND.(WAVLEN.LE.52000.)) ISETA=4
	    IF ((WAVLEN.GT.3700.).AND.(WAVLEN.LE.4150.)) ISETA=4

	    IF ((WAVLEN.GT.4150.).AND.(WAVLEN.LE.4300.)) THEN
	      JJ=1
	      WLOUT=4200.
	    ELSE IF ((WAVLEN.GT.4300.).AND.(WAVLEN.LE.4500.)) THEN
	      JJ=2
	      WLOUT=4400.
	    ELSE IF ((WAVLEN.GT.4500.).AND.(WAVLEN.LE.4700.)) THEN
	      JJ=3
	      WLOUT=4600.
	    ELSE IF ((WAVLEN.GT.4700.).AND.(WAVLEN.LE.4900.)) THEN
	      JJ=4
	      WLOUT=4800.
	    ELSE IF ((WAVLEN.GT.4900.).AND.(WAVLEN.LE.5100.)) THEN
	      JJ=5
	      WLOUT=5000.
	    ELSE IF ((WAVLEN.GT.5100.).AND.(WAVLEN.LE.5300.)) THEN
	      JJ=6
	      WLOUT=5200.
	    ELSE IF ((WAVLEN.GT.5300.).AND.(WAVLEN.LE.5500.)) THEN
	      JJ=7
	      WLOUT=5400.
	    ELSE IF ((WAVLEN.GT.5500.).AND.(WAVLEN.LE.5700.)) THEN
	      JJ=8
	      WLOUT=5600.
	    ELSE IF ((WAVLEN.GT.5700.).AND.(WAVLEN.LE.5900.)) THEN
	      JJ=9
	      WLOUT=5800.
	    ELSE IF ((WAVLEN.GT.5900.).AND.(WAVLEN.LE.6100.)) THEN
	      JJ=10
	      WLOUT=6000.
	    ELSE IF ((WAVLEN.GT.6100.).AND.(WAVLEN.LE.6300.)) THEN
	      JJ=11
	      WLOUT=6200.
	    ELSE IF ((WAVLEN.GT.6300.).AND.(WAVLEN.LE.6500.)) THEN
	      JJ=12
	      WLOUT=6400.
	    ELSE IF ((WAVLEN.GT.6500.).AND.(WAVLEN.LE.6700.)) THEN
	      JJ=13
	      WLOUT=6600.
	    ELSE IF ((WAVLEN.GT.6700.).AND.(WAVLEN.LE.6900.)) THEN
	      JJ=14
	      WLOUT=6800.
	    ELSE IF ((WAVLEN.GT.6900.).AND.(WAVLEN.LE.7100.)) THEN
	      JJ=15
	      WLOUT=7000.
	    ELSE IF ((WAVLEN.GT.7100.).AND.(WAVLEN.LE.7300.)) THEN
	      JJ=1
	      WLOUT=7200.
	    ELSE IF ((WAVLEN.GT.7300.).AND.(WAVLEN.LE.7500.)) THEN
	      JJ=2
	      WLOUT=7400.
	    ELSE IF ((WAVLEN.GT.7500.).AND.(WAVLEN.LE.7700.)) THEN
	      JJ=3
	      WLOUT=7600.
	    ELSE IF ((WAVLEN.GT.7700.).AND.(WAVLEN.LE.7900.)) THEN
	      JJ=4
	      WLOUT=7800.
	    ELSE IF ((WAVLEN.GT.7900.).AND.(WAVLEN.LE.8250.)) THEN
	      JJ=5
	      WLOUT=8000.
	    ELSE IF ((WAVLEN.GT.8250.).AND.(WAVLEN.LE.8750.)) THEN
	      JJ=6
	      WLOUT=8500.
	    ELSE IF ((WAVLEN.GT.8750.).AND.(WAVLEN.LE.9250.)) THEN
	      JJ=7
	      WLOUT=9000.
	    ELSE IF ((WAVLEN.GT.9250.).AND.(WAVLEN.LE.9750.)) THEN
	      JJ=8
	      WLOUT=9500.
	    ELSE IF ((WAVLEN.GT.9750.).AND.(WAVLEN.LE.10250.)) THEN
	      JJ=9
	      WLOUT=10000.
	    ELSE IF ((WAVLEN.GT.10250.).AND.(WAVLEN.LE.10750.)) THEN
	      JJ=10
	      WLOUT=10500.
	    ELSE IF ((WAVLEN.GT.10750.).AND.(WAVLEN.LE.11250.)) THEN
	      JJ=11
	      WLOUT=11000.
	    ELSE IF ((WAVLEN.GT.11250.).AND.(WAVLEN.LE.11750.)) THEN
	      JJ=12
	      WLOUT=11500.
	    ELSE IF ((WAVLEN.GT.11750.).AND.(WAVLEN.LE.12250.)) THEN
	      JJ=13
	      WLOUT=12000.
	    ELSE IF ((WAVLEN.GT.12250.).AND.(WAVLEN.LE.12750.)) THEN
	      JJ=14
	      WLOUT=12500.
	    ELSE IF ((WAVLEN.GT.12750.).AND.(WAVLEN.LE.13250.)) THEN
	      JJ=15
	      WLOUT=13000.
	    ELSE IF ((WAVLEN.GT.14500.).AND.(WAVLEN.LE.15500.)) THEN
	      JJ=1
	      WLOUT=15250.
	    ELSE IF ((WAVLEN.GT.15500.).AND.(WAVLEN.LE.16000.)) THEN
	      JJ=2
	      WLOUT=15750.
	    ELSE IF ((WAVLEN.GT.16000.).AND.(WAVLEN.LE.16500.)) THEN
	      JJ=3
	      WLOUT=16250.
	    ELSE IF ((WAVLEN.GT.16500.).AND.(WAVLEN.LE.17000.)) THEN
	      JJ=4
	      WLOUT=16750.
	    ELSE IF ((WAVLEN.GT.17000.).AND.(WAVLEN.LE.18000.)) THEN
	      JJ=5
	      WLOUT=17500.
	    ELSE IF ((WAVLEN.GT.19000.).AND.(WAVLEN.LE.21000.)) THEN
	      JJ=6
	      WLOUT=20000.
	    ELSE IF ((WAVLEN.GT.21000.).AND.(WAVLEN.LE.23000.)) THEN
	      JJ=7
	      WLOUT=22000.
	    ELSE IF ((WAVLEN.GT.23000.).AND.(WAVLEN.LE.25000.)) THEN
	      JJ=8
	      WLOUT=24000.
	    ELSE IF ((WAVLEN.GT.35000.).AND.(WAVLEN.LE.40000.)) THEN
	      JJ=9
	      WLOUT=37000.
	    ELSE IF ((WAVLEN.GT.44000.).AND.(WAVLEN.LE.48000.)) THEN
	      JJ=10
	      WLOUT=46000.
	    ELSE IF ((WAVLEN.GT.48000.).AND.(WAVLEN.LE.52000.)) THEN
	      JJ=11
	      WLOUT=50000.
	    ELSE IF ((WAVLEN.GT.3700.).AND.(WAVLEN.LE.3850.)) THEN
	      JJ=12
	      WLOUT=3800.
	    ELSE IF ((WAVLEN.GT.3850.).AND.(WAVLEN.LE.3950.)) THEN
	      JJ=13
	      WLOUT=3900.
	    ELSE IF ((WAVLEN.GT.3950.).AND.(WAVLEN.LE.4050.)) THEN
	      JJ=14
	      WLOUT=4000.
	    ELSE IF ((WAVLEN.GT.4050.).AND.(WAVLEN.LE.4150.)) THEN
	      JJ=15
	      WLOUT=4100.
	    ELSE
	      PRINT*,'SORRY, THE DESIRED WAVELENGTH IS CURRENTLY NOT ',
     1               'AVAILABLE IN MODCONINP'
	      PRINT*,'PLEASE CHANGE THAT FILE AND THE CODE, OR TRY ',
     1               'ANOTHER WAVELENGTH'
	      STOP
	    END IF
	  END IF
	  TESTWL=WLOUT
	  PEIN=1.
	  NT=1
	  IBAD=-1
	  CONANS='N'

	  DO 20 I=IBEG, IEND
C The following construct is only to keep the code from crashing whenever
C the gas pressure gets too small.
	    IF (I.EQ.IBAD+1) THEN
	      NPEL=NPKEEP
	      NKAPPA=NKKEEP
	      IBAD=-1
	    END IF
	    IF (PRESS(I).LT.0.01) THEN
	      IF (CONANS.NE.'Y') THEN
	        PRINT*,'LOW GAS PRESSURE', press(i), ' DO YOU WANT TO ', 
     1		' CONTINUE, BUT WITHOUT PE AND KAP CALC? (Y/N,DEF=N)'
                READ '(A)', CONANS
	      END IF
	      IF (CONANS.EQ.'Y') THEN
	        IF (NPEL.NE.99999) NPKEEP=NPEL
	        NPEL=99999
		PEL(I)=1.E-10
	        IF (NKAPPA.NE.99999) NKKEEP=NKAPPA
	        NKAPPA=99999
		KAPPA(I)=1.E-10
		IBAD=I
	      ELSE
	        PRINT*,'STOPPING DUE TO SMALL GAS PRESSURE (SUBR. PEABS)'
	        STOP
	      END IF
	    END IF

	    TIN=TEMP(I)

c	PRINT*,'PEABS:I,IBEG,IEND',I,IBEG,IEND
c	PRINT*,'PEABS:NPEL,NKAPPA',NPEL,NKAPPA
c	PRINT*,'PEABS:PEL',PEL(I)
c	PRINT*,'PEABS:TEMP,PRESS', TEMP(I), PRESS(I)
	    IF (NPEL.EQ.0) THEN
c	      print *, I, ':'
	      PG=PRESS(I)
	      CALL PEMAKE (TIN, PEIN, PG, PE)
	      PEIN=PE
	      PEL(I)=PE
	    ELSE
	      PE=PEL(I)
	    END IF
c	    print *,'tin,pein,pg,pe = ', tin,pein,pg,pe
c	PRINT*,'NEWT,NT,TIN,PE,ISETA,JJ,RO,ABSC,SCATC ',
c     1          NEWT,NT,TIN,PE,ISETA,JJ,RO,ABSC,SCATC 
c	    print *, 'NKAPPA,WAVANS', NKAPPA,WAVANS
	    IF ((NKAPPA.EQ.0).OR.(WAVANS.EQ.'Y')) THEN
c	      print *, 'ro in peabs =', ro
	      CALL ABSKO (NEWT, NT, TIN, PE, ISETA, JJ,RO, ABSC, SCATC)
c	      print *, 'ro in peabs =', ro
              KAPPA(I)=ABSC+SCATC	
c	      KAPPA(I)=ABSC		! Does it change ??  Yep.. but only for atmosph with T-inversion
	      RHO(I)=RO
	      NEWT=1
	    END IF
20      CONTINUE

	RETURN
	END


