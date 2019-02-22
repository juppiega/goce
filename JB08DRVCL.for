C     JB08DRVCL TEST CASE DRIVER FOR RUNNING JB2008
C
C     TEST PROGRAM DEVELOPED BY Bruce R. Bowman - June 2008
C
C     JB08DRVCL       Driver generates JB2008 density values
C                    *Inputs time and location data, outputs density and other parameters
C       JB2008        Density model in separate Fortran source file
C       SOLFSMY       Returns F10,S10,M10,Y10 and 81-centered values
C                    *Reads solar flux file
C       DTCVAL        Returns dTc value from input time
C                    *Reads dTc file and loads DTC common
C       SUNPOS        Computes sun position
C       THETA         Computes Greenwich hour angle
C
C     * Following files used:
C
C       DTCFILE.TXT input  unit 46 in DTCVAL    contains dTc values
C       DSTVAL.GAP  output unit 47 in DTCVAL    contains Dst gap time
C       SOLFSMY.TXT input  unit 24 in SOLFSMY   contains solar flux
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C
      DIMENSION SUN(2),SAT(3),TEMP(2)
      CHARACTER*80 INPFILE,OUTFILE
C
      PARAMETER (PI = 3.1415927D0)
      PARAMETER (DG2RD = PI / 180.D0)
C
      CHARACTER*5 IYRVAL,IDAYVAL,IHRVAL,IMINVAL
      CHARACTER*11 SECVAL
      CHARACTER*12 XHTVAL
      CHARACTER*11 XLATVAL,XLONVAL

      INTEGER IYR,IDAY,IHR,IMIN
      REAL SEC,XHT,XLAT,XLON

C     READ INPUT TEST DATA FROM COMMAND LINE

      CALL GETARG(1, IYRVAL)
      READ (IYRVAL,'(I5)') IYR
      CALL GETARG(2, IDAYVAL)
      READ (IDAYVAL,'(I5)') IDAY
      CALL GETARG(3, IHRVAL)
      READ (IHRVAL,'(I5)') IHR
      CALL GETARG(4, IMINVAL)
      READ (IMINVAL,'(I5)') IMIN
      CALL GETARG(5, SECVAL)
      READ (SECVAL,'(F7.3)') SEC
      CALL GETARG(6, XHTVAL)
      READ (XHTVAL,'(F9.2)') XHT
      CALL GETARG(7, XLATVAL)
      READ (XLATVAL,'(F8.2)') XLAT
      CALL GETARG(8, XLONVAL)
      READ (XLONVAL,'(F8.2)') XLON


C     DUMMY STATEMENT.
  10  WRITE(6,*)
C
C     CONVERT TIME TO DAYS SINCE 1950 JAN 0 AND MJD
      IF (IYR.GT.1000) IYR = IYR - 1900
      IF (IYR.LT.50) IYR = IYR + 100
      IYY = ((IYR-1)/4-12)
      IYY = (IYR-50)*365 + IYY
      D1950 = IYY*1.D0 + IDAY*1.D0+IHR/24.D0+IMIN/1440.D0+SEC/86400.D0
      AMJD = D1950 + 33281.0D0

C
C     READ SOLAR INDICES
C     USE 1 DAY LAG FOR F10 AND S10 FOR JB2008
      DLAG = 1.
      T1950 = D1950 - DLAG
      CALL SOLFSMY (T1950,F10,F10B,S10,S10B,XMXX,XMXXB,XYXX,XYXXB)
      IF ((F10.LT.40.).OR.(F10B.LT.40.))  GO TO 10
      IF ((S10.LT.40.).OR.(S10B.LT.40.))  GO TO 10

C      IYR = 2000
C
C     USE 2 DAY LAG FOR M10 FOR JB2008
      DLAG = 2.
      T1950 = D1950 - DLAG
      CALL SOLFSMY (T1950,XFXX,XFXXB,XSXX,XSXXB,XM10,XM10B,XYXX,XYXXB)
      IF ((XM10.LT.40.).OR.(XM10B.LT.40.))  GO TO 10


C
C     USE 5 DAY LAG FOR Y10 FOR JB2008
      DLAG = 5.
      T1950 = D1950 - DLAG
      CALL SOLFSMY (T1950,XFXX,XFXXB,XSXX,XSXXB,XMXX,XMXXB,Y10,Y10B)
      IF ((Y10.LT.40.).OR.(Y10B.LT.40.))  GO TO 10
C
C     READ GEOMAGNETIC STORM DTC VALUE
      CALL DTCVAL (D1950,IDTCVAL)
      DSTDTC = IDTCVAL
C
C     CONVERT POINT OF INTEREST LOCATION (RADIANS AND KM)
C     CONVERT LONGITUDE TO RA
      GWRAS  = THETA(D1950)
      SAT(1) = DMOD(GWRAS + XLON*DG2RD + 2.D0*PI, 2.D0*PI)
      SAT(2) = XLAT * DG2RD
      SAT(3) = XHT
C
C     GET SUN LOCATION (RAD)
C
      CALL SUNPOS(AMJD,SOLRAS,SOLDEC)
      SUN(1) = SOLRAS
      SUN(2) = SOLDEC
C
C     COMPUTE DENSITY KG/M3 RHO
C
      CALL JB2008 (AMJD,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,
     *             Y10,Y10B,DSTDTC,TEMP,RHO)
C
C      WRITE(10,101) IYR,IDAY,IHR,IMIN,SEC,XHT,XLAT,XLON,SUN,TEMP,RHO
C  101 FORMAT(4I5,F7.3,F9.2,2F8.2,2F10.4,2F8.1,D15.5)
C      WRITE(10,102) D1950,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC
C  102 FORMAT(F15.5,9F8.0/)

C     THE ABOVE WRITE STATEMENTS HAVE BEEN REPLACED WITH MULTIPLE WRITE STATEMENTS TO PRODUCE A COLUMN VECTOR .TXT FILE. 
C     FIRST LINE FROM ABOVE


C      UNCOMMENT THESE TO OUTPUT DATA TO THE COMMAND LINE.
        WRITE(6,'(I5)') IYR
        WRITE(6,'(I5)') IDAY
        WRITE(6,'(I5)') IHR
        WRITE(6,'(I5)') IMIN
        WRITE(6,'(F7.3)') SEC
        WRITE(6,'(F9.2)') XHT
        WRITE(6,'(F8.2)') XLAT
        WRITE(6,'(F8.2)') XLON
        WRITE(6,'(2F10.4)') SUN
        WRITE(6,'(2F8.1)') TEMP
        WRITE(6,'(D15.5)') RHO
C        WRITE(6,*)
C C C     SECOND LINE FROM ABOVE
        WRITE(6,'(F15.5)') D1950
        WRITE(6,'(F8.0)') F10
        WRITE(6,'(F8.0)') F10B
         WRITE(6,'(F8.0)') S10
         WRITE(6,'(F8.0)') S10B
         WRITE(6,'(F8.0)') XM10
         WRITE(6,'(F8.0)') XM10B
         WRITE(6,'(F8.0)') Y10
         WRITE(6,'(F8.0)') Y10B
         WRITE(6,'(F8.0)') DSTDTC
C     NOTE THAT RHO IS A DOUBLE-PRECISION OUTPUT, BUT MATLAB'S DEFAULT TYPE IS DOUBLE ANYWAY, SO IT WORKS JUST FINE! IF MATLAB CHANGES THIS, IT MAY NEED TO BE UPDATED.

C      WRITE(6,'(D15.5)') RHO

C       GO TO 10
C C
C    90 STOP
      END
C
C***********************************************************************
C
      SUBROUTINE DTCVAL (T1950,IDTCVAL)
C
C     INPUT T1950 AND READ FILE FOR VALUE FOR DTC
C     DTC ON 1 HR INCREMENT, INTERPOLATE FOR OUTPUT VALUE
C     PROGRAM WILL ABORT IF MISSING DATA FOUND
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*50 FDTC,FGAP
C
      COMMON /DTCDATA/ DT1950,DTDTC,IDTC(500000),IDTCRD,IDTCNM
      DIMENSION IDTCH(24)
C
C     SET UP I/O FILES
C
      KIUNIT = 46
      KGUNIT = 47
C
      FDTC = 'DTCFILE.TXT'
      FGAP = 'DTCVAL.GAP'
C
      IF (IDTCRD.NE.1) THEN
C
C       INITIALIZE DTC FILE AND DATA
        OPEN (KIUNIT,FILE=FDTC,ACCESS='SEQUENTIAL',STATUS='OLD')
C
C       INPUT DTC DATA
C       FORMAT EXAMPLE:
C   YYDDD DTC(DEG)
CTC 97125 010 015 025 035 060 066 078 105 120 135 150 145 141 132 122 104 086 076 045 040 034 045 079 096
C
        IDTCRD = 1
        IDTCNM = 0
        TDATA = 0.D0
        DTDIF = 1.01D0/24.D0
        DTDTC = 1.00D0/24.D0
        KNUM = 0
C
   10   READ (KIUNIT,101,END=30) NYR,NDAY,(IDTCH(I),I=1,24)
  101   FORMAT(3X,I3,I3,24I4)
C
        IF (NYR.LT.50) NYR = NYR + 100
        IYY = ((NYR-1)/4-12)
        IYY = (NYR-50)*365 + IYY
        D1950X = IYY + NDAY
        IF (KNUM.EQ.0) DT1950 = D1950X
C
        DO 20 I=1,24
          KNUM = KNUM + 1
          IHR  = I - 1
          D1950 = D1950X + IHR/24.D0
C
          DT = D1950 - TDATA
C
C         CHECK FOR DATA GAPS OR MISSING DATA
          IF (KNUM.GT.1.AND.(IDTCH(I).GT.2000.OR.DT.GT.DTDIF)) THEN
            OPEN (KGUNIT,FILE=FGAP,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
            CALL TMOUTD(D1950,NYR,DAY)
            WRITE (KGUNIT,104) D1950,NYR,DAY,IDTCH(I),DT1950
  104       FORMAT(F15.5,I5,F12.5,I10,F15.5)
            STOP
          ENDIF
C         LOAD DTC DATA INTO COMMON
          IDTC(KNUM)  = IDTCH(I)
          TDATA = D1950
   20   CONTINUE
        GO TO 10
C
   30   IDTCNM = KNUM
      ENDIF
C
      INDX = 1. + (T1950 - DT1950 + 0.0000001D0)*24.
C
      IF (INDX.LT.1.OR.INDX.GE.IDTCNM) THEN
C       ERROR IN INDEX
        WRITE(*,110)
  110   FORMAT('  DTCVAL ERROR - REQUESTED TIME OUTSIDE OF DTC VALUES ')
        STOP
      ENDIF
C
      IF (INDX.EQ.1) THEN
        IDTCVAL = IDTC(1)
        GO TO 50
      ENDIF
C
C     INTERPOLATE FOR TEMPERATURE VALUE
      IDTCVAL1 = IDTC(INDX)
      IDTCVAL2 = IDTC(INDX+1)
C
      FX = 1. + (T1950 - DT1950 + 0.0000001D0)*24.
      FAC = (FX - INDX)/1.D0
      IDTCVAL  = IDTCVAL1  + FAC*(IDTCVAL2-IDTCVAL1) + 0.5
C
   50 RETURN
      END SUBROUTINE DTCVAL
C
C***********************************************************************
C
      SUBROUTINE SOLFSMY (T1950,XF10,XF10B,XS10,XS10B,XM10,XM10B,
     *                          XY10,XY10B)
C
C     INPUT T1950 AND READ FILE FOR VALUES FOR F10, S10, M10, AND Y10
C      READ ONE TIME AND STORE ALL FILE VALUES IN COMMON FOR RETRIEVAL
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*50 FLUXUNIT
      CHARACTER*1  ADESC,A1SRC,A2SRC,A3SRC,A4SRC
C
      COMMON /CSOLFSMY/ FS1950,VF10(10000),VF10B(10000),
     *               VS10(10000),VS10B(10000),VM10(10000),VM10B(10000),
     *               VY10(10000),VY10B(10000),
     *               ISOLFSMRD,ISOLFSMTOT
C
C     INITIALIZE SOLAR FLUX VALUES
      KFUNIT = 24
      FLUXUNIT = 'SOLFSMY.TXT'
C
      IF (ISOLFSMRD.NE.1) THEN
C
C       INITIALIZE FLUX FILE AND DATA
        OPEN (KFUNIT,FILE=FLUXUNIT,ACCESS='SEQUENTIAL',STATUS='OLD')
C       INPUT 4 LINES OF DATA DESCRIPTION
        DO 5 I=1,4
    5   READ (KFUNIT,100,ERR=80,END=80) ADESC
  100   FORMAT(A1)
C
C       INPUT EXAMPLE:
C 2002  80   2452355.0 174.1 188.4 195.3 193.6 182.3 187.8 155.1 156.0  1F1D
C
        READ (KFUNIT,101,ERR=80,END=80) IY,IDAY,XF10,XF10B,XS10,
     *                                  XS10B,XM10,XM10B,XY10,XY10B,
     *                                  A1SRC,A2SRC,A3SRC,A4SRC
        IF (IY.LT.50) IY = IY + 100
        IYY = (IY-1)/4 - 12
        FS1950 = (IY-50)*365 + IYY + IDAY - 1
C
        BACKSPACE (KFUNIT)
        ISOLFSMRD = 1
        ISOLFSMTOT = 0
C
   10   READ (KFUNIT,101,ERR=80,END=20) IY,IDAY,XF10,XF10B,XS10,
     *                                  XS10B,XM10,XM10B,XY10,XY10B,
     *                                  A1SRC,A2SRC,A3SRC,A4SRC
  101   FORMAT(4X,I2,1X,I3,12X,8F6.1,2X,4A1)
        IF (IY.LT.50) IY = IY + 100
        IYY = (IY-1)/4 - 12
        D1950 = (IY-50)*365 + IYY + IDAY
        INDX = D1950 - FS1950
        ISOLFSMTOT = ISOLFSMTOT + 1
        IF (INDX.NE.ISOLFSMTOT) THEN
          WRITE(*,201) IY,IDAY
  201     FORMAT(' DATA GAP IN FLUX DATA',5X,I2,I4)
          STOP
        ENDIF
C
C       LOAD FLUX DATA INTO COMMON
        VF10(INDX)  = XF10
        VF10B(INDX) = XF10B
        VS10(INDX)  = XS10
        VS10B(INDX) = XS10B
        VM10(INDX)  = XM10
        VM10B(INDX) = XM10B
        VY10(INDX)  = XY10
        VY10B(INDX) = XY10B
C
        GO TO 10
C
      ENDIF
C
   20 INDX = T1950 - FS1950
C
      IF (INDX.LT.1.OR.INDX.EQ.ISOLFSMTOT) THEN
         WRITE(*,103)
  103    FORMAT(' TIME OUTSIDE FILE START - STOP TIME')
         STOP
      ENDIF
C
      XF10  = VF10(INDX)
      XF10B = VF10B(INDX)
      XS10  = VS10(INDX)
      XS10B = VS10B(INDX)
      XM10  = VM10(INDX)
      XM10B = VM10B(INDX)
      XY10  = VY10(INDX)
      XY10B = VY10B(INDX)
C
      IF (VF10(INDX).LT.40..OR.VF10B(INDX).LT.40.) THEN
        XF10  = 0.
        XF10B = 0.
      ENDIF
      IF (VS10(INDX).LT.40..OR.VS10B(INDX).LT.40.) THEN
        XS10  = 0.
        XS10B = 0.
      ENDIF
      IF (VM10(INDX).LT.40..OR.VM10B(INDX).LT.40.) THEN
        XM10  = 0.
        XM10B = 0.
      ENDIF
      IF (VY10(INDX).LT.40..OR.VY10B(INDX).LT.40.) THEN
        XY10  = 0.
        XY10B = 0.
      ENDIF
C
      RETURN
C
   80 WRITE(*,110)
  110 FORMAT(' ERROR IN SOLAR FLUX READ IN SOLFSMY - STOP ')
      STOP
C
      END SUBROUTINE SOLFSMY
C
C***********************************************************************
C
      SUBROUTINE SUNPOS (AMJD, SOLRAS, SOLDEC)
C
C     This subroutine returns the solar right ascension (SOLRAS) and
C     the solar declination (SOLDEC) as a function of the input
C     Modified Julian Date (AMJD).  The subroutine is accurate to
C     0.01 degrees between the years 1950 through 2050.
C
C     Reference:  The Astronomical Almanac, page C24.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER( HALF   = 1.D0 / 2.D0         )
      PARAMETER( PI     = 3.141592653589793D0 )
      PARAMETER( TWOPI  = 2.D0 * PI           )
      PARAMETER( DEGRAD = PI / 180.D0         )

C     Compute days since J2000.0

      D2000 = AMJD - 51544.5D0

C     Compute solar mean anomaly (SOLAN).

      SOLAN = 357.528D0 + 0.9856003 * D2000
      SOLAN = SOLAN * DEGRAD

C     Compute solar ecliptic longitude (ECLON) using
C     solar mean longitude (SOLON).

      SOLON = 280.460D0 + 0.9856474 * D2000
      SOLON = DMOD(SOLON,360.D0)
      IF (SOLON .LT. 0.D0) THEN
        SOLON = SOLON + 360.D0
      END IF
      ECLON = SOLON + 1.915D0 * DSIN(SOLAN) + 0.02D0 * DSIN(2.D0*SOLAN)
      ECLON = ECLON * DEGRAD

C     Compute obliquity of the ecliptic.

      EPS   = 23.439D0 - 0.0000004 * D2000
      EPS   = EPS * DEGRAD

C     Compute ecliptic longitude terms.

      SIN1L = DSIN(1.D0 * ECLON)
      SIN2L = DSIN(2.D0 * ECLON)
      SIN4L = DSIN(4.D0 * ECLON)

C     Compute obliquity terms.

      TANHALFEPS1 = DTAN(HALF * EPS)
      TANHALFEPS2 = TANHALFEPS1 * TANHALFEPS1
      TANHALFEPS4 = TANHALFEPS2 * TANHALFEPS2

C     Compute solar right ascension (SOLRAS) in radians and
C     compute solar declination (SOLDEC) in radians.

      SOLRAS = ECLON - TANHALFEPS2 * SIN2L + HALF * TANHALFEPS4 * SIN4L

      IF (SOLRAS .LT. 0.D0) THEN
        SOLRAS = SOLRAS + TWOPI
      ELSEIF (SOLRAS .GT. TWOPI) THEN
        SOLRAS = SOLRAS - TWOPI
      END IF

      SOLDEC = DASIN(DSIN(EPS) * SIN1L)

      RETURN
      END SUBROUTINE SUNPOS
C
C***********************************************************************
C
      REAL*8 FUNCTION THETA(T1950)
C
C     CALCULATES RIGHT ASCENSION OF GREENWICH AT T1950 (DAYS SINCE 1950)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DATA TWOPI/6.28318530717958648D0/
C
      NDAY  = T1950
      TFRAC = T1950 - NDAY
      IF (NDAY.LE.7305) THEN
C       COMPUTE THETA FROM 1950
        THETA = DMOD(1.7294446614D0 + 1.72027915246D-2*NDAY +
     *               6.3003880926D0*TFRAC, TWOPI)
      ELSE
C       COMPUTE THETA FROM 1970
C        7305.0 DAYS FROM JAN 0.0 1950 TO JAN 0.0 1970
C       25566.5 DAYS FROM JAN 0.5 1900 TO JAN 0.0 1970
C       18261.5 DAYS FROM JAN 0.5 1900 TO JAN 0.0 1950
C
        TS70  = T1950 - 7305.D0
        IDS70 = TS70
        TFRAC = TS70 - IDS70
        THETA = DMOD(1.73213438565D0 + 1.720279169407D-2*IDS70 +
     *               (1.720279169407D-2+TWOPI)*TFRAC +
     *               5.0755141943D-15*TS70**2, TWOPI)
      ENDIF
      IF (THETA.LT.0.D0) THETA = THETA + TWOPI
      RETURN
      END FUNCTION THETA
