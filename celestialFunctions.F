      MODULE celestialFunctions

      CONTAINS

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
C
C***********************************************************************
      END MODULE
