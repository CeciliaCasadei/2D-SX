SUBROUTINE calculate(radialTolerance, azimuthTolerance, pixelTolerance, &
                      orderedPeaksMatrix, &
                      possiblePredictedPattern, nMatches, &
                      nPeaks, nPredictedSpots)
IMPLICIT NONE

!f2py threadsafe
!f2py intent(in) nPeaks
!f2py intent(in) nPredictedSpots
!f2py intent(in) radialTolerance
!f2py intent(in) azimuthTolerance
!f2py intent(in) pixelTolerance
!f2py intent(in) orderedPeaksMatrix
!f2py intent(in) possiblePredictedPattern
!f2py intent(inplace) nMatches

! Arguments 
REAL(8) :: radialTolerance, azimuthTolerance, pixelTolerance
REAL(8), DIMENSION(0:nPeaks-1, 0:5) :: orderedPeaksMatrix
REAL(8), DIMENSION(0:nPredictedSpots-1, 0:12) :: possiblePredictedPattern
INTEGER, DIMENSION(1) :: nMatches
INTEGER :: nPeaks, nPredictedSpots

! Local variables
INTEGER :: i, j
REAL(8) :: expRadius, expAzimuth, &
           predictedRadius, predictedAzimuth, &
           deltaRadius, deltaAzimuth, phiTolerance
REAL(8), PARAMETER :: M_PI = 3.141592654D0

! Initializations
nMatches(1) = 0


LOOP_i : DO i=0, nPeaks-1
    IF (orderedPeaksMatrix(i,5) == 1) THEN
        expRadius  = orderedPeaksMatrix(i,3)
        expAzimuth = orderedPeaksMatrix(i,4)
        LOOP_j : DO j=0, nPredictedSpots-1
            predictedRadius  = possiblePredictedPattern(j,10) !Radius, pxls
            predictedAzimuth = possiblePredictedPattern(j,8)  !Azimuth in [0, 2pi]
            deltaRadius = predictedRadius - expRadius
            IF ((-radialTolerance <= deltaRadius) .AND. (deltaRadius <= radialTolerance) ) THEN
                phiTolerance = MIN( ((azimuthTolerance/180.D0) * M_PI), (pixelTolerance / predictedRadius) )
                deltaAzimuth = predictedAzimuth - expAzimuth
                IF (ABS(deltaAzimuth) <= phiTolerance) THEN
                    nMatches(1) = nMatches(1)+1
                ELSE IF ((predictedAzimuth <= phiTolerance/2.D0) .AND. (expAzimuth >= (2.D0 * M_PI - phiTolerance/2.D0) )) THEN
                    nMatches(1) = nMatches(1)+1
                ELSE IF ((expAzimuth <= phiTolerance/2.D0) .AND. (predictedAzimuth >= (2*M_PI - phiTolerance/2.D0) )) THEN
                    nMatches(1) = nMatches(1)+1
                END IF
            END IF
        END DO LOOP_j  
    END IF 
END DO LOOP_i

!PRINT *, "nMatches=", nMatches


END SUBROUTINE calculate

! >>> import numpy
! >>> import calculate_nMatches_F
! >>> print calculate_nMatches_F.calculate.__doc__
! COMPILE:
! f2py -c -m calculate_nMatches_F --fcompiler=gfortran --opt='-O3' calculate_nMatches_F.f90