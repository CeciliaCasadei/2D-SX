SUBROUTINE generate(radialTolerance, azimuthTolerance, pixelTolerance, &
                    indexedPeaksTable, orderedPeaksMatrix, possiblePredictedPattern, &
                    nMatches, nPeaks, nPredictedSpots)


IMPLICIT NONE

!f2py = threadsafe

!f2py intent(in) radialTolerance
!f2py intent(in) azimuthTolerance
!f2py intent(in) pixelTolerance
!f2py intent(inplace) indexedPeaksTable
!f2py intent(inplace) orderedPeaksMatrix
!f2py intent(in) possiblePredictedPattern
!f2py intent(in) nMatches
!f2py intent(in) nPeaks
!f2py intent(in) nPredictedSpots


! orderedPeaksMatrix needs to be modified
! indexedPeaksTable needs to be modified


! Arguments
REAL(8) :: radialTolerance, azimuthTolerance, pixelTolerance
REAL(8), DIMENSION(0:nMatches-1, 0:10) :: indexedPeaksTable   ! nMatches rows and 11 columns
REAL(8), DIMENSION(0:nPeaks-1, 0:5) :: orderedPeaksMatrix
REAL(8), DIMENSION(0:nPredictedSpots-1, 0:12) :: possiblePredictedPattern
INTEGER :: nMatches
INTEGER :: nPeaks, nPredictedSpots

! Local variables
INTEGER :: indexedPeaksTableRow
INTEGER :: i, j
REAL(8) :: expRadius, expAzimuth, &
           predictedRadius, predictedAzimuth, &
           deltaRadius, deltaAzimuth, phiTolerance
REAL(8), PARAMETER :: M_PI = 3.141592654D0

! initializations
indexedPeaksTableRow = 0


LOOP_i : DO i=0, nPeaks-1
    IF (orderedPeaksMatrix(i,5) == 1) THEN
        expRadius  = orderedPeaksMatrix(i,3)
        expAzimuth = orderedPeaksMatrix(i,4)
        LOOP_j : DO j=0, nPredictedSpots-1
            predictedRadius  = possiblePredictedPattern(j,10) !Radius, pxls
            predictedAzimuth = possiblePredictedPattern(j,8)  !Azimuth in [0, 2pi]
            deltaRadius = predictedRadius - expRadius
            IF ((ABS(deltaRadius) <= radialTolerance) .AND. (indexedPeaksTableRow < nMatches) ) THEN 
                phiTolerance = MIN( ((azimuthTolerance/180.D0) * M_PI), (pixelTolerance / predictedRadius) )
                indexedPeaksTable(indexedPeaksTableRow,0) = possiblePredictedPattern(j,0)                     ! h
                indexedPeaksTable(indexedPeaksTableRow,1) = possiblePredictedPattern(j,1)                     ! k
                indexedPeaksTable(indexedPeaksTableRow,2) = expRadius                                         ! experimental radius
                indexedPeaksTable(indexedPeaksTableRow,3) = expAzimuth                                        ! experimental azimuth
                indexedPeaksTable(indexedPeaksTableRow,4) = orderedPeaksMatrix(i,2)                 ! experimental intensity
                indexedPeaksTable(indexedPeaksTableRow,5) = ABS(deltaRadius)                  ! radial difference
                indexedPeaksTable(indexedPeaksTableRow,7) = i                                                 ! experimental peak n                                               
                indexedPeaksTable(indexedPeaksTableRow,8) = j                                                 ! predicted peak n
                indexedPeaksTable(indexedPeaksTableRow,9)  = predictedRadius                                  ! predicted radius
                indexedPeaksTable(indexedPeaksTableRow,10) = predictedAzimuth                                 ! predicted azimuth
                deltaAzimuth = predictedAzimuth - expAzimuth
                IF (ABS(deltaAzimuth) <= phiTolerance) THEN
                    orderedPeaksMatrix(i,5) = 0
                    indexedPeaksTable(indexedPeaksTableRow,6) = ABS(deltaAzimuth)            ! azimuth difference
                    indexedPeaksTableRow = indexedPeaksTableRow + 1
                ELSE IF ((predictedAzimuth <= phiTolerance/2.D0) .AND. (expAzimuth >= (2.D0 * M_PI - phiTolerance/2.D0) )) THEN
                    orderedPeaksMatrix(i,5) = 0
                    indexedPeaksTable(indexedPeaksTableRow,6) = predictedAzimuth + 2 * M_PI - expAzimuth
                    indexedPeaksTableRow = indexedPeaksTableRow + 1
                ELSE IF ((expAzimuth <= phiTolerance/2.D0) .AND. (predictedAzimuth >= (2 * M_PI - phiTolerance/2.D0) )) THEN
                    orderedPeaksMatrix(i,5) = 0
                    indexedPeaksTable(indexedPeaksTableRow,6) = expAzimuth + 2 * M_PI - predictedAzimuth
                    indexedPeaksTableRow = indexedPeaksTableRow + 1
                END IF
            END IF 
        END DO LOOP_j  
    END IF 
END DO LOOP_i
END SUBROUTINE generate

! f2py -c -m generate_indexedPeaksTable_F --fcompiler=gfortran --opt='-O3' generate_indexedPeaksTable_F.f90