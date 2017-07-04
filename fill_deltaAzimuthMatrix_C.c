#include <math.h>

void filling(int inPlaneAngle, int nPredictedSpots,
             double intenseExpPeakRadius, double intenseExpPeakAzimuth,
             double radialTolerance, double azimuthTolerance, double pixelTolerance,
             double predictionInOneOrientation[][13], double deltaAzimuthMatrix[][256]){

    int predictedSpot;
    double predictedRadius;
    double predictedAzimuth;
    double deltaRadius;
    double phiTolerance;
    double deltaAzimuth;

    for (predictedSpot=0; predictedSpot<nPredictedSpots; predictedSpot++){         
        predictedRadius  = predictionInOneOrientation[predictedSpot][10];                /* Radius, pxls */
        predictedAzimuth = predictionInOneOrientation[predictedSpot][8];                 /* Between 0 and 2pi */
        deltaRadius = predictedRadius-intenseExpPeakRadius;                             
        if ( (-radialTolerance <= deltaRadius) &&  (deltaRadius <= radialTolerance)  ) {   
            phiTolerance = azimuthTolerance/(180*M_PI);
            if (pixelTolerance/predictedRadius < phiTolerance) {phiTolerance = pixelTolerance/predictedRadius;}
            deltaAzimuth = predictedAzimuth-intenseExpPeakAzimuth;
            if ( (-phiTolerance <= deltaAzimuth ) && (deltaAzimuth <= phiTolerance) ) {
            deltaAzimuthMatrix[predictedSpot][inPlaneAngle-1] = fabs(deltaAzimuth);
            } 
            else if ( (predictedAzimuth <= phiTolerance/2) && (intenseExpPeakAzimuth >= (2*M_PI - phiTolerance/2) ) ) {
            deltaAzimuthMatrix[predictedSpot][inPlaneAngle-1] = predictedAzimuth + 2*M_PI - intenseExpPeakAzimuth;
            }
            else if ( (intenseExpPeakAzimuth <= phiTolerance/2) && (predictedAzimuth >= (2*M_PI - phiTolerance/2) ) ) {
            deltaAzimuthMatrix[predictedSpot][inPlaneAngle-1] = intenseExpPeakAzimuth + 2*M_PI - predictedAzimuth;
            }
            else {continue;}

} /* END IF */

} /* END FOR */

}

/*
gcc -O3 -std=gnu99 -c fill_deltaAzimuthMatrix_C.c
gcc -shared -o fill_deltaAzimuthMatrix_C.so fill_deltaAzimuthMatrix_C.o
*/