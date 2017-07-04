#include <math.h>
#include <stdio.h>

int calculate(int nPeaks, int nPredictedSpots,
              double radialTolerance, double azimuthTolerance, double pixelTolerance,
              double orderedPeaksMatrix[][6], 
              double possiblePredictedPattern[][13]){

int nMatches = 0;
int i;
int j;
double expRadius;
double expAzimuth;
double predictedRadius;
double predictedAzimuth;
double deltaRadius;
double deltaAzimuth;
double phiTolerance;

/*
printf ("%lf", orderedPeaksMatrix[2][3]);
printf ("%lf", orderedPeaksMatrix[nPeaks-1, 5]);
printf ("%lf", possiblePredictedPattern[7][6]);
printf ("%lf", possiblePredictedPattern[nPredictedSpots-1][12]);
*/

    for (i=0; i<nPeaks; i++){
        if (orderedPeaksMatrix[i][5] == 1){
            expRadius  = orderedPeaksMatrix[i][3];                      /* Radius, pxls */
            expAzimuth = orderedPeaksMatrix[i][4];                     /* Azimuth in [0, 2pi] */
            for (j=0; j<nPredictedSpots; j++){
                predictedRadius  = possiblePredictedPattern[j][10];     /* Radius, pxls */
                predictedAzimuth = possiblePredictedPattern[j][8];      /* Azimuth in [0, 2pi] */
                deltaRadius = predictedRadius-expRadius;
                if ( (-radialTolerance <= deltaRadius) &&  (deltaRadius <= radialTolerance)  ) { 
                    phiTolerance = (azimuthTolerance/180)*M_PI;
                    if (pixelTolerance/predictedRadius < phiTolerance) {phiTolerance = pixelTolerance/predictedRadius;}
                    deltaAzimuth = predictedAzimuth-expAzimuth;
                    if ( (-phiTolerance <= deltaAzimuth ) && (deltaAzimuth <= phiTolerance) ) {
                        nMatches = nMatches+1;
                        }
                    else if ( (predictedAzimuth <= phiTolerance/2) && (expAzimuth >= (2*M_PI - phiTolerance/2) ) ) {
                        nMatches = nMatches+1;
                        }
                    else if ( (expAzimuth <= phiTolerance/2) && (predictedAzimuth >= (2*M_PI - phiTolerance/2) ) ) {
                        nMatches = nMatches+1;
                        }
                }  else {continue;} 


} /* END FOR j */

}

} /* END FOR i */

/* printf ("%d", nMatches); */
return nMatches;

}

/*
gcc -O3 -std=gnu99 -c calculate_nMatches_C.c
gcc -shared -o calculate_nMatches_C.so calculate_nMatches_C.o
*/