#include <stdio.h>
#include <math.h>
#include "idl_export.h"


/* All arguments from IDL are passed by reference unless otherwise
   specified with the VALUE keyword to CALL_EXTERNAL */

//double CompareMFL4Alpha( IDL_LONG SemiNum, float *SemiCurveNum, float *SemiCurveX, float *SemiCurveY,
//     	                 IDL_LONG MFLNum, float *myMFLNum, float *myMFLX, float *myMFLY)

double CompareMFL4Alpha( IDL_LONG SemiNum, float *SemiCurveNum, float *SemiCurveX, float *SemiCurveY,
     	                 IDL_LONG MFLNum, float *myMFLNum, float *myMFLX, float *myMFLY,
						 IDL_STRING *CurrAlpha, IDL_STRING *myTraceIdx, IDL_STRING *myCurveIdx)

{
	unsigned int i, j, k, l;
	unsigned int semiStrtIdx, semiEndIdx;
	unsigned int MFLStrtIdx, MFLEndIdx;
	unsigned int MatchMFLStrtIdx, MatchMFLEndIdx;

    int ii, jj, kk;

	unsigned int boxStep;

	unsigned int boxCenCnt;
	float boxCenX[1000];
	float boxCenY[1000];
	
	unsigned int *closePtCnt;
	unsigned int maxClosePtCnt;
	unsigned int maxClosePtCntIdx;
	
	double distErr, myDistErr, minDistErr; 
	double avgDistErr=0.0; 
	unsigned int distErrCnt=0;

	double rmsErr, rmsDistErr=0.0;

	
	int SemiMaskArr[512][512];
	int MFLMaskArr[512][512];
	float MaskX, MaskY;
	unsigned int MaskCnt;


	int cntThreshold = 35;		// 1) 35,   2) 20,    3) 45

	int myBoxSz = 5;


	FILE *myFile, *myFile2;
	char fileName[200], fileName2[200];

	sprintf(fileName,"MFLMatchingIdx_%s_%s_%s.txt", IDL_STRING_STR(myTraceIdx), IDL_STRING_STR(CurrAlpha), IDL_STRING_STR(myCurveIdx));
	myFile = fopen(fileName, "w");


	sprintf(fileName2,"AvgRMSError_%s_%s.txt", IDL_STRING_STR(myTraceIdx), IDL_STRING_STR(myCurveIdx));
	myFile2 = fopen(fileName2, "a+");


	if ( ( closePtCnt = (unsigned int *)malloc(MFLNum*sizeof(unsigned int)) ) == NULL ){
		printf("\nError, memory not allocated.\n");
		exit(1);
	}


	semiStrtIdx = 0;
	semiEndIdx  = 0;


	for(i=0; i<(int)SemiNum; i++){	// for all SemiCurve
			
		semiEndIdx = semiStrtIdx + (int)SemiCurveNum[i] - 1;


		/*		
		// version 1
		if(SemiCurveNum[i] < 100.0){
		//if(SemiCurveNum[i] < 1000.0){
			//boxStep = 5;
			boxStep = (int)(SemiCurveNum[i]/5);
		}
		else{
			//boxStep = (int)(SemiCurveNum[i]/(int)(SemiCurveNum[i] / 15.0));
			//boxStep = (int)(SemiCurveNum[i]/(int)(SemiCurveNum[i] / 150.0));
			//boxStep = (int)(SemiCurveNum[i]/(int)(SemiCurveNum[i] / 100.0));
			boxStep = (int)(SemiCurveNum[i]/(int)(SemiCurveNum[i] / 50.0));
		}		
		*/

		// version 2 : used for real datasets
		//boxStep = (int)(SemiCurveNum[i] / ((int)(SemiCurveNum[i]/100.0)) );

		// version 3 : used for synthetic datasets (for short semi-curve)
		if(SemiCurveNum[i] <= 50.0){		
			boxStep = (int)(SemiCurveNum[i] / ((int)(SemiCurveNum[i]/5.0)) );
		}
		else if( (SemiCurveNum[i] > 50.0) && (SemiCurveNum[i] <= 100.0) ){		
			boxStep = (int)(SemiCurveNum[i] / ((int)(SemiCurveNum[i]/20.0)) );
		}
		else
		{
			boxStep = (int)(SemiCurveNum[i] / ((int)(SemiCurveNum[i]/50.0)) );
		}
		
		boxCenCnt = 0;
		// 1) get each SemiCurve's box Center Coordinates to compare		
		for(j=semiStrtIdx; j<=semiEndIdx; j=j+boxStep){

			boxCenX[boxCenCnt] = SemiCurveX[j];
			boxCenY[boxCenCnt] = SemiCurveY[j];			

			boxCenCnt++;

		}


tryLargeBox:

		// initial mask
		for(ii=0; ii<512; ii++){
			for(jj=0; jj<512; jj++){
				SemiMaskArr[ii][jj] = 0;				
			}
		}

		// Semi Mask
		for(j=0; j<boxCenCnt; j++){

			for(ii=-myBoxSz; ii<=myBoxSz; ii++){
				for(jj=-myBoxSz; jj<=myBoxSz; jj++){				

					MaskX = boxCenX[j]+(float)ii;
					MaskY = boxCenY[j]+(float)jj;

					if( (MaskX>=0.000000001) && (MaskX<511.000000001) && (MaskY>=0.000000001) && (MaskY<511.000000001) ){
						SemiMaskArr[(int)MaskX][(int)MaskY] = 1;
					}
				}
			}
		}

		MFLStrtIdx = 0;
		MFLEndIdx  = 0;
		// MFL Mask
		for(j=0; j<MFLNum; j++){

			MFLEndIdx = MFLStrtIdx + myMFLNum[j] - 1;				

			for(ii=0; ii<512; ii++){
				for(jj=0; jj<512; jj++){					
					MFLMaskArr[ii][jj]  = 0;
				}
			}
			
			// if there's matching loop (i.e., they are close to each other)
			for(k=MFLStrtIdx; k<=MFLEndIdx; k++){
				MFLMaskArr[ (int)myMFLX[k] ][ (int)myMFLY[k] ] = 1;
			}


			MaskCnt = 0;
			for(ii=0; ii<512; ii++){
				for(jj=0; jj<512; jj++){
					if( (SemiMaskArr[ii][jj]+MFLMaskArr[ii][jj]) == 2){
						MaskCnt++;
					}
				}
			}
			closePtCnt[j] = MaskCnt;

			MFLStrtIdx = MFLEndIdx + 1;
		}

		// counting matching mask points
		maxClosePtCnt    = 0;
		maxClosePtCntIdx = -1;
		for(j=0; j<MFLNum; j++){
			if( closePtCnt[j] > maxClosePtCnt){
				maxClosePtCnt = closePtCnt[j];
				maxClosePtCntIdx = j;
			}
		}


if(maxClosePtCnt == 0){
	myBoxSz += 5;
	goto tryLargeBox;
}



		fprintf(myFile, "%d %d\n", maxClosePtCntIdx, maxClosePtCnt);

		
		if(maxClosePtCnt > cntThreshold){

			distErrCnt++;

			distErr = 0.0;
			rmsErr  = 0.0;


			MatchMFLStrtIdx = 0;
			MatchMFLEndIdx  = 0;
			for(l=0; l<=maxClosePtCntIdx; l++){
				MatchMFLEndIdx = MatchMFLStrtIdx + myMFLNum[l] - 1;

				if(l!=maxClosePtCntIdx)
					MatchMFLStrtIdx = MatchMFLEndIdx + 1;
			}
			
			
			for(j=semiStrtIdx; j<=semiEndIdx; j++){		// compute distance error
				minDistErr = 99999999.0;

				for(k=MatchMFLStrtIdx; k<=MatchMFLEndIdx; k++){

					myDistErr = sqrt( (SemiCurveX[j]-myMFLX[k])*(SemiCurveX[j]-myMFLX[k])+
					                  (SemiCurveY[j]-myMFLY[k])*(SemiCurveY[j]-myMFLY[k]) );
					
					if( myDistErr <= minDistErr){
						minDistErr = myDistErr;
					}
				}

				distErr += minDistErr;

				rmsErr += (minDistErr*minDistErr);
			}

			avgDistErr += (distErr/(float)(semiEndIdx-semiStrtIdx+1));

			rmsDistErr += sqrt(rmsErr/(float)(semiEndIdx-semiStrtIdx+1));
		}

		semiStrtIdx = semiEndIdx + 1;

	}	// end of for(i=0; i<(int)SemiNum; i++){	// for all SemiCurve


	if(distErrCnt > 0)
	{
		avgDistErr = avgDistErr/(double)distErrCnt;
	
		rmsDistErr = rmsDistErr/(double)distErrCnt;
	}
	else
	{
		avgDistErr = 999.0;
	
		rmsDistErr = 999.0;
	}	


	fprintf(myFile2, "%f\n", rmsDistErr);
	fclose(myFile2);

	fprintf(myFile, "%d %d\n", -999, -999);
	fclose(myFile);


	return(avgDistErr);
}
