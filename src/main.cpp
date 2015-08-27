#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <cstdint>
#include "LDPC_Class.hpp"

#include <sstream>
#include <string> // this should be already included in <sstream>
#include <vector> // this should be already included in <sstream>
#include <map> // this should be already included in <sstream>

using namespace std;

#define MaxErrorFrames 200
#define Epsilon 0.00001
//#define Penalty 3.0
//#define InversePenalty (1.0/Penalty)
#define Maxiter 500
#define MaxCheckdegree 7
#define Population 30  //5*(3+1)  T=5
#define Generation 500
#define MutationNumber 6
#define DiffVar 0.5


//#define HFile "Margulis_2640_1320.txt"
#define HFile "matrix_2304_1152.txt"
//#define HFile "PEGirReg504x1008.txt"
//#define HFile "matrix_576_288.txt"



//PCHK
int **PCHK; //У������
int M; //��
int N; //��
int MaxRowDegree;
int MaxColDegree;
int *pRowDegree;
int *pColDegree;
int **A; //��Ϣ�ڵ�-
int **B; //У���ڵ�-

//used by ADMM LP decoding
int NumIter;
double **pReplicaZ;  //auxiliary replica variables for each check node j
double **pLagrangeMulti;  //Lagrangian multipliers for each check node j
double **pChecktoVariable;
double *pDecVec;
double *pObjCoef;
double *pRecvVec;

//BER/FER
__int64_t MaxTransBits = 100000000000000000L; //������������
__int64_t TotalBits; //�������ܱ�����
__int64_t ErrorBits; //�����ı�������
__int64_t TotalFrames; //�������ܵ�Block
__int64_t ErrorFrames; //������Block
double BER; //��������
double FER; //������

void ReadPCHK();
void ProjectOnPPd(double *pReplicaBeforeProj, int ReplicaLength, double* pReplicaZ);
void ProjectOnPPd_CSA(double *pReplicaBeforeProj, int ReplicaLength, double* pReplicaZ);
int SearchCut(double *pProjOneZero, int *pIndicatorVec, int ReplicaLength);
double CalculateOptimalPara(double *pReplicaBeforeProj,int *pIndicatorVec,int ReplicaLength);
int cmp(const void * a, const void * b);
double ADMMDecoder(double stdev,double *pBeta,int *pDegreeSize,int betaSize,double penalty);

int main(int argc, char* argv[])
{
	double gauss;
	double stdev;

	double SNR = 2.0;

	double dbtemp;
	double coderate;
	//__int64 i;
	int i,j,k;
	int count;
    double dbFER;

//	srand((unsigned)time(NULL));
	srand( 0 );

//	stdev = sqrt(1.0/(2*coderate*pow(10,SNR/10)));  //sqrt(-2.0*SIGMA*SIGMA*log(s)/s);

	int      degreeSize;
	float  **pBeta;
	float  **pPrevBeta;
	int     *pDegreeSize;
	double  *pPenalty;
	double  *pPrevPenalty;
	double  *pFER;
	double  *pPrevFER;
	int      BestPopulationIndex;
	double  *pBestBeta;
	double   bestBER;
	double   bestFER;
	int     *pMutationIndex; //MutationNumber

    int N                        = 2640;
    int NmK                      = 1320;
    int ADMMtype                 = 2;
    unsigned int fe_limit        = 25;
    unsigned int codeword_limit  = 100000000;
    unsigned int time_limit      = 0;
    unsigned int num_threads     = 8;
    unsigned int maxIts          = 200;



    for (int p = 1; p < argc; p++)
    {
        //
        //  SIMULATION PARAMETERS
        //
        if (strcmp(argv[p], "-N") == 0) {
            N = atoi(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-NmK") == 0) {
            NmK = atoi(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-fe_limit") == 0) {
            fe_limit = atoi(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-cw_limit") == 0) {
        	codeword_limit = atoi(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-iters") == 0) {
        	maxIts = atoi(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-snr") == 0) {
        	SNR = atof(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-threads") == 0) {
            num_threads = atoi(argv[p + 1]);
            p += 1;

        }else{
                exit( 0 );
        }
    }

    char ldpcFile[1024], codwFile[1024];
    sprintf(ldpcFile, "./mat/H_%dx%d.txt", N, NmK);
    sprintf(codwFile, "./cwd/codeword_%dx%d.txt", N, NmK);
    double rendement = (float) (N-NmK) / (float) (N);

    Simulator simu(N, NmK, ldpcFile, "AWGN", SNR, num_threads);

    map<int, int> histo;
    for(int i=0; i<N; i++)
    {
    	int deg = simu.mDecoder[0]->get_deg_vn( i );
        histo[ deg ] = deg;
    }

    //
    // ON RECHERCHE LES VALEURS MAXIMALES CONTENUE DANS LA MAP
    //
    std::vector<int> rHisto;   // populate this
    for (std::map<int, int>::iterator it=histo.begin(); it!=histo.end(); ++it){
        rHisto.push_back( it->second );
    }
    std::sort   (rHisto.begin(), rHisto.end());
    //std::reverse(rHisto.begin(),rHisto.end());

    printf("(II) VersioN: (%s, %s)\n", __DATE__, __TIME__);
    printf("(II) Code LDPC (N, K, N-K): (%d, %d, %d)\n", N, N-NmK, NmK);
    printf("(II) Rendement du code    : %.3f\n", rendement);
    printf("(II) # ITERATIONs du CODE : %d\n", maxIts);
    printf("(II) FER LIMIT FOR SIMU   : %d\n", fe_limit);
    printf("(II) SIMULATION  RANGE    : %.2f dB\n", SNR);
    printf("(II) LOADED H FILENAME    : %s\n", ldpcFile);
    printf("(II) codeword file        : %s\n", codwFile);

	degreeSize     = rHisto.size();
	pDegreeSize    = new int    [degreeSize];
	pBestBeta      = new double [degreeSize];
	pPenalty       = new double [Population];
	pPrevPenalty   = new double [Population];
	pFER           = new double [Population];
	pPrevFER       = new double [Population];
	pMutationIndex = new int [MutationNumber];

	pBeta     = new float*  [Population];
	pPrevBeta = new float*  [Population];
	for(i=0;i<Population;i++)
	{
		pBeta[i]     = new float  [degreeSize];
		pPrevBeta[i] = new float  [degreeSize];
	}
	//Population Generation
    for(int i=0; i<rHisto.size(); i++)
    {
    	printf("> deg[%d] = %d\n", i, rHisto.at(i));
    	pDegreeSize[i] = rHisto.at(i);
    }
/*
	pDegreeSize[0] = 2;
	pDegreeSize[1] = 3;
	pDegreeSize[2] = 6;
*/
	for(i=0;i<Population;i++)
	{
		for(j=0;j<degreeSize;j++)
		{
			pBeta[i][j] = 20*(rand()/(float)RAND_MAX);
		}
	}

	for(i=0;i<Population;i++)
	{
		pPenalty[i] = 15 * (rand()/(float)RAND_MAX);
	}

	bestFER = 2.0;


	for(i=0;i<Population;i++)
	{
		/*
		pFER[i] = ADMMDecoder(stdev, pBeta[i], pDegreeSize, degreeSize, pPenalty[i]);
		 */

	    Simulator ldpcsim(N, NmK, ldpcFile, "AWGN", SNR, num_threads);
		ldpcsim.SetCodeword(codwFile);
		ldpcsim.SetCommandLineOutput(1000000, -1);// create a command line output every 10000 simulations.
		int seed = rand();
		ldpcsim.SetTargets(fe_limit, codeword_limit, time_limit);// set target number of frame errors = 100
		ldpcsim.SetDecoder("ADMML2");
		ldpcsim.SetDecoderParameters(maxIts, 1e-5, pPenalty[i], pBeta[i], pDegreeSize, degreeSize);
		ldpcsim.SetChannelSeed(seed);
		ldpcsim.RunSim();
		pFER[i] = ldpcsim.GetLatestFerValue();

		printf("[%2d%%] - FER = %1.3e - PARAM [", (100*i) / Population, pFER[i]);
		for(j=0;j<degreeSize;j++) { printf("deg[%d]=%1.3f, ", pDegreeSize[j], pBeta[i][j]); }
		printf("mu=%1.3f]          \r", pPenalty[i]); fflush(stdout);

		if(pFER[i]<bestFER)
		{
			bestFER             = pFER[i];
			BestPopulationIndex = i;
		}
	}

//	printf("BestPopulationIndex: %d\n",BestPopulationIndex);

	for(int iGeneration=1;iGeneration<Generation;iGeneration++)
	{
		bestFER = 2.0;
		for(i=0;i<Population;i++)
		{
			if(pFER[i]<bestFER)
			{
				bestFER = pFER[i];
				BestPopulationIndex = i;
			}
		}

        printf("=> Gen: %d, bestFER: %1.2e, opt_paras: ",iGeneration, bestFER);
		for(j=0;j<degreeSize;j++) { printf("deg[%d]=%1.8f, ", pDegreeSize[j], pBeta[BestPopulationIndex][j]); }
		printf("mu=%1.8f]\n", pPenalty[BestPopulationIndex]);
//		for(i=0;i<degreeSize;i++)
//		{
//			printf("%1.8f ",pBeta[BestPopulationIndex][i]);
//		}
//		printf("mu=%1.8f\n",pPenalty[BestPopulationIndex]);


		for(i=0;i<Population;i++)
		{
			for(j=0;j<degreeSize;j++)
			{
				 pPrevBeta[i][j] = pBeta[i][j];
			}
			pPrevFER[i]     = pFER[i];
			pPrevPenalty[i] = pPenalty[i];
		}

		//2. mutation and test
        for(i=0;i<Population;i++)
		{
			//generate J=MutationNumber distinct random numbers   pMutationIndex
			for(j=0;j<MutationNumber;j++)
			{
				while(1)
				{
					int temp;
					temp = rand()%Population;
					if(temp!=i)
					{
						pMutationIndex[j] = temp;
						break;
					}
				}
			}

			double oddDiff,evenDiff;

			//DiffVar

			for(j=0;j<degreeSize;j++)
			{
				oddDiff = 0;
				evenDiff = 0;
				for(int k=0;k<MutationNumber;k=k+2)
				{
					oddDiff  += pPrevBeta[pMutationIndex[k]][j];
					evenDiff += pPrevBeta[pMutationIndex[k+1]][j];
				}
				pBeta[i][j] = pPrevBeta[BestPopulationIndex][j] + DiffVar*(oddDiff-evenDiff);
				if(pBeta[i][j]<0)
					pBeta[i][j] = 0.00001;
			}

			oddDiff = 0;
			evenDiff = 0;
			for(int k=0;k<MutationNumber;k=k+2)
			{
				oddDiff  += pPrevPenalty[pMutationIndex[k]];
				evenDiff += pPrevPenalty[pMutationIndex[k+1]];
			}

			pPenalty[i] = pPrevPenalty[BestPopulationIndex] + DiffVar*(oddDiff-evenDiff);
			if(pPenalty[i]<0)
				pPenalty[i] = 0.00001;
		}

		for(i=0;i<Population;i++)
		{
			/*
			pFER[i] = ADMMDecoder(stdev,pBeta[i],pDegreeSize,degreeSize,pPenalty[i]);
			*/
			Simulator ldpcsim(N, NmK, ldpcFile, "AWGN", SNR, num_threads);
			ldpcsim.SetCodeword(codwFile);
			ldpcsim.SetCommandLineOutput(10000, -1);// create a command line output every 10000 simulations.
			int seed = rand();
			ldpcsim.SetTargets(fe_limit, codeword_limit, time_limit);// set target number of frame errors = 100
			ldpcsim.SetDecoder("ADMML2");
			ldpcsim.SetDecoderParameters(maxIts, 1e-5, pPenalty[i], pBeta[i], pDegreeSize, degreeSize);
			ldpcsim.SetChannelSeed(seed);
			ldpcsim.RunSim();
			pFER[i] = ldpcsim.GetLatestFerValue();

			printf("[%2d%%] - FER = %1.2e - PARAM [", (100*i) / Population, pFER[i]);
			for(j=0;j<degreeSize;j++) { printf("deg[%d]=%1.3f, ", pDegreeSize[j], pBeta[i][j]); }
			printf("mu=%1.3f]          \r", pPenalty[i]); fflush(stdout);
		}

		//3. compare and update parameters
		for(i=0;i<Population;i++)
		{
			if(pFER[i]>pPrevFER[i])
			{
				for(j=0;j<degreeSize;j++)
				{
					pBeta[i][j] = pPrevBeta[i][j];
				}
				pPenalty[i] = pPrevPenalty[i];
				pFER[i] = pPrevFER[i];
			}
		}
	}


	delete [] pDegreeSize;
	delete [] pPenalty;
	delete [] pPrevPenalty;
	delete [] pFER;
	delete [] pPrevFER;
	delete [] pBestBeta;

	for(i=0;i<Population;i++)
		delete [] pBeta[i];
	delete [] pBeta;

	for(i=0;i<Population;i++)
		delete [] pPrevBeta[i];
	delete [] pPrevBeta;
	//realease memory
	delete [] pObjCoef;

	for(i=0;i<M;i++)
		delete [] pReplicaZ[i];
	delete [] pReplicaZ;

	for(i=0;i<M;i++)
		delete [] pLagrangeMulti[i];
	delete [] pLagrangeMulti;

	delete [] pRecvVec;
	delete [] pDecVec;

	delete [] pBeta;
	delete [] pDegreeSize;


	return 0;
}

//read parity-check matrix
void ReadPCHK()
{
	int i,j;
	int temp;
	FILE* fp;
	fp = fopen(HFile,"r");

	fscanf(fp,"%d",&N);
	fscanf(fp,"%d",&M);
	fscanf(fp,"%d",&MaxColDegree);
	fscanf(fp,"%d",&MaxRowDegree);

	//PHCK malloc memory
	PCHK = new int *[M];
	for(i=0;i<M;i++)
		PCHK[i] = new int [N];
	for(i=0;i<M;i++)
	{
		for(j=0;j<N;j++)
		{
			PCHK[i][j] = 0;
		}
	}

	//degrees of columns and rows malloc memory
	pRowDegree = new int [M];
	pColDegree = new int [N];

	//A,B malloc memory
	A = new int *[N];
	for(i=0;i<N;i++)
		A[i] = new int[MaxColDegree];
	B = new int *[M];
	for(i=0;i<M;i++)
		B[i] = new int[MaxRowDegree];

	for(i=0;i<N;i++)
		fscanf(fp,"%d",&pColDegree[i]);
	for(i=0;i<M;i++)
		fscanf(fp,"%d",&pRowDegree[i]);

	//read index by column
	for(i=0;i<N;i++)
	{
		for(j=0;j<MaxColDegree;j++) //pColDegree[i]
		{
			fscanf(fp,"%d",&temp);
			if(temp!=0)
			{
				A[i][j] = temp-1;
			}
			else
			{
				A[i][j] = 0;
			}
		}
	}

	//read index by row
	for(i=0;i<M;i++)
	{
		for(j=0;j<MaxRowDegree;j++)  //pRowDegree[i]
		{
			fscanf(fp,"%d",&temp);
			if(temp!=0)
			{
				PCHK[i][temp-1] = 1;
				B[i][j] = temp-1;
			}
			else
			{
				B[i][j] = 0;
			}
		}
	}

	fclose(fp);
}

// project pReplicaBeforeProj on PPd
void ProjectOnPPd(double *pReplicaBeforeProj, int ReplicaLength, double* pReplicaZ)
{
	int i,j,k;
	int sign;
	int ConstituentParity;
	int *pPermutationOrder;
	double dbtemp;
	double dbMaxBeta;
	double *pReplicaWithOrder;
	double *pProjectToZeroOne;
	double *pSetBeta;

	pPermutationOrder = new int[ReplicaLength];
	pReplicaWithOrder = new double[ReplicaLength];
	pProjectToZeroOne = new double[ReplicaLength];

	for(i=0;i<ReplicaLength;i++)
	{
		pPermutationOrder[i] = -1;
		pReplicaWithOrder[i] = 0;
		pProjectToZeroOne[i] = 0;
	}

	for(i=0;i<ReplicaLength;i++)
	{
		int index;
		double dbmax=-DBL_MAX;

		for(j=0;j<ReplicaLength;j++)
		{
			sign = 0;
			for(k=0;k<i;k++)
			{
				if(j==pPermutationOrder[k])
				{
					sign = 1;
					break;
				}
			}

			if(!sign)
			{
				if(pReplicaBeforeProj[j]>dbmax)
				{
					dbmax = pReplicaBeforeProj[j];
					index = j;
				}
			}
		}

		pReplicaWithOrder[i] = dbmax;
		pPermutationOrder[i] = index;
	}

	for(i=0;i<ReplicaLength;i++)
	{
		//pProjectToZeroOne
		if(pReplicaWithOrder[i]>1.0)
			pProjectToZeroOne[i] = 1.0;
		else if(pReplicaWithOrder[i]<0.0)
			pProjectToZeroOne[i] = 0.0;
		else
			pProjectToZeroOne[i] = pReplicaWithOrder[i];
	}

	//ConstituentParity = ?;
	dbtemp = 0.0;
	for(i=0;i<ReplicaLength;i++)
	{
		dbtemp = dbtemp + pProjectToZeroOne[i];
	}
	ConstituentParity = (int)dbtemp;
	if(ConstituentParity%2 != 0)
		ConstituentParity = ConstituentParity - 1;

	dbMaxBeta = 0.5*(pProjectToZeroOne[ConstituentParity] - pProjectToZeroOne[ConstituentParity+1]);

	//
	dbtemp = 0.0;
	for(i=0;i<ReplicaLength;i++)
	{
		if(i<=ConstituentParity)
			dbtemp = dbtemp + pProjectToZeroOne[i];
		else
			dbtemp = dbtemp - pProjectToZeroOne[i];
	}
	if(dbtemp <= ConstituentParity)
	{
		for(i=0;i<ReplicaLength;i++)
		{
			pReplicaZ[pPermutationOrder[i]] = pProjectToZeroOne[i]; //pPermutationOrder[i]
		}
	}
	else
	{
		int NumofBeta = 0;
		int a,b;
		for(i=0;i<ReplicaLength;i++)
		{
			if(i<=ConstituentParity)
			{
				if( (pReplicaWithOrder[i]-1)>0 && (pReplicaWithOrder[i]-1)<dbMaxBeta )
					NumofBeta++;
				if( (pReplicaWithOrder[i])>0 && (pReplicaWithOrder[i])<dbMaxBeta )
					NumofBeta++;
			}
			else
			{
				if( (-1*pReplicaWithOrder[i])>0 && (-1*pReplicaWithOrder[i])<dbMaxBeta )
					NumofBeta++;
				if( (-1*pReplicaWithOrder[i]+1)>0 && (-1*pReplicaWithOrder[i]+1)<dbMaxBeta )
					NumofBeta++;
			}
		}

		int totalBeta = NumofBeta+2;
		pSetBeta = new double [totalBeta];
		pSetBeta[0] = 0;
		pSetBeta[NumofBeta+1] = dbMaxBeta;
		NumofBeta = 1;
		for(i=0;i<ReplicaLength;i++)
		{
			if(i<=ConstituentParity)
			{
				if( (pReplicaWithOrder[i]-1)>0 && (pReplicaWithOrder[i]-1)<dbMaxBeta )
				{
					pSetBeta[NumofBeta++] = pReplicaWithOrder[i]-1;
				}
				if( (pReplicaWithOrder[i])>0 && (pReplicaWithOrder[i])<dbMaxBeta )
				{
					pSetBeta[NumofBeta++] = pReplicaWithOrder[i];
				}
			}
			else
			{
				if( (-1*pReplicaWithOrder[i])>0 && (-1*pReplicaWithOrder[i])<dbMaxBeta )
				{
					pSetBeta[NumofBeta++] = -1*pReplicaWithOrder[i];
				}
				if( (-1*pReplicaWithOrder[i]+1)>0 && (-1*pReplicaWithOrder[i]+1)<dbMaxBeta )
				{
					pSetBeta[NumofBeta++] = -1*pReplicaWithOrder[i]+1;
				}
			}
		}

		// sort pSetBeta quick sort
		qsort(pSetBeta,totalBeta,sizeof(pSetBeta[0]),cmp);

		int left = 0;
		int right = totalBeta-1;
		int mid;

		while(1)
		{
			mid = (left+right)/2;

			for(i=0;i<ReplicaLength;i++)
			{
				if(i<=ConstituentParity)
				{
					pProjectToZeroOne[i] = pReplicaWithOrder[i] - pSetBeta[mid];
				}
				else
				{
					pProjectToZeroOne[i] = pReplicaWithOrder[i] + pSetBeta[mid];
				}

				if(pProjectToZeroOne[i]>1)
				{
					pProjectToZeroOne[i] = 1;
				}
				else if (pProjectToZeroOne[i]<0)
				{
					pProjectToZeroOne[i] = 0;
				}
			}

			//dbtemp
			dbtemp = 0.0;
			for(i=0;i<ReplicaLength;i++)
			{
				if(i<=ConstituentParity)
				{
					dbtemp = dbtemp + pProjectToZeroOne[i];
				}
				else
				{
					dbtemp = dbtemp - pProjectToZeroOne[i];
				}
			}

			if(dbtemp>ConstituentParity)
			{
				if(right-left<=2)
				{
					left = mid;
					break;
				}
				left = mid;
			}
			else
			{
				if(right-left<=2)
				{
					right = mid;
					break;
				}
				right = mid;
			}
		}

		//
		double f_left, f_right,optima_beta;

		for(i=0;i<ReplicaLength;i++)
		{
			if(i<=ConstituentParity)
			{
				pProjectToZeroOne[i] = pReplicaWithOrder[i] - pSetBeta[left];
			}
			else
			{
				pProjectToZeroOne[i] = pReplicaWithOrder[i] + pSetBeta[left];
			}

			if(pProjectToZeroOne[i]>1)
			{
				pProjectToZeroOne[i] = 1;
			}
			else if (pProjectToZeroOne[i]<0)
			{
				pProjectToZeroOne[i] = 0;
			}
		}

		f_left = 0.0;
		for(i=0;i<ReplicaLength;i++)
		{
			if(i<=ConstituentParity)
			{
				f_left = f_left + pProjectToZeroOne[i];
			}
			else
			{
				f_left = f_left - pProjectToZeroOne[i];
			}
		}


		for(i=0;i<ReplicaLength;i++)
		{
			if(i<=ConstituentParity)
			{
				pProjectToZeroOne[i] = pReplicaWithOrder[i] - pSetBeta[right];
			}
			else
			{
				pProjectToZeroOne[i] = pReplicaWithOrder[i] + pSetBeta[right];
			}

			if(pProjectToZeroOne[i]>1)
			{
				pProjectToZeroOne[i] = 1;
			}
			else if (pProjectToZeroOne[i]<0)
			{
				pProjectToZeroOne[i] = 0;
			}
		}

		f_right = 0.0;
		for(i=0;i<ReplicaLength;i++)
		{
			if(i<=ConstituentParity)
			{
				f_right = f_right + pProjectToZeroOne[i];
			}
			else
			{
				f_right = f_right - pProjectToZeroOne[i];
			}
		}

		// optimal beta and projection of optimal beta into [0,1]

		optima_beta = (ConstituentParity-f_left)*(pSetBeta[right]-pSetBeta[left]);
		optima_beta = optima_beta/(f_right-f_left);
		optima_beta = optima_beta + pSetBeta[left];

		for(i=0;i<ReplicaLength;i++)
		{
			if(i<=ConstituentParity)
			{
				pProjectToZeroOne[i] = pReplicaWithOrder[i] - optima_beta;
			}
			else
			{
				pProjectToZeroOne[i] = pReplicaWithOrder[i] + optima_beta;
			}

			if(pProjectToZeroOne[i]>1)
			{
				pProjectToZeroOne[i] = 1;
			}
			else if (pProjectToZeroOne[i]<0)
			{
				pProjectToZeroOne[i] = 0;
			}
		}

		// inverse permutation  pPermutationOrder  pReplicaZ
		for(i=0;i<ReplicaLength;i++)
		{
			pReplicaZ[pPermutationOrder[i]] = pProjectToZeroOne[i];
		}

		delete [] pSetBeta;
	}

	delete [] pPermutationOrder;
	delete [] pReplicaWithOrder;
	delete [] pProjectToZeroOne;

}

void ProjectOnPPd_CSA(double *pReplicaBeforeProj, int ReplicaLength, double* pReplicaZ)
{
	int i,j,k;
	double *pProjOneZero;
	int *pIndicatorVec;
	int isBelongPPd;
	double optimalScalar;
	double dbtemp;

	pProjOneZero = new double [ReplicaLength];
	pIndicatorVec = new int [ReplicaLength];

	for(i=0;i<ReplicaLength;i++)
	{
		if(pReplicaBeforeProj[i]>1)
			pProjOneZero[i] = 1.0;
		else if(pReplicaBeforeProj[i]<0)
			pProjOneZero[i] = 0.0;
		else
			pProjOneZero[i] = pReplicaBeforeProj[i];
	}

	isBelongPPd = SearchCut(pProjOneZero,pIndicatorVec,ReplicaLength);

	if(isBelongPPd)
	{
		for(i=0;i<ReplicaLength;i++)
		{
			pReplicaZ[i] = pProjOneZero[i];
		}
	}
	else
	{
		optimalScalar = CalculateOptimalPara(pReplicaBeforeProj,pIndicatorVec,ReplicaLength);

		for(i=0;i<ReplicaLength;i++)
		{
			dbtemp = pReplicaBeforeProj[i] - optimalScalar * pIndicatorVec[i];
			if(dbtemp>1)
				pReplicaZ[i] = 1.0;
			else if(dbtemp<0)
				pReplicaZ[i] = 0.0;
			else
				pReplicaZ[i] = dbtemp;
		}
	}

	delete [] pProjOneZero;
	delete [] pIndicatorVec;
}

int SearchCut(double *pProjOneZero, int *pIndicatorVec, int ReplicaLength)
{
	int i,j;
	int NumIndicatorVec;
	int ClosetoHalf;
	int ret;
	double dbtemp,dbtemp1;

	NumIndicatorVec = 0;
	for(i=0;i<ReplicaLength;i++)
	{
		if(pProjOneZero[i]>0.5)
		{
			NumIndicatorVec++;
			pIndicatorVec[i] = 1;
		}
		else
			pIndicatorVec[i] = -1;
	}

	if(NumIndicatorVec % 2 == 0)
	{
		dbtemp = DBL_MAX;
		for(i=0;i<ReplicaLength;i++)
		{
			dbtemp1 = fabs(0.5 - pProjOneZero[i]);
			if(dbtemp1<dbtemp)
			{
				dbtemp = dbtemp1;
				ClosetoHalf = i;
			}
		}
		pIndicatorVec[ClosetoHalf] = -1*pIndicatorVec[ClosetoHalf];
	}

	int NumCuttingSet = 0;
	dbtemp = 0;
	for(i=0;i<ReplicaLength;i++)
	{
		dbtemp += pIndicatorVec[i] * pProjOneZero[i];
		if(pIndicatorVec[i]>0)
			NumCuttingSet++;
	}

	if(dbtemp > NumCuttingSet-1)
	{
		ret = 0;
	}
	else
	{
		ret = 1;
	}

	return ret;
}

double CalculateOptimalPara(double *pReplicaBeforeProj,int *pIndicatorVec,int ReplicaLength)
{
	int i,j;
	int NumCuttingSet;
	double delta,kesi;
	double dbtemp;
	double ret;

	NumCuttingSet = 0;
	dbtemp = 0;
    for(i=0;i<ReplicaLength;i++)
	{
		dbtemp += pIndicatorVec[i] * pReplicaBeforeProj[i];
		if(pIndicatorVec[i]>0)
			NumCuttingSet++;
	}

	delta = dbtemp - NumCuttingSet + 1;
	kesi = ReplicaLength;

	int itemp;
	int NumofT;
	double *T;

	NumofT = 0;
	for(i=0;i<ReplicaLength;i++)
	{
		if(pReplicaBeforeProj[i]>1 || pReplicaBeforeProj[i]<0)
			NumofT++;
	}

	if(NumofT>0)
	{
		T = new double [NumofT];
		itemp = 0;
		for(i=0;i<ReplicaLength;i++)
		{
			if(pReplicaBeforeProj[i]>1)
			{
				T[itemp++] = pReplicaBeforeProj[i] - 1;
			}
			if(pReplicaBeforeProj[i]<0)
			{
				T[itemp++] = -1*pReplicaBeforeProj[i];
			}
		}

		qsort(T,NumofT,sizeof(T[0]),cmp);

		for(i=0;i<NumofT;i++)
		{
			dbtemp = delta/kesi;
			if(dbtemp>T[i])
			{
				break;
			}
			else
			{
				delta = delta - T[i];
				kesi = kesi - 1;
			}
		}


		delete [] T;
	}

    ret = delta/kesi;
	return ret;
}

int cmp(const void * a, const void * b)
{
	return((*(double*)a-*(double*)b>0)? -1 : 1);
}

