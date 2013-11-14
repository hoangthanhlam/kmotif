/***********************************************************************************************************************************/
//				This C++ source code is written by Hoang Thanh Lam and Ninh Pham as a part of the NWO COMPASS project 
// Lam is a PhD student at Eindhoven University of Technology (TU/e), Eindhoven, the Netherlands
// Ninh is a graduate student at University of Technology in Ho Chi Minh city Vietnam
/***********************************************************************************************************************************/
// The code is written based on the original source code  by  Abdullah Al Mueen from University of California Riverside      
// We implement  the kmotif and the naive solutions for the top-k motif problem. For more information about the algorithm description
//	please refer to the technical report mentioned below.
/***********************************************************************************************************************************/
//										A small subproject in the hot FIFA World Cup summer 2010
// Feel free to re-use and re-distribute this source code for any purpose, please cite our following work when you re-use this code:                                 
/***********************************************************************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <deque>

using namespace std;

#define INF					9999999.0
#define EPS					0.00001

#define __DB_SIZE			5000
#define __MOTIF_LENGTH		128
#define __WINDOW_SIZE		1000
#define __K					10
#define __FILE_NAME			"C:\MATLAB7\work\Dataset_1\random_walk.dat"

//internal data structure
double **data;

//input variables

long long DB_SIZE;
long long MOTIF_LENGTH;
long long WINDOW_SIZE;
int K;
int ZZ = 0;

// Clock
double	TOTAL_TIME = 0;
double	DIST_TIME = 0;	
double	UPDATE_TIME = 0;

double  AVG_UPDATE_TIME = 0;
//double	AVG_DIST_TIME_PER_UPDATE = 0;
double	kMOTIF_SPACE_USAGE = 0;
bool	NEED_CALCULATE = false;


//Our data structure contains a nodes list, every node in this list is a node structure
struct node{
	double *data; //the vector representing the node
	list<struct point> sNN; //the list of the nearest neighbors in the suffix of the data structure starting from this node 
	int ID;
};

// every point is a part of the SNN list of a given node
struct point{
	double ed; //the Euclidean distance from the given node to the node with identifier toID  
	int fromID;
	int toID; 
};

//compare two point structures, comparison reference for sorting
bool pointCMP(const point &p1, const point &p2)
{
	if ( p1.ed != p2.ed)
		return p1.ed < p2.ed;
	else
		return p1.fromID > p2.fromID;
}

void error(int id)
{
    if ( id == 1 )
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
        printf("ERROR : Invalid Number of Arguments!!!\n\n");
    system("PAUSE");
    exit(1);

}

//Euclidean distance
double Euclidean_Dist(double *x, double *y, int p_iMotifLength , double p_dBSF )
{	
	//double dStartTime = clock();
    double dSum2 = 0;
    double dBSF2 = p_dBSF * p_dBSF;

    for ( int i = 0; (i < p_iMotifLength) && (dSum2 <= dBSF2); i++ )
        dSum2 += (x[i] - y[i]) * (x[i] - y[i]);

	//double dStopTime = clock();
	//DIST_TIME += ( dStopTime - dStartTime ) / CLOCKS_PER_SEC;

	//if ( NEED_CALCULATE )
	//	AVG_DIST_TIME_PER_UPDATE += ( dStopTime - dStartTime ) / CLOCKS_PER_SEC;

	return sqrt(dSum2);
}

/* This is a naive summary update function with space complexity O(kw)

	Input:	p_iK			: K in TopK Query
			p_iMotifLength	: The length of motif
*/

void naive_Update(deque<struct node> &nodes, int p_iK, int p_iMotifLength)
{
	deque<struct node>::iterator ni;
	list<struct point>::iterator low;
	struct point p;
	p.toID = nodes.back().ID;
	int iCount = 0;
	int mk=nodes.size()- p_iMotifLength;
	if(p_iK==1){
		for ( ni = nodes.begin(); iCount < mk ; ++ni )  // Eliminate trivial matching
		{
			p.fromID = ni->ID;
			iCount++;
		    if ( ni->sNN.size() < p_iK )
			{
				p.ed = Euclidean_Dist(nodes.back().data, ni->data, MOTIF_LENGTH, INF);
				ni->sNN.push_back(p);
			} 
			else 
			{
				p.ed = Euclidean_Dist(nodes.back().data, ni->data, MOTIF_LENGTH, ni->sNN.back().ed);
				if ( p.ed < ni->sNN.back().ed )
				{
					ni->sNN.pop_back();
					ni->sNN.push_back(p);
				}
			}
		}
	} else {
		for ( ni = nodes.begin(); iCount < mk; ++ni )  // Eliminate trivial matching
		{
			p.fromID = ni->ID;
			iCount++;
		    if ( ni->sNN.size() < p_iK )
			{
				p.ed = Euclidean_Dist(nodes.back().data, ni->data, MOTIF_LENGTH, INF);
				low = lower_bound(ni->sNN.begin(), ni->sNN.end(), p, pointCMP);
				ni->sNN.insert(low, p);
			} 
			else 
			{
				p.ed = Euclidean_Dist(nodes.back().data, ni->data, MOTIF_LENGTH, ni->sNN.back().ed);
				low = lower_bound(ni->sNN.begin(), ni->sNN.end(), p, pointCMP);
				if ( low != ni->sNN.end() )
				{
					ni->sNN.insert(low, p);
					ni->sNN.pop_back();
				}
			}
		}
	}
}


/* This is a another summary update function with provably prove average space complexity w+klogw 

	Input:	p_iK			: K in TopK Query
			p_iMotifLength	: The length of motif

*/
void kmotif_back_Update(deque<struct node> &nodes, int p_iK, int p_iMotifLength)
{
	deque<struct node>::iterator ni;
	deque<double> vf;
	struct point p;
	list<struct point>::iterator pi, low;
	deque<double>::iterator up;

	p.toID = nodes.back().ID;
	int iCount = 0;
	ni = nodes.end();
	ni--;

	while ( iCount < ((int)nodes.size() - 1) )
	{
		iCount++;
		if ( ni != nodes.begin() )
			ni--;
		p.fromID = ni->ID;
		for ( pi = ni->sNN.begin(); pi != ni->sNN.end(); pi++ )
		{
			up = upper_bound(vf.begin(), vf.end(), pi->ed);
			if ( up - vf.begin() <= p_iK )
				vf.insert(up, pi->ed);
			else
			{
				ni->sNN.erase(pi, ni->sNN.end());
				break;
			}
		}
		//avoid trivial matches 
		if ( iCount >= p_iMotifLength ) 
		{
			if ( vf.size() >= p_iK )
				p.ed = Euclidean_Dist(nodes.back().data, ni->data, MOTIF_LENGTH, vf[p_iK - 1]);
			else
				p.ed = Euclidean_Dist(nodes.back().data, ni->data, MOTIF_LENGTH, INF);

			up = upper_bound(vf.begin(), vf.end(), p.ed);
			if ( up - vf.begin() <= p_iK )
			{
				vf.insert(up, p.ed);
				low = lower_bound(ni->sNN.begin(), ni->sNN.end(), p, pointCMP);
				ni->sNN.insert(low, p);
			}
		}
	}

	if ( NEED_CALCULATE )
		kMOTIF_SPACE_USAGE += vf.size();	

}

//print top-k motif pairs
void query(deque<struct node> &nodes, int p_iK)
{
	list<struct point> output;
	deque<struct node>::iterator ni;
	list<struct point>::iterator pi;

	for ( ni = nodes.begin(); ni != nodes.end(); ni++)
	{
		for ( pi = ni->sNN.begin(); pi != ni->sNN.end(); pi++ )
		{
			output.push_back(*pi);
			//for(int it=0;it<=m;it++)
			//	cout <<ni->data[it]<< " a ";
			//cout<<l<<endl;
			//for( pi=ni->sNN.begin();pi!=ni->sNN.end();pi++)
			//	cout <<pi->ed<< " "<<pi->toID<< " ";
			//cout<<endl;
		}
	}

	// Not need to query
	output.sort(pointCMP);
	int iCount = 0;
	cout << "Top-k motif pairs are:" << endl;

	list<struct point> listMotif;
	list<struct point>::iterator iterlistMotif;

	for ( pi = output.begin(); iCount < p_iK && pi != output.end(); pi++)
	{
		int iFromID = pi->fromID;
		int iToID   = pi->toID;

		if ( listMotif.empty() )
		{
			listMotif.push_back(*pi);			
			cout << pi->ed << "  " << iFromID << " " << iToID << " | ";
			iCount++;
		}
		else
		{
			bool bIsTrueMotif = true;
			for ( iterlistMotif = listMotif.begin(); iterlistMotif != listMotif.end(); iterlistMotif++ )
			{
				if ( abs( ( iterlistMotif->fromID - iFromID) ) <= MOTIF_LENGTH || 
					 abs( ( iterlistMotif->toID   - iToID) ) <= MOTIF_LENGTH ||
					 abs( ( iterlistMotif->fromID - iToID) ) <= MOTIF_LENGTH || 
					 abs( ( iterlistMotif->toID   - iFromID) ) <= MOTIF_LENGTH )
				{
					bIsTrueMotif = false;					
					break;
				}
			}
			if ( bIsTrueMotif )
			{
				listMotif.push_back(*pi);
				cout << pi->ed << "  " << iFromID << " " << iToID << " | ";
				iCount++;
			}
		}
		
		
	}

	
	cout << endl;
	cout << "Size of the summary: " << output.size() << " elements" << endl;
	cout << endl;
}

int main(int argc, char *argv[])
{
	double dValue;
	double dStartTotalTime, dStopTotalTime, dStartUpdateTime, dStopUpdateTime;

	int iOption = 0; // 0 is heuristic algorithm and 1 is naive algorithm
	int iQueryPoint = 1; // The number of points for executing queries
    
    FILE * pFileName;
    
	deque<struct node> nodes;  
	struct node nd;

	// Initialize param
	DB_SIZE			= __DB_SIZE;
   	MOTIF_LENGTH	= __MOTIF_LENGTH;
    WINDOW_SIZE		= __WINDOW_SIZE;
	K				= __K;
    ZZ				= 1;

	// Initialize param with command line values
    if ( argc > iQueryPoint )
    {
		pFileName = fopen(argv[1], "r");
		if( pFileName == NULL ) 
			error(2);
		DB_SIZE = atol( argv[2] );
		WINDOW_SIZE = atol( argv[3] );
		MOTIF_LENGTH = atol( argv[4] );
		K = atol( argv[5] );
		iOption = atol( argv[6] );
		iQueryPoint = atol( argv[7] );
    }
    else
	{
		printf("Command file data_size window_length motif_length k naive_or_kmotif query_every_window_of_size_q\n");
		pFileName = fopen("EOG.txt", "r");
		if( pFileName == NULL ) 
			error(2);
	}
	
	dStartTotalTime = clock();
	srand(time(NULL));

	DB_SIZE = DB_SIZE - MOTIF_LENGTH + 1;

    data = (double **)malloc(sizeof(double *) * DB_SIZE);
    if( data == NULL  )
		error(1);

    //This is the loop that iterates over each clock tick
    int l = 1;
	int i = 0;
	int j = 0;
	int MOTIF_MAX_INDEX = MOTIF_LENGTH - 1; // Because index of an array start from 0
	double dSum , dSum2 , dMean, dSTD;
	
	dSum = dSum2 = 0;
	int cc=0;
	while( i < DB_SIZE  )
    {
		if ( fscanf(pFileName, "%lf", &dValue) ==  EOF )
			error(1);

		if ( i >= WINDOW_SIZE )
			NEED_CALCULATE = true;

		if ( i < DB_SIZE )
        {
            data[i] = (double *)malloc(sizeof(double) * MOTIF_LENGTH);
			if( data[i] == NULL )
                error(1);
        }

        dSum += dValue;
        dSum2 += dValue * dValue;
       
		// Problem - Does we need ( i >= DB_SIZE )
        for ( j = ( (i >= DB_SIZE) ? DB_SIZE - 1 : i ); j >= ( (i >= MOTIF_LENGTH) ? (i - MOTIF_MAX_INDEX) : 0 ); j--)
            data[j][i - j] = dValue;

        if ( i >= MOTIF_MAX_INDEX )
        {
            //z-Normalize the seq Node created LENGTH ticks before
            dMean = dSum / MOTIF_LENGTH;
            dSTD = dSum2 / MOTIF_LENGTH;
            dSTD = sqrt(dSTD - dMean * dMean);

			// Keep dSum and dSum2 for next interaction
           	dSum -= data[i - MOTIF_MAX_INDEX][0];
            dSum2 -= data[i - MOTIF_MAX_INDEX][0] * data[i - MOTIF_MAX_INDEX][0];

            if ( ZZ == 1 )
            {
				for ( long long u = 0 ; u < MOTIF_LENGTH ; u++ )
				{
					if ( dSTD > 0 )
						data[i - MOTIF_MAX_INDEX][u] = (data[i - MOTIF_MAX_INDEX][u] - dMean) / dSTD;
					else
						data[i - MOTIF_MAX_INDEX][u] = (data[i - MOTIF_MAX_INDEX][u] - dMean) / (dSTD + EPS);
				}
            }

			nd.data = data[i - MOTIF_MAX_INDEX];
			nd.ID = i - MOTIF_MAX_INDEX;
			nodes.push_back(nd);

			// Update structure
			//dStartUpdateTime = clock();
			
			if ( iOption == 1 ){
				naive_Update(nodes, K, MOTIF_LENGTH);				
			}else
				kmotif_back_Update(nodes, K, MOTIF_LENGTH);

			//dStopUpdateTime = clock();
			//UPDATE_TIME += ( dStopUpdateTime - dStartUpdateTime ) / CLOCKS_PER_SEC;
			
			//if ( NEED_CALCULATE )
				//AVG_UPDATE_TIME += ( dStopUpdateTime - dStartUpdateTime ) / CLOCKS_PER_SEC;

			if ( l % iQueryPoint == 0 )
			{
				/*if ( NEED_CALCULATE )
					printf("Average space usage / 2Klog(W) : %lf\n", kMOTIF_SPACE_USAGE / ( (l - WINDOW_SIZE) * 2 * K * ( log(double(WINDOW_SIZE)) / log(2.0) ) ) );*/
				cout << l << endl;
				query(nodes, K);
			}
			l++;
	     }//end else of if i >= k

	     if ( i >= WINDOW_SIZE )
         {
            nodes.pop_front();
            free(data[i - WINDOW_SIZE]);
		 }//if i>=w
	     i++;
    }//end while(fscanf////.....


	dStopTotalTime = clock();

	TOTAL_TIME = (dStopTotalTime - dStartTotalTime) / CLOCKS_PER_SEC;

	// Print Time
	printf("Total time : %lf\n", TOTAL_TIME);
	printf("Total update time : %lf\n", UPDATE_TIME);
	printf("Total distance time : %lf\n", DIST_TIME);

	// Average performance
	printf("Average time per update : %lf\n", AVG_UPDATE_TIME / (DB_SIZE - WINDOW_SIZE) );
	//printf("Average distance time per update : %lf\n", AVG_DIST_TIME_PER_UPDATE / (DB_SIZE - WINDOW_SIZE) );
	printf("Average space usage per update : %lf\n", kMOTIF_SPACE_USAGE / (DB_SIZE - WINDOW_SIZE) );	

	fclose(pFileName);
    system("PAUSE");
    return 1;   
}

