// ==============================================================================
//
//  Scatter3D.cpp
//  QTR
//
//  Created by Albert Lu on 8/6/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/10/18
//
//  Note:
//
// ==============================================================================

//#define BOOST_DISABLE_ASSERTS

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <omp.h>
#include <parallel/algorithm>
#include <vector>

#include "Containers.h"
#include "Error.h"
#include "Log.h"
#include "Parameters.h"
#include "Scatter3d.h"

#include "boost/multi_array.hpp"

using namespace QTR_NS;
using std::vector;
using std::cout;
using std::endl;

/* ------------------------------------------------------------------------------- */

Scatter3d::Scatter3d(class QTR *q)
{
    qtr = q;
    err = qtr->error;
    log = qtr->log;
    parameters = qtr->parameters;
    init();
}
/* ------------------------------------------------------------------------------- */

Scatter3d::~Scatter3d()
{     
    return;
}
/* ------------------------------------------------------------------------------- */

void Scatter3d::init()
{
    log->log("\n\n[Scatter3d] INIT starts ...\n");

    // General parameters
    I = {0,1}; // sqrt(-1)
    xZERO = {0,0}; // complex zero
    DIMENSIONS = parameters->scxd_dimensions;
    PERIOD = parameters->scxd_period;
    TIME = parameters->scxd_Tf;

    // Grid size
    H.resize(DIMENSIONS);
    Hi.resize(DIMENSIONS);
    Hisq.resize(DIMENSIONS);
    H[0] = parameters->scxd_h1;
    H[1] = parameters->scxd_h2;
    H[2] = parameters->scxd_h3;
    Hi[0] = 1 / H[0];
    Hi[1] = 1 / H[1];    
    Hi[2] = 1 / H[2];  
    Hisq[0] = 1 / pow(H[0],2);
    Hisq[1] = 1 / pow(H[1],2);    
    Hisq[2] = 1 / pow(H[2],2);  

    S.resize(DIMENSIONS);  
    kk = parameters->scxd_k;
    S << kk * Hisq[0], kk * Hisq[1], kk * Hisq[2]; 

    // Domain size and # grids
    Box.resize(DIMENSIONS * 2);
    Box[0] = parameters->scxd_xi1;
    Box[1] = parameters->scxd_xf1;
    Box[2] = parameters->scxd_xi2;
    Box[3] = parameters->scxd_xf2; 
    Box[4] = parameters->scxd_xi3;
    Box[5] = parameters->scxd_xf3;    

    BoxShape.resize(DIMENSIONS);
    BoxShape[0] = (int)( (Box[1] - Box[0]) / H[0] ) + 1;
    BoxShape[1] = (int)( (Box[3] - Box[2]) / H[1] ) + 1;
    BoxShape[2] = (int)( (Box[5] - Box[4]) / H[2] ) + 1;
	GRIDS_TOT = BoxShape[0] * BoxShape[1] * BoxShape[2];

    log->log("[Scatter3d] Number of grids = (%d,%d,%d)\n", BoxShape[0], BoxShape[1], BoxShape[2]);

    idx_x0 = (int) ( ( 0 - Box[0] ) / H[0] );
	log->log("[Scatter3d] x = 0 at idx = %d\n", idx_x0);

    // Potential: form
    Vmode.resize(DIMENSIONS);
    Vmode[0] = parameters->scxd_Vmode_1;
    Vmode[1] = parameters->scxd_Vmode_2;
    Vmode[2] = parameters->scxd_Vmode_3;
    hb = parameters->scxd_hb;
    m  = parameters->scxd_m;

    // Potential: HO specific
    w = parameters->scxd_w;
    aa = 0.5 * m * w / hb;
    log->log("[Scatter3d] aa = %lf\n", aa);

    // Potential: Eckart
    V0 = parameters->scxd_V0;
    alpha = parameters->scxd_alpha;
    ek2v = parameters->scxd_ek2v;
    Ek0 = V0 * ek2v;
    log->log("[Scatter3d] Ek0 = %lf\n", Ek0);

    // Potential: Related HO
    k0 = parameters->scxd_k0;  
    sig = parameters->scxd_sig;
    lan = parameters->scxd_lan;

    // Potential: Morse
    De = parameters->scxd_De;
    Da = parameters->scxd_Da;
    r0 = parameters->scxd_r0;
    Ld = sqrt(2 * m * De) / (Da * hb);  
    log->log("[Scatter3d] Ld = %lf\n", Ld);

    // Wavefunction parameters
    Wave0.resize(DIMENSIONS);
    Wave0[0] = parameters->scxd_x01;
    Wave0[1] = parameters->scxd_x02;
    Wave0[2] = parameters->scxd_x03;

    A.resize(DIMENSIONS);
    A[0] = parameters->scxd_a1;
    A[1] = parameters->scxd_a2;
    A[2] = parameters->scxd_a3;

    P.resize(DIMENSIONS);
    P[0] = parameters->scxd_p1;
    P[1] = parameters->scxd_p2;
    P[2] = parameters->scxd_p3;
    // Ｏverride
    P[0] = sqrt(2.0 * m * Ek0);
    P[1] = 0;
    P[2] = 0;                  
    log->log("[Scatter3d] P[0] = %lf\n", P[0]);

    // Truncate parameters
    isFullGrid = parameters->scxd_isFullGrid;
    TolH = parameters->scxd_TolH;       // Tolerance of probability density for Zero point Cutoff
    TolL = parameters->scxd_TolL;       // Tolerance of probability density for Edge point
    TolHd = parameters->scxd_TolHd;     // Tolerance of probability first diff for Zero point Cutoff
    TolLd = parameters->scxd_TolLd;     // Tolerance of probability density for Edge point
    ExReduce = parameters->scxd_ExReduce;   //Extrapolation reduce factor

    log->log("[Scatter3d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void Scatter3d::Evolve()
{
	#pragma omp declare reduction (merge : MeshIndex : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    log->log("[Scatter3d] Evolve starts ...\n");

    // Variables 
    int     index;
    bool    b1, b2, b3, b4, b5, b6;
    double  kh2m = kk * hb / 2 / m;
    double  k2hb = kk / hb;
    double  coeff;
    double  pftrans;
    double  t_0_begin, t_0_end, t_0_elapsed;
    double  t_1_begin, t_1_end, t_1_elapsed;
    MeshIndex tmpVec;  // temporary index container
     
    // normalization factor
    double  norm;
 
    // 3d Grid vector and indices
    Vector3i  grid;
    int     gx, gy, gz;
    double  xx, yy, zz;

    // Vector iterater
    std::vector<int>::iterator it;

    // Extrapolation 
    int count;
    std::vector<std::complex<double>> ExTBL;
    std::complex<double> val, val_min;
    std::complex<double> sum;
    double val_min_abs;

    // PF_trans
    std::vector<double> PF_trans;

    // Make DBi and DBi2
    if (isFullGrid)  {
    	DefineBoundary();
    }

    log->log("[Scatter3d] Initializing containers ...\n");

    t_1_begin = omp_get_wtime();

    // Initialize containers
    MeshCX3D F(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    F.reindex(0);

    MeshD3D PF(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    PF.reindex(0);
    
    MeshD3D PFdX(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    PFdX.reindex(0);  

    MeshD3D PFdY(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    PFdY.reindex(0);

    MeshD3D PFdZ(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    PFdZ.reindex(0); 

    MeshCX3D FF(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    FF.reindex(0);   

    MeshCX3D KK1(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    KK1.reindex(0);

    MeshCX3D KK2(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    KK2.reindex(0);  

    MeshCX3D KK3(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    KK3.reindex(0);

    MeshCX3D KK4(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]]);
    KK4.reindex(0);   

    #pragma omp parallel for

    for (int i = 0; i < BoxShape[0]; i ++)  {

        for (int j = 0; j < BoxShape[1]; j ++)  {

            for (int k = 0; k < BoxShape[2]; k++)  {
        
                F[i][j][k] = xZERO;
                PF[i][j][k] = 0.0;
                PFdX[i][j][k] = 0.0;
                PFdY[i][j][k] = 0.0;
                PFdZ[i][j][k] = 0.0;
            }
        }
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    log->log("[Scatter3d] Elapsed time (initializing containers) = %lf sec\n\n", t_1_elapsed); 

    // .........................................................................................

    log->log("[Scatter3d] Initializing wavefunction ...\n");  

    t_1_begin = omp_get_wtime();

    // Initialize wavefunction

    #pragma omp parallel for

    for (int i = 1; i < BoxShape[0] - 1 ; i ++)  {

        for (int j = 1; j < BoxShape[1] - 1 ; j ++)  {

            for (int k = 1; k < BoxShape[2] - 1 ; k++)  {
          
                F[i][j][k] = Wavefunction(Box[0] + i * H[0], Box[2] + j * H[1], Box[4] + k * H[2]);
            }
        }
    }

    // Normalization

    norm = 0.0;

    #pragma omp parallel for reduction (+:norm)

    for (int i = 0; i <  BoxShape[0]; i++)  {

        for (int j = 0; j < BoxShape[1]; j++)  {

            for (int k = 0; k < BoxShape[2]; k++)  {
                
				norm += std::abs(F[i][j][k] * std::conj(F[i][j][k]));
            }
        }
    }
    norm *= H[0] * H[1] * H[2];

    log->log("[Scatter3d] Normalization factor = %e\n",norm);
    norm = 1.0 / sqrt(norm);

    #pragma omp parallel for

    for (int i = 0; i < BoxShape[0]; i ++)  {

        for (int j = 0; j < BoxShape[1]; j ++)  {

            for (int k = 0; k < BoxShape[2]; k++)  {

                F[i][j][k] = norm * F[i][j][k];
            }
        }
    }     

    // Initial probability

    #pragma omp parallel for

    for (int i = 0; i < BoxShape[0]; i ++)  {

        for (int j = 0; j < BoxShape[1]; j ++)  {

            for (int k = 0; k < BoxShape[2]; k++)  {

                PF[i][j][k] = std::abs(F[i][j][k] * std::conj(F[i][j][k]));
            }
        }
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    log->log("[Scatter3d] Elapsed time (initializing wavefunction) = %lf sec\n\n", t_1_elapsed); 

    // .........................................................................................

    log->log("[Scatter3d] Computing finite differences ...\n");   

    t_1_begin = omp_get_wtime();

    // Pf 1st-order finite difference 

    // PFdX

    coeff = 1.0 / (2.0 * H[0]);

    #pragma omp parallel for

    for (int i = 1; i < BoxShape[0] - 1; i ++)  {

        for (int j = 1; j < BoxShape[1] - 1; j ++)  {

            for (int k = 1; k < BoxShape[2] - 1; k++)  {

                PFdX[i][j][k] = std::abs(F[i+1][j][k] - F[i-1][j][k]) * coeff;
            }
        }
    }    

    // PFdY    

    coeff = 1.0 / (2.0 * H[1]);

    #pragma omp parallel for

    for (int i = 1; i < BoxShape[0] - 1; i ++)  {

        for (int j = 1; j < BoxShape[1] - 1; j ++)  {

            for (int k = 1; k < BoxShape[2] - 1; k++)  {

                PFdY[i][j][k] = std::abs(F[i][j+1][k] - F[i][j-1][k]) * coeff;
            }
        }
    }  

    // PFdZ    

    coeff = 1.0 / (2.0 * H[2]);

    #pragma omp parallel for

    for (int i = 1; i < BoxShape[0] - 1; i ++)  {

        for (int j = 1; j < BoxShape[1] - 1; j ++)  {

            for (int k = 1; k < BoxShape[2] - 1; k++)  {

                PFdZ[i][j][k] = std::abs(F[i][j][k+1] - F[i][j][k-1]) * coeff;
            }
        }
    }  

    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    log->log("[Scatter3d] Elapsed time (finite differences) = %lf sec\n\n", t_1_elapsed); 

    // .........................................................................................

    // Truncate initial & edge point check

    log->log("[Scatter3d] Initial truncation ...\n");   

    t_1_begin = omp_get_wtime();

	if ( !isFullGrid )  // Truncate method
    {
        #pragma omp parallel for

        for (int i = 1; i < BoxShape[0] - 1; i++)  {

            for (int j = 1; j < BoxShape[1] - 1; j++)  {

                for (int k = 1; k < BoxShape[2] - 1; k++)  {
          
                    if ( ( PF[i][j][k] < TolH ) && ( PFdX[i][j][k] < TolHd ) && ( PFdY[i][j][k] < TolHd ) && ( PFdZ[i][j][k] < TolHd ) )
                        F[i][j][k] = xZERO;
                }
            }
        }

        // `````````````````````````````````````````````````````````````````
        
		#pragma omp parallel for reduction(merge: tmpVec)

        for (int i = 1; i < BoxShape[0] - 1; i++)  {

            for (int j = 1; j < BoxShape[1] - 1; j++)  {

                for (int k = 1; k < BoxShape[2] - 1; k++)  {
          
                    if ( F[i][j][k] != std::complex<double>(0,0) || ( PFdX[i][j][k] >= TolHd || PFdY[i][j][k] >= TolHd || PFdZ[i][j][k] >= TolHd ))  {

                        tmpVec.push_back(GridToIdx(i,j,k));                                   
                    }
                }
            }
        }
		tmpVec.swap(TA);
        tmpVec.clear();

        log->log("[Scatter3d] TA size = %d\n", TA.size());

        // `````````````````````````````````````````````````````````````````

        #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz, b1, b2, b3, b4, b5, b6)

        for (int i = 0; i < TA.size(); i++)
        {
            grid = IdxToGrid(TA[i]);
            gx = grid[0];
            gy = grid[1];
            gz = grid[2];

            b1 = F[ gx - 1 ][ gy  ][ gz   ] == xZERO;
            b2 = F[ gx + 1 ][ gy  ][ gz   ] == xZERO;
            b3 = F[ gx   ][ gy -1 ][ gz   ] == xZERO;
            b4 = F[ gx   ][ gy +1 ][ gz   ] == xZERO;
            b5 = F[ gx   ][ gy  ][ gz - 1 ] == xZERO;
            b6 = F[ gx   ][ gy  ][ gz + 1 ] == xZERO;
               
            if ( b1 || b2 || b3 || b4 || b5 || b6 )
            {
                b1 = ( PFdX[ gx - 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdY[ gx - 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdZ[ gx - 1 ][ gy   ][ gz   ] < TolHd );
                b2 = ( PFdX[ gx + 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdY[ gx + 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdZ[ gx + 1 ][ gy   ][ gz   ] < TolHd );
                b3 = ( PFdX[ gx   ][ gy - 1 ][ gz   ] < TolHd ) && ( PFdY[ gx   ][ gy - 1 ][ gz   ] < TolHd ) && ( PFdZ[ gx   ][ gy - 1 ][ gz   ] < TolHd );
                b4 = ( PFdX[ gx   ][ gy + 1 ][ gz   ] < TolHd ) && ( PFdY[ gx   ][ gy + 1 ][ gz   ] < TolHd ) && ( PFdZ[ gx   ][ gy + 1 ][ gz   ] < TolHd );
                b5 = ( PFdX[ gx   ][ gy   ][ gz - 1 ] < TolHd ) && ( PFdY[ gx   ][ gy   ][ gz - 1 ] < TolHd ) && ( PFdZ[ gx   ][ gy   ][ gz - 1 ] < TolHd );
                b6 = ( PFdX[ gx   ][ gy   ][ gz + 1 ] < TolHd ) && ( PFdY[ gx   ][ gy   ][ gz + 1 ] < TolHd ) && ( PFdZ[ gx   ][ gy   ][ gz + 1 ] < TolHd );

                if ( b1 || b2 || b3 || b4 || b5 || b6 ) {

                    tmpVec.push_back(TA[i]);  
                }
            }       
        }

        tmpVec.swap(TB);
        tmpVec.clear();

        log->log("[Scatter3d] TB size = %d\n", TB.size());


        // `````````````````````````````````````````````````````````````````

        #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz)

        for (int i = 0; i < TA.size(); i++)
        {
            grid = IdxToGrid(TA[i]);
            gx = grid[0];
            gy = grid[1];
            gz = grid[2];

            tmpVec.push_back(GridToIdx(gx + 1, gy  , gz  ));
            tmpVec.push_back(GridToIdx(gx - 1, gy  , gz  ));
            tmpVec.push_back(GridToIdx(gx  , gy + 1, gz  ));
            tmpVec.push_back(GridToIdx(gx  , gy - 1, gz  ));  
            tmpVec.push_back(GridToIdx(gx  , gy  , gz + 1));
            tmpVec.push_back(GridToIdx(gx  , gy  , gz - 1));            
        }

        // Combine TA and tmpVec
        TA.reserve(TA.size() + tmpVec.size());
        TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());
        tmpVec.clear();

        // Find unique elements
        __gnu_parallel::sort(TA.begin(),TA.end());
        it = std::unique (TA.begin(), TA.end()); 
        TA.resize(std::distance(TA.begin(),it));

        // Erase elements at the domain boundary (DBi)

		it = TA.begin();

		while ( it != TA.end() )
        {
			index = std::distance( TA.begin(), it );
            grid = IdxToGrid(TA[index]);
            gx = grid[0];
            gy = grid[1];
            gz = grid[2];

            if (gx == 0 || gy == 0 || gz == 0 || gx == BoxShape[0] - 1 || gy == BoxShape[1] - 1 || gz == BoxShape[2] - 1)  {
                it = TA.erase(it);
            }
			else  {
				++it;
			}
        }          
    }
    else  // Full grid approach
    {
        TB = DBi;
    }

	log->log("[Scatter3d] TA size = %d, TB size = %d\n", TA.size(), TB.size());
	log->log("[Scatter3d] DBi = %d DBi2 = %d\n", DBi.size(), DBi2.size());

    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    log->log("[Scatter3d] Elapsed time (initial truncation) = %lf sec\n\n", t_1_elapsed); 

    // .........................................................................................

    // Time iteration 

    log->log("=======================================================\n\n"); 
    log->log("[Scatter3d] Time interation starts ...\n"); 
    log->log("[Scatter3d] Number of steps = %d\n\n", (int)(TIME / kk)); 
    log->log("=======================================================\n\n"); 

    for (int tt = 0; tt < (int)(TIME / kk); tt ++)
    {
		t_0_begin = omp_get_wtime();
        t_1_begin = omp_get_wtime();

    	#pragma omp parallel for

    	for (int i = 1; i < BoxShape[0] - 1 ; i ++)  {

        	for (int j = 1; j < BoxShape[1] - 1 ; j ++)  {

            	for (int k = 1; k < BoxShape[2] - 1 ; k++)  {
          
                	FF[i][j][k] = xZERO;
                	KK1[i][j][k] = xZERO;
                	KK2[i][j][k] = xZERO;
                	KK3[i][j][k] = xZERO;
                	KK4[i][j][k] = xZERO;
            	}
        	}
    	}

		t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        //log->log("Elapsed time (omp-a-1: fill_n) = %lf sec\n", t_1_elapsed);   

        // Check if TB of f is higher than TolL
        
        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();

            TBL.clear();

            // omp-a-1
            #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz, b1, b2, b3, b4, b5, b6)

            for (int i = 0; i < TB.size(); i++)
            {
                grid = IdxToGrid(TB[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                b1 = PF[gx][gy][gz] >= TolL;
                b2 = PFdX[gx][gy][gz] >= TolLd;
                b3 = PFdY[gx][gy][gz] >= TolLd;
                b4 = PFdZ[gx][gy][gz] >= TolLd;
                // Not in DBi2
                b5 = gx > 1 && gy > 1 && gz > 1;
                b6 = gx < BoxShape[0] - 2 && gy < BoxShape[1] - 2 && gz < BoxShape[2] - 2;

                if ( (b1 || b2 || b3 || b4 ) && b5 && b6 )  {
                    tmpVec.push_back(TB[i]);
                }
            }
			tmpVec.swap(TBL);
            tmpVec.clear();
            TBL_P = TBL;

			t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-a-2: TBL, TBL_P) = %lf sec\n", t_1_elapsed);   
        }    

        isExtrapolate = false;

        // CASE 1: Truncating with extrapolation

        while ( TBL.size() != 0 && !isFullGrid )
        {
            isExtrapolate = true;

            // Extrapolation3D
            // .............................................................................................

            t_1_begin = omp_get_wtime();

            // Avoid unexpected arrangement of TBL
            __gnu_parallel::sort(TBL.begin(),TBL.end());
            it = std::unique (TBL.begin(), TBL.end()); 
            TBL.resize(std::distance(TBL.begin(),it)); 

            // Find extrapolation target
            ExFF.clear();

            for ( it = TBL.begin(); it != TBL.end(); it ++ )  {

                index = std::distance( TBL.begin(), it );
                grid = IdxToGrid(TBL[index]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];  

                if ( F[gx - 1][gy][gz] == xZERO )  {

                    ExFF.push_back(GridToIdx(gx - 1, gy, gz));
                }

                if ( F[gx + 1][gy][gz] == xZERO )  {

                    ExFF.push_back(GridToIdx(gx + 1, gy, gz));
                }

                if ( F[gx][gy - 1][gz] == xZERO )  {

                    ExFF.push_back(GridToIdx(gx, gy - 1, gz));
                }

                if ( F[gx][gy + 1][gz] == xZERO )  {

                    ExFF.push_back(GridToIdx(gx, gy + 1, gz));
                }

                if ( F[gx][gy][gz - 1] == xZERO )  {

                    ExFF.push_back(GridToIdx(gx, gy, gz - 1));
                }

                if ( F[gx][gy][gz + 1] == xZERO )  {

                    ExFF.push_back(GridToIdx(gx, gy, gz + 1));
                }
            }

            // ExFF & TBL set difference
            tmpVec.resize(ExFF.size() + TBL.size());
            __gnu_parallel::sort (TBL.begin(), TBL.end());
            __gnu_parallel::sort (ExFF.begin(), ExFF.end());
            it=std::set_difference( ExFF.begin(), ExFF.end(), TBL.begin(), TBL.end(), tmpVec.begin() );
            tmpVec.resize(it - tmpVec.begin()); 
            tmpVec.swap(ExFF);
            tmpVec.clear();

            // Find unique elements
            __gnu_parallel::sort(ExFF.begin(),ExFF.end());
            it = std::unique (ExFF.begin(), ExFF.end()); 
            ExFF.resize(std::distance(ExFF.begin(),it));

			t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-b-1: ExFF) = %lf sec\n", t_1_elapsed);   

            // .....................................................................

            // Find the direction of Outer to Edge points

            t_1_begin = omp_get_wtime();

            ExTBL.clear();

			it = ExFF.begin();
            
            while ( it != ExFF.end() )  {

                index = std::distance( ExFF.begin(), it );
                grid = IdxToGrid(ExFF[index]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                sum = xZERO;
                count = 0;

                isEmpty = true;
                val_min_abs = 100000000;

                if ( F[gx - 1][gy][gz] != xZERO )  {

					if ( std::abs(F[gx - 1][gy][gz]) < val_min_abs &&  F[gx - 2][gy][gz] != xZERO )  {
                    	val_min_abs = std::abs(F[gx - 1][gy][gz]);
                    	val_min = F[gx - 1][gy][gz];
					}

                    if (F[gx - 2][gy][gz] != xZERO)  {

                        val = exp( 2.0 * std::log(F[gx - 1][gy][gz]) - std::log(F[gx - 2][gy][gz]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }
                        isEmpty = false;
                    }
                }

                if ( F[gx + 1][gy][gz] != xZERO )  {

                    if ( std::abs(F[gx + 1][gy][gz]) < val_min_abs && F[gx + 2][gy][gz] != xZERO  )  {
                        val_min_abs = std::abs(F[gx + 1][gy][gz]);
                        val_min = F[gx + 1][gy][gz];
                    }

                    if ( F[gx + 2][gy][gz] != xZERO )  {

                        val = exp( 2.0 * std::log(F[gx + 1][gy][gz]) - std::log(F[gx + 2][gy][gz]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }

                        isEmpty = false;
                    }
                }

                if ( F[gx][gy - 1][gz] != xZERO )  {

                    if ( std::abs(F[gx][gy - 1][gz]) < val_min_abs && F[gx][gy - 2][gz] != xZERO  )  {

                        val_min_abs = std::abs(F[gx][gy - 1][gz]);
                        val_min = F[gx][gy - 1][gz];
                    }

                    if ( F[gx][gy - 2][gz] != xZERO )  {

                        val = exp( 2.0 * std::log(F[gx][gy - 1][gz]) - std::log(F[gx][gy - 2][gz]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }

                        isEmpty = false;
                    }
                }

                if ( F[gx][gy + 1][gz] != xZERO )  {

                    if ( std::abs(F[gx][gy + 1][gz]) < val_min_abs && F[gx][gy+2][gz] != xZERO )  {

                        val_min_abs = std::abs(F[gx][gy + 1][gz]);
                        val_min = F[gx][gy + 1][gz];
                    }

                    if ( F[gx][gy + 2][gz] != xZERO )  {

                        val = exp( 2.0 * std::log(F[gx][gy + 1][gz]) - std::log(F[gx][gy + 2][gz]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }

                        isEmpty = false;
                    }
                }

                if ( F[gx][gy][gz - 1] != xZERO )  {

                    if ( std::abs(F[gx][gy][gz - 1]) < val_min_abs && F[gx][gy][gz-2] != xZERO )  {

                        val_min_abs = std::abs(F[gx][gy][gz - 1]);
                        val_min = F[gx][gy][gz - 1];
                    }

                    if ( F[gx][gy][gz - 2] != xZERO )  {

                        val = exp( 2.0 * std::log(F[gx][gy][gz - 1]) - std::log(F[gx][gy][gz - 2]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }

                        isEmpty = false;
                    }
                }

                if ( F[gx][gy][gz + 1] != xZERO )  {

                    if ( std::abs(F[gx][gy][gz + 1]) < val_min_abs && F[gx][gy][gz+2] != xZERO )  {

                        val_min_abs = std::abs(F[gx][gy][gz + 1]);
                        val_min = F[gx][gy][gz + 1];
                    }

                    if ( F[gx][gy][gz + 2] != xZERO )  {

                        val = exp( 2.0 * std::log(F[gx][gy][gz + 1]) - std::log(F[gx][gy][gz + 2]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }

                        isEmpty = false;
                    }
                }

                if ( isEmpty )
                {
                    it = ExFF.erase(it);   

                }  else  {

                    // Assume the probability of outer points always smaller than edge point.
                    // if larger, then set the edge point with smallest P to outer
                    // point instead of using the extrapolation result.

                    if ( std::abs( std::complex<double> (std::real(sum) / count, std::imag(sum) / count) ) > val_min_abs * exp(-ExReduce * H[0]) )  {


                        ExTBL.push_back( val_min * exp(-ExReduce * H[0]) );
                    }
                    else  {

                        ExTBL.push_back( std::complex<double> (std::real(sum) / count, std::imag(sum) / count) );
                    }
					++it;
                }
            }

            #pragma omp parallel for private(grid, gx, gy, gz)

            for ( int i = 0; i < ExFF.size(); i++ )  {

                grid = IdxToGrid(ExFF[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                F[gx][gy][gz] = ExTBL[i];
            }

			t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-b-2: ExFF) = %lf sec\n", t_1_elapsed);  

            // ............................................................................................. Extrapolation

            // Check Extending nonzero Area

            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz)

            for (int i = 0; i < ExFF.size(); i++)
            {
                grid = IdxToGrid(ExFF[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                tmpVec.push_back(GridToIdx(gx  , gy  , gz  ));
                tmpVec.push_back(GridToIdx(gx + 1, gy  , gz  ));
                tmpVec.push_back(GridToIdx(gx - 1, gy  , gz  ));
                tmpVec.push_back(GridToIdx(gx  , gy + 1, gz  ));
                tmpVec.push_back(GridToIdx(gx  , gy - 1, gz  ));  
                tmpVec.push_back(GridToIdx(gx  , gy  , gz + 1));
                tmpVec.push_back(GridToIdx(gx  , gy  , gz - 1));            
            }

            // Combine TA and tmpVec
            TA.reserve(TA.size() + tmpVec.size());
            TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());
            tmpVec.clear();

            // Find unique elements
            __gnu_parallel::sort(TA.begin(),TA.end());
            it = std::unique (TA.begin(), TA.end()); 
            TA.resize(std::distance(TA.begin(),it));

			t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-c-1: CASE 1 TA) = %lf sec\n", t_1_elapsed); 
     
            // Runge–Kutta 4
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];

                KK1[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( F[gx - 1][gy][gz] - 2.0 * F[gx][gy][gz] + F[gx + 1][gy][gz] )
                                + Hisq[1] * ( F[gx][gy - 1][gz] - 2.0 * F[gx][gy][gz] + F[gx][gy + 1][gz] ) 
                                + Hisq[2] * ( F[gx][gy][gz - 1] - 2.0 * F[gx][gy][gz] + F[gx][gy][gz + 1] ) 
                               )  - I * k2hb * Potential(xx, yy, zz) * F[gx][gy][gz];
            }

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];

                KK2[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 0.5 * KK1[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 0.5 * KK1[gx + 1][gy][gz] ) )
                                + Hisq[1] * ( ( F[gx][gy - 1][gz] + 0.5 * KK1[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 0.5 * KK1[gx][gy + 1][gz] ) )
                                + Hisq[2] * ( ( F[gx][gy][gz - 1] + 0.5 * KK1[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 0.5 * KK1[gx][gy][gz + 1] ) )
                                ) - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] );
            }

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];
                
                KK3[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 0.5 * KK2[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 0.5 * KK2[gx + 1][gy][gz] ) )
                                + Hisq[1] * ( ( F[gx][gy - 1][gz] + 0.5 * KK2[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 0.5 * KK2[gx][gy + 1][gz] ) )
                                + Hisq[2] * ( ( F[gx][gy][gz - 1] + 0.5 * KK2[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 0.5 * KK2[gx][gy][gz + 1] ) )
                               )  - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] );
            }

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];

                KK4[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 1.0 * KK3[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 1.0 * KK3[gx + 1][gy][gz] ) )
                                + Hisq[1] * ( ( F[gx][gy - 1][gz] + 1.0 * KK3[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 1.0 * KK3[gx][gy + 1][gz] ) )
                                + Hisq[2] * ( ( F[gx][gy][gz - 1] + 1.0 * KK3[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 1.0 * KK3[gx][gy][gz + 1] ) )
                               )  - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] );
            }

            #pragma omp parallel for private(grid, gx, gy, gz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                FF[gx][gy][gz] = F[gx][gy][gz] + ( KK1[gx][gy][gz] + 2.0 * KK2[gx][gy][gz] + 2.0 * KK3[gx][gy][gz] + KK4[gx][gy][gz] ) / 6.0;
            }
            
            #pragma omp parallel for private(grid, gx, gy, gz)

            for (int i = 0; i < ExFF.size(); i++)  {

                grid = IdxToGrid(ExFF[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                PFdX[gx][gy][gz] = 0.5 * Hi[0] * std::abs( FF[gx + 1][gy][gz] - FF[gx - 1][gy][gz]);
                PFdY[gx][gy][gz] = 0.5 * Hi[1] * std::abs( FF[gx][gy + 1][gz] - FF[gx][gy - 1][gz]);
                PFdZ[gx][gy][gz] = 0.5 * Hi[2] * std::abs( FF[gx][gy][gz + 1] - FF[gx][gy][gz - 1]);
            }

			t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-c-2 CASE 1 RK4) = %lf sec\n", t_1_elapsed); 

            // check Multiple Expanding 
            // TBL = index of FF that FF(TBL) is higher than TolL

            TBL.clear();

            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz, b1, b2, b3, b4, b5, b6)

            for (int i = 0; i < ExFF.size(); i++)
            {
                grid = IdxToGrid(ExFF[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                b1 = std::abs(FF[gx][gy][gz] * std::conj(FF[gx][gy][gz])) >= TolH;
                b2 = PFdX[gx][gy][gz] >= TolHd;
                b3 = PFdY[gx][gy][gz] >= TolHd;
                b4 = PFdZ[gx][gy][gz] >= TolHd;
                b5 = gx > 1 && gy > 1 && gz > 1;
                b6 = gx < BoxShape[0] - 2 && gy < BoxShape[1] - 2 && gz < BoxShape[2] - 2;

                if (  ( b1 || b2 || b3 || b4 ) && b5 && b6 )  {
                    tmpVec.push_back(ExFF[i]);
                }
            }        
            tmpVec.swap(TBL);
            tmpVec.clear(); 

            // TBL & TBL_P set difference
            tmpVec.resize(TBL_P.size() + TBL.size());
            __gnu_parallel::sort (TBL.begin(), TBL.end());
            __gnu_parallel::sort (TBL_P.begin(), TBL_P.end());
            it=std::set_difference( TBL.begin(), TBL.end(), TBL_P.begin(), TBL_P.end(), tmpVec.begin() );
            tmpVec.resize(it - tmpVec.begin()); 
            tmpVec.swap(TBL);
            tmpVec.clear();

            // Combine TBL and TBL_P
            TBL_P.reserve(TBL_P.size() + TBL.size());
            TBL_P.insert(TBL_P.end(), TBL.begin(), TBL.end());

            // Find unique elements
            __gnu_parallel::sort(TBL_P.begin(),TBL_P.end());
            it = std::unique (TBL_P.begin(), TBL_P.end()); 
            TBL_P.resize(std::distance(TBL_P.begin(),it));

			t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-c-3 CASE 1 TBL TBL_P) = %lf sec\n", t_1_elapsed); 
        }

        // CASE 2: Truncating without extrapolation

        if ( !isExtrapolate && !isFullGrid )
        {
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];

                KK1[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( F[gx - 1][gy][gz] - 2.0 * F[gx][gy][gz] + F[gx + 1][gy][gz] )
                                + Hisq[1] * ( F[gx][gy - 1][gz] - 2.0 * F[gx][gy][gz] + F[gx][gy + 1][gz] ) 
                                + Hisq[2] * ( F[gx][gy][gz - 1] - 2.0 * F[gx][gy][gz] + F[gx][gy][gz + 1] ) 
                               )  - I * k2hb * Potential(xx, yy, zz) * F[gx][gy][gz];
            }

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];

                KK2[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 0.5 * KK1[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 0.5 * KK1[gx + 1][gy][gz] ) )
                                + Hisq[1] * ( ( F[gx][gy - 1][gz] + 0.5 * KK1[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 0.5 * KK1[gx][gy + 1][gz] ) )
                                + Hisq[2] * ( ( F[gx][gy][gz - 1] + 0.5 * KK1[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 0.5 * KK1[gx][gy][gz + 1] ) )
                                ) - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] );
            }

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];
                
                KK3[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 0.5 * KK2[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 0.5 * KK2[gx + 1][gy][gz] ) )
                                + Hisq[1] * ( ( F[gx][gy - 1][gz] + 0.5 * KK2[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 0.5 * KK2[gx][gy + 1][gz] ) )
                                + Hisq[2] * ( ( F[gx][gy][gz - 1] + 0.5 * KK2[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 0.5 * KK2[gx][gy][gz + 1] ) )
                               )  - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] );
            }

            #pragma omp parallel for private(grid, gx, gy, gz, xx, yy, zz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];
                xx = Box[0] + gx * H[0];
                yy = Box[2] + gy * H[1];
                zz = Box[4] + gz * H[2];

                KK4[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 1.0 * KK3[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 1.0 * KK3[gx + 1][gy][gz] ) )
                                + Hisq[1] * ( ( F[gx][gy - 1][gz] + 1.0 * KK3[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 1.0 * KK3[gx][gy + 1][gz] ) )
                                + Hisq[2] * ( ( F[gx][gy][gz - 1] + 1.0 * KK3[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 1.0 * KK3[gx][gy][gz + 1] ) )
                               )  - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] );
            }

            #pragma omp parallel for private(grid, gx, gy, gz)

            for (int i = 0; i < TA.size(); i++)  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                FF[gx][gy][gz] = F[gx][gy][gz] + ( KK1[gx][gy][gz] + 2.0 * KK2[gx][gy][gz] + 2.0 * KK3[gx][gy][gz] + KK4[gx][gy][gz] ) / 6.0;

            }

			t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-d-1 CASE 2 RK4) = %lf sec\n", t_1_elapsed); 
        }   
        else if ( !isExtrapolate && isFullGrid )
        {
            // CASE 3: Full grid

            #pragma omp parallel for private(gx, gy, gz, xx, yy, zz)

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        gx = i;
                        gy = j;
                        gz = k;
                		xx = Box[0] + gx * H[0];
                		yy = Box[2] + gy * H[1];
                 			zz = Box[4] + gz * H[2];

                        KK1[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( F[gx - 1][gy][gz] - 2.0 * F[gx][gy][gz] + F[gx + 1][gy][gz] )
                                        + Hisq[1] * ( F[gx][gy - 1][gz] - 2.0 * F[gx][gy][gz] + F[gx][gy + 1][gz] ) 
                                        + Hisq[2] * ( F[gx][gy][gz - 1] - 2.0 * F[gx][gy][gz] + F[gx][gy][gz + 1] ) 
                                    )  - I * k2hb * Potential(xx, yy, zz) * F[gx][gy][gz];
                    }
                }
            }

            #pragma omp parallel for private(gx, gy, gz, xx, yy, zz)

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        gx = i;
                        gy = j;
                        gz = k;
                		xx = Box[0] + gx * H[0];
                		yy = Box[2] + gy * H[1];
                		zz = Box[4] + gz * H[2];

                        KK2[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 0.5 * KK1[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 0.5 * KK1[gx + 1][gy][gz] ) )
                                        + Hisq[1] * ( ( F[gx][gy - 1][gz] + 0.5 * KK1[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 0.5 * KK1[gx][gy + 1][gz] ) )
                                        + Hisq[2] * ( ( F[gx][gy][gz - 1] + 0.5 * KK1[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 0.5 * KK1[gx][gy][gz + 1] ) )
                                    ) - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 0.5 * KK1[gx][gy][gz] );
                    }
                }
            }

            #pragma omp parallel for private(gx, gy, gz, xx, yy, zz)

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        gx = i;
                        gy = j;
                        gz = k;
                		xx = Box[0] + gx * H[0];
                		yy = Box[2] + gy * H[1];
                		zz = Box[4] + gz * H[2];
                
                        KK3[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 0.5 * KK2[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 0.5 * KK2[gx + 1][gy][gz] ) )
                                        + Hisq[1] * ( ( F[gx][gy - 1][gz] + 0.5 * KK2[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 0.5 * KK2[gx][gy + 1][gz] ) )
                                        + Hisq[2] * ( ( F[gx][gy][gz - 1] + 0.5 * KK2[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 0.5 * KK2[gx][gy][gz + 1] ) )
                                    )  - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 0.5 * KK2[gx][gy][gz] );
                    }
                }
            }

            #pragma omp parallel for private(gx, gy, gz, xx, yy, zz)

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        gx = i;
                        gy = j;
                        gz = k;
                		xx = Box[0] + gx * H[0];
                		yy = Box[2] + gy * H[1];
                		zz = Box[4] + gz * H[2];

                        KK4[gx][gy][gz] = I * kh2m * (  Hisq[0] * ( ( F[gx - 1][gy][gz] + 1.0 * KK3[gx - 1][gy][gz] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx + 1][gy][gz] + 1.0 * KK3[gx + 1][gy][gz] ) )
                                        + Hisq[1] * ( ( F[gx][gy - 1][gz] + 1.0 * KK3[gx][gy - 1][gz] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx][gy + 1][gz] + 1.0 * KK3[gx][gy + 1][gz] ) )
                                        + Hisq[2] * ( ( F[gx][gy][gz - 1] + 1.0 * KK3[gx][gy][gz - 1] ) - 2.0 * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] ) + ( F[gx][gy][gz + 1] + 1.0 * KK3[gx][gy][gz + 1] ) )
                                    )  - I * k2hb * Potential(xx, yy, zz) * ( F[gx][gy][gz] + 1.0 * KK3[gx][gy][gz] );
                    }
                }
            }

            #pragma omp parallel for private(gx, gy, gz)

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        gx = i;
                        gy = j;
                        gz = k;

                        FF[gx][gy][gz] = F[gx][gy][gz] + ( KK1[gx][gy][gz] + 2.0 * KK2[gx][gy][gz] + 2.0 * KK3[gx][gy][gz] + KK4[gx][gy][gz] ) / 6.0;
                    }
                }
            }
        }  

        t_1_begin = omp_get_wtime();

        // ff(t+1) Normailzed & go on
        norm = 0.0;

		if (!isFullGrid)  {

        	#pragma omp parallel for private(grid, gx, gy, gz) reduction (+:norm)

        	for (int i = 0; i < TA.size(); i++)  {

        		grid = IdxToGrid(TA[i]);
            	gx = grid[0];
            	gy = grid[1];
            	gz = grid[2];
            	norm += std::abs(FF[gx][gy][gz] * std::conj(FF[gx][gy][gz]));
        	}
		}  else  {

			#pragma omp parallel for private(grid, gx, gy, gz) reduction (+:norm)

        	for (int i = 0; i < BoxShape[0]; i ++)  {

            	for (int j = 0; j < BoxShape[1]; j ++)  {

                	for (int k = 0; k < BoxShape[2]; k++)  {

                		norm += std::abs(FF[i][j][k] * std::conj(FF[i][j][k]));
					}
				}
			}

		}

        norm *= H[0] * H[1] * H[2];
        norm = 1.0 / sqrt(norm);

        #pragma omp parallel for

        for (int i = 0; i < BoxShape[0]; i ++)  {

            for (int j = 0; j < BoxShape[1]; j ++)  {

                for (int k = 0; k < BoxShape[2]; k++)  {

                    FF[i][j][k] = norm * FF[i][j][k];
                    F[i][j][k] = FF[i][j][k];
					PF[i][j][k] = std::abs(F[i][j][k] * std::conj(F[i][j][k]));
                }
            }
        }     

		t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        //log->log("Elapsed time (omp-e-1 FF) = %lf sec\n", t_1_elapsed); 

        // Truncated_New Edge
        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();

            // PFdX

            coeff = 1.0 / (2.0 * H[0]);

            #pragma omp parallel for

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        PFdX[i][j][k] = std::abs(F[i+1][j][k] - F[i-1][j][k]) * coeff;
                    }
                }
            }

            // PFdY    

            coeff = 1.0 / (2.0 * H[1]);

            #pragma omp parallel for

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        PFdY[i][j][k] = std::abs(F[i][j+1][k] - F[i][j-1][k]) * coeff;
                    }
                }
            }  

            // PFdZ    

            coeff = 1.0 / (2.0 * H[2]);

            #pragma omp parallel for

            for (int i = 1; i < BoxShape[0] - 1; i ++)  {

                for (int j = 1; j < BoxShape[1] - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1; k++)  {

                        PFdZ[i][j][k] = std::abs(F[i][j][k+1] - F[i][j][k-1]) * coeff;
                    }
                }
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-2 PF) = %lf sec\n", t_1_elapsed); 

            // Truncate

            t_1_begin = omp_get_wtime();

            #pragma omp parallel for

            for (int i = 1; i < BoxShape[0] - 1 ; i ++)  {

                for (int j = 1; j < BoxShape[1]  - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1 ; k++)  {
            
                        if ( ( PF[i][j][k] < TolH ) && ( PFdX[i][j][k] < TolHd ) && ( PFdY[i][j][k] < TolHd ) && ( PFdZ[i][j][k] < TolHd ) )  {
                            F[i][j][k] = xZERO;
						}
                    }
                }
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-1 PF) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // TA

            #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz, b1, b2)

            for ( int i = 0; i < TA.size(); i++ )  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                b1 = F[gx][gy][gz] != xZERO;
                b2 = PFdX[gx][gy][gz] >= TolHd || PFdY[gx][gy][gz] >= TolHd || PFdZ[gx][gy][gz]>= TolHd;

                if ( b1 || b2 )  {

                    tmpVec.push_back(TA[i]);
                }
            }
            tmpVec.swap(TA);  
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-2 TA) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // TB

            #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz, b1, b2, b3, b4, b5, b6)

            for ( int i = 0; i < TA.size(); i++ )  {

                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

                b1 = F[ gx - 1 ][ gy  ][ gz   ] == xZERO;
                b2 = F[ gx + 1 ][ gy  ][ gz   ] == xZERO;
                b3 = F[ gx   ][ gy -1 ][ gz   ] == xZERO;
                b4 = F[ gx   ][ gy +1 ][ gz   ] == xZERO;
                b5 = F[ gx   ][ gy  ][ gz - 1 ] == xZERO;
                b6 = F[ gx   ][ gy  ][ gz + 1 ] == xZERO;
               
                if ( b1 || b2 || b3 || b4 || b5 || b6 )
                {
                    b1 = ( PFdX[ gx - 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdY[ gx - 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdZ[ gx - 1 ][ gy   ][ gz   ] < TolHd );
                    b2 = ( PFdX[ gx + 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdY[ gx + 1 ][ gy   ][ gz   ] < TolHd ) && ( PFdZ[ gx + 1 ][ gy   ][ gz   ] < TolHd );
                    b3 = ( PFdX[ gx   ][ gy - 1 ][ gz   ] < TolHd ) && ( PFdY[ gx   ][ gy - 1 ][ gz   ] < TolHd ) && ( PFdZ[ gx   ][ gy - 1 ][ gz   ] < TolHd );
                    b4 = ( PFdX[ gx   ][ gy + 1 ][ gz   ] < TolHd ) && ( PFdY[ gx   ][ gy + 1 ][ gz   ] < TolHd ) && ( PFdZ[ gx   ][ gy + 1 ][ gz   ] < TolHd );
                    b5 = ( PFdX[ gx   ][ gy   ][ gz - 1 ] < TolHd ) && ( PFdY[ gx   ][ gy   ][ gz - 1 ] < TolHd ) && ( PFdZ[ gx   ][ gy   ][ gz - 1 ] < TolHd );
                    b6 = ( PFdX[ gx   ][ gy   ][ gz + 1 ] < TolHd ) && ( PFdY[ gx   ][ gy   ][ gz + 1 ] < TolHd ) && ( PFdZ[ gx   ][ gy   ][ gz + 1 ] < TolHd );

                    if (b1 || b2 || b3 || b4 || b5 || b6) {

                        tmpVec.push_back(TA[i]);  
                    }
                }
            }
            tmpVec.swap(TB);
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-3 TB) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // TA expansion
            #pragma omp parallel for reduction(merge: tmpVec) private(grid, gx, gy, gz)                 
            for (int i = 0; i < TA.size(); i++)
            {
                grid = IdxToGrid(TA[i]);
                gx = grid[0];
                gy = grid[1];
                gz = grid[2];

				if (gx + 1 != BoxShape[0] - 1)
                	tmpVec.push_back(GridToIdx(gx + 1, gy  , gz  ));

				if (gx - 1 != 0)
                	tmpVec.push_back(GridToIdx(gx - 1, gy  , gz  ));

				if (gy + 1 != BoxShape[1] - 1)
                	tmpVec.push_back(GridToIdx(gx  , gy + 1, gz  ));

				if (gy - 1 != 0)
                	tmpVec.push_back(GridToIdx(gx  , gy - 1, gz  ));

				if (gz + 1 != BoxShape[2] - 1)
                	tmpVec.push_back(GridToIdx(gx  , gy  , gz + 1));

				if (gz - 1 != 0)
                	tmpVec.push_back(GridToIdx(gx  , gy  , gz - 1));            
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-4-1 push_back) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // Combine TA and tmpVec
            TA.reserve(TA.size() + tmpVec.size());
            TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-4-2 combine) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // Find unique elements
            __gnu_parallel::sort(TA.begin(),TA.end());

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-4-2-1 sort) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            it = std::unique (TA.begin(), TA.end()); 

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-4-2-2 unique) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            TA.resize(std::distance(TA.begin(),it));

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-3-4-2-3 resize) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();
            
			#pragma omp parallel for

            for (int i = 1; i < BoxShape[0] - 1 ; i ++)  {

                for (int j = 1; j < BoxShape[1]  - 1; j ++)  {

                    for (int k = 1; k < BoxShape[2] - 1 ; k++)  {

						PF[i][j][k] = std::abs( F[i][j][k] * std::conj(F[i][j][k]) );

                    }
                }
            }
        }

        if ( (tt + 1) % PERIOD == 0 )
        {
			// Compute transmittance

            t_1_begin = omp_get_wtime();

            pftrans = 0.0;

            #pragma omp parallel for reduction (+:pftrans)

            for (int i = idx_x0; i < BoxShape[0]; i ++)  {

                for (int j = 0; j < BoxShape[1]; j ++)  {

                    for (int k = 0; k < BoxShape[2]; k++)  {

                        pftrans+=PF[i][j][k];
                    }
                }
            }
            pftrans *= H[0] * H[1] * H[2];

            PF_trans.push_back(pftrans);

            log->log("[Scatter3d] Time %lf, Trans = %e\n", ( tt + 1 ) * kk, pftrans);

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            //log->log("Elapsed time (omp-e-2 trans) = %lf sec\n", t_1_elapsed); 

            t_0_end = omp_get_wtime();
            t_0_elapsed = t_0_end - t_0_begin;
 
            log->log("[Scatter3d] Step: %d, Elapsed time: %lf sec\n", tt + 1, t_0_elapsed);

			if (!isFullGrid)  {

            	log->log("[Scatter3d] TA size = %d, TB size = %d\n", TA.size(), TB.size());
            	log->log("[Scatter3d] TA / total grids = %lf\n", ( TA.size() * 1.0 ) / GRIDS_TOT);
			}
        	log->log("\n........................................................\n\n");
        }           

    } // Time iteration 

    log->log("[Scatter3d] Evolve done.\n");
}
/* ------------------------------------------------------------------------------- */

inline std::complex<double> Scatter3d::Wavefunction(double x, double y, double z)
{
    return exp( -A[0] * pow(x - Wave0[0], 2) + I / hb * P[0] * (x - Wave0[0])) * exp( (Ld - 0.5) * (std::log(2 * Ld) - Da * (y - r0)) - Ld * exp(-Da * (y - r0))) * exp( (Ld - 0.5) * (std::log(2 * Ld) - Da * (z - r0)) - Ld * exp(-Da * (z - r0)));
}

/* ------------------------------------------------------------------------------- */

inline double Scatter3d::Potential(double x, double y, double z)
{
    return V0 * pow(cosh(alpha * x), -2.0) + De * pow(1.0 - exp( - Da * ( y - r0 )), 2.0) - De + De * pow(1.0 - exp( - Da * ( z - r0 )), 2.0) - De;
}

/* ------------------------------------------------------------------------------- */

Vector3i Scatter3d::IdxToGrid(int idx)
{
    int z = idx % BoxShape[2];
    int y = (idx % (BoxShape[1] * BoxShape[2]))/ BoxShape[2];
    int x = idx / (BoxShape[1] * BoxShape[2]);

    Vector3i grid;
    grid << x, y, z;

    return grid;
}
/* ------------------------------------------------------------------------------- */

inline int Scatter3d::GridToIdx(int x, int y, int z)
{
    return x * BoxShape[1] * BoxShape[2] + y * BoxShape[2] + z;
}
/* ------------------------------------------------------------------------------- */

inline void Scatter3d::DefineBoundary()
{
    // Find elements of DBi and DBi2

    for (int j = 0; j < BoxShape[1]; j++)  {
        for (int k = 0; k < BoxShape[2]; k++)  {

            DBi.push_back(GridToIdx(0,j,k));
            DBi.push_back(GridToIdx(BoxShape[0]-1,j,k));

            DBi2.push_back(GridToIdx(0,j,k));
            DBi2.push_back(GridToIdx(1,j,k));
            DBi2.push_back(GridToIdx(2,j,k));
            DBi2.push_back(GridToIdx(BoxShape[0]-1,j,k));
            DBi2.push_back(GridToIdx(BoxShape[0]-2,j,k));
            DBi2.push_back(GridToIdx(BoxShape[0]-3,j,k));
        }
    }

    for (int i = 0; i < BoxShape[0]; i++)  {
        for (int k = 0; k < BoxShape[2]; k++)  {

            DBi.push_back(GridToIdx(i,0,k));
            DBi.push_back(GridToIdx(i,BoxShape[1]-1,k));

            DBi2.push_back(GridToIdx(i,0,k));
            DBi2.push_back(GridToIdx(i,1,k));
            DBi2.push_back(GridToIdx(i,2,k));
            DBi2.push_back(GridToIdx(i,BoxShape[1]-1,k));
            DBi2.push_back(GridToIdx(i,BoxShape[1]-2,k));
            DBi2.push_back(GridToIdx(i,BoxShape[1]-3,k));
        }
    }  

    for (int i = 0; i < BoxShape[0]; i++)  {
        for (int j = 0; j < BoxShape[1]; j++)  {

            DBi.push_back(GridToIdx(i,j,0));
            DBi.push_back(GridToIdx(i,j,BoxShape[2]-1));

            DBi2.push_back(GridToIdx(i,j,0));
            DBi2.push_back(GridToIdx(i,j,1));
            DBi2.push_back(GridToIdx(i,j,2));
            DBi2.push_back(GridToIdx(i,j,BoxShape[2]-1));
            DBi2.push_back(GridToIdx(i,j,BoxShape[2]-2));
            DBi2.push_back(GridToIdx(i,j,BoxShape[2]-3));
        }
    }  

    // Find unique elements

    std::vector<int>::iterator it;
    __gnu_parallel::sort(DBi.begin(),DBi.end());
    it = std::unique (DBi.begin(), DBi.end()); 
    DBi.resize(std::distance(DBi.begin(),it)); 

    __gnu_parallel::sort(DBi2.begin(),DBi2.end());
    it = std::unique (DBi2.begin(), DBi2.end()); 
    __gnu_parallel::sort(DBi2.begin(),DBi2.end());
    DBi2.resize(std::distance(DBi2.begin(),it));
}
/* ------------------------------------------------------------------------------- */

