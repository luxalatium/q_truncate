// ==============================================================================
//
//  Scatter2d.cpp
//  QTR
//
//  Created by Albert Lu on 10/27/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 1/22/19
//
//  Note:
//
// ==============================================================================

#define BOOST_DISABLE_ASSERTS

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <omp.h>
#include <parallel/algorithm>
#include <vector>

#include "Constants.h"
#include "Containers.h"
#include "Error.h"
#include "Log.h"
#include "Parameters.h"
#include "Scatter2d.h"

#include "boost/multi_array.hpp"

using namespace QTR_NS;
using std::vector;
using std::cout;
using std::endl;

/* ------------------------------------------------------------------------------- */

// DEFINE POTENTIAL TYPE

#define POTMIN (-1000.0)

#if defined POT_ECKMO

#define WAVEFUNCTION(x1,x2) Wavefunction_EckMO(x1,x2)
#define POTENTIAL(x1,x2) Potential_EckMO(x1,x2)
#define POTNAME "EckMO"

#elif defined POT_ECKHO

#define WAVEFUNCTION(x1,x2) Wavefunction_EckHO(x1,x2)
#define POTENTIAL(x1,x2) Potential_EckHO(x1,x2)
#define POTNAME "EckHO"

#elif defined POT_GAUMO

#define WAVEFUNCTION(x1,x2) Wavefunction_GauMO(x1,x2)
#define POTENTIAL(x1,x2) Potential_GauMO(x1,x2)
#define POTNAME "GauMO"

#elif defined POT_GAUHO

#define WAVEFUNCTION(x1,x2) Wavefunction_GauHO(x1,x2)
#define POTENTIAL(x1,x2) Potential_GauHO(x1,x2)
#define POTNAME "GauHO"

#else

#define WAVEFUNCTION(x1,x2) Wavefunction_HH(x1,x2)
#define POTENTIAL(x1,x2) Potential_HH(x1,x2)
#define POTNAME "HH"

#endif
/* ------------------------------------------------------------------------------- */

Scatter2d::Scatter2d(class QTR *q)
{
    qtr = q;
    err = qtr->error;
    log = qtr->log;
    parameters = qtr->parameters;
    init();
} 
/* ------------------------------------------------------------------------------- */

Scatter2d::~Scatter2d()
{     
    return;
}
/* ------------------------------------------------------------------------------- */

void Scatter2d::init()
{
    log->log("\n\n[Scatter2d] INIT starts ...\n");
    log->log("\n\n[Scatter2d] Potential type: %s\n", POTNAME);

    // General parameters
    I = {0,1}; // sqrt(-1)
    PI_INV = 1.0 / PI; // 1/PI
    xZERO = {0,0}; // complex zero
    DIMENSIONS = parameters->scxd_dimensions;
    PERIOD = parameters->scxd_period;
    SORT_PERIOD = parameters->scxd_sortperiod;
    PRINT_PERIOD = parameters->scxd_printperiod;
    TIME = parameters->scxd_Tf;
    QUIET = parameters->quiet;
    TIMING = parameters->timing;
    isTrans = parameters->scxd_isTrans;
    isAcf = parameters->scxd_isAcf;
    isPrintEdge = parameters->scxd_isPrintEdge;
    isPrintDensity = parameters->scxd_isPrintDensity;

    // Grid size
    H.resize(DIMENSIONS);
    Hi.resize(DIMENSIONS);
    Hisq.resize(DIMENSIONS);
    S.resize(DIMENSIONS);  
    kk = parameters->scxd_k;
    H[0] = parameters->scxd_h1;
    H[1] = parameters->scxd_h2;

    for (unsigned int i = 0; i < DIMENSIONS; i ++)  {

        Hi[i] = 1 / H[i];
        Hisq[i] = 1 / pow(H[i],2);
        S[i] = kk * Hisq[i];
    }

    // Domain size and # grids
    Box.resize(DIMENSIONS * 2);
    Box[0] = parameters->scxd_xi1;
    Box[1] = parameters->scxd_xf1;
    Box[2] = parameters->scxd_xi2;
    Box[3] = parameters->scxd_xf2; 
    BoxShape.resize(DIMENSIONS);
    GRIDS_TOT = 1;
    log->log("[Scatter2d] Number of grids = (");

    for (unsigned int i = 0; i < DIMENSIONS; i ++)  {

        BoxShape[i] = (int)( (Box[2 * i + 1] - Box[2 * i]) / H[i] ) + 1;
        GRIDS_TOT *= BoxShape[i];

        if ( i < DIMENSIONS - 1 )
            log->log("%d, ", BoxShape[i]);
        else
            log->log("%d)\n", BoxShape[i]);
    }

    idx_x0 = (int) ( ( 0 - Box[0] ) / H[0] );
    M1 = BoxShape[1];
    log->log("[Scatter2d] x = 0 at idx = %d\n", idx_x0);

    // Parameters
    hb = parameters->scxd_hb;
    m  = parameters->scxd_m;

    // Potential: HO specific
    w = parameters->scxd_w;
    aa = 0.5 * m * w / hb;
    log->log("[Scatter2d] aa = %lf\n", aa);

    // Potential: Eckart
    V0 = parameters->scxd_V0;
    alpha = parameters->scxd_alpha;
    ek2v = parameters->scxd_ek2v;
    Ek0 = V0 * ek2v;
    log->log("[Scatter2d] Ek0 = %lf\n", Ek0);

    // Potential: Morse
    De = parameters->scxd_De;
    Da = parameters->scxd_Da;
    r0 = parameters->scxd_r0;
    Ld = sqrt(2 * m * De) / (Da * hb);  
    log->log("[Scatter2d] Ld = %lf\n", Ld);

    // Potential: Henon-Heiles
    lambda = parameters->scxd_lambda;

    // Potential: HO
    sigma = parameters->scxd_sigma;
    k0 = parameters->scxd_k0;

    // Potential: Gaussian
    beta = parameters->scxd_beta;

    // Wavefunction parameters
    Wave0.resize(DIMENSIONS);
    Wave0[0] = parameters->scxd_x01;
    Wave0[1] = parameters->scxd_x02;

    A.resize(DIMENSIONS);
    A[0] = parameters->scxd_a1;
    A[1] = parameters->scxd_a2;

    // Ｏverride for HO
    A[1] = 0.5 * sqrt(m * k0) / hb;
    log->log("[Scatter2d] A[1] = %lf\n", A[1]);

    P.resize(DIMENSIONS);
    P[0] = parameters->scxd_p1;
    P[1] = parameters->scxd_p2;
    // Ｏverride
    P[0] = sqrt(2.0 * m * Ek0);
    P[1] = 0;
    log->log("[Scatter2d] P[0] = %lf\n", P[0]);

    // Truncate parameters
    isFullGrid = parameters->scxd_isFullGrid;
    TolH = parameters->scxd_TolH;    // Tolerance of probability density for Zero point Cutoff
    TolL = parameters->scxd_TolL;    // Tolerance of probability density for Edge point
    TolHd = parameters->scxd_TolHd;  // Tolerance of probability first diff for Zero point Cutoff
    TolLd = parameters->scxd_TolLd;  // Tolerance of probability density for Edge point
    ExReduce = parameters->scxd_ExReduce; //Extrapolation reduce factor

    // Spectrum
    dk = parameters->scxd_dk;
    kMax = parameters->scxd_kmax;

    log->log("[Scatter2d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void Scatter2d::Evolve()
{
    #pragma omp declare reduction (merge : MeshIndex : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    log->log("[Scatter2d] Evolve starts ...\n");

    // Files
    FILE *pfile;

    // Variables 
    int index;
    int n1, n2;
    bool b1, b2, b3, b4, b5;
    double kh2m = kk * hb / 2 / m;
    double k2hb = kk / hb;
    double coeff1 = 0.0;
    double coeff2 = 0.0;
    double pftrans;
    double t_0_begin, t_0_end;
    double t_1_begin, t_1_end;
    double t_0_elapsed = 0.0;
    double t_1_elapsed = 0.0;
    std::complex<double> num_cpx; // complex number holder

    // Core computation time (RK4, normalization, initialization, etc)
    double t_full = 0.0;
    double t_truncate = 0.0;

    // Overhead time (truncate)
    double t_overhead = 0.0;

    // temporary index container
    MeshIndex tmpVec; 

    // Boundary layer container for extrapolation loop
    MeshIndex ExBD;     
     
    // normalization factor
    double norm;
    double norm_re;
    double norm_im;

    // potential value
    double pot_val;
 
    //  2d Grid vector and indices
    VectorXi grid;
    int g1, g2;
    double xx1, xx2;

    // Vector iterater
    vector<int>::iterator it;

    // Extrapolation 
    int count;
    vector<std::complex<double>> ExTBL;
    std::complex<double> val, val_min;
    std::complex<double> sum;
    double val_min_abs;
    bool isFirstExtrp;

    // Neighborlist
    int nneigh = 0;
    vector<vector<int>> neighlist;
    vector<int> neighs(DIMENSIONS);

    // PF_trans
    vector<double> PF_trans;
    PF_trans.push_back(0.0);

    // Auto-correlation
    vector<vector<double>> ACFunc;
    vector<double> acf(2);
    acf = {1.0, 0.0};
    ACFunc.push_back(acf);

    // Spectrum
    vector<double> Spectrum;

    // Make DBi and DBi2
    if ( isFullGrid )  
    {
      t_1_begin = omp_get_wtime();
      DefineBoundary();
      t_1_end = omp_get_wtime();
      t_1_elapsed = t_1_end - t_1_begin;
      t_full += t_1_elapsed;
      if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (define boundary) = %lf sec\n\n", t_1_elapsed);
    }

    log->log("[Scatter2d] Initializing containers ...\n");

    // Initialize containers

    t_0_begin = omp_get_wtime();

    MeshCX2D F(boost::extents[BoxShape[0]][BoxShape[1]]);
    F.reindex(0);
    MeshCX2D F0_STAR(boost::extents[BoxShape[0]][BoxShape[1]]);
    F0_STAR.reindex(0);
    MeshD2D PF(boost::extents[BoxShape[0]][BoxShape[1]]);
    PF.reindex(0);
    MeshD2D POT(boost::extents[BoxShape[0]][BoxShape[1]]);
    POT.reindex(0);
    MeshCX2D FF(boost::extents[BoxShape[0]][BoxShape[1]]);
    FF.reindex(0);   
    MeshCX2D KK1(boost::extents[BoxShape[0]][BoxShape[1]]);
    KK1.reindex(0);
    MeshCX2D KK2(boost::extents[BoxShape[0]][BoxShape[1]]);
    KK2.reindex(0);  
    MeshCX2D KK3(boost::extents[BoxShape[0]][BoxShape[1]]);
    KK3.reindex(0); 
    MeshCX2D KK4(boost::extents[BoxShape[0]][BoxShape[1]]);
    KK4.reindex(0);
    MeshD2D PFdX1(boost::extents[BoxShape[0]][BoxShape[1]]);
    PFdX1.reindex(0);  
    MeshD2D PFdX2(boost::extents[BoxShape[0]][BoxShape[1]]);
    PFdX2.reindex(0);
    Mask2D TAMask(boost::extents[BoxShape[0]][BoxShape[1]]);
    TAMask.reindex(0);

    #pragma omp parallel for
    for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
        for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

            F[i1][i2] = xZERO;
            PF[i1][i2] = 0.0;
            POT[i1][i2] = POTMIN - 1.0;
            FF[i1][i2] = xZERO;
            KK1[i1][i2] = xZERO;
            KK2[i1][i2] = xZERO;
            KK3[i1][i2] = xZERO;
            KK4[i1][i2] = xZERO;
        }
    }

    if ( isAcf )  {

        #pragma omp parallel for
        for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                F0_STAR[i1][i2] = xZERO;
            }
        }
    }

    if ( !isFullGrid )  {

        t_1_begin = omp_get_wtime();
        
        #pragma omp parallel for
        for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                PFdX1[i1][i2] = 0.0;
                PFdX2[i1][i2] = 0.0;
                TAMask[i1][i2] = 0;
            }
        }

        nneigh = 0;

        for (int d = 1; d <= 4; d ++)  {
            for (n1 = -d; n1 <= d; n1 ++)  {

                n2 = d - abs(n1);

                if (n2 != 0)  {
                    neighs = {n1,n2};
                    neighlist.push_back(neighs);
                    neighs = {n1,-n2};
                    neighlist.push_back(neighs);
                    nneigh += 2;
                }
                else  {
                    neighs = {n1,0};
                    neighlist.push_back(neighs);
                    nneigh += 1;
                }
            }
        }
        log->log("[Scatter2d] nneigh = %d\n", nneigh); 
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_overhead += t_1_elapsed;
    }
    t_0_end = omp_get_wtime();
    t_0_elapsed = t_0_end - t_0_begin;
    t_full += t_0_elapsed;

    if ( !isFullGrid )
        t_truncate += t_0_elapsed - t_1_elapsed; // subtract overhead

    if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing containers) = %lf sec\n\n", t_0_elapsed); 

    // .........................................................................................

    log->log("[Scatter2d] Initializing wavefunction ...\n");  

    t_1_begin = omp_get_wtime();

    // Initialize wavefunction
    #pragma omp parallel for
    for (unsigned int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {
        for (unsigned int i2 = 1; i2 < BoxShape[1] - 1 ; i2 ++)  {
       
            F[i1][i2] = WAVEFUNCTION(Box[0] + i1 * H[0], Box[2] + i2 * H[1]);
        }
    }

    // Normalization
    norm = 0.0;

    #pragma omp parallel for reduction (+:norm)
    for (unsigned int i1 = 0; i1 <  BoxShape[0]; i1 ++)  {
        for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                
            norm += std::abs(F[i1][i2] * std::conj(F[i1][i2]));
        }
    }
    norm *= H[0] * H[1];
    log->log("[Scatter2d] Normalization factor = %e\n",norm);
    norm = 1.0 / sqrt(norm);

    #pragma omp parallel for
    for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
        for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

            F[i1][i2] = norm * F[i1][i2];
            PF[i1][i2] = std::abs(F[i1][i2] * std::conj(F[i1][i2]));
        }
    }

    if ( isAcf )  {

        #pragma omp parallel for
        for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                F0_STAR[i1][i2] = std::conj(F[i1][i2]);
            }
        }  
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    t_full += t_1_elapsed;
    t_truncate += t_1_elapsed;
    if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (initializing wavefunction) = %lf sec\n\n", t_1_elapsed); 

    // .........................................................................................

    // PF 1st-order finite difference 

    if ( !isFullGrid )
    {
        t_1_begin = omp_get_wtime();

        log->log("[Scatter2d] Computing finite differences ...\n");   

        coeff1 = 1.0 / (2.0 * H[0]);
        coeff2 = 1.0 / (2.0 * H[1]);

        #pragma omp parallel for
        for (unsigned int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
            for (unsigned int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                PFdX1[i1][i2] = std::abs(F[i1+1][i2] - F[i1-1][i2]) * coeff1;
                PFdX2[i1][i2] = std::abs(F[i1][i2+1] - F[i1][i2-1]) * coeff2;
            }
        }    
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_overhead += t_1_elapsed;
        if (!QUIET && TIMING) log->log("[Scatter2d] Elapsed time (finite differences) = %lf sec\n\n", t_1_elapsed); 
    }
    // .........................................................................................

    // Initial truncation & edge point check

    if ( !isFullGrid )
    {
        t_1_begin = omp_get_wtime();

        log->log("[Scatter2d] Initial truncation ...\n");

        // Truncation
        #pragma omp parallel for private(b1, b2, b3)
        for (unsigned int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
            for (unsigned int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                b1 = PF[i1][i2] < TolH;
                b2 = PFdX1[i1][i2] < TolHd;
                b3 = PFdX2[i1][i2] < TolHd;
        
                if ( b1 && b2 && b3 )
                {
                    F[i1][i2] = xZERO;
                } 
                else if (PF[i1][i2] <  TolH)
                {
                    log->log("%d %d %d %d %d %e %e %e %d \n", i1,i2,b1,b2,b3, PF[i1][i2], PFdX1[i1][i2], PFdX2[i1][i2], F[i1][i2] != xZERO );

                }
           }
        }   
        // TA 
        #pragma omp parallel for reduction(merge: tmpVec) private(b1, b2, b3)
        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
            for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                b1 = F[i1][i2] != xZERO;
                b2 = PFdX1[i1][i2] >= TolHd;
                b3 = PFdX2[i1][i2] >= TolHd;
            
                if ( b1 || ( b2 || b3 ))  {

                    log->log("%d %d %d %d %d %e %e %e\n",i1,i2, b1,b2,b3,PF[i1][i2], PFdX1[i1][i2], PFdX2[i1][i2] );
                    tmpVec.push_back(GridToIdx(i1,i2));   
                    TAMask[i1][i2] = 1;                                
                }
            }
        }
        tmpVec.swap(TA);
        tmpVec.clear();
        log->log("[Scatter2d] TA size = %d\n", TA.size());

        // `````````````````````````````````````````````````````````````````

        // TB
        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, b1, b2, b3, b4)
        for (int i = 0; i < TA.size(); i++)
        {
            g2 = TA[i] % M1;
            g1 = TA[i] / M1;

            b1 = F[g1-1][g2] == xZERO;
            b2 = F[g1+1][g2] == xZERO;
            b3 = F[g1][g2-1] == xZERO;
            b4 = F[g1][g2+1] == xZERO;
               
            if ( b1 || b2 || b3 || b4 )
            {
                b1 = ( PFdX1[g1-1][g2] < TolHd ) && ( PFdX2[g1-1][g2] < TolHd ) ;
                b2 = ( PFdX1[g1+1][g2] < TolHd ) && ( PFdX2[g1+1][g2] < TolHd ) ;
                b3 = ( PFdX1[g1][g2-1] < TolHd ) && ( PFdX2[g1][g2-1] < TolHd ) ;
                b4 = ( PFdX1[g1][g2+1] < TolHd ) && ( PFdX2[g1][g2+1] < TolHd ) ;
 
                if ( b1 || b2 || b3 || b4 ) {

                    tmpVec.push_back(TA[i]);  
                }
            }       
        }
        tmpVec.swap(TB);
        tmpVec.clear();
        log->log("[Scatter2d] TB size = %d\n", TB.size());

        // `````````````````````````````````````````````````````````````````

        // TA expansion
        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2)
        for (int i = 0; i < TA.size(); i++)
        {
            g2 = TA[i] % M1;
            g1 = TA[i] / M1;

            if ( g1 + 1 != BoxShape[0] - 1 && TAMask[g1+1][g2] == 0 )  {
                tmpVec.push_back(GridToIdx(g1+1,g2));
            }
            if ( g1 - 1 != 0 && TAMask[g1-1][g2] == 0 )  {
                tmpVec.push_back(GridToIdx(g1-1,g2));
            }
            if ( g2 + 1 != BoxShape[1] - 1 && TAMask[g1][g2+1] == 0 )  {
                tmpVec.push_back(GridToIdx(g1,g2+1));
            }
            if ( g2 - 1 != 0 && TAMask[g1][g2-1] == 0 )  {
                tmpVec.push_back(GridToIdx(g1,g2-1));
            }          
        }
        // Find unique elements
        __gnu_parallel::sort(tmpVec.begin(),tmpVec.end());
        it = std::unique (tmpVec.begin(), tmpVec.end()); 
        tmpVec.resize(std::distance(tmpVec.begin(),it)); 

        // Combine TA and tmpVec
        TA.reserve(TA.size() + tmpVec.size());
        TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());

        // Update TA Mask
        #pragma omp parallel for private(g1, g2)
        for (int i = 0; i < tmpVec.size(); i++)
        {
            g2 = tmpVec[i] % M1;
            g1 = tmpVec[i] / M1;
            TAMask[g1][g2] = 1;
        }      
        tmpVec.clear();

        // Sort TA
        __gnu_parallel::sort (TA.begin(), TA.end());

        log->log("[Scatter2d] TA size = %d, TB size = %d\n", TA.size(), TB.size());

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_overhead += t_1_elapsed;

        if (!QUIET && TIMING)  {
            log->log("[Scatter2d] Elapsed time (initial truncation) = %lf sec\n\n", t_1_elapsed);
            log->log("[Scatter2d] Initialization core computation time: %lf sec\n", t_truncate);            
            log->log("[Scatter2d] Initialization overhead: %lf sec\n", t_overhead); 
        }
    }
    else  // Full grid approach
    {
        TB = DBi;
        log->log("[Scatter2d] DBi = %d DBi2 = %d\n\n", DBi.size(), DBi2.size());
        log->log("[Scatter2d] Initialization core computation time: %lf sec\n", t_full);   
    }
    // .........................................................................................

    // Time iteration 

    log->log("=======================================================\n\n"); 
    log->log("[Scatter2d] Time interation starts ...\n"); 
    log->log("[Scatter2d] Number of steps = %d\n\n", (int)(TIME / kk)); 
    log->log("=======================================================\n\n"); 

    for (int tt = 0; tt < (int)(TIME / kk); tt ++)
    {
        t_0_begin = omp_get_wtime(); 

        if ( tt % PRINT_PERIOD == 0 )
        {
            if ( isPrintEdge  && !isFullGrid )  {

                pfile = fopen ("edge.dat","a");
                fprintf(pfile, "%d %lf %lu\n", tt, tt * kk, TB.size());

                for (int i = 0; i < TB.size(); i++)
                {
                    g2 = TB[i] % M1;
                    g1 = TB[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    fprintf(pfile, "%d %d %lf %lf\n", g1, g2, xx1, xx2);          
                }
                fclose(pfile);
            }
            if ( isPrintDensity && !isFullGrid )  {

                pfile = fopen ("density.dat","a");
                fprintf(pfile, "%d %lf %lu\n", tt, tt * kk, TA.size());

                for (int i = 0; i < TA.size(); i++)
                {
                    g2 = TA[i] % M1;
                    g1 = TA[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    fprintf(pfile, "%d %d %lf %lf %.16e\n", g1, g2, xx1, xx2, PF[g1][g2]);          
                }
                fclose(pfile);
            }
            if ( isPrintDensity && isFullGrid )  {

                pfile = fopen ("density.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, BoxShape[0] * BoxShape[1] );

                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                        g1 = i1;
                        g2 = i2;
                        xx1 = Box[0] + g1 * H[0];
                        xx2 = Box[2] + g2 * H[1];
                        fprintf(pfile, "%d %d %lf %lf %.16e\n", g1, g2, xx1, xx2, PF[g1][g2]);
                    }
                }
                fclose(pfile);      
            }
        }

        // Check if TB of f is higher than TolL
        
        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();
            t_truncate = 0.0;
            t_overhead = 0.0;

            TBL.clear();

            // TBL an TBL_P
            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, b1, b2, b3, b4, b5)
            for (int i = 0; i < TB.size(); i++)
            {
                g2 = TB[i] % M1;
                g1 = TB[i] / M1;

                b1 = PF[g1][g2] >= TolL;
                b2 = PFdX1[g1][g2] >= TolLd;
                b3 = PFdX2[g1][g2] >= TolLd;

                // Not in DBi2
                b4 = g1 > 2 && g2 > 2;
                b5 = g1 < BoxShape[0] - 3 && g2 < BoxShape[1] - 3;

                if ( (b1 || b2 || b3 ) && b4 && b5 )  {
                    tmpVec.push_back(TB[i]);
                }
            }
            tmpVec.swap(TBL);
            tmpVec.clear();
            TBL_P = TBL;

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-a-2: TBL, TBL_P) = %lf sec\n", t_1_elapsed);   
            //if (!QUIET) log->log("TBL size = %d TBL_P size = %d\n", TBL.size(), TBL_P.size());
        }
        else  
        {
            t_full = 0.0;
        }
        isExtrapolate = false;
        isFirstExtrp = true;
        // .........................................................................................

        // CASE 1: Truncating with extrapolation

        while ( TBL.size() != 0 && !isFullGrid )
        {
            isExtrapolate = true;

            // Extrapolation
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
                g2 = TBL[index] % M1;
                g1 = TBL[index] / M1;

                if ( F[g1-1][g2] == xZERO )  {
                    ExFF.push_back(GridToIdx(g1-1,g2));
                }
                if ( F[g1+1][g2] == xZERO )  {
                    ExFF.push_back(GridToIdx(g1+1,g2));
                }
                if ( F[g1][g2-1] == xZERO )  {
                    ExFF.push_back(GridToIdx(g1,g2-1));
                }
                if ( F[g1][g2+1] == xZERO )  {
                    ExFF.push_back(GridToIdx(g1,g2+1));
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
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-b-1: ExFF) = %lf sec\n", t_1_elapsed);   

            // .....................................................................

            // Find the direction of Outer to Edge points

            t_1_begin = omp_get_wtime();
            ExTBL.clear();
            it = ExFF.begin();
            
            while ( it != ExFF.end() )  {

                index = std::distance( ExFF.begin(), it );
                g2 = ExFF[index] % M1;
                g1 = ExFF[index] / M1;
                sum = xZERO;
                count = 0;
                isEmpty = true;
                val_min_abs = 100000000;

                if ( F[g1 - 1][g2] != xZERO )  {

                    if ( std::abs(F[g1 - 1][g2]) < val_min_abs &&  F[g1 - 2][g2] != xZERO )  {
                        val_min_abs = std::abs(F[g1 - 1][g2]);
                        val_min = F[g1 - 1][g2];
                    }
                    if ( F[g1 - 2][g2] != xZERO )  {

                        val = exp( 2.0 * std::log(F[g1 - 1][g2]) - std::log(F[g1 - 2][g2]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }
                        isEmpty = false;
                    }
                }

                if ( F[g1 + 1][g2] != xZERO )  {

                    if ( std::abs(F[g1 + 1][g2]) < val_min_abs && F[g1 + 2][g2] != xZERO  )  {
                        val_min_abs = std::abs(F[g1 + 1][g2]);
                        val_min = F[g1 + 1][g2];
                    }
                    if ( F[g1 + 2][g2] != xZERO )  {

                        val = exp( 2.0 * std::log(F[g1 + 1][g2]) - std::log(F[g1 + 2][g2]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }
                        isEmpty = false;
                    }
                }

                if ( F[g1][g2 - 1] != xZERO )  {

                    if ( std::abs(F[g1][g2 - 1]) < val_min_abs && F[g1][g2 - 2] != xZERO  )  {

                        val_min_abs = std::abs(F[g1][g2 - 1]);
                        val_min = F[g1][g2 - 1];
                    }
                    if ( F[g1][g2 - 2] != xZERO )  {

                        val = exp( 2.0 * std::log(F[g1][g2 - 1]) - std::log(F[g1][g2 - 2]) );

                        if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                        {
                            sum += val;
                            count += 1;
                        }
                        isEmpty = false;
                    }
                }

                if ( F[g1][g2 + 1] != xZERO )  {

                    if ( std::abs(F[g1][g2 + 1]) < val_min_abs && F[g1][g2 + 2] != xZERO )  {

                        val_min_abs = std::abs(F[g1][g2 + 1]);
                        val_min = F[g1][g2 + 1];
                    }
                    if ( F[g1][g2 + 2] != xZERO )  {

                        val = exp( 2.0 * std::log(F[g1][g2 + 1]) - std::log(F[g1][g2 + 2]) );

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

            #pragma omp parallel for private(g1, g2)
            for ( int i = 0; i < ExFF.size(); i++ )  {
                g2 = ExFF[i] % M1;
                g1 = ExFF[i] / M1;
                F[g1][g2] = ExTBL[i];
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-b-2: ExFF) = %lf sec\n", t_1_elapsed);  

            // ............................................................................................. Extrapolation

            if ( isFirstExtrp )  {

                // Check Extending nonzero Area

                t_1_begin = omp_get_wtime();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2)
                for (int i = 0; i < ExFF.size(); i++)
                {
                    g2 = ExFF[i] % M1;
                    g1 = ExFF[i] / M1;

                    if (TAMask[g1][g2] == 0)
                        tmpVec.push_back(GridToIdx(g1,g2));
                    if (TAMask[g1+1][g2] == 0)
                        tmpVec.push_back(GridToIdx(g1+1,g2));
                    if (TAMask[g1-1][g2] == 0)
                        tmpVec.push_back(GridToIdx(g1-1,g2));
                    if (TAMask[g1][g2+1] == 0)
                        tmpVec.push_back(GridToIdx(g1,g2+1));
                    if (TAMask[g1][g2-1] == 0)
                        tmpVec.push_back(GridToIdx(g1,g2-1));        
                }

                // Find unique elements
                __gnu_parallel::sort(tmpVec.begin(),tmpVec.end());
                it = std::unique (tmpVec.begin(), tmpVec.end()); 
                tmpVec.resize(std::distance(tmpVec.begin(),it));

                // Combine TA and tmpVec
                TA.reserve(TA.size() + tmpVec.size());
                TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());

                // Update TA Mask
                #pragma omp parallel for private(g1, g2)
                for (int i = 0; i < tmpVec.size(); i++)
                {
                    g2 = tmpVec[i] % M1;
                    g1 = tmpVec[i] / M1;
                    TAMask[g1][g2] = 1;
                }      
                tmpVec.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-c-1: CASE 1 TA) = %lf sec\n", t_1_elapsed); 
        
                // Runge–Kutta 4
                t_1_begin = omp_get_wtime();

                // RK4-1
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < TA.size(); i++)  {

                    g2 = TA[i] % M1;
                    g1 = TA[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK1[g1][g2] = I * kh2m * (  Hisq[0] * ( F[g1 - 1][g2] - 2.0 * F[g1][g2] + F[g1 + 1][g2] )
                                    + Hisq[1] * ( F[g1][g2 - 1] - 2.0 * F[g1][g2] + F[g1][g2 + 1] ) 
                                )  - I * k2hb * pot_val * F[g1][g2];
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_truncate += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-1: CASE 1 KK1) = %lf sec\n", t_1_elapsed);

                // RK4-2
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < TA.size(); i++)  {

                    g2 = TA[i] % M1;
                    g1 = TA[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK2[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK1[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK1[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK1[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK1[g1][g2 + 1] ) )
                                    ) - I * k2hb * pot_val * ( F[g1][g2] + 0.5 * KK1[g1][g2] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_truncate += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-2: CASE 1 KK2) = %lf sec\n", t_1_elapsed);

                // RK4-3
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < TA.size(); i++)  {

                    g2 = TA[i] % M1;
                    g1 = TA[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK3[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK2[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK2[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK2[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK2[g1][g2 + 1] ) )
                                    ) - I * k2hb * pot_val * ( F[g1][g2] + 0.5 * KK2[g1][g2] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_truncate += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-3: CASE 1 KK3) = %lf sec\n", t_1_elapsed);

                // RK4-4
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < TA.size(); i++)  {

                    g2 = TA[i] % M1;
                    g1 = TA[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK4[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 1.0 * KK3[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1 + 1][g2] + 1.0 * KK3[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 1.0 * KK3[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1][g2 + 1] + 1.0 * KK3[g1][g2 + 1] ) )
                                    ) - I * k2hb * pot_val * ( F[g1][g2] + 1.0 * KK3[g1][g2] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_truncate += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-4: CASE 1 KK4) = %lf sec\n", t_1_elapsed);

                // RK4-5
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2)
                for (int i = 0; i < TA.size(); i++)  {

                    g2 = TA[i] % M1;
                    g1 = TA[i] / M1;
                    FF[g1][g2] = F[g1][g2] + ( KK1[g1][g2] + 2.0 * KK2[g1][g2] + 2.0 * KK3[g1][g2] + KK4[g1][g2] ) / 6.0;
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_truncate += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-5: CASE 1 FF) = %lf sec\n", t_1_elapsed);

                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2)
                for (int  i = 0; i < ExFF.size(); i++)  {

                    g2 = ExFF[i] % M1;
                    g1 = ExFF[i] / M1;
                    PFdX1[g1][g2] = 0.5 * Hi[0] * std::abs( FF[g1+1][g2] - FF[g1-1][g2] );
                    PFdX2[g1][g2] = 0.5 * Hi[1] * std::abs( FF[g1][g2+1] - FF[g1][g2-1] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-6: CASE 1 PF) = %lf sec\n", t_1_elapsed);
            }
            else  
            {
                // Extrapolation loop when multiple expanding occured

                t_1_begin = omp_get_wtime();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2)
                for (int i = 0; i < ExFF.size(); i++)
                {
                    g2 = ExFF[i] % M1;
                    g1 = ExFF[i] / M1;

                    if (TAMask[g1][g2] == 0)
                        tmpVec.push_back(GridToIdx(g1,g2));
                    if (TAMask[g1+1][g2] == 0)
                        tmpVec.push_back(GridToIdx(g1+1,g2));
                    if (TAMask[g1-1][g2] == 0)
                        tmpVec.push_back(GridToIdx(g1-1,g2));
                    if (TAMask[g1][g2+1] == 0)
                        tmpVec.push_back(GridToIdx(g1,g2+1));
                    if (TAMask[g1][g2-1] == 0)
                        tmpVec.push_back(GridToIdx(g1,g2-1));        
                }

                // Find unique elements tmpVec
                __gnu_parallel::sort(tmpVec.begin(),tmpVec.end());
                it = std::unique (tmpVec.begin(), tmpVec.end()); 
                tmpVec.resize(std::distance(tmpVec.begin(),it));

                // Combine TA and tmpVec
                TA.reserve(TA.size() + tmpVec.size());
                TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());

                // Update TA Mask
                #pragma omp parallel for private(g1, g2)
                for (int i = 0; i < tmpVec.size(); i++)
                {
                    g2 = tmpVec[i] % M1;
                    g1 = tmpVec[i] / M1;
                    TAMask[g1][g2] = 1;
                }      
                tmpVec.clear();

                #pragma omp parallel for reduction(merge: ExBD) private(g1, g2, n1, n2)
                for (int i = 0; i < ExFF.size(); i++)
                {
                    g2 = ExFF[i] % M1;
                    g1 = ExFF[i] / M1;

                    ExBD.push_back(ExFF[i]);

                    for (int j = 0; j < nneigh; j ++)  {

                        n1 = neighlist[j][0];
                        n2 = neighlist[j][1];

                        if (TAMask[g1+n1][g2+n2] == 1)
                            ExBD.push_back(GridToIdx(g1+n1,g2+n2));
                    }
                }

                // Find unique elements (ExBD)
                __gnu_parallel::sort(ExBD.begin(),ExBD.end());
                it = std::unique (ExBD.begin(), ExBD.end()); 
                ExBD.resize(std::distance(ExBD.begin(),it));

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET) log->log("ExBD size = %d\n", ExBD.size()); 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-cx-1: CASE 1 ExBD) = %lf sec\n", t_1_elapsed); 

                // Runge–Kutta 4
                t_1_begin = omp_get_wtime();

                // RK4-1
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g2 = ExBD[i] % M1;
                    g1 = ExBD[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK1[g1][g2] = I * kh2m * (  Hisq[0] * ( F[g1 - 1][g2] - 2.0 * F[g1][g2] + F[g1 + 1][g2] )
                                    + Hisq[1] * ( F[g1][g2 - 1] - 2.0 * F[g1][g2] + F[g1][g2 + 1] ) 
                                )  - I * k2hb * pot_val * F[g1][g2];
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-1: CASE 1 KK1) = %lf sec\n", t_1_elapsed);

                // RK4-2
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g2 = ExBD[i] % M1;
                    g1 = ExBD[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK2[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK1[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK1[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK1[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK1[g1][g2 + 1] ) )
                                    ) - I * k2hb * pot_val * ( F[g1][g2] + 0.5 * KK1[g1][g2] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-2: CASE 1 KK2) = %lf sec\n", t_1_elapsed);

                // RK4-3
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g2 = ExBD[i] % M1;
                    g1 = ExBD[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK3[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK2[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK2[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK2[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK2[g1][g2 + 1] ) )
                                    ) - I * k2hb * pot_val * ( F[g1][g2] + 0.5 * KK2[g1][g2] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-3: CASE 1 KK3) = %lf sec\n", t_1_elapsed);

                // RK4-4
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2, xx1, xx2, pot_val)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g2 = ExBD[i] % M1;
                    g1 = ExBD[i] / M1;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    if ( POT[g1][g2] < POTMIN )  {
                        pot_val = POTENTIAL(xx1, xx2);
                        POT[g1][g2] = pot_val;
                    }
                    else  {
                        pot_val = POT[g1][g2];
                    }

                    KK4[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 1.0 * KK3[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1 + 1][g2] + 1.0 * KK3[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 1.0 * KK3[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1][g2 + 1] + 1.0 * KK3[g1][g2 + 1] ) )
                                    ) - I * k2hb * pot_val * ( F[g1][g2] + 1.0 * KK3[g1][g2] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-4: CASE 1 KK4) = %lf sec\n", t_1_elapsed);

                // RK4-5
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g2 = ExBD[i] % M1;
                    g1 = ExBD[i] / M1;
                    FF[g1][g2] = F[g1][g2] + ( KK1[g1][g2] + 2.0 * KK2[g1][g2] + 2.0 * KK3[g1][g2] + KK4[g1][g2] ) / 6.0;
                }
                ExBD.clear();
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-5: CASE 1 FF) = %lf sec\n", t_1_elapsed);

                t_1_begin = omp_get_wtime();
                #pragma omp parallel for private(g1, g2)
                for (int  i = 0; i < ExFF.size(); i++)  {

                    g2 = ExFF[i] % M1;
                    g1 = ExFF[i] / M1;
                    PFdX1[g1][g2] = 0.5 * Hi[0] * std::abs( FF[g1+1][g2] - FF[g1-1][g2] );
                    PFdX2[g1][g2] = 0.5 * Hi[1] * std::abs( FF[g1][g2+1] - FF[g1][g2-1] );
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-6: CASE 1 PF) = %lf sec\n", t_1_elapsed);
            }

            // Check Multiple Expanding 
            // TBL = index of FF that FF(TBL) is higher than TolL

            t_1_begin = omp_get_wtime();
            TBL.clear();

            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, b1, b2, b3, b4, b5)
            for (int i = 0; i < ExFF.size(); i++)
            {
                g2 = ExFF[i] % M1;
                g1 = ExFF[i] / M1;
                b1 = std::abs(FF[g1][g2] * std::conj(FF[g1][g2])) >= TolH;
                b2 = PFdX1[g1][g2] >= TolHd;
                b3 = PFdX2[g1][g2] >= TolHd;
                b4 = g1 > 2 && g2 > 2;
                b5 = g1 < BoxShape[0] - 3 && g2 < BoxShape[1] - 3;

                if (  ( b1 || b2 || b3 ) && b4 && b5 )  {
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
            t_overhead += t_1_elapsed;
            if (!QUIET) log->log("TBL size = %d TBL_P size = %d\n", TBL.size(), TBL_P.size()); 
            if (!QUIET && TIMING) log->log("Elapsed time (omp-c-3 CASE 1 TBL TBL_P) = %lf sec\n", t_1_elapsed); 

            // Update isFirstExtrp
            isFirstExtrp = (TBL.size() == 0) ? 1 : 0;
        }
        // .........................................................................................

        // CASE 2: Truncating without extrapolation

        if ( !isExtrapolate && !isFullGrid )
        {
            t_1_begin = omp_get_wtime();

            // RK4-1
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i = 0; i < TA.size(); i++)  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;
                xx1 = Box[0] + g1 * H[0];
                xx2 = Box[2] + g2 * H[1];

                KK1[g1][g2] = I * kh2m * (  Hisq[0] * ( F[g1 - 1][g2] - 2.0 * F[g1][g2] + F[g1 + 1][g2] )
                                  + Hisq[1] * ( F[g1][g2 - 1] - 2.0 * F[g1][g2] + F[g1][g2 + 1] ) 
                               )  - I * k2hb * POTENTIAL(xx1, xx2) * F[g1][g2];
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_truncate += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-11: CASE 2 KK1) = %lf sec\n", t_1_elapsed);

            // RK4-2
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i = 0; i < TA.size(); i++)  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;
                xx1 = Box[0] + g1 * H[0];
                xx2 = Box[2] + g2 * H[1];

                KK2[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK1[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK1[g1 + 1][g2] ) )
                                  + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK1[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK1[g1][g2 + 1] ) )
                                ) - I * k2hb * POTENTIAL(xx1, xx2) * ( F[g1][g2] + 0.5 * KK1[g1][g2] );
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_truncate += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-12: CASE 2 KK2) = %lf sec\n", t_1_elapsed);

            // RK4-3
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i = 0; i < TA.size(); i++)  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;
                xx1 = Box[0] + g1 * H[0];
                xx2 = Box[2] + g2 * H[1];

                KK3[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK2[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK2[g1 + 1][g2] ) )
                                  + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK2[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK2[g1][g2 + 1] ) )
                                ) - I * k2hb * POTENTIAL(xx1, xx2) * ( F[g1][g2] + 0.5 * KK2[g1][g2] );
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_truncate += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-13: CASE 2 KK3) = %lf sec\n", t_1_elapsed);

            // RK4-4
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i = 0; i < TA.size(); i++)  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;
                xx1 = Box[0] + g1 * H[0];
                xx2 = Box[2] + g2 * H[1];

                KK4[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 1.0 * KK3[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1 + 1][g2] + 1.0 * KK3[g1 + 1][g2] ) )
                                  + Hisq[1] * ( ( F[g1][g2 - 1] + 1.0 * KK3[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1][g2 + 1] + 1.0 * KK3[g1][g2 + 1] ) )
                                ) - I * k2hb * POTENTIAL(xx1, xx2) * ( F[g1][g2] + 1.0 * KK3[g1][g2] );
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_truncate += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-14: CASE 2 KK4) = %lf sec\n", t_1_elapsed);

            // RK4-5
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2)
            for (int i = 0; i < TA.size(); i++)  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;

                FF[g1][g2] = F[g1][g2] + ( KK1[g1][g2] + 2.0 * KK2[g1][g2] + 2.0 * KK3[g1][g2] + KK4[g1][g2] ) / 6.0;
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_truncate += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-15: CASE 2 FF) = %lf sec\n", t_1_elapsed);
        } 
        else if ( !isExtrapolate && isFullGrid )
        {
            // .........................................................................................

            // CASE 3: Full grid

            // RK4-1
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                    g1 = i1;
                    g2 = i2;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    KK1[g1][g2] = I * kh2m * (  Hisq[0] * ( F[g1 - 1][g2] - 2.0 * F[g1][g2] + F[g1 + 1][g2] )
                                    + Hisq[1] * ( F[g1][g2 - 1] - 2.0 * F[g1][g2] + F[g1][g2 + 1] ) 
                                    )  - I * k2hb * POTENTIAL(xx1, xx2) * F[g1][g2];
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_full += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-21: CASE 3 KK1) = %lf sec\n", t_1_elapsed);

            // RK4-2
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                    g1 = i1;
                    g2 = i2;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    KK2[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK1[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK1[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK1[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK1[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK1[g1][g2 + 1] ) )
                                    ) - I * k2hb * POTENTIAL(xx1, xx2) * ( F[g1][g2] + 0.5 * KK1[g1][g2] );
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_full += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-22: CASE 3 KK2) = %lf sec\n", t_1_elapsed);

            // RK4-3
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                    g1 = i1;
                    g2 = i2;
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];

                    KK3[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 0.5 * KK2[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1 + 1][g2] + 0.5 * KK2[g1 + 1][g2] ) )
                                    + Hisq[1] * ( ( F[g1][g2 - 1] + 0.5 * KK2[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 0.5 * KK2[g1][g2] ) + ( F[g1][g2 + 1] + 0.5 * KK2[g1][g2 + 1] ) )
                                    ) - I * k2hb * POTENTIAL(xx1, xx2) * ( F[g1][g2] + 0.5 * KK2[g1][g2] );
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_full += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-23: CASE 3 KK3) = %lf sec\n", t_1_elapsed);

            // RK4-4
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2, xx1, xx2)
            for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                        g1 = i1;
                        g2 = i2;
                        xx1 = Box[0] + g1 * H[0];
                        xx2 = Box[2] + g2 * H[1];

                        KK4[g1][g2] = I * kh2m * (  Hisq[0] * ( ( F[g1 - 1][g2] + 1.0 * KK3[g1 - 1][g2] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1 + 1][g2] + 1.0 * KK3[g1 + 1][g2] ) )
                                        + Hisq[1] * ( ( F[g1][g2 - 1] + 1.0 * KK3[g1][g2 - 1] ) - 2.0 * ( F[g1][g2] + 1.0 * KK3[g1][g2] ) + ( F[g1][g2 + 1] + 1.0 * KK3[g1][g2 + 1] ) )
                                        ) - I * k2hb * POTENTIAL(xx1, xx2) * ( F[g1][g2] + 1.0 * KK3[g1][g2] );
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_full += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-24: CASE 3 KK4) = %lf sec\n", t_1_elapsed);

            // RK4-5
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2)
            for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {
                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                    g1 = i1;
                    g2 = i2;
                    FF[g1][g2] = F[g1][g2] + ( KK1[g1][g2] + 2.0 * KK2[g1][g2] + 2.0 * KK3[g1][g2] + KK4[g1][g2] ) / 6.0;
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_full += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-5: CASE 3 FF) = %lf sec\n", t_1_elapsed); 
        }
        // .........................................................................................

        // FF(t+1) Normailzed & go on

        t_1_begin = omp_get_wtime();
        norm = 0.0;

        if (!isFullGrid)  {

            #pragma omp parallel for private(g1, g2) reduction (+:norm)
            for (int i = 0; i < TA.size(); i++)  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;
                norm += std::abs(FF[g1][g2] * std::conj(FF[g1][g2]));
            }
        }  
        else  {

            #pragma omp parallel for reduction (+:norm)
            for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                    norm += std::abs(FF[i1][i2] * std::conj(FF[i1][i2]));
                }
            }
        }
        log->log("norm = %lf\n", norm);
        norm *= H[0] * H[1];
        norm = 1.0 / sqrt(norm);

        #pragma omp parallel for
        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                FF[i1][i2] = norm * FF[i1][i2];
                F[i1][i2] = FF[i1][i2];
                PF[i1][i2] = std::abs(F[i1][i2] * std::conj(F[i1][i2]));
            }
        }     
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1 FF) = %lf sec\n", t_1_elapsed); 

        // Truncated_New Edge
        if ( !isFullGrid )
        {
            // PFdX
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for
            for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                    PFdX1[i1][i2] = std::abs(F[i1+1][i2] - F[i1-1][i2]) * coeff1;
                    PFdX2[i1][i2] = std::abs(F[i1][i2+1] - F[i1][i2-1]) * coeff2;
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-2 PF) = %lf sec\n", t_1_elapsed); 

            // Truncate
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for
            for (int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {
                for (int i2 = 1; i2 < BoxShape[1]  - 1; i2 ++)  {
            
                    if ( ( PF[i1][i2] < TolH ) && ( PFdX1[i1][i2] < TolHd ) && ( PFdX2[i1][i2] < TolHd ) )  {
                        F[i1][i2] = xZERO;
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-1 PF) = %lf sec\n", t_1_elapsed);
            
            // TA
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, b1, b2)
            for ( int i = 0; i < TA.size(); i++ )  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;
                b1 = F[g1][g2] != xZERO;
                b2 = PFdX1[g1][g2] >= TolHd || PFdX2[g1][g2] >= TolHd;
                TAMask[g1][g2] = 0;

                if ( b1 || b2 )  {
                    tmpVec.push_back(TA[i]);
                    TAMask[g1][g2] = 1;
                }
            }
            tmpVec.swap(TA);  
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-2 TA) = %lf sec\n", t_1_elapsed);
            
            // TB
            t_1_begin = omp_get_wtime(); 
            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, b1, b2, b3, b4)
            for ( int i = 0; i < TA.size(); i++ )  {

                g2 = TA[i] % M1;
                g1 = TA[i] / M1;

                b1 = F[g1-1][g2] == xZERO;
                b2 = F[g1+1][g2] == xZERO;
                b3 = F[g1][g2-1] == xZERO;
                b4 = F[g1][g2+1] == xZERO;

                if ( b1 || b2 || b3 || b4 )
                {
                    b1 = ( PFdX1[g1-1][g2] < TolHd ) && ( PFdX2[g1-1][g2] < TolHd )  ;
                    b2 = ( PFdX1[g1+1][g2] < TolHd ) && ( PFdX2[g1+1][g2] < TolHd )  ;
                    b3 = ( PFdX1[g1][g2-1] < TolHd ) && ( PFdX2[g1][g2-1] < TolHd )  ;
                    b4 = ( PFdX1[g1][g2+1] < TolHd ) && ( PFdX2[g1][g2+1] < TolHd )  ;
     
                    if ( b1 || b2 || b3 || b4 ) {

                        tmpVec.push_back(TA[i]);  
                    }
                }
            }
            tmpVec.swap(TB);
            tmpVec.clear();
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-3 TB) = %lf sec\n", t_1_elapsed);
        
            // TA expansion
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2)                 
            for (int i = 0; i < TA.size(); i++)
            {
                g2 = TA[i] % M1;
                g1 = TA[i] / M1;

                if (g1 + 1 != BoxShape[0] - 1 && TAMask[g1+1][g2] == 0)  {
                    tmpVec.push_back(GridToIdx(g1+1,g2));
                }
                if (g1 - 1 != 0 && TAMask[g1-1][g2] == 0)  {
                    tmpVec.push_back(GridToIdx(g1-1,g2));
                }
                if (g2 + 1 != BoxShape[1] - 1 && TAMask[g1][g2+1] == 0)  {
                    tmpVec.push_back(GridToIdx(g1,g2+1));
                }
                if (g2 - 1 != 0 && TAMask[g1][g2-1] == 0)  {
                    tmpVec.push_back(GridToIdx(g1,g2-1));
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4-1 push_back) = %lf sec\n", t_1_elapsed);
            
            // Find unique elements
            t_1_begin = omp_get_wtime();
            __gnu_parallel::sort(tmpVec.begin(),tmpVec.end());
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4-2-1 sort) = %lf sec\n", t_1_elapsed);
            
            t_1_begin = omp_get_wtime();
            it = std::unique (tmpVec.begin(), tmpVec.end()); 
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4-2-2 unique) = %lf sec\n", t_1_elapsed);
            
            t_1_begin = omp_get_wtime();
            tmpVec.resize(std::distance(tmpVec.begin(),it));
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4-2-3 resize) = %lf sec\n", t_1_elapsed);
            
            // Combine TA and tmpVec
            t_1_begin = omp_get_wtime();
            TA.reserve(TA.size() + tmpVec.size());
            TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4-2-4 combine) = %lf sec\n", t_1_elapsed);
            
            // Update TA Mask
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for private(g1, g2)
            for (int i = 0; i < tmpVec.size(); i++)
            {
                g2 = tmpVec[i] % M1;
                g1 = tmpVec[i] / M1;
                TAMask[g1][g2] = 1;
            }      
            tmpVec.clear();
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4-2-5 update mask) = %lf sec\n", t_1_elapsed);
            
            t_1_begin = omp_get_wtime();
            #pragma omp parallel for
            for (int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {
                for (int i2 = 1; i2 < BoxShape[1]  - 1; i2 ++)  {

                    PF[i1][i2] = std::abs( F[i1][i2] * std::conj(F[i1][i2]) );
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4-2-6 PF) = %lf sec\n", t_1_elapsed);
        }

        if ( (tt + 1) % PERIOD == 0 )
        {
            // REPORT MEASUREMENTS

            // ----------------------------------------------------------------------------
            // Compute Transmittance
            // ----------------------------------------------------------------------------

            if (isTrans)  {
            
                t_1_begin = omp_get_wtime();
                pftrans = 0.0;

                #pragma omp parallel for reduction (+:pftrans)
                for (int i1 = idx_x0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        pftrans+=PF[i1][i2];
                    }
                }
                pftrans *= H[0] * H[1];
                PF_trans.push_back(pftrans);
                log->log("[Scatter2d] Time %lf, Trans = %e\n", ( tt + 1 ) * kk, pftrans);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-x-2 trans) = %lf sec\n", t_1_elapsed); 
            }
            // ----------------------------------------------------------------------------
            // Compute auto-correlation function
            // ----------------------------------------------------------------------------

            if (isAcf)  {

                t_1_begin = omp_get_wtime();
                norm = 0.0;

                #pragma omp parallel for reduction (+:norm)
                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        norm += std::real(F0_STAR[i1][i2] * FF[i1][i2]);
                    }
                }
                norm_re = norm * H[0] * H[1];
                norm = 0.0;

                #pragma omp parallel for reduction (+:norm)
                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        norm += std::imag(F0_STAR[i1][i2] * FF[i1][i2]);
                    }
                }
                norm_im = norm * H[0] * H[1];
                acf = {norm_re, norm_im};
                ACFunc.push_back(acf);

                //log->log("[Scatter2d] Step: %d, Time = %f ACF(Re) = %.16e \n", tt + 1, kk * (tt + 1), norm_re);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if ( !QUIET && TIMING ) log->log("Elapsed time (omp-x-3 trans) = %lf sec\n", t_1_elapsed); 
             }
        }

        // Reset

        t_1_begin = omp_get_wtime();

        #pragma omp parallel for
        for (int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {
            for (int i2 = 1; i2 < BoxShape[1] - 1 ; i2 ++)  {
          
                FF[i1][i2] = xZERO;
                KK1[i1][i2] = xZERO;
                KK2[i1][i2] = xZERO;
                KK3[i1][i2] = xZERO;
                KK4[i1][i2] = xZERO;
            }
        }
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if ( !QUIET && TIMING ) log->log("Elapsed time (omp-e-4: reset) = %lf sec\n", t_1_elapsed);  

        if ( (tt + 1) % SORT_PERIOD == 0 && !isFullGrid )
        {
            t_1_begin = omp_get_wtime();
            __gnu_parallel::sort (TA.begin(), TA.end());
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_truncate += t_1_elapsed;
            if ( !QUIET && TIMING ) log->log("Elapsed time (omp-e-5: sort) = %lf sec\n", t_1_elapsed);
        }

        if ( (tt + 1) % PERIOD == 0 )
        {   
            t_0_end = omp_get_wtime();
            t_0_elapsed = t_0_end - t_0_begin;

            if ( !QUIET ) log->log("[Scatter2d] Step: %d, Elapsed time: %lf sec\n", tt + 1, t_0_elapsed);

            if ( !isFullGrid && !QUIET )  {

                log->log("[Scatter2d] Core computation time = %lf\n", t_truncate);
                log->log("[Scatter2d] Overhead time = %lf\n", t_overhead);
                log->log("[Scatter2d] TA size = %d, TB size = %d\n", TA.size(), TB.size());
                log->log("[Scatter2d] TA / total grids = %lf\n", ( TA.size() * 1.0 ) / GRIDS_TOT);
            }
            else if ( isFullGrid && !QUIET )  {

                log->log("[Scatter2d] Core computation time = %lf\n", t_full);
            }
            if ( !QUIET ) log->log("\n........................................................\n\n");
        }         
    } // Time iteration 

    // ACF
    if ( isAcf )
    {
         for ( int i = 0; i < ACFunc.size(); i ++ )
         {
             log->log("[Scatter2d] Time = %f ACF = ( %.16e , %.16e ) \n", i * PERIOD * kk, ACFunc[i][0], ACFunc[i][1]);
         }
    }

    // Compute spectrum from ACFunc
    if ( isAcf )
    {
        double tau_window = 15.0;

        t_1_begin = omp_get_wtime();

        norm = 0.0;

        for ( int i = 0; i < int(kMax / dk) + 1; i ++)  {
            for ( int j = 0; j < ACFunc.size(); j ++)  {

                norm += ACFunc[j][0] * cos(i * dk * j * PERIOD * kk) - ACFunc[j][1] * sin(i * dk * j * PERIOD * kk) * exp(-pow(j * PERIOD * kk / tau_window, 2));
            }
            norm *= i * dk / PI * (PERIOD * kk);
            Spectrum.push_back(norm);
            log->log("[Scatter2d] k = %f Spectrum = %.16e \n", i * dk, norm);
        }

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-x-3 spectrum) = %lf sec\n", t_1_elapsed);
    }

    log->log("[Scatter2d] Evolve done.\n");
}
/* ------------------------------------------------------------------------------- */

// Translational component (1): Eckart
// Vibrational component (1): Morse

inline std::complex<double> Scatter2d::Wavefunction_EckMO(double x1, double x2)
{
    return exp( -A[0] * pow(x1 - Wave0[0], 2) + I / hb * P[0] * x1) * exp( - Da * (Ld - 0.5) * (x2 - r0)  - Ld * exp(- Da * (x2 - r0)) );
}
/* ------------------------------------------------------------------------------- */

inline double Scatter2d::Potential_EckMO(double x1, double x2)
{
    return V0 * pow(cosh(alpha * x1), -2.0) + (1.0 - sigma * exp(- lambda * x1 * x1)) * De * pow(1.0 - exp(- Da * x2), 2.0);
}
/* ------------------------------------------------------------------------------- */

// Translational component (1): Gaussian
// Vibrational component (1): Morse

inline std::complex<double> Scatter2d::Wavefunction_GauMO(double x1, double x2)
{
    return exp( -A[0] * pow(x1 - Wave0[0], 2) + I / hb * P[0] * x1) * exp( - Da * (Ld - 0.5) * (x2 - r0) - Ld * exp(- Da * (x2 - r0)) );
}
/* ------------------------------------------------------------------------------- */

inline double Scatter2d::Potential_GauMO(double x1, double x2)
{
    return V0 * exp(-beta * x1 * x1) + (1.0 - sigma * exp(- lambda * x1 * x1)) * De * pow(1.0 - exp( - Da * x2), 2.0);
}
/* ------------------------------------------------------------------------------- */

// Translational component (1): Eckart
// Vibrational component (1): Harmonic

inline std::complex<double> Scatter2d::Wavefunction_EckHO(double x1, double x2)
{
    return exp( -A[0] * pow(x1 - Wave0[0], 2) + I / hb * P[0] * x1 - A[1] * pow(x2, 2) );
}
/* ------------------------------------------------------------------------------- */

inline double Scatter2d::Potential_EckHO(double x1, double x2)
{     
    return V0 * pow(cosh(alpha * x1), -2.0) + 0.5 * k0 * (1.0 - sigma * exp(- lambda * x1 * x1)) * x2 * x2;
}
/* ------------------------------------------------------------------------------- */

// Translational component (1): Gaussian
// Vibrational component (1): Harmonic

inline std::complex<double> Scatter2d::Wavefunction_GauHO(double x1, double x2)
{
    return exp( -A[0] * pow(x1 - Wave0[0], 2) + I / hb * P[0] * x1 - A[1] * pow(x2, 2) );
}
/* ------------------------------------------------------------------------------- */

inline double Scatter2d::Potential_GauHO(double x1, double x2)
{     
    return V0 * exp(-beta * x1 * x1) + 0.5 * k0 * (1.0 - sigma * exp(- lambda * x1 * x1)) * x2 * x2;
}
/* ------------------------------------------------------------------------------- */

// Henon-Heiles Potential

inline std::complex<double> Scatter2d::Wavefunction_HH(double x1, double x2)
{
    return pow(PI_INV, 0.25 * DIMENSIONS) * exp( - 0.5 * ((x1 - 2.0) * (x1 - 2.0) + (x2 - 2.0) * (x2 - 2.0)) );
}
/* ------------------------------------------------------------------------------- */

inline double Scatter2d::Potential_HH(double x1, double x2)
{     
    return 0.5 * (x1 * x1 + x2 * x2) + lambda * (x1 * x1 * x2 - (x2 * x2 * x2) / 3.0);
}
/* ------------------------------------------------------------------------------- */

VectorXi Scatter2d::IdxToGrid(int idx)
{
    int x1 = idx / M1;
    int x2 = idx % M1;

    VectorXi grid;
    grid.resize(DIMENSIONS);
    grid << x1, x2;

    return grid;
}
/* ------------------------------------------------------------------------------- */

inline int Scatter2d::GridToIdx(int x1, int x2)
{
    return x1 * M1 + x2;
}
/* ------------------------------------------------------------------------------- */

inline void Scatter2d::DefineBoundary()
{
    // Find elements of DBi and DBi2

    // i1 
    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

        DBi.push_back(GridToIdx(0, i2));
        DBi.push_back(GridToIdx(BoxShape[0]-1, i2));

        DBi2.push_back(GridToIdx(0, i2));
        DBi2.push_back(GridToIdx(1, i2));
        DBi2.push_back(GridToIdx(2, i2));
        DBi2.push_back(GridToIdx(BoxShape[0]-1, i2));
        DBi2.push_back(GridToIdx(BoxShape[0]-2, i2));
        DBi2.push_back(GridToIdx(BoxShape[0]-3, i2));
    }

    // i2

    for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

        DBi.push_back(GridToIdx(i1, 0));
        DBi.push_back(GridToIdx(i1, BoxShape[1]-1));

        DBi2.push_back(GridToIdx(i1, 0));
        DBi2.push_back(GridToIdx(i1, 1));
        DBi2.push_back(GridToIdx(i1, 2));
        DBi2.push_back(GridToIdx(i1, BoxShape[1]-1));
        DBi2.push_back(GridToIdx(i1, BoxShape[1]-2));
        DBi2.push_back(GridToIdx(i1, BoxShape[1]-3));
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

