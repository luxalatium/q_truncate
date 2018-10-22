// ==============================================================================
//
//  Imt4d.cpp
//  QTR
//
//  Created by Albert Lu on 10/3/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/5/18
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

#include "Containers.h"
#include "Error.h"
#include "Log.h"
#include "Parameters.h"
#include "Imt4d.h"

#include "boost/multi_array.hpp"

using namespace QTR_NS;
using std::vector;
using std::cout;
using std::endl;

/* ------------------------------------------------------------------------------- */

Imt4d::Imt4d(class QTR *q)
{
        qtr = q;
        err = qtr->error;
        log = qtr->log;
        parameters = qtr->parameters;
        init();
} 
/* ------------------------------------------------------------------------------- */

Imt4d::~Imt4d()
{       
        return;
}
/* ------------------------------------------------------------------------------- */

void Imt4d::init()
{
        log->log("\n\n[Imt4d] INIT starts ...\n");

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
        S.resize(DIMENSIONS);  
        kk = parameters->scxd_k;

        H[0] = parameters->scxd_h1;
        H[1] = parameters->scxd_h2;
        H[2] = parameters->scxd_h3;
        H[3] = parameters->scxd_h4;

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
        Box[4] = parameters->scxd_xi3;
        Box[5] = parameters->scxd_xf3; 
        Box[6] = parameters->scxd_xi4;
        Box[7] = parameters->scxd_xf4; 

        BoxShape.resize(DIMENSIONS);
        GRIDS_TOT = 1;

        log->log("[Imt4d] Number of grids = (");

        for (unsigned int i = 0; i < DIMENSIONS; i ++)  {

                BoxShape[i] = (int)( (Box[2 * i + 1] - Box[2 * i]) / H[i] ) + 1;
                GRIDS_TOT *= BoxShape[i];

                if (i < DIMENSIONS - 1)
                        log->log("%d, ", BoxShape[i]);
                else
                        log->log("%d)\n", BoxShape[i]);
        }

        M1 = BoxShape[1] * BoxShape[2] * BoxShape[3];
        M2 = BoxShape[2] * BoxShape[3];
        M3 = BoxShape[3];

        // Potential: form

        //Vmode.resize(DIMENSIONS);
        //Vmode[0] = parameters->scxd_Vmode_1;
        //Vmode[1] = parameters->scxd_Vmode_2;

        hb = parameters->scxd_hb;
        m  = parameters->scxd_m;

        // Potential: HO specific
        w = parameters->scxd_w;
        mw2 = m * w * w;
        aa = 0.5 * m * w / hb;
        log->log("[Imt4d] aa = %lf\n", aa);

        // Potential: Eckart
        V0 = parameters->scxd_V0;
        alpha = parameters->scxd_alpha;
        ek2v = parameters->scxd_ek2v;
        Ek0 = V0 * ek2v;
        log->log("[Imt4d] Ek0 = %lf\n", Ek0);

        // Potential: Related HO
        k0 = parameters->scxd_k0;  
        sig = parameters->scxd_sig;
        lan = parameters->scxd_lan;

        // Potential: Morse
        De = parameters->scxd_De;
        Da = parameters->scxd_Da;
        r0 = parameters->scxd_r0;
        Ld = sqrt(2 * m * De) / (Da * hb);    
        log->log("[Imt4d] Ld = %lf\n", Ld);

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
        //P[0] = sqrt(2.0 * m * Ek0);
        //P[1] = 0;                                
        //log->log("[Imt4d] P[0] = %lf\n", P[0]);

        // Truncate parameters

        isFullGrid = parameters->scxd_isFullGrid;
        TolH = parameters->scxd_TolH;           // Tolerance of probability density for Zero point Cutoff
        TolL = parameters->scxd_TolL;           // Tolerance of probability density for Edge point
        TolHd = parameters->scxd_TolHd;         // Tolerance of probability first diff for Zero point Cutoff
        TolLd = parameters->scxd_TolLd;         // Tolerance of probability density for Edge point
        ExReduce = parameters->scxd_ExReduce;   //Extrapolation reduce factor

        log->log("[Imt4d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void Imt4d::Evolve()
{
	#pragma omp declare reduction (merge : MeshIndex : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

        log->log("[Imt4d] Evolve starts ...\n");

        // Variables 
        int       index;
        bool      b1, b2, b3, b4, b5, b6, b7, b8;
        double    k2hb = kk / hb;
        double    h2m = hb / 2 / m;
	double    hsq2m = hb * hb / 2 / m;
        double    coeff1, coeff2, coeff3, coeff4;
        double    t_0_begin, t_0_end, t_0_elapsed;
        double    t_1_begin, t_1_end, t_1_elapsed;
        MeshIndex tmpVec;  // temporary index container
       
        // normalization factor
        double    norm;
 
	// Expectation value of energy
	double energy;
 
        // 3d Grid vector and indices
        VectorXi  grid;
        int       g1, g2, g3, g4;
        double    xx1, xx2, xx3, xx4;

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

        log->log("[Imt4d] Initializing containers ...\n");

        t_1_begin = omp_get_wtime();

        // Initialize containers
        MeshCX4D F(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]][BoxShape[3]]);
        F.reindex(0);

        MeshD4D PF(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]][BoxShape[3]]);
        PF.reindex(0);
        
        MeshD4D PFdX1(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]][BoxShape[3]]);
        PFdX1.reindex(0);  

        MeshD4D PFdX2(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]][BoxShape[3]]);
        PFdX2.reindex(0);

        MeshD4D PFdX3(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]][BoxShape[3]]);
        PFdX3.reindex(0);

        MeshD4D PFdX4(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]][BoxShape[3]]);
        PFdX4.reindex(0);

        MeshCX4D FF(boost::extents[BoxShape[0]][BoxShape[1]][BoxShape[2]][BoxShape[3]]);
        FF.reindex(0);

        #pragma omp parallel for

        for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                        for (unsigned int i3 = 0; i3 < BoxShape[2]; i3 ++)  {

                                for (unsigned int i4 = 0; i4 < BoxShape[3]; i4 ++)  {

                                        F[i1][i2][i3][i4] = xZERO;
                                        PF[i1][i2][i3][i4] = 0.0;
                                        PFdX1[i1][i2][i3][i4] = 0.0;
                                        PFdX2[i1][i2][i3][i4] = 0.0;
                                        PFdX3[i1][i2][i3][i4] = 0.0;
                                        PFdX4[i1][i2][i3][i4] = 0.0;
                                }
                        }
                }
        }
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        log->log("[Imt4d] Elapsed time (initializing containers) = %lf sec\n\n", t_1_elapsed); 

        // .........................................................................................

        log->log("[Imt4d] Initializing wavefunction ...\n");    

        t_1_begin = omp_get_wtime();

        // Initialize wavefunction

        #pragma omp parallel for

        for (unsigned int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {

                for (unsigned int i2 = 1; i2 < BoxShape[1] - 1 ; i2 ++)  {

                        for (unsigned int i3 = 1; i3 < BoxShape[2] - 1 ; i3 ++)  {

                                for (unsigned int i4 = 1; i4 < BoxShape[3] - 1 ; i4 ++)  {
                  
                                        F[i1][i2][i3][i4] = Wavefunction(Box[0] + i1 * H[0], Box[2] + i2 * H[1], Box[4] + i3 * H[2], Box[6] + i4 * H[3]);
                                }
                        }
                }
        }

        // Normalization

        norm = 0.0;

        #pragma omp parallel for reduction (+:norm)

        for (unsigned int i1 = 0; i1 <  BoxShape[0]; i1 ++)  {

                for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                        for (unsigned int i3 = 0; i3 < BoxShape[2]; i3 ++)  {

                                for (unsigned int i4 = 0; i4 < BoxShape[3]; i4 ++)  {
                                
			                norm += std::abs(F[i1][i2][i3][i4] * std::conj(F[i1][i2][i3][i4]));
                                }
                        }
                }
        }
        norm *= H[0] * H[1] * H[2] * H[3];

        log->log("[Imt4d] Normalization factor = %e\n",norm);
        norm = 1.0 / sqrt(norm);

        #pragma omp parallel for

        for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                        for (unsigned int i3 = 0; i3 < BoxShape[2]; i3 ++)  {

                                for (unsigned int i4 = 0; i4 < BoxShape[3]; i4 ++)  {

                                        F[i1][i2][i3][i4] = norm * F[i1][i2][i3][i4];
					PF[i1][i2][i3][i4] = std::abs(F[i1][i2][i3][i4] * std::conj(F[i1][i2][i3][i4]));
                                }
                        }
                }
        }       

        // Initial probability

        //#pragma omp parallel for

        //for (unsigned int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

          //      for (unsigned int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

            //            for (unsigned int i3 = 0; i3 < BoxShape[2]; i3 ++)  {

              //                  for (unsigned int i4 = 0; i4 < BoxShape[3]; i4 ++)  {

                //                        PF[i1][i2][i3][i4] = std::abs(F[i1][i2][i3][i4] * std::conj(F[i1][i2][i3][i4]));
                  //              }
                    //    }
               // }
        //}
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        log->log("[Imt4d] Elapsed time (initializing wavefunction) = %lf sec\n\n", t_1_elapsed); 

        // .........................................................................................

        log->log("[Imt4d] Computing finite differences ...\n");   

        t_1_begin = omp_get_wtime();

        // Pf 1st-order finite difference 

        // PFdX1

        coeff1 = 1.0 / (2.0 * H[0]);
        coeff2 = 1.0 / (2.0 * H[1]);
        coeff3 = 1.0 / (2.0 * H[2]);
        coeff4 = 1.0 / (2.0 * H[3]);

        #pragma omp parallel for

        for (unsigned int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                for (unsigned int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                        for (unsigned int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                for (unsigned int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                        PFdX1[i1][i2][i3][i4] = std::abs(F[i1+1][i2][i3][i4] - F[i1-1][i2][i3][i4]) * coeff1;
                                        PFdX2[i1][i2][i3][i4] = std::abs(F[i1][i2+1][i3][i4] - F[i1][i2-1][i3][i4]) * coeff2;
                                        PFdX3[i1][i2][i3][i4] = std::abs(F[i1][i2][i3+1][i4] - F[i1][i2][i3-1][i4]) * coeff3;
                                        PFdX4[i1][i2][i3][i4] = std::abs(F[i1][i2][i3][i4+1] - F[i1][i2][i3][i4-1]) * coeff4;
                                }
                        }
                }
        }        

        // PFdX2       

        /*coeff = 1.0 / (2.0 * H[1]);

        #pragma omp parallel for

        for (unsigned int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                for (unsigned int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                        for (unsigned int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                for (unsigned int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                        PFdX2[i1][i2][i3][i4] = std::abs(F[i1][i2+1][i3][i4] - F[i1][i2-1][i3][i4]) * coeff;
                                }
                        }
                }
        } */ 

        // PFdX3     

        /*coeff = 1.0 / (2.0 * H[2]);

        #pragma omp parallel for

        for (unsigned int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                for (unsigned int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                        for (unsigned int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                for (unsigned int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                        PFdX3[i1][i2][i3][i4] = std::abs(F[i1][i2][i3+1][i4] - F[i1][i2][i3-1][i4]) * coeff;
                                }
                        }
                }
        } */ 

        // PFdX4   

        /*coeff = 1.0 / (2.0 * H[3]);

        #pragma omp parallel for

        for (unsigned int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                for (unsigned int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                        for (unsigned int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                for (unsigned int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                        PFdX4[i1][i2][i3][i4] = std::abs(F[i1][i2][i3][i4+1] - F[i1][i2][i3][i4-1]) * coeff;
                                }
                        }
                }
        } */ 

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        log->log("[Imt4d] Elapsed time (finite differences) = %lf sec\n\n", t_1_elapsed); 

        // .........................................................................................

        // Truncate initial & edge point check

        log->log("[Imt4d] Initial truncation ...\n");   

        t_1_begin = omp_get_wtime();

	if ( !isFullGrid )  // Truncate method
        {
                #pragma omp parallel for private(b1, b2, b3, b4, b5)

                for (unsigned int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                        for (unsigned int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                for (unsigned int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                        for (unsigned int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                b1 = PF[i1][i2][i3][i4] < TolH;
                                                b2 = PFdX1[i1][i2][i3][i4] < TolHd;
                                                b3 = PFdX2[i1][i2][i3][i4] < TolHd;
                                                b4 = PFdX3[i1][i2][i3][i4] < TolHd;
                                                b5 = PFdX4[i1][i2][i3][i4] < TolHd;

                                                if ( b1 && b2 && b3 && b4 && b5 )
                                                        F[i1][i2][i3][i4] = xZERO;
                                        }
                                }
                        }
                }

                // `````````````````````````````````````````````````````````````````
                
		#pragma omp parallel for reduction(merge: tmpVec) private(b1, b2, b3, b4, b5)

                for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                        for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                        for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                b1 = F[i1][i2][i3][i4] != std::complex<double>(0,0);
                                                b2 = PFdX1[i1][i2][i3][i4] >= TolHd;
                                                b3 = PFdX2[i1][i2][i3][i4] >= TolHd;
                                                b4 = PFdX3[i1][i2][i3][i4] >= TolHd;
                                                b5 = PFdX4[i1][i2][i3][i4] >= TolHd;
                        
                                                if ( b1 || ( b2 || b3 || b4 || b5 ))  {

                                                        tmpVec.push_back(GridToIdx(i1,i2,i3,i4));                                                                     
                                                }
                                        }
                                }
                        }
                }
		tmpVec.swap(TA);
                tmpVec.clear();

                log->log("[Imt4d] TA size = %d\n", TA.size());

                // `````````````````````````````````````````````````````````````````

                #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4, b1, b2, b3, b4, b5, b6, b7, b8)

                for (int i = 0; i < TA.size(); i++)
                {
                        g1 = TA[i] / M1;
                        g2 = (TA[i] % M1) / M2;
                        g3 = (TA[i] % M2) / M3;
                        g4 = TA[i] % M3;

                        b1 = F[ g1 - 1 ][ g2     ][ g3     ][ g4     ] == xZERO;
                        b2 = F[ g1 + 1 ][ g2     ][ g3     ][ g4     ] == xZERO;
                        b3 = F[ g1     ][ g2 - 1 ][ g3     ][ g4     ] == xZERO;
                        b4 = F[ g1     ][ g2 + 1 ][ g3     ][ g4     ] == xZERO;
                        b5 = F[ g1     ][ g2     ][ g3 - 1 ][ g4     ] == xZERO;
                        b6 = F[ g1     ][ g2     ][ g3 + 1 ][ g4     ] == xZERO;
                        b7 = F[ g1     ][ g2     ][ g3     ][ g4 - 1 ] == xZERO;
                        b8 = F[ g1     ][ g2     ][ g3     ][ g4 + 1 ] == xZERO;
                           
                        if ( b1 || b2 || b3 || b4 || b5 || b6 || b7 || b8 )
                        {
                                b1 = ( PFdX1[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd );
                                b2 = ( PFdX1[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd );
                                b3 = ( PFdX1[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd );
                                b4 = ( PFdX1[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd );
                                b5 = ( PFdX1[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd );
                                b6 = ( PFdX1[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd );
                                b7 = ( PFdX1[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd );
                                b8 = ( PFdX1[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd );

                                if ( b1 || b2 || b3 || b4 || b5 || b6 || b7 || b8 ) {

                                        tmpVec.push_back(TA[i]);  
                                }
                        }             
                }

                tmpVec.swap(TB);
                tmpVec.clear();

                log->log("[Imt4d] TB size = %d\n", TB.size());


                // `````````````````````````````````````````````````````````````````

                #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4)

                for (int i = 0; i < TA.size(); i++)
                {
                        g1 = TA[i] / M1;
                        g2 = (TA[i] % M1) / M2;
                        g3 = (TA[i] % M2) / M3;
                        g4 = TA[i] % M3;

                        if (g1 + 1 != BoxShape[0] - 1)
                                tmpVec.push_back(GridToIdx(g1 + 1, g2, g3, g4));

                        if (g1 - 1 != 0)
                                tmpVec.push_back(GridToIdx(g1 - 1, g2, g3, g4));

                        if (g2 + 1 != BoxShape[1] - 1)
                                tmpVec.push_back(GridToIdx(g1, g2 + 1, g3, g4));

                        if (g2 - 1 != 0)
                                tmpVec.push_back(GridToIdx(g1, g2 - 1, g3, g4)); 

                        if (g3 + 1 != BoxShape[2] - 1)
                                tmpVec.push_back(GridToIdx(g1, g2, g3 + 1, g4));

                        if (g3 - 1 != 0)
                                tmpVec.push_back(GridToIdx(g1, g2, g3 - 1, g4)); 

                        if (g4 + 1 != BoxShape[3] - 1)
                                tmpVec.push_back(GridToIdx(g1, g2, g3, g4 + 1));

                        if (g4 - 1 != 0)
                                tmpVec.push_back(GridToIdx(g1, g2, g3, g4 - 1)); 
                }

                // Combine TA and tmpVec
                TA.reserve(TA.size() + tmpVec.size());
                TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());
                tmpVec.clear();

                // Find unique elements
                __gnu_parallel::sort(TA.begin(),TA.end());
                it = std::unique (TA.begin(), TA.end()); 
                TA.resize(std::distance(TA.begin(),it));              
        }
        else  // Full grid approach
        {
                TB = DBi;
        }

	log->log("[Imt4d] TA size = %d, TB size = %d\n", TA.size(), TB.size());
	log->log("[Imt4d] DBi = %d DBi2 = %d\n", DBi.size(), DBi2.size());

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        log->log("[Imt4d] Elapsed time (initial truncation) = %lf sec\n\n", t_1_elapsed); 

        // .........................................................................................

        // Time iteration 

        log->log("=======================================================\n\n"); 
        log->log("[Imt4d] Time interation starts ...\n"); 
        log->log("[Imt4d] Number of steps = %d\n\n", (int)(TIME / kk)); 
        log->log("=======================================================\n\n"); 

        for (int tt = 0; tt < (int)(TIME / kk); tt ++)
        {
		t_0_begin = omp_get_wtime();
                t_1_begin = omp_get_wtime();

        	#pragma omp parallel for

        	for (int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {

                	for (int i2 = 1; i2 < BoxShape[1] - 1 ; i2 ++)  {

                                for (int i3 = 1; i3 < BoxShape[2] - 1 ; i3 ++)  {

                                        for (int i4 = 1; i4 < BoxShape[3] - 1 ; i4 ++)  {
                  
                                                FF[i1][i2][i3][i4] = xZERO;
                                        }
                                }
                	}
        	}

		t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                log->log("Elapsed time (omp-a-1: fill_n) = %lf sec\n", t_1_elapsed);   

                // Check if TB of f is higher than TolL
                
                if ( !isFullGrid )
                {
                        t_1_begin = omp_get_wtime();

                        TBL.clear();

                        // omp-a-1
                        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4, b1, b2, b3, b4, b5, b6, b7)

                        for (int i = 0; i < TB.size(); i++)
                        {
                        	g1 = TB[i] / M1;
                        	g2 = (TB[i] % M1) / M2;
                        	g3 = (TB[i] % M2) / M3;
                        	g4 = TB[i] % M3;

                                b1 =    PF[g1][g2][g3][g4] >= TolL;
                                b2 = PFdX1[g1][g2][g3][g4] >= TolLd;
                                b3 = PFdX2[g1][g2][g3][g4] >= TolLd;
                                b4 = PFdX3[g1][g2][g3][g4] >= TolLd;
                                b5 = PFdX4[g1][g2][g3][g4] >= TolLd;

                                // Not in DBi2
                                b6 = g1 > 1 && g2 > 1 && g3 > 1 && g4 > 1;
                                b7 = g1 < BoxShape[0] - 2 && g2 < BoxShape[1] - 2 && g3 < BoxShape[2] - 2 && g4 < BoxShape[3] - 2;

                                if ( (b1 || b2 || b3 || b4 || b5) && b6 && b7 )  {
                                        tmpVec.push_back(TB[i]);
                                }
                        }
			tmpVec.swap(TBL);
                        tmpVec.clear();
                        TBL_P = TBL;

			t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-a-2: TBL, TBL_P) = %lf sec\n", t_1_elapsed);   
                }      

                isExtrapolate = false;

                // CASE 1: Truncating with extrapolation

                while ( TBL.size() != 0 && !isFullGrid )
                {
                        isExtrapolate = true;

                        // Extrapolation2D
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
                        	g1 = TBL[index] / M1;
                        	g2 = (TBL[index] % M1) / M2;
                        	g3 = (TBL[index] % M2) / M3;
                        	g4 = TBL[index] % M3;

                                if ( F[g1 - 1][g2][g3][g4] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1 - 1, g2, g3, g4));
                                }

                                if ( F[g1 + 1][g2][g3][g4] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1 + 1, g2, g3, g4));
                                }

                                if ( F[g1][g2 - 1][g3][g4] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1, g2 - 1, g3, g4));
                                }

                                if ( F[g1][g2 + 1][g3][g4] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1, g2 + 1, g3, g4));
                                }

                                if ( F[g1][g2][g3 - 1][g4] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1, g2, g3 - 1, g4));
                                }

                                if ( F[g1][g2][g3 + 1][g4] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1, g2, g3 + 1, g4));
                                }

                                if ( F[g1][g2][g3][g4 - 1] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1, g2, g3, g4 - 1));
                                }

                                if ( F[g1][g2][g3][g4 + 1] == xZERO )  {

                                        ExFF.push_back(GridToIdx(g1, g2, g3, g4 + 1));
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
                        log->log("Elapsed time (omp-b-1: ExFF) = %lf sec\n", t_1_elapsed);   

                        // .....................................................................

                        // Find the direction of Outer to Edge points

                        t_1_begin = omp_get_wtime();

                        ExTBL.clear();

			it = ExFF.begin();
                        
                        while ( it != ExFF.end() )  {

                                index = std::distance( ExFF.begin(), it );

                        	g1 = ExFF[index] / M1;
                       		g2 = (ExFF[index] % M1) / M2;
                        	g3 = (ExFF[index] % M2) / M3;
                        	g4 = ExFF[index] % M3;

                                sum = xZERO;
                                count = 0;

                                isEmpty = true;
                                val_min_abs = 100000000;

                                if ( F[g1 - 1][g2][g3][g4] != xZERO )  {

					if ( std::abs(F[g1 - 1][g2][g3][g4]) < val_min_abs &&  F[g1 - 2][g2][g3][g4] != xZERO )  {
                                        	val_min_abs = std::abs(F[g1 - 1][g2][g3][g4]);
                                        	val_min = F[g1 - 1][g2][g3][g4];
					}

                                        if (F[g1 - 2][g2][g3][g4] != xZERO)  {

                                                val = exp( 2.0 * std::log(F[g1 - 1][g2][g3][g4]) - std::log(F[g1 - 2][g2][g3][g4]) );

                                                if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                                                {
                                                        sum += val;
                                                        count += 1;
                                                }
                                                isEmpty = false;
                                        }
                                }

                                if ( F[g1 + 1][g2][g3][g4] != xZERO )  {

                                        if ( std::abs(F[g1 + 1][g2][g3][g4]) < val_min_abs && F[g1 + 2][g2][g3][g4] != xZERO  )  {
                                                val_min_abs = std::abs(F[g1 + 1][g2][g3][g4]);
                                                val_min = F[g1 + 1][g2][g3][g4];
                                        }

                                        if ( F[g1 + 2][g2][g3][g4] != xZERO )  {

                                                val = exp( 2.0 * std::log(F[g1 + 1][g2][g3][g4]) - std::log(F[g1 + 2][g2][g3][g4]) );

                                                if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                                                {
                                                        sum += val;
                                                        count += 1;
                                                }

                                                isEmpty = false;
                                        }
                                }

                                if ( F[g1][g2 - 1][g3][g4] != xZERO )  {

                                        if ( std::abs(F[g1][g2 - 1][g3][g4]) < val_min_abs && F[g1][g2 - 2][g3][g4] != xZERO  )  {

                                                val_min_abs = std::abs(F[g1][g2 - 1][g3][g4]);
                                                val_min = F[g1][g2 - 1][g3][g4];
                                        }

                                        if ( F[g1][g2 - 2][g3][g4] != xZERO )  {

                                                val = exp( 2.0 * std::log(F[g1][g2 - 1][g3][g4]) - std::log(F[g1][g2 - 2][g3][g4]) );

                                                if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                                                {
                                                        sum += val;
                                                        count += 1;
                                                }

                                                isEmpty = false;
                                        }
                                }

                                if ( F[g1][g2 + 1][g3][g4] != xZERO )  {

                                        if ( std::abs(F[g1][g2 + 1][g3][g4]) < val_min_abs && F[g1][g2 + 2][g3][g4] != xZERO )  {

                                                val_min_abs = std::abs(F[g1][g2 + 1][g3][g4]);
                                                val_min = F[g1][g2 + 1][g3][g4];
                                        }

                                        if ( F[g1][g2 + 2][g3][g4] != xZERO )  {

                                                val = exp( 2.0 * std::log(F[g1][g2 + 1][g3][g4]) - std::log(F[g1][g2 + 2][g3][g4]) );

                                                if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                                                {
                                                        sum += val;
                                                        count += 1;
                                                }

                                                isEmpty = false;
                                        }
                                }

                                if ( F[g1][g2][g3 - 1][g4] != xZERO )  {

                                        if ( std::abs(F[g1][g2][g3 - 1][g4]) < val_min_abs && F[g1][g2][g3 - 2][g4] != xZERO  )  {

                                                val_min_abs = std::abs(F[g1][g2][g3 - 1][g4]);
                                                val_min = F[g1][g2][g3 - 1][g4];
                                        }

                                        if ( F[g1][g2][g3 - 2][g4] != xZERO )  {

                                                val = exp( 2.0 * std::log(F[g1][g2][g3 - 1][g4]) - std::log(F[g1][g2][g3 - 2][g4]) );

                                                if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                                                {
                                                        sum += val;
                                                        count += 1;
                                                }

                                                isEmpty = false;
                                        }
                                }

                                if ( F[g1][g2][g3 + 1][g4] != xZERO )  {

                                        if ( std::abs(F[g1][g2][g3 + 1][g4]) < val_min_abs && F[g1][g2][g3 + 2][g4] != xZERO )  {

                                                val_min_abs = std::abs(F[g1][g2][g3 + 1][g4]);
                                                val_min = F[g1][g2][g3 + 1][g4];
                                        }

                                        if ( F[g1][g2][g3 + 2][g4] != xZERO )  {

                                                val = exp( 2.0 * std::log(F[g1][g2][g3 + 1][g4]) - std::log(F[g1][g2][g3 + 2][g4]) );

                                                if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                                                {
                                                        sum += val;
                                                        count += 1;
                                                }

                                                isEmpty = false;
                                        }
                                }

                                if ( F[g1][g2][g3][g4 - 1] != xZERO )  {

                                        if ( std::abs(F[g1][g2][g3][g4 - 1]) < val_min_abs && F[g1][g2][g3][g4 - 2] != xZERO  )  {

                                                val_min_abs = std::abs(F[g1][g2][g3][g4 - 1]);
                                                val_min = F[g1][g2][g3][g4 - 1];
                                        }

                                        if ( F[g1][g2][g3][g4 - 2] != xZERO )  {

                                                val = exp( 2.0 * std::log(F[g1][g2][g3][g4 - 1]) - std::log(F[g1][g2][g3][g4 - 2]) );

                                                if ( !isnan( real(val) ) && !isnan( imag(val) ) && !isinf( real(val) ) && !isinf( imag(val) ) )  
                                                {
                                                        sum += val;
                                                        count += 1;
                                                }

                                                isEmpty = false;
                                        }
                                }

                                if ( F[g1][g2][g3][g4 + 1] != xZERO )  {

                                        if ( std::abs(F[g1][g2][g3][g4 + 1]) < val_min_abs && F[g1][g2][g3][g4 + 2] != xZERO )  {

                                                val_min_abs = std::abs(F[g1][g2][g3][g4 + 1]);
                                                val_min = F[g1][g2][g3][g4 + 1];
                                        }

                                        if ( F[g1][g2][g3][g4 + 2] != xZERO )  {

                                                val = exp( 2.0 * std::log(F[g1][g2][g3][g4 + 1]) - std::log(F[g1][g2][g3][g4 + 2]) );

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

                        #pragma omp parallel for private(g1, g2, g3, g4)

                        for ( int i = 0; i < ExFF.size(); i++ )  {

                        	g1 = ExFF[i] / M1;
                        	g2 = (ExFF[i] % M1) / M2;
                        	g3 = (ExFF[i] % M2) / M3;
                        	g4 = ExFF[i] % M3;
                                F[g1][g2][g3][g4] = ExTBL[i];
                        }

			t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-b-2: ExFF) = %lf sec\n", t_1_elapsed);  

                        // ............................................................................................. Extrapolation

                        // Check Extending nonzero Area

                        t_1_begin = omp_get_wtime();

                        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4)

                        for (int i = 0; i < ExFF.size(); i++)
                        {
                        	g1 = ExFF[i] / M1;
                        	g2 = (ExFF[i] % M1) / M2;
                        	g3 = (ExFF[i] % M2) / M3;
                        	g4 = ExFF[i] % M3;

                                tmpVec.push_back(GridToIdx( g1    , g2    , g3     , g4     ));
                                tmpVec.push_back(GridToIdx( g1 + 1, g2    , g3     , g4     ));
                                tmpVec.push_back(GridToIdx( g1 - 1, g2    , g3     , g4     ));
                                tmpVec.push_back(GridToIdx( g1,     g2 + 1, g3     , g4     ));
                                tmpVec.push_back(GridToIdx( g1,     g2 - 1, g3     , g4     )); 
                                tmpVec.push_back(GridToIdx( g1,     g2    , g3 + 1 , g4     ));
                                tmpVec.push_back(GridToIdx( g1,     g2    , g3 - 1 , g4     ));     
                                tmpVec.push_back(GridToIdx( g1,     g2    , g3     , g4 + 1 ));
                                tmpVec.push_back(GridToIdx( g1,     g2    , g3     , g4 - 1 ));                        
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
                        log->log("Elapsed time (omp-c-1: CASE 1 TA) = %lf sec\n", t_1_elapsed); 
         
                        // Main iteration
                        t_1_begin = omp_get_wtime();

                        #pragma omp parallel for private(g1, g2, g3, g4, xx1, xx2, xx3, xx4)

                        for (int i = 0; i < TA.size(); i++)  {

	                        g1 = TA[i] / M1;
        	                g2 = (TA[i] % M1) / M2;
                	        g3 = (TA[i] % M2) / M3;
                        	g4 = TA[i] % M3;
                                xx1 = Box[0] + g1 * H[0];
                                xx2 = Box[2] + g2 * H[1];
                                xx3 = Box[4] + g3 * H[2];
                                xx4 = Box[6] + g4 * H[3];

                                FF[g1][g2][g3][g4] = h2m * ( 
                                                     S[0] * ( F[g1 - 1][g2][g3][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1 + 1][g2][g3][g4] ) 
                                                   + S[1] * ( F[g1][g2 - 1][g3][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2 + 1][g3][g4] ) 
                                                   + S[2] * ( F[g1][g2][g3 - 1][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2][g3 + 1][g4] ) 
                                                   + S[3] * ( F[g1][g2][g3][g4 - 1] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2][g3][g4 + 1] ) 
                                                ) 
                                                + F[g1][g2][g3][g4] * ( 1.0 - k2hb * Potential(xx1, xx2, xx3, xx4));
                        }

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-d-1: CASE 1 FF) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                        
                        #pragma omp parallel for private(g1, g2, g3, g4)

                        for (int i = 0; i < ExFF.size(); i++)  {

	                        g1 = ExFF[i] / M1;
        	                g2 = (ExFF[i] % M1) / M2;
                	        g3 = (ExFF[i] % M2) / M3;
                        	g4 = ExFF[i] % M3;

                                PFdX1[g1][g2][g3][g4] = 0.5 * Hi[0] * std::abs( FF[g1 + 1][g2][g3][g4] - FF[g1 - 1][g2][g3][g4]);
                                PFdX2[g1][g2][g3][g4] = 0.5 * Hi[1] * std::abs( FF[g1][g2 + 1][g3][g4] - FF[g1][g2 - 1][g3][g4]);
                                PFdX3[g1][g2][g3][g4] = 0.5 * Hi[2] * std::abs( FF[g1][g2][g3 + 1][g4] - FF[g1][g2][g3 - 1][g4]);
                                PFdX4[g1][g2][g3][g4] = 0.5 * Hi[3] * std::abs( FF[g1][g2][g3][g4 + 1] - FF[g1][g2][g3][g4 - 1]);
			}

			t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-d-2: CASE 1 PF) = %lf sec\n", t_1_elapsed);

                        // check Multiple Expanding 
                        // TBL = index of FF that FF(TBL) is higher than TolL

                        TBL.clear();

                        t_1_begin = omp_get_wtime();

                        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4, b1, b2, b3, b4, b5, b6, b7)

                        for (int i = 0; i < ExFF.size(); i++)
                        {
	                        g1 = ExFF[i] / M1;
        	                g2 = (ExFF[i] % M1) / M2;
                	        g3 = (ExFF[i] % M2) / M3;
                        	g4 = ExFF[i] % M3;

                                b1 = std::abs(FF[g1][g2][g3][g4] * std::conj(FF[g1][g2][g3][g4])) >= TolH;
                                b2 = PFdX1[g1][g2][g3][g4] >= TolHd;
                                b3 = PFdX2[g1][g2][g3][g4] >= TolHd;
                                b4 = PFdX3[g1][g2][g3][g4] >= TolHd;
                                b5 = PFdX4[g1][g2][g3][g4] >= TolHd;
                                b6 = g1 > 1 && g2 > 1 && g3 > 1 && g4 > 1;
                                b7 = g1 < BoxShape[0] - 2 && g2 < BoxShape[1] - 2 && g3 < BoxShape[2] - 2 && g4 < BoxShape[3] - 2;

                                if (  ( b1 || b2 || b3 || b4 || b5 ) && b6 && b7 )  {
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
                        log->log("Elapsed time (omp-d-3 CASE 1 TBL TBL_P) = %lf sec\n", t_1_elapsed); 
                }

                // CASE 2: Truncating without extrapolation

                if ( !isExtrapolate && !isFullGrid )
                {
                        t_1_begin = omp_get_wtime();

                        #pragma omp parallel for private(g1, g2, g3, g4, xx1, xx2, xx3, xx4)

                        for (int i = 0; i < TA.size(); i++)  {

	                        g1 = TA[i] / M1;
        	                g2 = (TA[i] % M1) / M2;
                	        g3 = (TA[i] % M2) / M3;
                        	g4 = TA[i] % M3;

                                xx1 = Box[0] + g1 * H[0];
                                xx2 = Box[2] + g2 * H[1];
                                xx3 = Box[4] + g3 * H[2];
                                xx4 = Box[6] + g4 * H[3];

                                FF[g1][g2][g3][g4] = h2m * ( 
                                                         S[0] * ( F[g1 - 1][g2][g3][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1 + 1][g2][g3][g4] ) 
                                                       + S[1] * ( F[g1][g2 - 1][g3][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2 + 1][g3][g4] ) 
                                                       + S[2] * ( F[g1][g2][g3 - 1][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2][g3 + 1][g4] ) 
                                                       + S[3] * ( F[g1][g2][g3][g4 - 1] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2][g3][g4 + 1] )
                                                       ) 
                                                       + F[g1][g2][g3][g4] * ( 1.0 - k2hb * Potential(xx1, xx2, xx3, xx4));
                        }

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-d-1 FF) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                }   
                else if ( !isExtrapolate && isFullGrid )
                {
                        t_1_begin = omp_get_wtime();
                        // CASE 3: Full grid

                        #pragma omp parallel for private(g1, g2, g3, g4, xx1, xx2, xx3, xx4)

                        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                        g1 = i1;
                                                        g2 = i2;
                                                        g3 = i3;
                                                        g4 = i4;
                                                        xx1 = Box[0] + g1 * H[0];
                                                        xx2 = Box[2] + g2 * H[1];
                                                        xx3 = Box[4] + g3 * H[2];
                                                        xx4 = Box[6] + g4 * H[3];

                                                        FF[g1][g2][g3][g4] = h2m * ( 
                                                                                  S[0] * ( F[g1 - 1][g2][g3][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1 + 1][g2][g3][g4] ) 
                                                                                + S[1] * ( F[g1][g2 - 1][g3][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2 + 1][g3][g4] ) 
                                                                                + S[2] * ( F[g1][g2][g3 - 1][g4] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2][g3 + 1][g4] ) 
                                                                                + S[3] * ( F[g1][g2][g3][g4 - 1] - 2.0 * F[g1][g2][g3][g4] + F[g1][g2][g3][g4 + 1] )
                                                                                ) 
                                                                                + F[g1][g2][g3][g4] * ( 1.0 - k2hb * Potential(xx1, xx2, xx3, xx4));
                                                }
                                        }
                                }
                        }
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-d-1: CASE 3 FF) = %lf sec\n", t_1_elapsed);
                }

                t_1_begin = omp_get_wtime();

                // ff(t+1) Normailzed & go on
                norm = 0.0;

		if (!isFullGrid)  {

                	#pragma omp parallel for private(g1, g2, g3, g4) reduction (+:norm)

                	for (int i = 0; i < TA.size(); i++)  {

                        	g1 = TA[i] / M1;
                        	g2 = (TA[i] % M1) / M2;
                        	g3 = (TA[i] % M2) / M3;
                        	g4 = TA[i] % M3;
                        	norm += std::abs(FF[g1][g2][g3][g4] * std::conj(FF[g1][g2][g3][g4]));
                	}
		}  else  {

			#pragma omp parallel for reduction (+:norm)

                	for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                        	for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                                        for (int i3 = 0; i3 < BoxShape[2]; i3 ++)  {

                                                for (int i4 = 0; i4 < BoxShape[3]; i4 ++)  {

                                	                norm += std::abs(FF[i1][i2][i3][i4] * std::conj(FF[i1][i2][i3][i4]));
                                                }
                                        }
				}
			}
		}

                norm *= H[0] * H[1] * H[2] * H[3];
                norm = 1.0 / sqrt(norm);

                #pragma omp parallel for

                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                        for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                                for (int i3 = 0; i3 < BoxShape[2]; i3 ++)  {

                                        for (int i4 = 0; i4 < BoxShape[3]; i4 ++)  {

                                                FF[i1][i2][i3][i4] = norm * FF[i1][i2][i3][i4];
                                                F[i1][i2][i3][i4] = FF[i1][i2][i3][i4];
                                                PF[i1][i2][i3][i4] = std::abs(F[i1][i2][i3][i4] * std::conj(F[i1][i2][i3][i4]));
                                        }
                                }
                        }
                }       

		t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                log->log("Elapsed time (omp-e-1 FF) = %lf sec\n", t_1_elapsed); 

                // Truncated_New Edge
                if ( !isFullGrid )
                {
                        t_1_begin = omp_get_wtime();

                        // PFdX1

                        coeff1 = 1.0 / (2.0 * H[0]);
                        coeff2 = 1.0 / (2.0 * H[1]);
                        coeff3 = 1.0 / (2.0 * H[2]);
                        coeff4 = 1.0 / (2.0 * H[3]);

                        #pragma omp parallel for

                        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                        PFdX1[i1][i2][i3][i4] = std::abs(F[i1+1][i2][i3][i4] - F[i1-1][i2][i3][i4]) * coeff1;
							PFdX2[i1][i2][i3][i4] = std::abs(F[i1][i2+1][i3][i4] - F[i1][i2-1][i3][i4]) * coeff2;
							PFdX3[i1][i2][i3][i4] = std::abs(F[i1][i2][i3+1][i4] - F[i1][i2][i3-1][i4]) * coeff3;
							PFdX4[i1][i2][i3][i4] = std::abs(F[i1][i2][i3][i4+1] - F[i1][i2][i3][i4-1]) * coeff4;
                                                }
                                        }
                                }
                        }

                        // PFdX2

                        /*coeff = 1.0 / (2.0 * H[1]);

                        #pragma omp parallel for

                        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                        PFdX2[i1][i2][i3][i4] = std::abs(F[i1][i2+1][i3][i4] - F[i1][i2-1][i3][i4]) * coeff;
                                                }
                                        }
                                }
                        }*/

                        // PFdX3

                        /*coeff = 1.0 / (2.0 * H[2]);

                        #pragma omp parallel for

                        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                        PFdX3[i1][i2][i3][i4] = std::abs(F[i1][i2][i3+1][i4] - F[i1][i2][i3-1][i4]) * coeff;
                                                }
                                        }
                                }
                        }*/

                        // PFdX4

                        /*coeff = 1.0 / (2.0 * H[3]);

                        #pragma omp parallel for

                        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                        PFdX4[i1][i2][i3][i4] = std::abs(F[i1][i2][i3][i4+1] - F[i1][i2][i3][i4-1]) * coeff;
                                                }
                                        }
                                }
                        }*/

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-2 PF) = %lf sec\n", t_1_elapsed); 

                        // Truncate

                        t_1_begin = omp_get_wtime();

                        #pragma omp parallel for

                        for (int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1]  - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2]  - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3]  - 1; i4 ++)  {
                        
                                                        if ( ( PF[i1][i2][i3][i4] < TolH ) && ( PFdX1[i1][i2][i3][i4] < TolHd ) && ( PFdX2[i1][i2][i3][i4] < TolHd ) && ( PFdX3[i1][i2][i3][i4] < TolHd ) && ( PFdX4[i1][i2][i3][i4] < TolHd ) )  {
                                                                F[i1][i2][i3][i4] = xZERO;
                                                        }
                                                }
                                        }
                                }
                        }

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-1 PF) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();

                        // TA

                        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4, b1, b2)

                        for ( int i = 0; i < TA.size(); i++ )  {

	                        g1 = TA[i] / M1;
        	                g2 = (TA[i] % M1) / M2;
                	        g3 = (TA[i] % M2) / M3;
                        	g4 = TA[i] % M3;

                                b1 = F[g1][g2][g3][g4] != xZERO;
                                b2 = PFdX1[g1][g2][g3][g4] >= TolHd || PFdX2[g1][g2][g3][g4] >= TolHd || PFdX3[g1][g2][g3][g4] >= TolHd || PFdX4[g1][g2][g3][g4] >= TolHd;

                                if ( b1 || b2 )  {

                                        tmpVec.push_back(TA[i]);
                                }
                        }
                        tmpVec.swap(TA);  
                        tmpVec.clear();

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-2 TA) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();

                        // TB

                        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4, b1, b2, b3, b4, b5, b6, b7, b8)

                        for ( int i = 0; i < TA.size(); i++ )  {

	                        g1 = TA[i] / M1;
        	                g2 = (TA[i] % M1) / M2;
                	        g3 = (TA[i] % M2) / M3;
                        	g4 = TA[i] % M3;

                                b1 = F[ g1 - 1 ][ g2     ][ g3     ][ g4     ] == xZERO;
                                b2 = F[ g1 + 1 ][ g2     ][ g3     ][ g4     ] == xZERO;
                                b3 = F[ g1     ][ g2 - 1 ][ g3     ][ g4     ] == xZERO;
                                b4 = F[ g1     ][ g2 + 1 ][ g3     ][ g4     ] == xZERO;
                                b5 = F[ g1     ][ g2     ][ g3 - 1 ][ g4     ] == xZERO;
                                b6 = F[ g1     ][ g2     ][ g3 + 1 ][ g4     ] == xZERO;
                                b7 = F[ g1     ][ g2     ][ g3     ][ g4 - 1 ] == xZERO;
                                b8 = F[ g1     ][ g2     ][ g3     ][ g4 + 1 ] == xZERO;

                                if ( b1 || b2 || b3 || b4 || b5 || b6 || b7 || b8 )
                                {
                                        b1 = ( PFdX1[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1 - 1 ][ g2     ][ g3     ][ g4     ] < TolHd );
                                        b2 = ( PFdX1[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1 + 1 ][ g2     ][ g3     ][ g4     ] < TolHd );
                                        b3 = ( PFdX1[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2 - 1 ][ g3     ][ g4     ] < TolHd );
                                        b4 = ( PFdX1[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2 + 1 ][ g3     ][ g4     ] < TolHd );
                                        b5 = ( PFdX1[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3 - 1 ][ g4     ] < TolHd );
                                        b6 = ( PFdX1[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3 + 1 ][ g4     ] < TolHd );
                                        b7 = ( PFdX1[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3     ][ g4 - 1 ] < TolHd );
                                        b8 = ( PFdX1[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd ) && ( PFdX2[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd ) && ( PFdX3[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd ) && ( PFdX4[ g1     ][ g2     ][ g3     ][ g4 + 1 ] < TolHd );

                                        if ( b1 || b2 || b3 || b4 || b5 || b6 || b7 || b8 ) {

                                                tmpVec.push_back(TA[i]);  
                                        }
                                }
                        }
                        tmpVec.swap(TB);
                        tmpVec.clear();

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-3 TB) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();

                        // TA expansion
                        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4)                               
                        for (int i = 0; i < TA.size(); i++)
                        {
                        	g1 = TA[i] / M1;
                        	g2 = (TA[i] % M1) / M2;
                        	g3 = (TA[i] % M2) / M3;
                       		g4 = TA[i] % M3;

				if (g1 + 1 != BoxShape[0] - 1)
                                	tmpVec.push_back(GridToIdx(g1 + 1, g2, g3, g4));

				if (g1 - 1 != 0)
                                	tmpVec.push_back(GridToIdx(g1 - 1, g2, g3, g4));

				if (g2 + 1 != BoxShape[1] - 1)
                                	tmpVec.push_back(GridToIdx(g1, g2 + 1, g3, g4));

				if (g2 - 1 != 0)
                                	tmpVec.push_back(GridToIdx(g1, g2 - 1, g3, g4)); 

				if (g3 + 1 != BoxShape[2] - 1)
                                	tmpVec.push_back(GridToIdx(g1, g2, g3 + 1, g4));

				if (g3 - 1 != 0)
                                	tmpVec.push_back(GridToIdx(g1, g2, g3 - 1, g4));  

				if (g4 + 1 != BoxShape[3] - 1)
                                	tmpVec.push_back(GridToIdx(g1, g2, g3, g4 + 1));

				if (g4 - 1 != 0)
                                	tmpVec.push_back(GridToIdx(g1, g2, g3, g4 - 1));                     
                        }

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-4-1 push_back) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();

                        // Combine TA and tmpVec
                        TA.reserve(TA.size() + tmpVec.size());
                        TA.insert(TA.end(), tmpVec.begin(), tmpVec.end());
                        tmpVec.clear();

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-4-2 combine) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();

                        // Find unique elements
                        __gnu_parallel::sort(TA.begin(),TA.end());

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-4-2-1 sort) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();

                        it = std::unique (TA.begin(), TA.end()); 

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-4-2-2 unique) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();

                        TA.resize(std::distance(TA.begin(),it));

                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        log->log("Elapsed time (omp-e-3-4-2-3 resize) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                        
			#pragma omp parallel for

                        for (int i1 = 1; i1 < BoxShape[0] - 1 ; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1]  - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2]  - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3]  - 1; i4 ++)  {

					                PF[i1][i2][i3][i4] = std::abs( F[i1][i2][i3][i4] * std::conj(F[i1][i2][i3][i4]) );
                                                }
                                        }
                                }
                        }
                }

                if ( (tt + 1) % PERIOD == 0 )
                {
			// Compute <E>
			// --------------------------------------------------------------
			
                        t_1_begin = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;

			energy = 0.0;

                        #pragma omp parallel for reduction (+:energy)

                        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                        energy += std::real( std::conj(F[i1][i2][i3][i4]) * ( 

                                                                        -hsq2m * ( 
                                                                                Hisq[0] * ( F[i1-1][i2][i3][i4] - 2.0 * F[i1][i2][i3][i4] + F[i1+1][i2][i3][i4] ) 
                                                                                + Hisq[1] * ( F[i1][i2-1][i3][i4] - 2.0 * F[i1][i2][i3][i4] + F[i1][i2+1][i3][i4] ) 
                                                                                + Hisq[2] * ( F[i1][i2][i3-1][i4] - 2.0 * F[i1][i2][i3][i4] + F[i1][i2][i3+1][i4] ) 
                                                                                + Hisq[3] * ( F[i1][i2][i3][i4-1] - 2.0 * F[i1][i2][i3][i4] + F[i1][i2][i3][i4+1] )
                                                                        ) + Potential(Box[0] + i1 * H[0], Box[2] + i2 * H[1], Box[4] + i3 * H[2], Box[6] + i4 * H[3]) * F[i1][i2][i3][i4])       
                                                        );
                                                }
                                        }
                                }
                        }

			t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;

			log->log("Elapsed time (omp-f-1 resize) = %lf sec\n", t_1_elapsed);

                        t_1_begin = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;

                	norm = 0.0;

                        #pragma omp parallel for reduction (+:norm)

                        for (int i1 = 1; i1 < BoxShape[0] - 1; i1 ++)  {

                                for (int i2 = 1; i2 < BoxShape[1] - 1; i2 ++)  {

                                        for (int i3 = 1; i3 < BoxShape[2] - 1; i3 ++)  {

                                                for (int i4 = 1; i4 < BoxShape[3] - 1; i4 ++)  {

                                                        norm += std::abs(F[i1][i2][i3][i4] * std::conj(F[i1][i2][i3][i4]));
                                                }
                                        }
                                }
                        }

			t_1_end = omp_get_wtime();
			t_1_elapsed = t_1_end - t_1_begin;

			log->log("Elapsed time (omp-f-2 resize) = %lf sec\n", t_1_elapsed);

			energy /= norm;

			log->log("[Imt4d] Time %lf, <E> = %e, norm = %e\n", ( tt + 1 ) * kk, energy, norm);

			// --------------------------------------------------------------
			
                        t_0_end = omp_get_wtime();
                        t_0_elapsed = t_0_end - t_0_begin;
 
                        log->log("[Imt4d] Step: %d, Elapsed time: %lf sec\n", tt + 1, t_0_elapsed);

			if (!isFullGrid)  {

                        	log->log("[Imt4d] TA size = %d, TB size = %d\n", TA.size(), TB.size());
                        	log->log("[Imt4d] TA / total grids = %lf\n", ( TA.size() * 1.0 ) / GRIDS_TOT);
			}
                	log->log("\n........................................................\n\n");
                }                     

        } // Time iteration 

        log->log("[Imt4d] Evolve done.\n");
}
/* ------------------------------------------------------------------------------- */

inline std::complex<double> Imt4d::Wavefunction(double x1, double x2, double x3, double x4)
{       
        // HO-MO
        return exp( -A[0] * ( pow(x1 - Wave0[0], 2) + pow(x2 - Wave0[1], 2) + pow(x3 - Wave0[2], 2) + pow(x4 - Wave0[3], 2) ) );
}

/* ------------------------------------------------------------------------------- */

inline double Imt4d::Potential(double x1, double x2, double x3, double x4)
{       
	return 0.5 * mw2 * (x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4);
        // HO-MO
        //return 0.5 * mw2 * x1 * x1 + De * ( exp(-2.0 * Da * x2) - 2.0 * exp(-Da * x2) ) + De * ( exp(-2.0 * Da * x3) - 2.0 * exp(-Da * x3) ) + De * ( exp(-2.0 * Da * x4) - 2.0 * exp(-Da * x4) );
}

/* ------------------------------------------------------------------------------- */

VectorXi Imt4d::IdxToGrid(int idx)
{
        int x4 = idx % M3;
        int x3 = (idx % M2) / M3;
        int x2 = (idx % M1) / M2;
        int x1 = idx / M1;

        VectorXi grid;
        grid.resize(DIMENSIONS);
        grid << x1, x2, x3, x4;

        return grid;
}
/* ------------------------------------------------------------------------------- */

inline int Imt4d::GridToIdx(int x1, int x2, int x3, int x4)
{
        return x1 * M1 + x2 * M2 + x3 * M3 + x4;
}
/* ------------------------------------------------------------------------------- */

inline void Imt4d::DefineBoundary()
{
        // Find elements of DBi and DBi2

        // i1 
        for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                for (int i3 = 0; i3 < BoxShape[2]; i3 ++) { 

                        for (int i4 = 0; i4 < BoxShape[3]; i4 ++) { 

                                DBi.push_back(GridToIdx(0, i2, i3, i4));
                                DBi.push_back(GridToIdx(BoxShape[0]-1, i2, i3, i4));

                                DBi2.push_back(GridToIdx(0, i2, i3, i4));
                                DBi2.push_back(GridToIdx(1, i2, i3, i4));
                                DBi2.push_back(GridToIdx(2, i2, i3, i4));
                                DBi2.push_back(GridToIdx(BoxShape[0]-1, i2, i3, i4));
                                DBi2.push_back(GridToIdx(BoxShape[0]-2, i2, i3, i4));
                                DBi2.push_back(GridToIdx(BoxShape[0]-3, i2, i3, i4));
                        }
                }
        }

        // i2

        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                for (int i3 = 0; i3 < BoxShape[2]; i3 ++)  {

                        for (int i4 = 0; i4 < BoxShape[3]; i4 ++) { 

                                DBi.push_back(GridToIdx(i1, 0, i3, i4));
                                DBi.push_back(GridToIdx(i1, BoxShape[1]-1, i3, i4));

                                DBi2.push_back(GridToIdx(i1, 0, i3, i4));
                                DBi2.push_back(GridToIdx(i1, 1, i3, i4));
                                DBi2.push_back(GridToIdx(i1, 2, i3, i4));
                                DBi2.push_back(GridToIdx(i1, BoxShape[1]-1, i3, i4));
                                DBi2.push_back(GridToIdx(i1, BoxShape[1]-2, i3, i4));
                                DBi2.push_back(GridToIdx(i1, BoxShape[1]-3, i3, i4));
                        }
                }
        }

        // i3

        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                        for (int i4 = 0; i4 < BoxShape[3]; i4 ++) { 

                                DBi.push_back(GridToIdx(i1, i2, 0, i4));
                                DBi.push_back(GridToIdx(i1, i2, BoxShape[2]-1, i4));

                                DBi2.push_back(GridToIdx(i1, i2, 0, i4));
                                DBi2.push_back(GridToIdx(i1, i2, 1, i4));
                                DBi2.push_back(GridToIdx(i1, i2, 2, i4));
                                DBi2.push_back(GridToIdx(i1, i2, BoxShape[2]-1, i4));
                                DBi2.push_back(GridToIdx(i1, i2, BoxShape[2]-2, i4));
                                DBi2.push_back(GridToIdx(i1, i2, BoxShape[2]-3, i4));
                        }
                }
        }

        // i4

        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {

                        for (int i3 = 0; i3 < BoxShape[2]; i3 ++) { 

                                DBi.push_back(GridToIdx(i1, i2, i3, 0));
                                DBi.push_back(GridToIdx(i1, i2, i3, BoxShape[3]-1));

                                DBi2.push_back(GridToIdx(i1, i2, i3, 0));
                                DBi2.push_back(GridToIdx(i1, i2, i3, 1));
                                DBi2.push_back(GridToIdx(i1, i2, i3, 2));
                                DBi2.push_back(GridToIdx(i1, i2, i3, BoxShape[3]-1));
                                DBi2.push_back(GridToIdx(i1, i2, i3, BoxShape[3]-2));
                                DBi2.push_back(GridToIdx(i1, i2, i3, BoxShape[3]-3));
                        }
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

