// ==============================================================================
//
//  Imt3d.h
//  QTR
//
//  Created by Albert Lu on 10/3/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/3/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_IMT3D_H
#define QTR_IMT3D_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
        
        class Imt3d {
                
        public:
                Imt3d(class QTR *q);
                ~Imt3d();
  
                void                          Evolve();
                VectorXi                      IdxToGrid(int idx);
                inline int                    GridToIdx(int x1, int x2, int x3);
                inline std::complex<double>   Wavefunction(double x1, double x2, double x3);
                inline void                   DefineBoundary();
                inline double                 Potential(double x1, double x2, double x3);
                
        private:

                void                  init();
                
                QTR                   *qtr;
                Error                 *err;
                Log                   *log;
                Parameters            *parameters;

                // General parameters
                std::complex<double>  I;        // sqrt(-1)
                std::complex<double>  xZERO;    // complex zero
                int                   DIMENSIONS;
                int                   PERIOD; 
		int                   GRIDS_TOT;
                double                TIME;   

                // Grid size
                double                kk;    // time resolution
                VectorXd              H;     // grid size
                VectorXd              Hi;    // inverse grid size
                VectorXd              Hisq;  // inverse grid size square         
                VectorXd              S;    

                // Domain size
                VectorXd              Box;
                VectorXi              BoxShape;
                int                   M1;
                int                   M2;

                // Potential parameters
                //int                   idx_x0;
                //VectorXi              Vmode;                
                double                hb;
                double                m;
                double                w;      // HO specific
                double                mw2;
                double                aa;
                double                V0;     // Eckart potential 
                double                ek2v;
                double                alpha;        
                double                Ek0;          
                double                k0;     // Related HO
                double                sig;              
                double                lan;                 
                double                De;     // Morse
                double                Da;
                double                r0;
                double                Ld;

                // Wavefunction
                VectorXd              Wave0;
                VectorXd              A;
                VectorXd              P;

                // Truncate parameters
                bool                  isEmpty;
                bool                  isFullGrid; 
                bool                  isExtrapolate;               
                double                TolH;
                double                TolL;
                double                TolHd;
                double                TolLd;
                double                ExReduce;

                // Domains
                MeshIndex             TA;
                MeshIndex             TB;      // Truncation boundary
                MeshIndex             TBL;
                MeshIndex             TBL_P;                  
                MeshIndex             DBi;     // Grid boundary
                MeshIndex             DBi2;    // Extrapolation-restricted area
                MeshIndex             ExFF;
                MeshIndex             ExFF2;
        };
}

#endif /* QTR_IMT3D_H */
