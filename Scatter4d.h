// ==============================================================================
//
//  Scatter4d.h
//  QTR
//
//  Created by Albert Lu on 9/9/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/27/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_SCATTER4D_H
#define QTR_SCATTER4D_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
    
    class Scatter4d {
        
    public:
        Scatter4d(class QTR *q);
        ~Scatter4d();
  
        void                          Evolve();
        VectorXi                      IdxToGrid(int idx);
        inline int                    GridToIdx(int x1, int x2, int x3, int x4);
        inline std::complex<double>   Wavefunction_MO(double x1, double x2, double x3, double x4);
        inline std::complex<double>   Wavefunction_HH(double x1, double x2, double x3, double x4);
        inline void                   DefineBoundary();
        inline double                 Potential_MO(double x1, double x2, double x3, double x4);
        inline double                 Potential_HH(double x1, double x2, double x3, double x4);
        
    private:

        void            init();
        QTR             *qtr;
        Error           *err;
        Log             *log;
        Parameters      *parameters;

        // General parameters
        std::complex<double>  I;      // sqrt(-1)
        std::complex<double>  xZERO;  // complex zero
        int                   DIMENSIONS;
        int                   PERIOD; 
	    int                   GRIDS_TOT;
        double                TIME;
        double                PI_INV;  // 1/pi

        // Grid size
        double          kk;  // time resolution
        VectorXd        H;   // grid size
        VectorXd        Hi;  // inverse grid size
        VectorXd        Hisq;  // inverse grid size square     
        VectorXd        S;  

        // Domain size
        VectorXd        Box;
        VectorXi        BoxShape;
        int             M1;
        int             M2;
        int             M3;

        // Potential parameters
        int             idx_x0;
        VectorXi        Vmode;        
        double          hb;
        double          m;
        double          w;    // HO specific
        double          aa;
        double          V0;   // Eckart potential 
        double          ek2v;
        double          alpha;    
        double          Ek0;      
        double          k0;   // Related HO
        double          sig;        
        double          lan;         
        double          De;   // Morse
        double          Da;
        double          r0;
        double          Ld;
        double          lambda; // HH

        // Wavefunction
        VectorXd        Wave0;
        VectorXd        A;
        VectorXd        P;

        // Truncate parameters
        bool            isEmpty;
        bool            isFullGrid; 
        bool            isExtrapolate;         
        double          TolH;
        double          TolL;
        double          TolHd;
        double          TolLd;
        double          ExReduce;

        // Domains
        MeshIndex       TA;
        MeshIndex       TB;    // Truncation boundary
        MeshIndex       TBL;
        MeshIndex       TBL_P;          
        MeshIndex       DBi;   // Grid boundary
        MeshIndex       DBi2;  // Extrapolation-restricted area
        MeshIndex       ExFF;
        MeshIndex       ExFF2;
    };
}

#endif /* QTR_SCATTER4D_H */
