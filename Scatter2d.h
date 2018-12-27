// ==============================================================================
//
//  Scatter2d.h
//  QTR
//
//  Created by Albert Lu on 10/27/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 12/26/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_SCATTER2D_H
#define QTR_SCATTER2D_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
    
    class Scatter2d {
        
    public:
        Scatter2d(class QTR *q);
        ~Scatter2d();
  
        void                          Evolve();
        VectorXi                      IdxToGrid(int idx);
        inline int                    GridToIdx(int x1, int x2);
        inline std::complex<double>   WAVEFUNCTION(double x1, double x2);
        inline std::complex<double>   Wavefunction_EckMO(double x1, double x2);
        inline std::complex<double>   Wavefunction_EckHO(double x1, double x2);
        inline std::complex<double>   Wavefunction_GauMO(double x1, double x2);
        inline std::complex<double>   Wavefunction_GauHO(double x1, double x2);
        inline std::complex<double>   Wavefunction_HH(double x1, double x2);
        inline void                   DefineBoundary();
        inline double                 POTENTIAL(double x1, double x2);
        inline double                 Potential_EckMO(double x1, double x2);
        inline double                 Potential_EckHO(double x1, double x2);
        inline double                 Potential_GauMO(double x1, double x2);
        inline double                 Potential_GauHO(double x1, double x2);
        inline double                 Potential_HH(double x1, double x2);

    private:

        void            init();
        QTR             *qtr;
        Error           *err;
        Log             *log;
        Parameters      *parameters;

        // General parameters
        std::complex<double>  I;      // sqrt(-1)
        std::complex<double>  xZERO;  // complex zero
        int             DIMENSIONS;
        int             PERIOD; 
        int             GRIDS_TOT;
        bool            QUIET;
        bool            TIMING;
        double          TIME;   
        double          PI_INV;  // 1/pi

        // Grid size
        double          kk;    // time resolution
        VectorXd        H;     // grid size
        VectorXd        Hi;    // inverse grid size
        VectorXd        Hisq;  // inverse grid size square     
        VectorXd        S;  

        // Domain size
        VectorXd        Box;
        VectorXi        BoxShape;
        int             M1;

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
        double          sigma;  // EckHO
        double          beta;   // GauHO

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

        // Output
        bool            isTrans;
        bool            isAcf;
    };
}

#endif /* QTR_SCATTER2D_H */
