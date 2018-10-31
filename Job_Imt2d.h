// ==============================================================================
//
//  Job_Imt2d.h
//  QTR
//
//  Created by Albert Lu on 9/24/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/24/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_IMT2D_H
#define QTR_JOB_IMT2D_H

#include "Job.h"
#include "Imt2d.h"

namespace QTR_NS  {

    class JobImt2d: public Job    {
        
    public:
        JobImt2d(class QTR *);
        ~JobImt2d(void);
        
        void run(class QTR *);
        
    private:
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        Imt2d      *imt;
    };
}
#endif /* QTR_JOB_IMT2D_H */
