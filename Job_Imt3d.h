// ==============================================================================
//
//  Job_Imt3d.h
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

#ifndef QTR_JOB_IMT3D_H
#define QTR_JOB_IMT3D_H

#include "Job.h"
#include "Imt3d.h"

namespace QTR_NS  {

    class JobImt3d: public Job    {
        
    public:
        JobImt3d(class QTR *);
        ~JobImt3d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        Imt3d      *imt;
    };
}
#endif /* QTR_JOB_IMT3D_H */
