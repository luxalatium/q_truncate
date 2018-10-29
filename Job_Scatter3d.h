// ==============================================================================
//
//  Job_Scatter3d.h
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 8/6/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_SCATTER3D_H
#define QTR_JOB_SCATTER3D_H

#include "Job.h"
#include "Scatter3d.h"

namespace QTR_NS  {

    class JobScatter3d: public Job    {
        
    public:
        JobScatter3d(class QTR *);
        ~JobScatter3d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        Scatter3d  *scat;
    };
}
#endif /* QTR_JOB_SCATTER3D_H */
