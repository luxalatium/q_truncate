// ==============================================================================
//
//  Job.h
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/27/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_H
#define QTR_JOB_H

#include "Qtr.h"

namespace QTR_NS {
    
    class Job {
        
    public:
        
        virtual ~Job() {};
        
        virtual void run(class QTR *) = 0;
        
        static Job *getJob(class QTR *);

        static const char IMT2D[];
        static const char IMT3D[];
        static const char IMT4D[];
        static const char SCATTER2D[];
        static const char SCATTER3D[];
        static const char SCATTER4D[];
        static const char TEST[];
    };
}
#endif /* QTR_JOB_H */
