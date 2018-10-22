// ==============================================================================
//
//  Job_Scatter4d.h
//  QTR
//
//  Created by Albert Lu on 9/10/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/10/18
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_SCATTER4D_H
#define QTR_JOB_SCATTER4D_H

#include "Job.h"
#include "Scatter4d.h"

namespace QTR_NS  {

        class JobScatter4d: public Job      {
                
        public:
                JobScatter4d(class QTR *);
                ~JobScatter4d(void);
                
                void run(class QTR *);
                
        private:
                
                Log        *log;
                Parameters *parameters;
                QTR        *qtr;
                Scatter4d  *scat;
        };
}
#endif /* QTR_JOB_SCATTER4D_H */
