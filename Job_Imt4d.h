// ==============================================================================
//
//  Job_Imt4d.h
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

#ifndef QTR_JOB_IMT4D_H
#define QTR_JOB_IMT4D_H

#include "Job.h"
#include "Imt4d.h"

namespace QTR_NS  {

        class JobImt4d: public Job      {
                
        public:
                JobImt4d(class QTR *);
                ~JobImt4d(void);
                
                void run(class QTR *);
                
        private:
                
                Log        *log;
                Parameters *parameters;
                QTR        *qtr;
                Imt4d      *imt;
        };
}
#endif /* QTR_JOB_IMT4D_H */
