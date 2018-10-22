// ==============================================================================
//
//  Job_Test.h
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

#ifndef QTR_JOB_TEST_H
#define QTR_JOB_TEST_H

#include "Job.h"

namespace QTR_NS  {

        class JobTest: public Job      {
                
        public:
                JobTest(class QTR *);
                ~JobTest(void);
                
                void run(class QTR *);   
        private:
                Log        *log;
                Parameters *parameters;
        };
}
#endif /* QTR_JOB_TEST_H */
