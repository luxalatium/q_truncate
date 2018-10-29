// ==============================================================================
//
//  Job_Test.cpp
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

#include <string>

#include "Job_Test.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobTest::JobTest(class QTR *qtr)
{
    log = qtr->log;
    parameters = qtr->parameters;
}
/* ------------------------------------------------------------------------------- */

JobTest::~JobTest(){ }

/* ------------------------------------------------------------------------------- */

void JobTest::run(class QTR *qtr)
{     
    log->log("[Job_Test] test done - %d.", 1);
}
/* ------------------------------------------------------------------------------- */

