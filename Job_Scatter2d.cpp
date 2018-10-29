// ==============================================================================
//
//  Job_Scatter2d.cpp
//  QTR
//
//  Created by Albert Lu on 10/27.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/27/18
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Scatter2d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Scatter2d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobScatter2d::JobScatter2d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    scat = new Scatter2d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobScatter2d::~JobScatter2d()
{
    delete scat;
}

/* ------------------------------------------------------------------------------- */

void JobScatter2d::run(class QTR *qtr)
{     
    scat->Evolve();
    log->log("[Job_Scatter2d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

