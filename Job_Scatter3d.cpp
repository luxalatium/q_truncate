// ==============================================================================
//
//  Job_Scatter3d.cpp
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 8/4/18
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Scatter3d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Scatter3d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobScatter3d::JobScatter3d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    scat = new Scatter3d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobScatter3d::~JobScatter3d()
{
    delete scat;
}

/* ------------------------------------------------------------------------------- */

void JobScatter3d::run(class QTR *qtr)
{     
    scat->Evolve();
    log->log("[Job_Scatter3d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

