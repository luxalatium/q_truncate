// ==============================================================================
//
//  Job_Imt2d.cpp
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

#include <string>

#include "Job_Imt2d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Imt2d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobImt2d::JobImt2d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    imt = new Imt2d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobImt2d::~JobImt2d()
{
    delete imt;
}

/* ------------------------------------------------------------------------------- */

void JobImt2d::run(class QTR *qtr)
{     
    imt->Evolve();
    log->log("[Job_Imt2d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

