// ==============================================================================
//
//  Job_Imt4d.cpp
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

#include <string>

#include "Job_Imt4d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Imt4d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobImt4d::JobImt4d(class QTR *q)
{
        qtr = q;
        log = qtr->log;
        parameters = qtr->parameters;
        imt = new Imt4d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobImt4d::~JobImt4d()
{
        delete imt;
}

/* ------------------------------------------------------------------------------- */

void JobImt4d::run(class QTR *qtr)
{       
        imt->Evolve();
        log->log("[Job_Imt4d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

