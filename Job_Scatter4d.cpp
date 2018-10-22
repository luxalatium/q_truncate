// ==============================================================================
//
//  Job_Scatter4d.cpp
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

#include <string>

#include "Job_Scatter4d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Scatter4d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobScatter4d::JobScatter4d(class QTR *q)
{
        qtr = q;
        log = qtr->log;
        parameters = qtr->parameters;
        scat = new Scatter4d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobScatter4d::~JobScatter4d()
{
        delete scat;
}

/* ------------------------------------------------------------------------------- */

void JobScatter4d::run(class QTR *qtr)
{       
        scat->Evolve();
        log->log("[Job_Scatter4d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

