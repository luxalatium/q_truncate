// ==============================================================================
//
//  Job_Imt3d.cpp
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

#include "Job_Imt3d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Imt3d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobImt3d::JobImt3d(class QTR *q)
{
        qtr = q;
        log = qtr->log;
        parameters = qtr->parameters;
        imt = new Imt3d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobImt3d::~JobImt3d()
{
        delete imt;
}

/* ------------------------------------------------------------------------------- */

void JobImt3d::run(class QTR *qtr)
{       
        imt->Evolve();
        log->log("[Job_Imt3d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

