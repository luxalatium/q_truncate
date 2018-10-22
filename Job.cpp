// ==============================================================================
//
//  Job.cpp
//  QTR
//
//  Created by Albert Lu on 8/4/18.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/3/18
//
//  Note:
//
// ==============================================================================

# include "Job.h"
# include "Parameters.h"

# include "Job_Test.h"
# include "Job_Imt2d.h"
# include "Job_Imt3d.h"
# include "Job_Imt4d.h"
# include "Job_Scatter3d.h"
# include "Job_Scatter4d.h"

using namespace QTR_NS;

const char Job::TEST[] = "test";
const char Job::IMT2D[] = "imt2d";
const char Job::IMT3D[] = "imt3d";
const char Job::IMT4D[] = "imt4d";
const char Job::SCATTER3D[] = "scatter3d";
const char Job::SCATTER4D[] = "scatter4d";

/* ------------------------------------------------------------------------------- */

Job *Job::getJob(class QTR *qtr) {
        
        Job *job = NULL;
        
        if (qtr->parameters->job == TEST)
        {
                job = new JobTest(qtr);
        }
        else if (qtr->parameters->job == IMT2D)
        {
                job = new JobImt2d(qtr);
        }
        else if (qtr->parameters->job == IMT3D)
        {
                job = new JobImt3d(qtr);
        }
        else if (qtr->parameters->job == IMT4D)
        {
                job = new JobImt4d(qtr);
        }
        else if (qtr->parameters->job == SCATTER3D)
        {
                job = new JobScatter3d(qtr);
        }
        else if (qtr->parameters->job == SCATTER4D)
        {
                job = new JobScatter4d(qtr);
        }
        return job;
}
/* ----------------------------------------------------------------- */
