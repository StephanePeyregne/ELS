#include <cstdlib>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

#include "obsData.h"
#include "modelprob.h"
#include "hmmresults.h"


// These should be parameters:
//#define EXTERNAL_ULINDI_DERIVED_PROB 0.01

using namespace std ;

class hmm
{
	public:
                hmm( obsSequence *seq) ;

                void computeFwdBwd(modelProb *prob, hmmResults *result);
		void computeFwdBwd3states(modelProb *prob, hmmResults *result);
		void writeLogFile(const char* fileName, string logInfo);
		void writeLogFile3states(const char* fileName, string logInfo);
                void writeOutputFile(void);
		void writeOutputFile3states(void);

                double prob_obs_internal( obsSite dat, probSegregDInternal prob );
                double prob_obs_external( obsSite dat );
		
	protected:

                obsSequence *mObsSequence;
                modelProb   *mModelProb;
                hmmResults  *mHmmResults;

                //temporary hmm results
//                vector<double> mfwd_stateE_scaled;
//                vector<double> mfwd_stateI_scaled;
//                vector<double> mbwd_stateI_scaled;
//                vector<double> mbwd_stateE_scaled;
//                vector<bool> mExternal;

//                double mlogLikelihood;

} ;


