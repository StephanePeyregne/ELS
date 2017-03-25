#ifndef BAUMWELCH_H
#define BAUMWELCH_H

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

using namespace std ;

class baumWelch
{
    public:
        baumWelch(obsSequence *data);
        void computeBaumWelch(hmmResults *results, modelProb *model );
	void computeBaumWelch3states(hmmResults *results, modelProb *model );

    protected:
            obsSequence *mObsSequence;
            hmmResults  *mHmmResults;
            modelProb *mModelProb;
};

#endif // BAUMWELCH_H

