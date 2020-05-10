#ifndef BAUMWELCH_H
#define BAUMWELCH_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "hmmresults.h"
#include "modelprob.h"
#include "obsData.h"

using namespace std;

class baumWelch
{
  public:
    baumWelch(obsSequence *data);
    void computeBaumWelch(hmmResults *results, modelProb *model);
    void computeBaumWelch3states(hmmResults *results, modelProb *model);

  protected:
    obsSequence *mObsSequence;
    hmmResults *mHmmResults;
    modelProb *mModelProb;
};

#endif // BAUMWELCH_H
