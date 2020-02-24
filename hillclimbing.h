#ifndef HILLCLIMBING_H
#define HILLCLIMBING_H

#include "baumwelch.h"
#include "hmm.h"
#include <nlopt.hpp>

// Random search default values
#define MAX_ITERATIONS 10000
#define WINDOW_INT 4000
#define WINDOW_EXT 2000
#define WINDOW_LONGEXT 4000
#define MAX_WINDOW_SHRINKAGE 500
#define WINDOW_LERATE 0.2

#define DEFAULT_OPT_TOLERANCE 0.0001

enum
{
    IDX_STAY_INTERNAL = 0,
    IDX_STAY_EXTERNAL,
    IDX_STAY_LONG_INTERNAL,
    IDX_LERATE
};

class hillClimbing
{
  private:
    nlopt::opt *mOpt;

    obsSequence *mObsSequence;

    modelProb *mModelProb;
    modelProb *mModelProbTemp;

    hmmResults *mHmmResults;
    hmmResults *mHmmResultsTemp;

    hmm *mHmm;
    baumWelch *mBaumWelch;

    double mConverThrld;
    double mOptTolerance; // Tolerance de fin d'optimisation
    const char *logFileName;

    void optimisationBaumWelch(modelProb *model, hmmResults *results, int max_iteration);

  public:
    // hillClimbing();
    hillClimbing(const char *sequenceFile,
                 const char *configFile,
                 double lengthInternal,
                 double lengthExternal,
                 double probFixDInternal,
                 double probFixDExternal,
                 double externalUlindiDerivedProb);
    hillClimbing(const char *sequenceFile,
                 const char *configFile,
                 double lengthInternal,
                 double lengthExternal,
                 double lengthLongExternal,
                 double LongExternalrate,
                 double probFixDInternal,
                 double probFixDExternal,
                 double externalUlindiDerivedProb);
    ~hillClimbing();

    void setLogFileName(const char *);

    void hmmOnly();
    void hmmOnly3states();

    // void randomSearch(int maxIteration);
    // void randomSearch(int maxIteration, double windowInt, double windowExt);
    void randomSearch(
        double converThrld, int maxIteration = MAX_ITERATIONS, double windowInt = WINDOW_INT, double windowExt = WINDOW_EXT, double maxWindowShrinkage = MAX_WINDOW_SHRINKAGE);

    void randomSearch3states(double converThrld,
                             int maxIteration = MAX_ITERATIONS,
                             double windowInt = WINDOW_INT,
                             double windowExt = WINDOW_EXT,
                             double windowLongExt = WINDOW_LONGEXT,
                             double windowLErate = WINDOW_LERATE,
                             double maxWindowShrinkage = MAX_WINDOW_SHRINKAGE);


    double nelderMeadSimplex(double converThrld);
    double likelihoodFunction(const std::vector<double> &x, std::vector<double> &grad, void *f_data);
};

#endif // HILLCLIMBING_H
