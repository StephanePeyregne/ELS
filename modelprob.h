#ifndef MODELPROB_H
#define MODELPROB_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


class probSegregDInternal
{
  public:
    probSegregDInternal();
    void load(istream &in);
    int getSize(void)
    {
        return probs.size();
    };
    void resize(int value)
    {
        probs.resize(value);
    };
    double getProb(int num_derived);
    void setProbsAt(int index, double value);
    string printAll(void);

  private:
    vector<double> probs;  // probabilities: 0 is not used, $1 .. coverage-1$ is for derived alleles 1 to coverage-1.
    size_t fixed_coverage; // used to check whether input-lines are wrong
};


class modelProb
{
  public:
    modelProb();

    char loadProb(const char *fileName);

    // clang-format off
    void setStayInternal(double value)		{mStayInternal = value;};
    void setStayExternal(double value)		{mStayExternal = value;};
    void setStayLongExternal(double value)      {mStayLongExternal = value;};
    void setLErate(double value)                {mLErate = value;};
    void setProbFixDInternal(double value)	{mProbFixDInternal = value;};
    void setProbFixDExternal(double value)	{mProbFixDExternal = value;};
    void setProbSegregDExternal(double value)	{mProbSegregDExternal = value;};
    void setProbSegregDInternalAt(int i, double value) {obs.setProbsAt(i, value);};
    void resizeProbSegregDInternal(int size) {obs.resize(size);};
    int getSizeProbSegregDInternal(void) {return obs.getSize();};
    probSegregDInternal getProbSegregDInternal(void) {return obs;};//read only access

    double getStayInternal(void)		{return mStayInternal;};
    double getStayExternal(void)		{return mStayExternal;};
    double getStayLongExternal(void)            {return mStayLongExternal;};
    double getLErate(void)                      {return mLErate;};
    double getProbFixDInternal(void)		{return mProbFixDInternal;};
    double getProbFixDExternal(void)		{return mProbFixDExternal;};
    double getProbSegregDExternal(void)         {return mProbSegregDExternal;};
    // clang-format on


  protected:
    probSegregDInternal obs;
    // probSegregDInternal backup; //used only in the optimization step to get back to old values if the new ones are not better

    double mStayInternal;
    double mStayExternal;
    double mStayLongExternal;
    double mLErate;
    double mProbFixDInternal;
    double mProbFixDExternal;
    double mProbSegregDExternal;
};


#endif // MODELPROB_H
