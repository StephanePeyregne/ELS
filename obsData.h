#ifndef OBSDATA_H
#define OBSDATA_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


class obsSite
{
    friend std::istream &operator>>(std::istream &, obsSite &);

  public:
    obsSite();

    // clang-format off
		string getChr(void)			{return mChr;};
		long getLocation(void)			{return mLocation;};
		int getCoverage(void)			{return mCoverage;};
		int getDerived(void)			{return mDerived;};
		char getClint(void)			{return mClint;};
		char getUlindi(void)			{return mUlindi;};
		double getDist(void)			{return mDist;};
		double getDistBin(void);
		
		void setChr(string value)		{mChr = value;};
		void setLocation(long value)		{mLocation = value;};
		void setCoverage(int value)		{mCoverage = value;};
		void setDerived(int value)		{mDerived = value;};
		void setClint(char value)		{mClint = value;};
		void setUlindi(char value)		{mUlindi = value;};
		void setDist(double value)		{mDist = value;};
    // clang-format on

    string print(void);

  private:
    string mChr;
    long mLocation;
    int mCoverage;
    int mDerived;
    char mClint;
    char mUlindi;
    double mDist;
};

class obsSequence
{
  public:
    obsSequence();
    char loadSequence(const char *fileName);

    int size(void);
    obsSite at(int i);

  protected:
    vector<obsSite> mSequence;
};


#endif
