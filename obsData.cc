
#include "obsData.h"


obsSite::obsSite()
{
    mChr = "";
    mLocation = 0;
    mCoverage = 0;
    mDerived = 0;
    mClint = 0;
    mUlindi = 0;
    mDist = 0.0;
}


// redefine >> operator
std::istream &operator>>(std::istream &Stream, obsSite &Obj)
{
    Stream >> Obj.mChr >> Obj.mLocation >> Obj.mClint >> Obj.mCoverage >> Obj.mDerived >> Obj.mUlindi >> Obj.mDist;
    return Stream;
}


obsSequence::obsSequence()
{
}

int obsSequence::size(void)
{
    return mSequence.size();
}

obsSite obsSequence::at(int i)
{
    return mSequence[i];
}

char obsSequence::loadSequence(const char *fileName)
{
    ifstream dataStream(fileName, ifstream::in);
    if (!dataStream)
    {
        cerr << "Error: Can't open inputfile `" << fileName << "'" << endl;
        return 1;
    }
    while (dataStream)
    {
        string s;
        getline(dataStream, s);
        istringstream ss(s);
        obsSite tmp;
        ss >> tmp;

        if (tmp.getChr() != "")
            mSequence.push_back(tmp);
    }
}
