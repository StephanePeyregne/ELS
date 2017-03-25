#ifndef HMMRESULTS_H
#define HMMRESULTS_H

#include <vector>

using namespace std;

class hmmResults
{
public:
    hmmResults();
    vector<double> mfwd_stateE_scaled;
    vector<double> mfwd_stateLE_scaled;
    vector<double> mfwd_stateI_scaled;
    vector<double> mbwd_stateI_scaled;
    vector<double> mbwd_stateE_scaled;
    vector<double> mbwd_stateLE_scaled;
    vector<int> mExternal;

    double mlogLikelihood;
};

#endif // HMMRESULTS_H

