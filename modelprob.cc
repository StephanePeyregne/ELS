#include "modelprob.h"

probSegregDInternal::probSegregDInternal()
{
    probs.resize(1024, 0.); // assuming this is the limit.
}

void probSegregDInternal::load(istream &in)
{
    int max_num_derived = 0;

    // probs.resize( 1024, 0. ) ; // assuming this is the limit.
    string s;
    while (in)
    {
        getline(in, s);
        stringstream ss(s, stringstream::in);
        int num_derived;
        ss >> num_derived;
        double prob_derived;
        ss >> prob_derived;
        if (num_derived > 1023)
        {
            cerr << "Can't deal with coverage > 1023. (Adjust in code observation.cc)" << endl;
            exit(1);
        }
        probs[num_derived] = prob_derived;
        if (num_derived > max_num_derived)
        {
            max_num_derived = num_derived;
        }
    }
    probs.resize(max_num_derived + 1);
}

double probSegregDInternal::getProb(int num_derived)
{
    return probs[num_derived];
}

void probSegregDInternal::setProbsAt(int index, double value)
{
    if (index <= probs.size())
    {
        probs[index] = value;
        //	cerr << "probs[index]" << probs[index] << "index:" << index << endl;
    }
    else
    {
        cerr << "Can't setProbsAt index > vector size (observation.cc)   size: " << probs.size() << "   index: " << index << endl;
        exit(1);
    }
}

string probSegregDInternal::printAll(void)
{
    ostringstream oss;
    string myStr;

    for (int i = 1; i < probs.size() - 1; i++) // elements-1 with '\t'
    {
        oss << probs[i];
        oss << "\t";
    }
    oss << probs[probs.size() - 1]; // last element without <<'\t'
    myStr = oss.str();
    return myStr;
}


modelProb::modelProb()
{
    mStayInternal = 0.0;
    mStayExternal = 0.0;
    mStayLongExternal = 0.0;
    mLErate = 0.0;
    mProbFixDInternal = 0.0;
    mProbFixDExternal = 0.0;
    mProbSegregDExternal = 0.0;
}

char modelProb::loadProb(const char *fileName)
{
    // Initialize observation_prob using emission probability config file
    ifstream dataStream(fileName, ifstream::in);
    if (!dataStream)
    {
        cerr << "Error: Can't open internal state emission probabilities config file `" << fileName << "'" << endl;
        return 1;
    }
    obs.load(dataStream);
    // dataStream.clear();
    // dataStream.seekg(0, std::ios::beg);
    // backup.load( dataStream ) ;
}
