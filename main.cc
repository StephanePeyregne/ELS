#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <popt.h>
#include <cstring>
#include <time.h>

using namespace std ;

#include "hillclimbing.h"

int main( int argc, const char *argv[] ) 
{

    // Commandline parameters
    int length_internal = 0 ;
    int length_external = 0 ;
    int length_long_external = 0 ;
    double LErate = 0. ;
    char *configfile = 0 ;
    double prob_fixD_internal = 0. ;
    double prob_fixD_external = 0. ;
    double converThrld = 0.1 ;
    bool baum_welch = false ;
    bool nelderMead = false;
    int additionalParameters = 0 ;
    char *logFileName = 0;

    double EXTERNAL_ULINDI_DERIVED_PROB = 0.01;// This should be a parameter

    struct poptOption optionsTable[] = {
            { "length_internal",                    'L', POPT_ARG_INT, &length_internal, 0, "Average length of internal region", "length" },
            { "length_external",                    'l', POPT_ARG_INT, &length_external, 0, "Average length of external region", "length" },
            { "emission_probabilities",             'e', POPT_ARG_STRING, &configfile, 0, "Emission probabilities for internal state by frequency", "file" },
            { "Baum-Welch",                         'B', POPT_ARG_DOUBLE, &converThrld, 0, "Baum-Welch algorithm to re-estimate parameters", "value" },
            { "NelderMead",                         'N', POPT_ARG_DOUBLE, &converThrld, 0, "NelderMead algorithm", "converThrld" },
            { "prob_fixD_internal",                 'F', POPT_ARG_DOUBLE, &prob_fixD_internal, 0, "Probability for sharing fixed derived when in an internal region", "prob" },
            { "prob_fixD_external",                 'f', POPT_ARG_DOUBLE, &prob_fixD_external, 0, "Probability for sharing fixed derived when in an external region", "prob" },
            { "length_long_external",               'S', POPT_ARG_INT, &length_long_external, 0, "Average length of long external regions", "length"},
            { "long_external_rate",                 'r', POPT_ARG_DOUBLE, &LErate, 0, "Proportion of external regions produced by positive selection", "rate"},
            { "logfile",                            'o', POPT_ARG_STRING, &logFileName, 0, "Set log file name", "path/name"},
            POPT_AUTOHELP { NULL, 0, 0, NULL, 0 } //Should we change something here?
    };

    poptContext optCon;
    optCon = poptGetContext( "sweep_hmm", argc, argv, optionsTable, 0 ) ;
    poptSetOtherOptionHelp( optCon, "[OPTIONS]* <inputfile>" ) ;
    int rc = poptGetNextOpt( optCon ) ;
    if (rc != -1) {
            poptPrintUsage( optCon, stderr, 0 ) ;
            return 1;
    }


    // inputfile present?
    const char* infile = poptGetArg( optCon ) ;
    if (!infile) {
            cerr << "Error: need inputfile as argument." << endl;
            poptPrintUsage( optCon, stderr, 0 ) ;
            return 1;
    }

    // all other parameters present?
    if ( length_internal <= 0 || length_external <= 0 || !configfile || prob_fixD_internal <= 0. || prob_fixD_external <= 0. ||prob_fixD_internal >= 1. || prob_fixD_external >= 1. )
    {
            cerr << "Error: length_internal, length_external, configfile, prob_fixD_internal and prob_fixD_external must be set." << endl ;
            poptPrintUsage( optCon, stderr, 0 ) ;
            return 1;
    }

    // Baum-Welch and/or 3-state algorithm turned on?
    for(int i=0;i<argc;i++)
    {
        if(strcmp(argv[i], "-B")==0)
        {
            baum_welch = true;
        }

	if((strcmp(argv[i], "-S")==0) || (strcmp(argv[i], "-r")==0))
	{
	    additionalParameters += 1;
	}

        if(strcmp(argv[i], "-N")==0)
        {
            nelderMead = true;
        }
    }
    if(additionalParameters==1)
    {
	    cerr << "Error: length_long_external and long_external_rate must be both set when using the 3-state model." << endl ;
            poptPrintUsage( optCon, stderr, 0 ) ;
            return 1;
    }

	if(additionalParameters == 0)
	{
		hillClimbing mAlgorithm(infile, configfile, length_internal, length_external, prob_fixD_internal, prob_fixD_external, EXTERNAL_ULINDI_DERIVED_PROB);

        //Si param nom fichier log
        if(logFileName != 0)
	    {
            mAlgorithm.setLogFileName(logFileName);
	    }

	if(baum_welch == false)
		{
		mAlgorithm.hmmOnly();
		}
	else
		{
		mAlgorithm.randomSearch(converThrld);
		}
	}
	else if(additionalParameters == 2)
	{
		hillClimbing mAlgorithm(infile, configfile, length_internal, length_external, length_long_external, LErate, prob_fixD_internal, prob_fixD_external, EXTERNAL_ULINDI_DERIVED_PROB);

        //Si param nom fichier log
        if(logFileName != 0)
	    {
            mAlgorithm.setLogFileName(logFileName);
	    }

        if(baum_welch == false && nelderMead == false)
		{
			mAlgorithm.hmmOnly3states();
		}
        else if(nelderMead == true)
        {
            mAlgorithm.nelderMeadSimplex(converThrld);
        }
	else
        {
            mAlgorithm.randomSearch3states(converThrld);
        }
	}
}


