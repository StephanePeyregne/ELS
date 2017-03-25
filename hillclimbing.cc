#include "hillclimbing.h"
#include <ctime>

double wrapObjFunc(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
	hillClimbing *obj = static_cast<hillClimbing*>(f_data);
	return obj->likelihoodFunction(x, grad, f_data);
}


hillClimbing::hillClimbing(const char* sequenceFile, const char* configFile,
			   double lengthInternal, double lengthExternal,
                           double probFixDInternal, double probFixDExternal,
                           double externalUlindiDerivedProb)
{
    // convergence_20160330_222635.log
//    time_t t = time(0);   // get time now
//    struct tm * now = localtime( & t );
//    logFileName = "./convergence_"
//         + (now->tm_year + 1900).str()
//         + (now->tm_mon + 1).str()
//         +  now->tm_mday + '_'
//         +  now->tm_hour
//         +  now->tm_min
//         +  now->tm_sec
//        + ".log";

    logFileName = "./convergence.log";

    mOptTolerance = DEFAULT_OPT_TOLERANCE;

    mObsSequence = new obsSequence() ;

    mModelProb = new modelProb();
    mModelProbTemp = NULL;

    mHmmResults = new hmmResults();
    mHmmResultsTemp = NULL;

    mObsSequence->loadSequence(sequenceFile);
    mModelProb->loadProb(configFile);

    //New mHmm class
    mHmm = new hmm(mObsSequence);
    mBaumWelch = new baumWelch(mObsSequence);

    //Init mHmm members
    mModelProb->setStayInternal(lengthInternal);
    mModelProb->setStayExternal(lengthExternal);
    mModelProb->setProbFixDInternal(probFixDInternal);
    mModelProb->setProbFixDExternal(probFixDExternal);
    mModelProb->setProbSegregDExternal(externalUlindiDerivedProb);
}

hillClimbing::hillClimbing(const char* sequenceFile, const char* configFile,
                                  double lengthInternal, double lengthExternal,
			          double lengthLongExternal, double LongExternalrate,
                                  double probFixDInternal, double probFixDExternal,
                                  double externalUlindiDerivedProb)
{
    // convergence_20160330_222635.log
    //time_t t = time(0);   // get time now
    //struct tm * now = localtime( & t );
    //logFileName = "./convergence_"
    //     + (now->tm_year + 1900).str()
    //     + (now->tm_mon + 1).str()
    //     +  now->tm_mday + '_'
    //     +  now->tm_hour
    //     +  now->tm_min
    //     +  now->tm_sec
    //     + ".log";



    logFileName = "./convergence.log";
    mOptTolerance = DEFAULT_OPT_TOLERANCE;

    mObsSequence = new obsSequence() ;

    mModelProb = new modelProb();
    mModelProbTemp = NULL;

    mHmmResults = new hmmResults();
    mHmmResultsTemp = NULL;

    mObsSequence->loadSequence(sequenceFile);
    mModelProb->loadProb(configFile);

    //New mHmm class
    mHmm = new hmm(mObsSequence);
    mBaumWelch = new baumWelch(mObsSequence);

    //Init mHmm members
    mModelProb->setStayInternal(lengthInternal);
    mModelProb->setStayExternal(lengthExternal);
    mModelProb->setStayLongExternal(lengthLongExternal);
    mModelProb->setLErate(LongExternalrate);
    mModelProb->setProbFixDInternal(probFixDInternal);
    mModelProb->setProbFixDExternal(probFixDExternal);
    mModelProb->setProbSegregDExternal(externalUlindiDerivedProb);
}

hillClimbing::~hillClimbing()
{
    delete mObsSequence;
    delete mModelProb;
    delete mHmmResults;
    delete mHmm;
    delete mBaumWelch;
}

void hillClimbing::setLogFileName(const char* name)
{
    logFileName = name;
}

void hillClimbing::hmmOnly()
{
    mHmm->computeFwdBwd(mModelProb, mHmmResults);
    mHmm->writeOutputFile();
}

void hillClimbing::hmmOnly3states()
{
    mHmm->computeFwdBwd3states(mModelProb, mHmmResults);
    mHmm->writeOutputFile3states();
}
/*
void hillClimbing::randomSearch(int maxIteration)
{
    randomSearch(maxIteration, WINDOW_INT, WINDOW_EXT, MAX_WINDOW_SHRINKAGE);
}

void hillClimbing::randomSearch(int maxIteration, double windowInt, double windowExt)
{
    randomSearch(maxIteration, windowInt, windowExt, MAX_WINDOW_SHRINKAGE);
}
*/
void hillClimbing::randomSearch(double converThrld, int maxIteration, double windowInt, double windowExt, double maxWindowShrinkage)
{
    mHmm->computeFwdBwd(mModelProb, mHmmResults);

    double windowShrinkage  = 1;
    srand(time(NULL));

    // Run the random search

    double stayInternal;
    double stayExternal;
    map<int,double> proba_obsD_int ;

    stayInternal = mModelProb->getStayInternal();
    stayExternal = mModelProb->getStayExternal();

    for (int n=1; n<=maxIteration; n++)
    {
        int step=1 ;
        double baum_welch_LogLikelihood=2*mHmmResults->mlogLikelihood;
        int max_iteration = 40;//Change this too

        mModelProbTemp = new modelProb();
        mHmmResultsTemp = new hmmResults();

        mModelProbTemp->resizeProbSegregDInternal(mModelProb->getSizeProbSegregDInternal());

        // Propose new values for the parameters of the transition probabilities
        mModelProbTemp->setStayInternal(mModelProb->getStayInternal() + (rand() % (int)(windowInt/windowShrinkage)) - (int)((windowInt/windowShrinkage)/2));
        mModelProbTemp->setStayExternal(mModelProb->getStayExternal() + (rand() % (int)(windowExt/windowShrinkage)) - (int)((windowExt/windowShrinkage)/2));
        if (mModelProbTemp->getStayInternal() < 0.0)
                mModelProbTemp->setStayInternal(-mModelProbTemp->getStayInternal());
        else if(mModelProbTemp->getStayInternal() == 0.0)
		cerr << "stayInt = 0" << endl;
        if (mModelProbTemp->getStayExternal() < 0.0)
                mModelProbTemp->setStayExternal(-mModelProbTemp->getStayExternal());
        else if(mModelProbTemp->getStayExternal() == 0.0)
                cerr << "stayExt = 0" << endl;


        // We keep this parameter throughout the algorithm because somehow we cannot re-estimate it
        mModelProbTemp->setProbSegregDExternal(mModelProb->getProbSegregDExternal());

            // Run the Baum-Welch algorithm with the new values

        mHmmResultsTemp->mlogLikelihood = mHmmResults->mlogLikelihood;

        // La boucle while correspond a une methode convergenceBaumWelch autonome
        // qui serait appellee par les differents algo qui proposent des valeurs stayInternal et External (convergence des transitions)
        while((step<max_iteration)&&((mHmmResultsTemp->mlogLikelihood - baum_welch_LogLikelihood)>converThrld))
        {
            //ostringstream oss;
            //oss << step;
            //string logInfo = oss.str();
            //mHmm.writeLogFile("./convergence.log",logInfo);

            step++;
            //keep previous likelihood to detect convergence
            baum_welch_LogLikelihood = mHmmResultsTemp->mlogLikelihood;
            mBaumWelch->computeBaumWelch(mHmmResults, mModelProbTemp);
            mHmm->computeFwdBwd(mModelProbTemp, mHmmResultsTemp);
        }

        // Accept or reject the new values.
        if (mHmmResults->mlogLikelihood < mHmmResultsTemp->mlogLikelihood)
        {
            //Keep new model
            delete mModelProb;
            mModelProb = mModelProbTemp;

            //Update Results
            delete mHmmResults;
            mHmmResults = mHmmResultsTemp;

            if (windowShrinkage>1)
                    windowShrinkage--;

            ostringstream oss;
            oss << step;
            string logInfo = oss.str();
            mHmm->writeLogFile(logFileName,logInfo);
        }
        else
        {
            //Reject new model
            delete mModelProbTemp;
            //Reject results
            delete mHmmResultsTemp;
            mModelProbTemp = NULL;
            mHmmResultsTemp = NULL;
            windowShrinkage++;
        }

        // Convergence?
        if (windowShrinkage==maxWindowShrinkage)
            break;
    }
    //write output file with last parameters
    mHmm->computeFwdBwd(mModelProb, mHmmResults);//TODO methode pour set HmmResults afin de changer le pointeur pour ecrire OutputFile
    mHmm->writeOutputFile();
 
}

void hillClimbing::randomSearch3states(double converThrld, int maxIteration, double windowInt, double windowExt, double windowLongExt, double windowLErate, double maxWindowShrinkage)
{
    //cout << "Here0" << flush;
    mHmm->computeFwdBwd3states(mModelProb, mHmmResults);

    double windowShrinkage  = 1;
    srand(time(NULL));

    // Run the random search

    double stayInternal;
    double stayExternal;
    double stayLongExternal;
    double LErate;
    map<int,double> proba_obsD_int ;

    stayInternal = mModelProb->getStayInternal();
    stayExternal = mModelProb->getStayExternal();
    stayLongExternal = mModelProb->getStayLongExternal();
    LErate = mModelProb->getLErate();

    for (int n=1; n<=maxIteration; n++)
    {
        //cout << "Here" << flush;
        int step=1 ;
        double baum_welch_LogLikelihood=2*mHmmResults->mlogLikelihood;
        int max_iteration = 40;//Change this too

        mModelProbTemp = new modelProb();
        mHmmResultsTemp = new hmmResults();

        mModelProbTemp->resizeProbSegregDInternal(mModelProb->getSizeProbSegregDInternal());

        // Propose new values for the parameters of the transition probabilities
        mModelProbTemp->setStayInternal(mModelProb->getStayInternal() + (rand() % (int)(windowInt/windowShrinkage)) - (int)((windowInt/windowShrinkage)/2));
        mModelProbTemp->setStayExternal(mModelProb->getStayExternal() + (rand() % (int)(windowExt/windowShrinkage)) - (int)((windowExt/windowShrinkage)/2));
        mModelProbTemp->setStayLongExternal(mModelProb->getStayLongExternal() + (rand() % (int)(windowLongExt/windowShrinkage)) - (int)((windowLongExt/windowShrinkage)/2));
        mModelProbTemp->setLErate(mModelProb->getLErate() + (((double)rand() / RAND_MAX) * (double)(windowLErate/windowShrinkage)) - (double)((windowLErate/windowShrinkage)/2));

        if (mModelProbTemp->getStayInternal() < 0.0)
                mModelProbTemp->setStayInternal(-mModelProbTemp->getStayInternal());
        else if(mModelProbTemp->getStayInternal() == 0.0)
                cerr << "stayInt = 0" << endl;

        if (mModelProbTemp->getStayExternal() < 0.0)
                mModelProbTemp->setStayExternal(-mModelProbTemp->getStayExternal());
        else if(mModelProbTemp->getStayExternal() == 0.0)
                cerr << "stayExt = 0" << endl;

        if (mModelProbTemp->getStayLongExternal() < mModelProbTemp->getStayExternal())
            mModelProbTemp->setStayLongExternal(2*mModelProbTemp->getStayExternal() - mModelProbTemp->getStayLongExternal());
        else if(mModelProbTemp->getStayLongExternal() == 0.0)
                    cerr << "stayLongExt = 0" << endl;

        if (mModelProbTemp->getLErate() < 0.0)
            mModelProbTemp->setLErate(-mModelProbTemp->getLErate());
        else if (mModelProbTemp->getLErate() > 1.0)
            mModelProbTemp->setLErate(2*mModelProbTemp->getLErate() - 1);
	

        // We keep this parameter throughout the algorithm because somehow we cannot re-estimate it 
        mModelProbTemp->setProbSegregDExternal(mModelProb->getProbSegregDExternal());

        // Run the Baum-Welch algorithm with the new values

        mHmmResultsTemp->mlogLikelihood = mHmmResults->mlogLikelihood;

        // La boucle while correspond a une methode convergenceBaumWelch autonome
        // qui serait appellee par les differents algo qui proposent des valeurs stayInternal et External (convergence des transitions)
        while((step<max_iteration)&&((mHmmResultsTemp->mlogLikelihood - baum_welch_LogLikelihood)>converThrld))
        {
            //ostringstream oss;
            //oss << step;
            //string logInfo = oss.str();
            //mHmm.writeLogFile("./convergence.log",logInfo);
            step++;
            //keep previous likelihood to detect convergence
            baum_welch_LogLikelihood = mHmmResultsTemp->mlogLikelihood;
            mBaumWelch->computeBaumWelch3states(mHmmResults, mModelProbTemp);
            mHmm->computeFwdBwd3states(mModelProbTemp, mHmmResultsTemp);
        }
	
        // Accept or reject the new values.
        if (mHmmResults->mlogLikelihood < mHmmResultsTemp->mlogLikelihood)
        {
            //Keep new model
            delete mModelProb;
            mModelProb = mModelProbTemp;

            //Update Results
            delete mHmmResults;
            mHmmResults = mHmmResultsTemp;

            if (windowShrinkage>1)
                    windowShrinkage--;

            ostringstream oss;
            oss << step;
            string logInfo = oss.str();
            mHmm->writeLogFile3states(logFileName,logInfo);
        }
        else
        {
            //Reject new model
            delete mModelProbTemp;
            //Reject results
            delete mHmmResultsTemp;
            mModelProbTemp = NULL;
            mHmmResultsTemp = NULL;
            windowShrinkage++;
        }

        // Convergence?
        if (windowShrinkage==maxWindowShrinkage)
            break;
    }
    //write output file with last parameters
    mHmm->computeFwdBwd3states(mModelProb, mHmmResults);//TODO methode pour set HmmResults afin de changer le pointeur pour ecrire OutputFile
    mHmm->writeOutputFile3states();

}





//optimisation d'un modele
void hillClimbing::optimisationBaumWelch(modelProb *model, hmmResults *results, int max_iteration)
{
    int step = 1;
    double baum_welch_LogLikelihood;

    //Start from base model
    baum_welch_LogLikelihood = mHmmResults->mlogLikelihood;
    mBaumWelch->computeBaumWelch3states(mHmmResults, model);
    mHmm->computeFwdBwd3states(model, results);

    while((step<max_iteration)&&(results->mlogLikelihood - baum_welch_LogLikelihood)>mConverThrld)
    {
        step++;
        //keep previous likelihood to detect convergence
        baum_welch_LogLikelihood = results->mlogLikelihood;
        mBaumWelch->computeBaumWelch3states(results, model);
        mHmm->computeFwdBwd3states(model, results);
    }
}

//Function pour nlopt

//double f(const std::vector<double> &x, std::vector<double> &grad, void* f_data);

double hillClimbing::likelihoodFunction(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
    double likelihood;

    mModelProbTemp = new modelProb();
    mHmmResultsTemp = new hmmResults();

    mModelProbTemp->resizeProbSegregDInternal(mModelProb->getSizeProbSegregDInternal());

    // Propose new values for the parameters of the transition probabilities
    mModelProbTemp->setStayInternal(x.at(IDX_STAY_INTERNAL));
    mModelProbTemp->setStayExternal(x.at(IDX_STAY_EXTERNAL));
    mModelProbTemp->setStayLongExternal(x.at(IDX_STAY_LONG_INTERNAL));
    mModelProbTemp->setLErate(x.at(IDX_LERATE));

    // We keep this parameter throughout the algorithm because somehow we cannot re-estimate it
    mModelProbTemp->setProbSegregDExternal(mModelProb->getProbSegregDExternal());

    // Run the Baum-Welch algorithm with the new values

    mHmmResultsTemp->mlogLikelihood = mHmmResults->mlogLikelihood;

    optimisationBaumWelch(mModelProbTemp, mHmmResultsTemp, 40);

    likelihood = mHmmResultsTemp->mlogLikelihood;

    //Write log file
    ostringstream oss;
    //oss << step;
    string logInfo = oss.str();
    mHmm->writeLogFile3states(logFileName,logInfo);

    //Reject new model
    delete mModelProbTemp;
    //Reject results
    delete mHmmResultsTemp;
    mModelProbTemp = NULL;
    mHmmResultsTemp = NULL;

    return likelihood;
}

double constraint(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
    return x.at(IDX_STAY_EXTERNAL) - x.at(IDX_STAY_LONG_INTERNAL); //stayExternal - stayLongExternal <= 0
}

//Algo de recherche avec nlopt
double hillClimbing::nelderMeadSimplex(double converThrld)
{
    mConverThrld = converThrld;

    mHmm->computeFwdBwd3states(mModelProb, mHmmResults);

    //Optimisation du modèle de base
    optimisationBaumWelch(mModelProb, mHmmResults, 40);

    //A partir du modèle précédent on utilise nlopt pour rechercher le maximum en fonctions des 4 paramètres :
        //StayInternal
        //StayExternal
        //StayLongExternal
        //LErate

    //Init nlopt
//    mOpt = new nlopt::opt(nlopt::LN_NELDERMEAD, 4);
    mOpt = new nlopt::opt(nlopt::LN_COBYLA, 4);

    //Borner limites nlopt
    mOpt->set_lower_bounds(0.0);

    double mybounds[] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, 1.0}; //HUGE_VAL inclus dans math.h
    vector<double> values(mybounds, mybounds + sizeof(mybounds) / sizeof(double) );
    mOpt->set_upper_bounds(values);
    //mOpt->set_upper_bounds(vector<double>({HUGE_VAL, HUGE_VAL, HUGE_VAL, 1.0}, 4)); //HUGE_VAL math.h


    //Contrainte stayLongExtrnal > stayExternal
    mOpt->add_inequality_constraint(constraint, NULL, 0.0);

    //Criteria
    mOpt->set_ftol_abs(mOptTolerance);
    mOpt->set_maxeval(1000);//Max iteration
    mOpt->set_max_objective(wrapObjFunc, this);

    //String de description de l'algorithme obtenue avec
    //const char *nlopt::opt::get_algorithm_name() const;



    //Initialise parametres de recherche
    vector<double> x(4,0.0);
    double opt_f;
    x.at(IDX_STAY_INTERNAL) = mModelProb->getStayInternal();
    x.at(IDX_STAY_EXTERNAL) = mModelProb->getStayExternal();
    x.at(IDX_STAY_LONG_INTERNAL) = mModelProb->getStayLongExternal();
    x.at(IDX_LERATE) = mModelProb->getLErate();

    //appel nlopt pour faire la recherche
    char result; //1 octet : va de -128 jusqu'à +127
    result = mOpt->optimize(x, opt_f);

    if(result > 0)
    {
        //Opt ok
        //write output file with best parameters

	mModelProb->setStayInternal(x.at(IDX_STAY_INTERNAL));
        mModelProb->setStayExternal(x.at(IDX_STAY_EXTERNAL));
        mModelProb->setStayLongExternal(x.at(IDX_STAY_LONG_INTERNAL));
        mModelProb->setLErate(x.at(IDX_LERATE));

        optimisationBaumWelch(mModelProb, mHmmResults, 40);
        //likelihoodFunction(x, grad,  NULL);
        mHmm->writeOutputFile3states();

    }
    else
    {
        cout << "Nlopt optimisation failed: " << result << "\t";//erreur
    }
    delete mOpt;
}

/*
 *
 * class YourClass
{
public:
    double f(unsigned int i, const double*, double*, void*);
};

double wrapper(unsigned int i, const double* a, double* b, void* o) {
    // not sure what you'd need the last param for, but you say you have it...
    reinterpet_cast<YourClass*>(o)->f(i, a, b, o);
}

and then:

YourClass myP;
opt.set_min_objective(wrapper, &myP);
*/
