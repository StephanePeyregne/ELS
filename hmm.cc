#include "hmm.h"
#include <iomanip>

hmm::hmm( obsSequence *seq)
{
        mObsSequence = seq;
        mModelProb = NULL;
        mHmmResults = NULL;
}

double hmm::prob_obs_internal( obsSite dat, probSegregDInternal prob )
{
	if ( dat.getDerived() < dat.getCoverage() ) 
		if ( dat.getUlindi() == 'D' ) 
                        return prob.getProb( dat.getDerived() ) ;
                else return 1.-prob.getProb( dat.getDerived() ) ;
	else 
		if ( dat.getUlindi() == 'D' ) 
                        return mModelProb->getProbFixDInternal() ;
                else return 1.-mModelProb->getProbFixDInternal() ;
}

double hmm::prob_obs_external( obsSite dat )
{
	if ( dat.getDerived() < dat.getCoverage() ) 
		if ( dat.getUlindi() == 'D' )
                        return mModelProb->getProbSegregDExternal() ;
                else return 1.-mModelProb->getProbSegregDExternal() ;
	else 
		if ( dat.getUlindi() == 'D' )
                        return mModelProb->getProbFixDExternal() ;
                else return 1.-mModelProb->getProbFixDExternal() ;
}

void hmm::computeFwdBwd(modelProb *model, hmmResults *result)
{

    mModelProb = model;
    mHmmResults = result;

	//////////
	// forward probabilities 
        mHmmResults->mfwd_stateE_scaled.resize( mObsSequence->size(), 0.0 ) ; // external
        mHmmResults->mfwd_stateI_scaled.resize( mObsSequence->size(), 0.0 ) ; // internal
	// first one by hand
        vector<double> fwd_scaling( mObsSequence->size(), 0.0 ) ;
	// Initial state = internal
        fwd_scaling[0] = prob_obs_internal( mObsSequence->at(0), mModelProb->getProbSegregDInternal() ) ;
	mHmmResults->mfwd_stateI_scaled[0] = 1. ; // == fwd_stateI[0]/fwd_scaling[0] ;
	// rest
        for ( int i = 1 ; i < mObsSequence->size() ; ++i ) {
		
            double I_tmp = prob_obs_internal( mObsSequence->at(i), mModelProb->getProbSegregDInternal() )
                                                                        *mHmmResults->mfwd_stateI_scaled[i-1]*1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayInternal()) +
                                  prob_obs_internal( mObsSequence->at(i), mModelProb->getProbSegregDInternal() )
                                                                  *mHmmResults->mfwd_stateE_scaled[i-1]*(1.-1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayExternal())) ;
                double E_tmp =	( prob_obs_external( mObsSequence->at(i) )
                                * mHmmResults->mfwd_stateI_scaled[i-1]*(1.-1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayInternal())) +
                                  prob_obs_external( mObsSequence->at(i) )
                                  * mHmmResults->mfwd_stateE_scaled[i-1]*1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayExternal()) ) ;

		fwd_scaling[i] = I_tmp + E_tmp ;

                mHmmResults->mfwd_stateI_scaled[i] = I_tmp/fwd_scaling[i] ;
                mHmmResults->mfwd_stateE_scaled[i] = E_tmp/fwd_scaling[i] ;
	}

	//////////////
	// Log-likelihood
	mHmmResults->mlogLikelihood = 0.0;
        for ( int i = 1 ; i < mObsSequence->size() ; ++i ) {
                mHmmResults->mlogLikelihood += log(fwd_scaling[i]);
	}

	//////////////////
	// backward probabilities 
        mHmmResults->mbwd_stateI_scaled.resize( mObsSequence->size(), 0.0 ) ;
        mHmmResults->mbwd_stateE_scaled.resize( mObsSequence->size(), 0.0 ) ;
        mHmmResults->mbwd_stateI_scaled[mObsSequence->size()-1] = 1. ;
        mHmmResults->mbwd_stateE_scaled[mObsSequence->size()-1] = 1. ;
        for ( int i = mObsSequence->size()-2 ; i >= 0 ; i-- ) {
		
		double scale = fwd_scaling[i+1] ;
                mHmmResults->mbwd_stateI_scaled[i] = (exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayInternal())
                * prob_obs_internal( mObsSequence->at(i+1), mModelProb->getProbSegregDInternal() )*mHmmResults->mbwd_stateI_scaled[i+1] +
                                 (1.-exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayInternal()))
                                 * prob_obs_external( mObsSequence->at(i+1) )*mHmmResults->mbwd_stateE_scaled[i+1])/scale ;
				 
                mHmmResults->mbwd_stateE_scaled[i] = ((1.-exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayExternal()))
                * prob_obs_internal( mObsSequence->at(i+1), mModelProb->getProbSegregDInternal() )*mHmmResults->mbwd_stateI_scaled[i+1] +
                                 exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayExternal()) * prob_obs_external( mObsSequence->at(i+1) )
                                 *mHmmResults->mbwd_stateE_scaled[i+1])/scale ;

        }

	//////////////////////
	// posterior decoding
        mHmmResults->mExternal.resize( mObsSequence->size(), 0 ) ; // == 1 if external
        for ( size_t i = 0 ; i < mObsSequence->size() ; ++i ) {
                double val = (mHmmResults->mfwd_stateI_scaled[i]*mHmmResults->mbwd_stateI_scaled[i]) ;
                if ( val < 0.2 ) mHmmResults->mExternal[i] = 1 ;
                else mHmmResults->mExternal[i] = 0 ;
        }
}

void hmm::computeFwdBwd3states(modelProb *model, hmmResults *result)
{

    mModelProb = model;
    mHmmResults = result;

        //////////
        // forward probabilities 
        mHmmResults->mfwd_stateE_scaled.resize( mObsSequence->size(), 0.0 ) ; // external
	mHmmResults->mfwd_stateLE_scaled.resize( mObsSequence->size(), 0.0 ) ; // long external
        mHmmResults->mfwd_stateI_scaled.resize( mObsSequence->size(), 0.0 ) ; // internal
        // first one by hand
        vector<double> fwd_scaling( mObsSequence->size(), 0.0 ) ;
        // Initial state = internal
        fwd_scaling[0] = prob_obs_internal( mObsSequence->at(0), mModelProb->getProbSegregDInternal() ) ;
        mHmmResults->mfwd_stateI_scaled[0] = 1. ; // == fwd_stateI[0]/fwd_scaling[0] ;
        // rest
        for ( int i = 1 ; i < mObsSequence->size() ; ++i ) {

            double I_tmp = prob_obs_internal( mObsSequence->at(i), mModelProb->getProbSegregDInternal() )
                                                *mHmmResults->mfwd_stateI_scaled[i-1]*1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayInternal()) +
                                  prob_obs_internal( mObsSequence->at(i), mModelProb->getProbSegregDInternal() )
                                                *mHmmResults->mfwd_stateE_scaled[i-1]*(1.-1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayExternal())) +
				  prob_obs_internal( mObsSequence->at(i), mModelProb->getProbSegregDInternal() )
                                                *mHmmResults->mfwd_stateLE_scaled[i-1]*(1.-1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayLongExternal())) ;

            double E_tmp =  ( prob_obs_external( mObsSequence->at(i) )
                                		* mHmmResults->mfwd_stateI_scaled[i-1]*(1-mModelProb->getLErate())*(1.-1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayInternal())) +
                                  prob_obs_external( mObsSequence->at(i) )
                                  		* mHmmResults->mfwd_stateE_scaled[i-1]*1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayExternal()) ) ;

            double LE_tmp = ( prob_obs_external( mObsSequence->at(i) )
                                		* mHmmResults->mfwd_stateI_scaled[i-1]*mModelProb->getLErate()*(1.-1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayInternal())) +
                                  prob_obs_external( mObsSequence->at(i) )
                                  		* mHmmResults->mfwd_stateLE_scaled[i-1]*1./exp(mObsSequence->at(i).getDist()/mModelProb->getStayLongExternal()) ) ;

            fwd_scaling[i] = I_tmp + E_tmp + LE_tmp ;

            mHmmResults->mfwd_stateI_scaled[i] = I_tmp/fwd_scaling[i] ;
            mHmmResults->mfwd_stateE_scaled[i] = E_tmp/fwd_scaling[i] ;
	    mHmmResults->mfwd_stateLE_scaled[i] = LE_tmp/fwd_scaling[i] ;
        }

        //////////////
        // Log-likelihood
        mHmmResults->mlogLikelihood = 0.0;
        for ( int i = 1 ; i < mObsSequence->size() ; ++i ) {
                mHmmResults->mlogLikelihood += log(fwd_scaling[i]);
        }
	//////////////////
        // backward probabilities 
        mHmmResults->mbwd_stateI_scaled.resize( mObsSequence->size(), 0.0 ) ;
        mHmmResults->mbwd_stateE_scaled.resize( mObsSequence->size(), 0.0 ) ;
	mHmmResults->mbwd_stateLE_scaled.resize( mObsSequence->size(), 0.0 ) ;
        mHmmResults->mbwd_stateI_scaled[mObsSequence->size()-1] = 1. ;
        mHmmResults->mbwd_stateE_scaled[mObsSequence->size()-1] = 1. ;
	mHmmResults->mbwd_stateLE_scaled[mObsSequence->size()-1] = 1. ;
        for ( int i = mObsSequence->size()-2 ; i >= 0 ; i-- ) {

                double scale = fwd_scaling[i+1] ;
                mHmmResults->mbwd_stateI_scaled[i] = (exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayInternal())
                				* prob_obs_internal( mObsSequence->at(i+1), mModelProb->getProbSegregDInternal() )*mHmmResults->mbwd_stateI_scaled[i+1] +
                                 (1-mModelProb->getLErate())*(1.-exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayInternal()))
                                 		* prob_obs_external( mObsSequence->at(i+1) )*mHmmResults->mbwd_stateE_scaled[i+1] +
				 mModelProb->getLErate()*(1.-exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayInternal()))
                                                * prob_obs_external( mObsSequence->at(i+1) )*mHmmResults->mbwd_stateLE_scaled[i+1])/scale ;

                mHmmResults->mbwd_stateE_scaled[i] = ((1.-exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayExternal()))
                				* prob_obs_internal( mObsSequence->at(i+1), mModelProb->getProbSegregDInternal() )*mHmmResults->mbwd_stateI_scaled[i+1] +
                                 exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayExternal()) * prob_obs_external( mObsSequence->at(i+1) )
                                 		*mHmmResults->mbwd_stateE_scaled[i+1])/scale ;

		mHmmResults->mbwd_stateLE_scaled[i] = ((1.-exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayLongExternal()))
                                                * prob_obs_internal( mObsSequence->at(i+1), mModelProb->getProbSegregDInternal() )*mHmmResults->mbwd_stateI_scaled[i+1] +
                                 exp(-1.*mObsSequence->at(i+1).getDist()/mModelProb->getStayLongExternal()) * prob_obs_external( mObsSequence->at(i+1) )
                                                *mHmmResults->mbwd_stateLE_scaled[i+1])/scale ;

        }

        //////////////////////
        // posterior decoding
        mHmmResults->mExternal.resize( mObsSequence->size(), 0 ) ; // == 1 if external
        for ( size_t i = 0 ; i < mObsSequence->size() ; ++i ) {
                double val = (mHmmResults->mfwd_stateE_scaled[i]*mHmmResults->mbwd_stateE_scaled[i]) ;
		double val2 = (mHmmResults->mfwd_stateLE_scaled[i]*mHmmResults->mbwd_stateLE_scaled[i]) ;
                if ( val > 0.8 ) mHmmResults->mExternal[i] = 1 ;
		else if ( val2 > 0.8 ) mHmmResults->mExternal[i] = 2 ;
                else mHmmResults->mExternal[i] = 0 ;
        }
}


void hmm::writeLogFile(const char* fileName, string logInfo)
{
	ofstream myfile;
	myfile.open(fileName, std::ios_base::app);
	myfile 	<< std::setprecision (15) << logInfo << "\t"
                << mHmmResults->mlogLikelihood << "\t"
                << mModelProb->getStayInternal() << "\t"
                << mModelProb->getStayExternal() << "\t"
                << mModelProb->getProbFixDInternal() << "\t"
                << mModelProb->getProbFixDExternal() << "\t"
                << mModelProb->getProbSegregDExternal() << "\t"
                << mModelProb->getProbSegregDInternal().printAll()
		<< "\n";
	myfile.close();
}

void hmm::writeLogFile3states(const char* fileName, string logInfo)
{
        ofstream myfile;
        myfile.open(fileName, std::ios_base::app);
        myfile  << std::setprecision (15) << logInfo << "\t"
                << mHmmResults->mlogLikelihood << "\t"
		<< mModelProb->getLErate() << "\t"
                << mModelProb->getStayInternal() << "\t"
                << mModelProb->getStayExternal() << "\t"
		<< mModelProb->getStayLongExternal() << "\t"
                << mModelProb->getProbFixDInternal() << "\t"
                << mModelProb->getProbFixDExternal() << "\t"
                << mModelProb->getProbSegregDExternal() << "\t"
                << mModelProb->getProbSegregDInternal().printAll()
                << "\n";
        myfile.close();
}

void hmm::writeOutputFile(void) 
{
                for(size_t i=1; i<mObsSequence->size(); ++i) {
                cout << mObsSequence->at(i).getChr() << "\t"
                                << mObsSequence->at(i).getLocation() << "\t"
                                << mObsSequence->at(i).getClint() << "\t"
                                << mObsSequence->at(i).getCoverage() << "\t"
                                << mObsSequence->at(i).getDerived() << "\t"
                                << mObsSequence->at(i).getUlindi() << "\t"
                                << mObsSequence->at(i).getDist() << "\t"
                                << mHmmResults->mExternal[i] << "\t"
                                << mHmmResults->mfwd_stateE_scaled[i]*mHmmResults->mbwd_stateE_scaled[i] << "\n" ;
	}
}

void hmm::writeOutputFile3states(void)
{
                for(size_t i=1; i<mObsSequence->size(); ++i) {
                cout << mObsSequence->at(i).getChr() << "\t"
                                << mObsSequence->at(i).getLocation() << "\t"
                                << mObsSequence->at(i).getClint() << "\t"
                                << mObsSequence->at(i).getCoverage() << "\t"
                                << mObsSequence->at(i).getDerived() << "\t"
                                << mObsSequence->at(i).getUlindi() << "\t"
                                << mObsSequence->at(i).getDist() << "\t"
                                << mHmmResults->mExternal[i] << "\t"
                                << mHmmResults->mfwd_stateE_scaled[i]*mHmmResults->mbwd_stateE_scaled[i] << "\t"
				<< mHmmResults->mfwd_stateLE_scaled[i]*mHmmResults->mbwd_stateLE_scaled[i] << "\n" ;
        }
}
