#include "baumwelch.h"

baumWelch::baumWelch(obsSequence *data)
{
    mObsSequence = data;
}

void baumWelch::computeBaumWelch(hmmResults *results ,modelProb *model)
{
    mHmmResults = results;
    mModelProb = model;

    //////////////////////
    // maximum likelihood estimation of the transition and emission probabilities

    map<int,double> pseudo_count_obsD_int ;
    map<int,double> pseudo_count_obsD_ext ;
    map<int,double> pseudo_count_obsA_int ;
    map<int,double> pseudo_count_obsA_ext ;
    map<int,double> total_pseudo_count_obs_int ;
    map<const char*,double> total_pseudo_count_obs_ext ;
    //map<int,double> total_pseudo_count_obs_ext ;

    map<int,double> proba_obsD_int ;
    map<const char*,double> proba_obsD_ext ;

    for ( int i = 1 ; i < mObsSequence->size()-1 ; ++i ) {

            if (mObsSequence->at(i+1).getUlindi() == 'D'){
                    pseudo_count_obsD_int[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateI_scaled[i+1] *  mHmmResults->mbwd_stateI_scaled[i+1] ;
                    pseudo_count_obsD_ext[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateE_scaled[i+1] *  mHmmResults->mbwd_stateE_scaled[i+1] ;
            }
            else {
                    pseudo_count_obsA_int[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateI_scaled[i+1] *  mHmmResults->mbwd_stateI_scaled[i+1] ;
                    pseudo_count_obsA_ext[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateE_scaled[i+1] *  mHmmResults->mbwd_stateE_scaled[i+1] ;
            }
    }

    for ( std::map<int,double>::iterator it = pseudo_count_obsD_int.begin(); it != pseudo_count_obsD_int.end(); ++it ) {
            int derived = it->first ;
            total_pseudo_count_obs_int[derived] += pseudo_count_obsD_int[derived] ;
            //total_pseudo_count_obs_ext[derived] += pseudo_count_obsD_ext[derived] ;
            if (derived != pseudo_count_obsD_int.rbegin()->first ){
                    total_pseudo_count_obs_ext["SNPs"] += pseudo_count_obsD_ext[derived];
            }
            else {
                    total_pseudo_count_obs_ext["fixed"] += pseudo_count_obsD_ext[derived];
            }
    }

    for ( std::map<int,double>::iterator it = pseudo_count_obsA_int.begin(); it != pseudo_count_obsA_int.end(); ++it ) {
            int derived = it->first ;
            total_pseudo_count_obs_int[derived] += pseudo_count_obsA_int[derived] ;
            //total_pseudo_count_obs_ext[derived] += pseudo_count_obsA_ext[derived] ;
            if (derived != pseudo_count_obsD_int.rbegin()->first ){
                    total_pseudo_count_obs_ext["SNPs"] += pseudo_count_obsA_ext[derived];
            }
            else {
                    total_pseudo_count_obs_ext["fixed"] += pseudo_count_obsA_ext[derived];
            }
    }

    for ( std::map<int,double>::iterator it = total_pseudo_count_obs_int.begin(); it != total_pseudo_count_obs_int.end(); ++it ) {
            int derived = it->first ;
            proba_obsD_int[derived] = pseudo_count_obsD_int[derived] / total_pseudo_count_obs_int[derived] ;
            if (derived != total_pseudo_count_obs_int.rbegin()->first ) {
                    proba_obsD_ext["SNPs"] += pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext["SNPs"] ;
            //	proba_obsD_ext["SNPs"] =  proba_obsD_ext["SNPs"] * ( pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext[derived] ) ;
            }
            else {
                    proba_obsD_ext["fixed"] = pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext["fixed"] ;
                    //proba_obsD_ext["fixed"] = pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext[derived] ;
            }
    }

    // Update the emission probabilities
    mModelProb->setProbFixDInternal(proba_obsD_int.rbegin()->second);
    mModelProb->setProbFixDExternal(proba_obsD_ext["fixed"]);

    //for ( int i=1; i<proba_obsD_int.size(); ++i)
    for (std::map<int,double>::iterator it = proba_obsD_int.begin(); it != proba_obsD_int.end(); ++it)
    {
            int i = it->first ;
            mModelProb->setProbSegregDInternalAt(i,proba_obsD_int[i]);
    }
    //mModelProb->setProbSegregDExternal(proba_obsD_ext["SNPs"]);

}

void baumWelch::computeBaumWelch3states(hmmResults *results ,modelProb *model)
{
    mHmmResults = results;
    mModelProb = model;

    //////////////////////
    // maximum likelihood estimation of the transition and emission probabilities

    map<int,double> pseudo_count_obsD_int ;
    map<int,double> pseudo_count_obsD_ext ;
    map<int,double> pseudo_count_obsA_int ;
    map<int,double> pseudo_count_obsA_ext ;
    map<int,double> total_pseudo_count_obs_int ;
    map<const char*,double> total_pseudo_count_obs_ext ;
    //map<int,double> total_pseudo_count_obs_ext ;

    map<int,double> proba_obsD_int ;
    map<const char*,double> proba_obsD_ext ;

    for ( int i = 1 ; i < mObsSequence->size()-1 ; ++i ) {

            if (mObsSequence->at(i+1).getUlindi() == 'D'){
                    pseudo_count_obsD_int[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateI_scaled[i+1] *  mHmmResults->mbwd_stateI_scaled[i+1] ;
                    pseudo_count_obsD_ext[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateE_scaled[i+1] *  mHmmResults->mbwd_stateE_scaled[i+1] + mHmmResults->mfwd_stateLE_scaled[i+1] *  mHmmResults->mbwd_stateLE_scaled[i+1];
            }
            else {
                    pseudo_count_obsA_int[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateI_scaled[i+1] *  mHmmResults->mbwd_stateI_scaled[i+1] ;
                    pseudo_count_obsA_ext[mObsSequence->at(i+1).getDerived()] += mHmmResults->mfwd_stateE_scaled[i+1] *  mHmmResults->mbwd_stateE_scaled[i+1] + mHmmResults->mfwd_stateLE_scaled[i+1] *  mHmmResults->mbwd_stateLE_scaled[i+1];
            }
    }

    for ( std::map<int,double>::iterator it = pseudo_count_obsD_int.begin(); it != pseudo_count_obsD_int.end(); ++it ) {
            int derived = it->first ;
            total_pseudo_count_obs_int[derived] += pseudo_count_obsD_int[derived] ;
            //total_pseudo_count_obs_ext[derived] += pseudo_count_obsD_ext[derived] ;
            if (derived != pseudo_count_obsD_int.rbegin()->first ){
                    total_pseudo_count_obs_ext["SNPs"] += pseudo_count_obsD_ext[derived];
            }
            else {
                    total_pseudo_count_obs_ext["fixed"] += pseudo_count_obsD_ext[derived];
            }
    }

    for ( std::map<int,double>::iterator it = pseudo_count_obsA_int.begin(); it != pseudo_count_obsA_int.end(); ++it ) {
            int derived = it->first ;
            total_pseudo_count_obs_int[derived] += pseudo_count_obsA_int[derived] ;
            //total_pseudo_count_obs_ext[derived] += pseudo_count_obsA_ext[derived] ;
            if (derived != pseudo_count_obsD_int.rbegin()->first ){
                    total_pseudo_count_obs_ext["SNPs"] += pseudo_count_obsA_ext[derived];
            }
            else {
                    total_pseudo_count_obs_ext["fixed"] += pseudo_count_obsA_ext[derived];
            }
    }

    for ( std::map<int,double>::iterator it = total_pseudo_count_obs_int.begin(); it != total_pseudo_count_obs_int.end(); ++it ) {
            int derived = it->first ;
            proba_obsD_int[derived] = pseudo_count_obsD_int[derived] / total_pseudo_count_obs_int[derived] ;
            if (derived != total_pseudo_count_obs_int.rbegin()->first ) {
                    proba_obsD_ext["SNPs"] += pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext["SNPs"] ;
            //  proba_obsD_ext["SNPs"] =  proba_obsD_ext["SNPs"] * ( pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext[derived] ) ;
            }
            else {
                    proba_obsD_ext["fixed"] = pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext["fixed"] ;
                    //proba_obsD_ext["fixed"] = pseudo_count_obsD_ext[derived] / total_pseudo_count_obs_ext[derived] ;
            }
    }

    // Update the emission probabilities
    mModelProb->setProbFixDInternal(proba_obsD_int.rbegin()->second);
    mModelProb->setProbFixDExternal(proba_obsD_ext["fixed"]);

    //for ( int i=1; i<proba_obsD_int.size(); ++i)
    for (std::map<int,double>::iterator it = proba_obsD_int.begin(); it != proba_obsD_int.end(); ++it)
    {
            int i = it->first ;
            mModelProb->setProbSegregDInternalAt(i,proba_obsD_int[i]);
    }
    //mModelProb->setProbSegregDExternal(proba_obsD_ext["SNPs"]);

}
