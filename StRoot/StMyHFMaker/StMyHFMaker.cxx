/*
**************************************************
* A Maker to Analyze Heavy Flavour Physics
* And Reconstruct D0 and JPSI Masons Based on 
* StPicoDstarMixedMaker class.
* Authors: Yousif El-Feky <yousifjoo@aucegypt.edu> 
**************************************************
*/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh" 
// #include "phys_constants.h" // commented for now
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StMemStat.h"

#include "StMyCuts.h"
#include "StMyHFMaker.h"


ClassImp(StMyHFMaker)
//______________________________________________________________
//  Constructor
StMyHFMaker::StMyHFMaker(
                const char* name,TString const inputFilesList,
                TString const outBaseName, StPicoDstMaker* picoDstMaker):
                StMaker(name), mPicoDstMaker(picoDstMaker),
                mInputFilesList(inputFilesList), 
                mOutFileBaseName(outBaseName)
            {}
//______________________________________________________________
// Destructor
StMyHFMaker::~StMyHFMaker()
{}

//______________________________________________________________
Int_t StMyHFMaker::Init()
{
  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  oFile = new TFile(mOutFileBaseName+".root", "RECREATE");
  initHistograms();
  initNTuples();
  
  return kStOK;
}
//______________________________________________________________
Int_t StMyHFMaker::Finish(){
    oFile->cd();
    /*Write Histograms/NTuples*/
    writeHistograms();
    oFile->Close();
    return kStOK;
}

//______________________________________________________________
Int_t StMyHFMaker::Make()
{


    return kStOK;
}




//______________________________________________________________
bool StMyHFMaker::isGoodEvent(StPicoEvent const* const picoEvent)const{
    TVector3 pVer = picoEvent->primaryVertex();
    return (pVer.z() < EventCuts::vZ_max && pVer.z() > EventCuts::vZ_min) &&
            (fabs(pVer.z() - picoEvent->vzVpd()) < EventCuts::vZVpdvZ) &&
            !((fabs(pVer.x()) < EventCuts::vError) &&
              (fabs(pVer.y()) < EventCuts::vError) &&
              (fabs(pVer.z()) < EventCuts::vError))&&
            (sqrt(TMath::Power(pVer.x(),2) + TMath::Power(pVer.y(),2))<=EventCuts::vR);
}//Check StMyCuts.h
//______________________________________________________________
bool StMyHFMaker::isGoodTrigger(StPicoEvent const* const picoEvent)const{
  for (auto trg : EventCuts::triggers)
  {
    if (picoEvent->isTrigger(trg)) return true;
    else return false;
  }
}//Check StMyCuts.h
//______________________________________________________________
bool StMyHFMaker::isGoodTrack(StPicoTrack const* trk)const{
    return ((trk->gPt() > TrackCuts::gPt)&&
            ((trk->gMom().Eta()) < TrackCuts::Eta)&&
            (trk->nHitsFit() > TrackCuts::nHitsFit)&&
            (trk->nHitsDedx() > TrackCuts::nHitsDedx)&&
            (((trk->nHitsFit())/(trk->nHitsDedx())) >= 
                TrackCuts::nHitsFit2Dedx)
        );
}//Check StMyCuts.h
//______________________________________________________________
bool StMyHFMaker::isElectron(StPicoTrack const* trk)const{}
//______________________________________________________________
bool StMyHFMaker::isPion(StPicoTrack const* trk)const{}
//______________________________________________________________
bool StMyHFMaker::isKaon(StPicoTrack const* trk)const{}
//______________________________________________________________
void StMyHFMaker::initHistograms(){
    hNevent = new TH1D("hNevnet","hNevnet",1000,0,1000);
}
//______________________________________________________________
void StMyHFMaker::initNTuples(){}
//______________________________________________________________
void StMyHFMaker::writeHistograms(){
    hNevent->Write();
}