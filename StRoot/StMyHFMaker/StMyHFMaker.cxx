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
#include "phys_constants.h"
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
//Constructor
StMyHFMaker::StMyHFMaker(
                const char* name,TString const inputFilesList,
                TString const outBaseName, StPicoDstMaker* picoDstMaker):
                StMaker(name), mPicoDstMaker(picoDstMaker),
                mInputFilesList(inputFilesList), 
                mOutFileBaseName(outBaseName)
            {}
//______________________________________________________________
//Destructor
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
    oFile->Close();
}

//______________________________________________________________
Int_t StMyHFMaker::Make()
{

}





















//______________________________________________________________
bool StMyHFMaker::isGoodEvent(StPicoEvent const* const picoEvent)const{}
//______________________________________________________________
bool StMyHFMaker::isGoodTrigger(StPicoEvent const* const picoEvent)const{}
//______________________________________________________________
bool StMyHFMaker::isGoodTrack(StPicoTrack const* trk)const{}
//______________________________________________________________
bool StMyHFMaker::isElectron(StPicoTrack const* trk)const{}
//______________________________________________________________
bool StMyHFMaker::isPion(StPicoTrack const* trk)const{}
//______________________________________________________________
bool StMyHFMaker::isKaon(StPicoTrack const* trk)const{}
//______________________________________________________________
void StMyHFMaker::initHistograms(){}
//______________________________________________________________
void StMyHFMaker::initNTuples(){}