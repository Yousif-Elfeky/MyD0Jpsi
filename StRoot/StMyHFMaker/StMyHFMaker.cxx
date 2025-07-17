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
//Constructor
StMyHFMaker::StMyHFMaker(
                const char* name,TString const inputFilesList,
                TString const outBaseName, StPicoDstMaker* picoDstMaker):
                StMaker(name), mPicoDstMaker(picoDstMaker),
                mInputFilesList(inputFilesList), 
                mOutFileBaseName(outBaseName)
            {}
//Destructor
StMyHFMaker::~StMyHFMaker()
{}

Int_t StMyHFMaker::Init()
{
  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  oFile = new TFile(mOutFileBaseName+".root", "RECREATE");
  initHistograms();
  initNTuples();
  
  return kStOK;
}
Int_t StMyHFMaker::Finish(){
    oFile->cd();
    /*Write Histograms/NTuples*/
    oFile->Close();
}

Int_t StMyHFMaker::Make()
{

}





















bool StMyHFMaker::isGoodEvent(StPicoEvent const* const picoEvent)const{}
bool StMyHFMaker::isGoodTrigger(StPicoEvent const* const picoEvent)const{}
bool StMyHFMaker::isGoodTrack(StPicoTrack const* trk)const{}
bool StMyHFMaker::isElectron(StPicoTrack const* trk)const{}
bool StMyHFMaker::isPion(StPicoTrack const* trk)const{}
bool StMyHFMaker::isKaon(StPicoTrack const* trk)const{}
void StMyHFMaker::initHistograms(){}
void StMyHFMaker::initNTuples(){}