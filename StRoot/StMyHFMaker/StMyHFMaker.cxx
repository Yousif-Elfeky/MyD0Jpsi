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
#include "../StPicoCharmContainers/StKaonPion.h"

#ifndef C_C_LIGHT
#define C_C_LIGHT 299792458
#endif // !C_C_LIGHT

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

  // Write Histograms/NTuples
  writeHistograms();
  
  oFile->Close();
  
  deleteHistograms();

  return kStOK;
}

//______________________________________________________________
Int_t StMyHFMaker::Make()
{
  if (!mPicoDstMaker)
  {
  LOG_WARN << "No PicoDstMaker! Skip! " << endm;
  return kStWarn;
  }
  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
  LOG_WARN << "No PicoDst! Skip! " << endm;
  return kStWarn;
  }  
  StPicoEvent const * picoEvent = picoDst->event();
  if (!picoEvent) {
    LOG_WARN << "No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  electroninfo.clear();
  positroninfo.clear();  

  mRunId = picoEvent->runId();
  hNevent->Fill(runnum[mRunId]);
  if(!isGoodEvent(picoEvent))return kStOK;
  
  TPCVer = picoEvent->primaryVertex();hVzTPC->Fill(TPCVer.z());
  float vz = picoEvent->primaryVertex().z();
  VPDvz = picoEvent->vzVpd();hVzVPD->Fill(VPDvz);
  Vr = std::sqrt(TMath::Power(TPCVer.x(),2)+TMath::Power(TPCVer.y(),2));hVr->Fill(Vr);
  
  mCentralityBin = getCentralityBin(picoEvent->grefMult());
    if (mCentralityBin < 0) {
      return kStOK;
  }

  nTracks = picoDst->numberOfTracks();
  std::vector<unsigned int> idxPicoPions;
  std::vector<unsigned int> idxPicoKaons;

  mEventPlane1 = calcEventPlane(picoDst, picoEvent, 1); // For v1
  mEventPlane2 = calcEventPlane(picoDst, picoEvent, 2); // For v2
  float eventPlane = calcEventPlane(picoDst, picoEvent, 2);
  hEventPlane->Fill(eventPlane);
  // Track Loop
  for (UInt_t itrack=0;itrack<nTracks;itrack++){
    StPicoTrack* trk = picoDst->track(itrack);
    mom = trk->pMom();
    dca_to_pv = trk->gDCA(TPCVer.x(), TPCVer.y(), TPCVer.z());
    if(!isGoodTrack(trk))continue;
    
    beta = getTofBeta(trk);
    tofmatch = (beta!=std::numeric_limits<float>::quiet_NaN()) && beta>0;
    

    if(isElectron(trk,tofmatch,beta,dca_to_pv)){
      hNsigmaElectron->Fill(trk->nSigmaElectron());
      pairElectrons(trk);
      
    }
    if(isPion(trk,tofmatch,beta,dca_to_pv)){
      idxPicoPions.push_back(itrack);
    }
    if(isKaon(trk,tofmatch,beta,dca_to_pv)){
      idxPicoKaons.push_back(itrack);
    }
  }
  makeJPSI();
  makeD0(picoDst,TPCVer,idxPicoKaons,idxPicoPions);

  return kStOK;
}
//______________________________________________________________
bool StMyHFMaker::isElectron(StPicoTrack const* trk, bool tofMatch, float beta, float DCA) const {
    bool isTPCElectron = false;
    if (trk->pMom().Mag() > 0.8) {
        isTPCElectron = trk->nSigmaElectron() < JPSI_Cuts::nSigmaElectron_max &&
                        trk->nSigmaElectron() > JPSI_Cuts::nSigmaElectron_min;
    } else {
        isTPCElectron = trk->nSigmaElectron() < JPSI_Cuts::nSigmaElectron_max &&
                        trk->nSigmaElectron() > (3 * trk->pMom().Mag() + JPSI_Cuts::nSigmaElectron_lowhmom);
    }

    bool isTOFElectron = tofMatch ? std::abs(1. / beta - 1.) < JPSI_Cuts::oneOverBetaElectron : false;

    return isTPCElectron && isTOFElectron && DCA < 3.0;
}
//______________________________________________________________
bool StMyHFMaker::isPion(StPicoTrack const* trk, bool tofMatch, float beta, float DCA) const {
  
  bool tpcPion = std::abs(trk->nSigmaPion()) < D0_Cuts::nSigmaPion;
  if (!tpcPion) return false;
  if(tofMatch)
  {
    float p = trk->gMom().Mag();
    if (p < 1.e-9) return false;
    float pion_beta_expected = p / sqrt(p * p + M_PION * M_PION);
    return std::abs(1. / beta - 1. / pion_beta_expected) < D0_Cuts::oneOverBetaPion &&
            DCA > D0_Cuts::DCA_pi;
  }
  return tpcPion && tofMatch;
}
//______________________________________________________________
bool StMyHFMaker::isKaon(StPicoTrack const* trk, bool tofMatch, float beta, float DCA) const {
  
  bool tpcKaon = std::abs(trk->nSigmaKaon()) < D0_Cuts::nSigmaKaon;
  if (!tpcKaon) return false;
  if(tofMatch)
  {
    float p = trk->gMom().Mag();
    if (p < 1.e-9) return false;
    float kaon_beta_expected = p / sqrt(p * p + M_KAON * M_KAON);
    return std::abs(1. / beta - 1. / kaon_beta_expected) < D0_Cuts::oneOverBetaKaon &&
            DCA > D0_Cuts::DCA_k;
  } 
  return tpcKaon && tofMatch;
}
//______________________________________________________________
void StMyHFMaker::pairElectrons(StPicoTrack const* trk){
  particleinfo.charge = trk->charge();
  particleinfo.px = mom.Px();
  particleinfo.py = mom.Py();
  particleinfo.pz = mom.Pz();
  particleinfo.Energy = sqrt(pow(M_ELECTRON,2.0)+pow(mom.Mag(),2.0));
  particleinfo.Eta = mom.Eta();
  particleinfo.Phi = mom.Phi();
  particleinfo.Pt = mom.Perp();
  particleinfo.trackId = trk->id();
  if(trk->charge()<0) electroninfo.push_back(particleinfo);
  else if(trk->charge()>0) positroninfo.push_back(particleinfo);
}
//______________________________________________________________
void StMyHFMaker::pairPions(StPicoTrack const* trk){
  particleinfo.charge = trk->charge();
  particleinfo.px = mom.Px();
  particleinfo.py = mom.Py();
  particleinfo.pz = mom.Pz();
  particleinfo.Energy = sqrt(pow(M_PION,2.0)+pow(mom.Mag(),2.0));
  particleinfo.Eta = mom.Eta();
  particleinfo.Phi = mom.Phi();
  particleinfo.Pt = mom.Perp();
  particleinfo.trackId = trk->id();
  if(trk->charge()>0) pionplusinfo.push_back(particleinfo);
  else if(trk->charge()<0) pionminusinfo.push_back(particleinfo);
}
//______________________________________________________________
void StMyHFMaker::pairKaons(StPicoTrack const* trk){
  particleinfo.charge = trk->charge();
  particleinfo.px = mom.Px();
  particleinfo.py = mom.Py();
  particleinfo.pz = mom.Pz();
  particleinfo.Energy = sqrt(pow(M_KAON,2.0)+pow(mom.Mag(),2.0));
  particleinfo.Eta = mom.Eta();
  particleinfo.Phi = mom.Phi();
  particleinfo.Pt = mom.Perp();
  particleinfo.trackId = trk->id();
  if(trk->charge()>0) kaonplusinfo.push_back(particleinfo);
  else if(trk->charge()<0) kaonminusinfo.push_back(particleinfo);
}
//______________________________________________________________
void StMyHFMaker::makeJPSI(){

  TLorentzVector pair_4v(0,0,0,0), p1_4v(0,0,0,0), p2_4v(0,0,0,0);

  // 1. Unlike-Sign Signal Pairs (e+ e-)
  for(const auto& electron : electroninfo) 
  {
    p1_4v.SetPxPyPzE(electron.px, electron.py, electron.pz, electron.Energy);
    for(const auto& positron : positroninfo) 
    {
      p2_4v.SetPxPyPzE(positron.px, positron.py, positron.pz, positron.Energy);
      pair_4v = p1_4v + p2_4v;
      hMee_ULike->Fill(pair_4v.M());
      if (pair_4v.M() > 2.9 && pair_4v.M() < 3.2) {
      float jpsi_pt = pair_4v.Perp();
      float jpsi_phi = pair_4v.Phi();
      hMassPt_Jpsi[mCentralityBin]->Fill(pair_4v.M(), jpsi_pt);
      pV1vsPt_Jpsi[mCentralityBin]->Fill(jpsi_pt, TMath::Cos(1. * (jpsi_phi - mEventPlane1)));
      pV2vsPt_Jpsi[mCentralityBin]->Fill(jpsi_pt, TMath::Cos(2. * (jpsi_phi - mEventPlane2)));
  }
    }
  }
  // 2. Like-Sign Background: Positron Pairs (e+ e+)
  for(uint i = 0; i < positroninfo.size(); ++i) 
  {
    p1_4v.SetPxPyPzE(positroninfo[i].px, positroninfo[i].py, positroninfo[i].pz, positroninfo[i].Energy);
    for(uint j = i + 1; j < positroninfo.size(); ++j) 
    { 
      p2_4v.SetPxPyPzE(positroninfo[j].px, positroninfo[j].py, positroninfo[j].pz, positroninfo[j].Energy);
      pair_4v = p1_4v + p2_4v;
      hMee_Like1->Fill(pair_4v.M());
    }
  }
  // 3. Like-Sign Background: Electron Pairs (e- e-)
  for(uint i = 0; i < electroninfo.size(); ++i) {
    p1_4v.SetPxPyPzE(electroninfo[i].px, electroninfo[i].py, electroninfo[i].pz, electroninfo[i].Energy);
    for(uint j = i + 1; j < electroninfo.size(); ++j) 
    { 
      p2_4v.SetPxPyPzE(electroninfo[j].px, electroninfo[j].py, electroninfo[j].pz, electroninfo[j].Energy);
      pair_4v = p1_4v + p2_4v;
      hMee_Like2->Fill(pair_4v.M());
    }
  }
}
//______________________________________________________________
void StMyHFMaker::makeD0(StPicoDst const* picoDst, TVector3 const& pVtx,
                         std::vector<unsigned int> const& kaonIndices,
                         std::vector<unsigned int> const& pionIndices) {

float const bField = picoDst->event()->bField();

for (unsigned int iK = 0; iK < kaonIndices.size(); ++iK) {
    StPicoTrack const* kaonTrack = picoDst->track(kaonIndices[iK]);
    if (!kaonTrack) continue;

    for (unsigned int iPi = 0; iPi < pionIndices.size(); ++iPi) {
      if (kaonIndices[iK] == pionIndices[iPi]) continue;

      StPicoTrack const* pionTrack = picoDst->track(pionIndices[iPi]);
      if (!pionTrack) continue;

      StKaonPion kaonPion(*kaonTrack, *pionTrack, pVtx, bField);

      if ((kaonPion.dcaDaughters() < D0_Cuts::DCA_12) &&
          (kaonPion.decayLength() > D0_Cuts::DecayLength) &&
          (cos(kaonPion.pointingAngle()) > D0_Cuts::cos_theta) &&
          (kaonPion.kaonDca() > D0_Cuts::DCA_k) &&
          (kaonPion.pionDca() > D0_Cuts::DCA_pi))
      {
        int k_charge = kaonTrack->charge();
        int pi_charge = pionTrack->charge();

        // Unlike-Sign: K- pi+
        if (k_charge < 0 && pi_charge > 0) {
          hMpik_ULike1->Fill(kaonPion.m());
          if (kaonPion.m() > 1.82 && kaonPion.m() < 1.91) {
          float d0_pt = kaonPion.pt();
          float d0_phi = kaonPion.phi();
          hMassPt_D0[mCentralityBin]->Fill(kaonPion.m(), d0_pt);
          pV1vsPt_D0[mCentralityBin]->Fill(d0_pt, TMath::Cos(1. * (d0_phi - mEventPlane1)));
          pV2vsPt_D0[mCentralityBin]->Fill(d0_pt, TMath::Cos(2. * (d0_phi - mEventPlane2)));
          }
        }
        // Unlike-Sign: K+ pi-
        else if (k_charge > 0 && pi_charge < 0) {
          hMpik_ULike2->Fill(kaonPion.m());
          if (kaonPion.m() > 1.82 && kaonPion.m() < 1.91) {
          float d0_pt = kaonPion.pt();
          float d0_phi = kaonPion.phi();
          hMassPt_D0[mCentralityBin]->Fill(kaonPion.m(), d0_pt);
          pV1vsPt_D0[mCentralityBin]->Fill(d0_pt, TMath::Cos(1. * (d0_phi - mEventPlane1)));
          pV2vsPt_D0[mCentralityBin]->Fill(d0_pt, TMath::Cos(2. * (d0_phi - mEventPlane2)));
          }
        }
        // Like-Sign: K+ pi+
        else if (k_charge > 0 && pi_charge > 0) {
          hMpik_Like1->Fill(kaonPion.m());
        }
        // Like-Sign: K- pi-
        else if (k_charge < 0 && pi_charge < 0) {
          hMpik_Like2->Fill(kaonPion.m());
        }
      }
    }
  }
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
  }
  return false;
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
void StMyHFMaker::initNTuples(){}
//______________________________________________________________
int StMyHFMaker::getTotalNRuns(){
    
  ifstream readnum;
  readnum.open(mRunNumList);

  if (!readnum.is_open()) {cout << "Error: Could not open run number list file: " << mRunNumList << endl;return EXIT_FAILURE; }
  int tmpRunNum; 
  int index = 0; 

  while (readnum >> tmpRunNum) {
    runnum.insert(pair<int, int>(tmpRunNum, index));
    if (DEBUG) cout << "Read run number: " << tmpRunNum << " -> assigned id: " << index << endl;
    index++;
  }
  readnum.close();

  return runnum.size();
}
//______________________________________________________________
double StMyHFMaker::getTofBeta(StPicoTrack const* const trk) const
{
  int index2tof = trk->bTofPidTraitsIndex();
  double beta = std::numeric_limits<float>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        TVector3 const vtx3 = mPicoDstMaker->picoDst()->event()->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
        StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  } 
  return beta;
}
//______________________________________________________________
float StMyHFMaker::calcEventPlane(StPicoDst const* const picoDst, StPicoEvent const* picoEvent, const int n) const
{
	float cossum_nocorrection = 0.;
	float sinsum_nocorrection = 0.;
	TVector3 pVtx = picoEvent->primaryVertex();
	const int nTrack = picoDst->numberOfTracks();
	for (int iTrack = 0; iTrack < nTrack; iTrack++) {
		StPicoTrack* mTrack = (StPicoTrack*)picoDst->track(iTrack);
		if (!mTrack) continue;
		float dca = mTrack->gDCA(pVtx.x(), pVtx.y(), pVtx.z());
    if(!isGoodTrack(mTrack))continue;
		const float pt = mTrack->pMom().Perp();
		const float phi = mTrack->pMom().Phi();
		const float cos_part_nocorrection = pt * cos(n*phi);
		const float sin_part_nocorrection = pt * sin(n*phi);
		cossum_nocorrection += cos_part_nocorrection;
		sinsum_nocorrection += sin_part_nocorrection;
	}
	TVector2 Q_nocorrection(cossum_nocorrection, sinsum_nocorrection);
	float eventPlane_nocorrection = Q_nocorrection.Phi();
	if (eventPlane_nocorrection < 0) eventPlane_nocorrection += 2.0*TMath::Pi();
	return eventPlane_nocorrection / n;
}
//______________________________________________________________
int StMyHFMaker::getCentralityBin(int grefmult) const {
  int grefmult_cuts_10_percent[] = {256, 185, 132, 93, 63, 40, 24, 14};
    if (grefmult > grefmult_cuts_10_percent[0]) return 0;
  for (int i = 0; i < 7; ++i) { 
    if (grefmult > grefmult_cuts_10_percent[i+1] && grefmult <= grefmult_cuts_10_percent[i]) {
      return i + 1;
    }
  }    
  if (grefmult > grefmult_cuts_10_percent[7] && grefmult <= grefmult_cuts_10_percent[6]) {
    return 7;
  }
  return -1;
}
//______________________________________________________________
void StMyHFMaker::initHistograms(){
  
  int nRuns = getTotalNRuns();

  hNevent = new TH1D("hNevnet","Number of Events",nRuns,0,nRuns);                      //TODO: Declare the limits of 
  hVzTPC = new TH1D("hVzTPC","TPC_{Vz}",200,EventCuts::vZ_min-20,EventCuts::vZ_max+20);//the histograms for better 
  hVzVPD = new TH1D("hVzVPD","VPD_{Vz}",200,EventCuts::vZ_min-20,EventCuts::vZ_max+20);//readability 
  hEventPlane = new TH1D("hEventPlane", "hEventPlane", 200, -0.1, 3.3);
  hVr = new TH1D("hVr","V_{r}",200,0,EventCuts::vR+1);
  hNsigmaElectron = new TH1D("hNsigmaElectron", "n^{#sigma}_{e}",100,-5,5);
  hMee_ULike = new TH1D("hMee_ULike","e^{+} e^{-} + Background",xBins,0,4);
  hMee_Like1 = new TH1D("hMee_Like1","e^{+} e^{+}",xBins,0,4);
  hMee_Like2 = new TH1D("hMee_Like2","e^{-} e^{-}",xBins,0,4);
  hMpik_ULike1 = new TH1D("hMpik_ULike1","K^{-} #pi^{+}",xBins,0,3);
  hMpik_ULike2 = new TH1D("hMpik_ULike2","K^{+} #pi^{-}",xBins,0,3);
  hMpik_Like1 = new TH1D("hMpik_Like1","K^{+} #pi^{+}",xBins,0,3);
  hMpik_Like2 = new TH1D("hMpik_Like2","K^{-} #pi^{-}",xBins,0,3);
}
//______________________________________________________________
void StMyHFMaker::writeHistograms(){
  hNevent->Write();
  hVzTPC->Write();
  hEventPlane->Write();
  hVzVPD->Write();
  hVr->Write();
  hNsigmaElectron->Write();
  hMee_ULike->Write();
  hMee_Like1->Write();
  hMee_Like2->Write();
  hMpik_ULike1->Write();
  hMpik_ULike2->Write();
  hMpik_Like1->Write();
  hMpik_Like2->Write();
  for (int i = 0; i < N_CENT_BINS; ++i) {
    hMassPt_Jpsi[i]->Write();
    pV1vsPt_Jpsi[i]->Write();
    pV2vsPt_Jpsi[i]->Write();
    hMassPt_D0[i]->Write();
    pV1vsPt_D0[i]->Write();
    pV2vsPt_D0[i]->Write();
  }
}
//______________________________________________________________
void StMyHFMaker::deleteHistograms(){
  delete hNevent;
  delete hVzTPC;
  delete hEventPlane;
  delete hVzVPD;
  delete hVr;
  delete hNsigmaElectron;
  delete hMee_ULike;
  delete hMee_Like1;
  delete hMee_Like2;
  delete hMpik_ULike1;
  delete hMpik_ULike2;
  delete hMpik_Like1;
  delete hMpik_Like2;
  for (int i = 0; i < N_CENT_BINS; ++i) {
    delete hMassPt_Jpsi[i];
    delete pV1vsPt_Jpsi[i];
    delete pV2vsPt_Jpsi[i];
    delete hMassPt_D0[i];
    delete pV1vsPt_D0[i];
    delete pV2vsPt_D0[i];
  }
}
