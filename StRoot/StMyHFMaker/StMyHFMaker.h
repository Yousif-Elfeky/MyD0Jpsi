#ifndef StMyHFMaker_h
#define StMyHFMaker_h
/*
**************************************************
* A Maker to Analyze Heavy Flavour Physics
* And Reconstruct D0 and JPSI Masons Based on 
* StPicoDstarMixedMaker class.
* Authors: Yousif El-Feky <yousifjoo@aucegypt.edu> 
**************************************************
*/
#include "TChain.h"
#include "TH1.h"
#include "TString.h"
#include "TVector3.h"
#include "StMaker.h"

class TString;
class TFile;
class TNtuple;
class TVector3;
class StPicoTrack;
class StPicoDstMaker;
class StPicoEvent;

struct Particle{
    Short_t charge;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t Energy;
    Float_t Eta;
    Float_t Phi;
};

class StMyHFMaker : public StMaker
{
    public:
        StMyHFMaker(const char* name,TString const inputFilesList,
            TString const outBaseName, StPicoDstMaker* picoDstMaker);    
        virtual ~StMyHFMaker();
        
        virtual Int_t Init();
        virtual Int_t Make();
        virtual Int_t Finish();
        void setRunNumList(string list){mRunNumList = list;}
        void getBadruns(string inputFileName);
        void setDebug(bool d){DEBUG = d;}
    
    private:
        StMyHFMaker(){}
        void initHistograms();
        void writeHistograms();
        void deleteHistograms();
        void initNTuples();
        int getTotalNRuns();
        // void writeNTuples();
        bool isGoodEvent(StPicoEvent const* const picoEvent)const;
        bool isGoodTrigger(StPicoEvent const* const picoEvent)const;
        bool isGoodTrack(StPicoTrack const* trk)const;
        bool isPion(StPicoTrack const* trk)const;
        bool isKaon(StPicoTrack const* trk)const;
        bool isElectron(StPicoTrack const* trk)const;
        double getTofBeta(StPicoTrack const* const trk) const;
        bool isBadrun(Int_t runId);
        
        StPicoDstMaker* mPicoDstMaker;
        TString mInputFilesList;
        TString mOutFileBaseName;
        TFile* oFile;
        std::map<int,int> runnum;
        string mRunNumList;
        vector<int> mBadRun;
        bool DEBUG;
        
        //Physics
        TVector3 TPCVer;
        Float_t VPDvz;
        Float_t Vr;
        UInt_t nTracks;
        Double_t beta;
        bool tofmatch;
        bool isTPCElectron;
        bool isTOFElectron;
        // NTuples Here.
        TNtuple* mPion;
        TNtuple* mKaon;
        TNtuple* mElectron;
        TNtuple* mPositron;
        TNtuple* mD0;
        TNtuple* mJPSI;
        
        // Histograms Here.
        Int_t xBins = 10000;
        Int_t yBins = 10000;
        TH1D* hNevent;
        TH1D* hVzTPC;
        TH1D* hVzVPD;
        TH1D* hVr;
        // Things
        int  mRunId;

ClassDef(StMyHFMaker, 1)
};

inline void StMyHFMaker::getBadruns(string inputFileName){
    ifstream fin(inputFileName.c_str());
    if(!fin){
      cout <<"no Bad runs list" << endl;
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t runId = 0 ;
    while( fin >> runId ) {
      mBadRun.push_back(runId);
    }
    cout << "get Bad runs list [OK]" << endl;
}
inline  bool StMyHFMaker::isBadrun(Int_t runId){
    vector<Int_t>::iterator iter = std::find(mBadRun.begin(), mBadRun.end(), runId);
    return ( iter != mBadRun.end() ) ;
}

#endif // !StMyHFMaker_h
