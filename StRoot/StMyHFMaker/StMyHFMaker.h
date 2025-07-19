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
#include "TH2.h"
#include "TProfile.h"
#include "TString.h"
#include "TVector3.h"
#include "StMaker.h"
#include <vector>

class TString;
class TFile;
class TNtuple;
class TVector3;
class StPicoDst;
class StPicoTrack;
class StPicoDstMaker;
class StPicoEvent;
class TH1D;
class TH2D;
class TProfile;

struct Particle{
    Short_t charge;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t Energy;
    Float_t Eta;
    Float_t Phi;
    Float_t Pt;
    UInt_t trackId;
};

const int N_CENT_BINS = 9; // For 0-10, 10-20, ..., 70-80%
const int N_PT_BINS_JPSI = 10;
const int N_PT_BINS_D0 = 10;

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
        // void writeNTuples();
        void pairElectrons(StPicoTrack const* trk);
        void pairPions(StPicoTrack const* trk);
        void pairKaons(StPicoTrack const* trk);
        void makeJPSI();
        void makeD0(StPicoDst const* picoDst, TVector3 const& pVtx,
                         std::vector<unsigned int> const& kaonIndices,
                         std::vector<unsigned int> const& pionIndices);
        int getTotalNRuns();
        double getTofBeta(StPicoTrack const* const trk) const;
        bool isGoodEvent(StPicoEvent const* const picoEvent)const;
        bool isGoodTrigger(StPicoEvent const* const picoEvent)const;
        bool isGoodTrack(StPicoTrack const* trk)const;
        bool isElectron(StPicoTrack const* trk, bool tofMatch, float beta, float DCA)const;
        bool isPion(StPicoTrack const* trk, bool tofMatch, float beta, float DCA)const;
        bool isKaon(StPicoTrack const* trk, bool tofMatch, float beta, float DCA)const;
        bool isBadrun(Int_t runId);
        float calcEventPlane(StPicoDst const* const picoDst, StPicoEvent const* picoEvent, const int n) const;
        int mCentralityBin;
        float mEventPlane1;
        float mEventPlane2;

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
        TVector3 mom;
        float M_ELECTRON=0.000511;//GeV
        float M_PION=0.139570;//GeV
        float M_KAON=0.493677;//GeV
        float dca_to_pv;
        // NTuples Here.
        TNtuple* mPion;
        TNtuple* mKaon;
        TNtuple* mElectron;
        TNtuple* mPositron;
        TNtuple* mD0;
        TNtuple* mJPSI;
        
        // Histograms Here.
        Int_t xBins = 40000;
        Int_t yBins = 40000;
        TH1D* hNevent;
        TH1D* hEventPlane;
        TH1D* hVzTPC;
        TH1D* hVzVPD;
        TH1D* hVr;
        TH1D* hNsigmaElectron;
        TH1D* hMee_ULike;
        TH1D* hMee_Like1;
        TH1D* hMee_Like2;
        TH1D* hMpik_ULike1;
        TH1D* hMpik_ULike2;
        TH1D* hMpik_Like1;
        TH1D* hMpik_Like2;
        TH2D* hMassPt_Jpsi[N_CENT_BINS];
        TProfile* pV1vsPt_Jpsi[N_CENT_BINS];
        TProfile* pV2vsPt_Jpsi[N_CENT_BINS];
        TH2D* hMassPt_D0[N_CENT_BINS];
        TProfile* pV1vsPt_D0[N_CENT_BINS];
        TProfile* pV2vsPt_D0[N_CENT_BINS];
        // Things
        int  mRunId;
        Particle particleinfo;
        vector<Particle> electroninfo;
        vector<Particle> positroninfo;
        vector<Particle> pionplusinfo;
        vector<Particle> pionminusinfo;
        vector<Particle> kaonplusinfo;
        vector<Particle> kaonminusinfo;


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
