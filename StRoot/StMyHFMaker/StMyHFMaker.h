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
#include "TString.h"
#include "StMaker.h"

class TString;
class TFile;
class TNtuple;
class StPicoTrack;
class StPicoDstMaker;
class StPicoEvent;

class StMyHFMaker : public StMaker
{
    public:
        StMyHFMaker(const char* name,TString const inputFilesList,
            TString const outBaseName, StPicoDstMaker* picoDstMaker);    
        virtual ~StMyHFMaker();
        
        virtual Int_t Init();
        virtual Int_t Make();
        virtual Int_t Finish();
        void getBadruns(string inputFileName);
    
    private:
        StMyHFMaker(){}
        void initHistograms();
        void initNTuples();
        bool isGoodEvent(StPicoEvent const* const picoEvent)const;
        bool isGoodTrigger(StPicoEvent const* const picoEvent)const;
        bool isGoodTrack(StPicoTrack const* trk)const;
        bool isPion(StPicoTrack const* trk)const;
        bool isKaon(StPicoTrack const* trk)const;
        bool isElectron(StPicoTrack const* trk)const;
        float getTofBeta(StPicoTrack const* const trk) const;
        bool isBadrun(Int_t runId);
        
        StPicoDstMaker* mPicoDstMaker;
        TString mInputFilesList;
        TString mOutFileBaseName;
        TFile* oFile;
        std::map<int,int> runnum;
        string mRunNumList;
        vector<int> mBadRun;

        // Put Histograms and NTuples Here.
        TNtuple* mPion;
        TNtuple* mKaon;
        TNtuple* mElectron;
        TNtuple* mPositron;
        TNtuple* mD0;
        TNtuple* mJPSI;

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
