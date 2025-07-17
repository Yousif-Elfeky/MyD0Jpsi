#include <string>
#include <fstream>
void load();
void runStMyHFMaker(TString picolist="input.list" , TString outFileName="test")
{
  TStopwatch*   stopWatch = new TStopwatch();
  stopWatch->Start();

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  // load();
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StMyHFMaker");

  chain = new StChain();
  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(2, picolist, "picoDstMaker");
  StMyHFMaker*  picoHFMaker = new StMyHFMaker("picoHFMaker",picolist , outFileName.Data(), picoDstMaker);
  picoHFMaker->getBadruns("StRoot/macros/badRun.list");
  picoHFMaker->setDebug(true);
  // -------------- USER variables -------------------------
  chain->Init();
  int nEntries = picoDstMaker->chain()->GetEntries();
  cout<<"Processing "<<nEntries<<" events...\n";
  for (int iEvent = 0; iEvent < nEntries; ++iEvent)
  {
    chain->Clear();
    if(iEvent && iEvent%1000 == 0) cout<<"... finished processing "<<iEvent<<" events.\n";

    int iret = chain->Make();
    if (iret)
    {
      cout << "Bad return code!" << iret << '\n';
      break;
    }
  }
  cout<<"Finished processing "<<nEntries<<" events.\n";

  chain->Finish();
  delete chain;

  stopWatch->Stop();   
  stopWatch->Print();
}
void load(){
  gSystem->Load("StTpcDb");
  gSystem->Load("StEvent");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("libgen_Tables");
  gSystem->Load("libsim_Tables");
  gSystem->Load("libglobal_Tables");
  gSystem->Load("StMagF");

  gSystem->Load("St_g2t.so");
  gSystem->Load("St_geant_Maker.so");
  gSystem->Load("StAssociationMaker");
  gSystem->Load("StMcAnalysisMaker");
  gSystem->Load("libgeometry_Tables");   
  gSystem->Load("StTriggerUtilities");

  gSystem->Load("StEmcUtil");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcSimulatorMaker");

  gSystem->Load("StDbLib");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("St_db_Maker");

  gSystem->Load("StMtdHitMaker");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StMtdCalibMaker");
  gSystem->Load("StBTofUtil");
  gSystem->Load("StVpdCalibMaker");

}
