#ifndef StMyCuts_hh
#define StMyCuts_hh
#include "Rtypes.h"

namespace EventCuts{
    const int vZ_min = -30; //cm
    const int vZ_max = 30; //cm
    const int vZVpdvZ = 2; //|TPC_vZ - VPD_Vz|(cm)
    const int vR = 2; //cm
    const float vError = 1.0e-5; //cm
    std::array<unsigned int, 3> const triggers = {1,2,3}; //PlaceHolder Triggers
}//EventCuts
namespace TrackCuts{
    const float gPt = 0.2; // Global pT > 0.2
    const int nHitsFit = 15;  // nHitsFit > 15
    const int nHitsDedx = 5; // nHitsDedx > 5
    const float nHitsFit2Dedx = 0.51; // ((nHitsFit)/(nHitsDedx)) >= 0.51
    const int DCA = 1; // DCA < 1
    const int Eta = 1; // |η| < 1
}//TrackCuts
namespace D0_Cuts{

}//D0_Cuts
namespace JPSI_Cuts{
    float nSigmaElectron_min = -0.75; //nSigmaElectron > -0.75 if mom>0.8 GeV
    float nSigmaElectron_max = 2;  //nSigmaElectron < 2 for all mom
    float nSigmaElectron_lowhmom = -3.15; // nSigmaElectron > 3*mom.Mag()-3.15 if mom<0.8 GeV
    float oneOverBetaElectron = 0.025; // (1/(β-1)) < 0.025
}//JPSI_Cuts
#endif // !StMyCuts_hh
