#ifndef StMyCuts_hh
#define StMyCuts_hh
#include "Rtypes.h"

namespace EventCuts{
    const int vZ_min = -100; //cm
    const int vZ_max = 100; //cm
    const int vZVpdvZ = 6; //|TPC_vZ - VPD_Vz|(cm)
    const int vR = 2; //cm
    const float vError = 1.0e-5; //cm
    std::array<unsigned int, 3> const triggers = {1,2,3}; //PlaceHolder Triggers
}//EventCuts
namespace TrackCuts{
    const float gPt = 0.5; // Global pT > 0.5
    const int nHitsFit = 20;  // nHitsFit > 15
    const int nHitsDedx = 5; // nHitsDedx > 5
    const float nHitsFit2Dedx = 0.51; // ((nHitsFit)/(nHitsDedx)) >= 0.51
    const float DCA = 0.005; // DCA > 0.005 cm
    const int Eta = 1; // |η| < 1
}//TrackCuts
namespace D0_Cuts{
    int nSigmaPion = 2;
    int nSigmaKaon = 2;
    float cos_theta = 0.95; //cos(the pointing angle θ) >0.95
    float DecayLength = 0.025; //Decay length >0.025cm.
    float DCA_12 = 0.02; // DCA Between two particles < 0.02
    float DCA_pi = 0.008; // DCA Pion > 0.008
    float DCA_k = 0.008; //  DCA Kaon > 0.008
    float oneOverBetaPion = 0.025; // (1/(β-1)) < 0.025
    float oneOverBetaKaon = 0.025; // (1/(β-1)) < 0.025
}//D0_Cuts
namespace JPSI_Cuts{
    float nSigmaElectron_min = -0.75; //nSigmaElectron > -0.75 if mom>0.8 GeV
    float nSigmaElectron_max = 2;  //nSigmaElectron < 2 for all mom
    float nSigmaElectron_lowhmom = -3.15; // nSigmaElectron > 3*mom.Mag()-3.15 if mom<0.8 GeV
    float oneOverBetaElectron = 0.025; // (1/(β-1)) < 0.025
}//JPSI_Cuts
#endif // !StMyCuts_hh
