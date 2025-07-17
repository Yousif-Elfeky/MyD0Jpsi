#ifndef StMyCuts_hh
#define StMyCuts_hh
#include "Rtypes.h"

namespace EventCuts{
    const int vZ_min = -30; //cm
    const int vZ_max = 30; //cm
    const int vZVpdvZ = 2; //|TPC_vZ - VPD_Vz|(cm)
    const int vR = 2; //cm
    const float vError = 1.0e-5; //cm
    std::array<unsigned int, 3> const triggers = {1,2,3}; //PlaceHolder
}//EventCuts
namespace TrackCuts{
    const int gPt = 1; // Global pT > 1
    const int nHitsFit = 1;  // nHitsFit > 1
    const int nHitsDedx = 1; // nHitsDedx > 1
    const float nHitsFit2Dedx = 0.51; // ((nHitsFit)/(nHitsDedx)) >= 0.51
    // const int Dca = 1; // DCA < 1
    const int Eta = 1; // |η| < 1
}//TrackCuts
namespace D0_Cuts{

}//D0_Cuts
namespace JPSI_Cuts{

}//JPSI_Cuts
#endif // !StMyCuts_hh
