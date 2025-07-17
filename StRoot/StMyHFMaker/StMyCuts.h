#ifndef StMyCuts_hh
#define StMyCuts_hh
#include "Rtypes.h"

namespace EventCuts{
    int vZ_min = -30; //cm
    int vZ_max = 30; //cm
    int vZVpdvZ = 2; //|TPC_vZ - VPD_Vz|(cm)
    int vR = 2; //cm
    float vError = 1.0e-5; //cm
    std::array<unsigned int, 3> const triggers = {1,2,3}; //PlaceHolder
}//EventCuts
namespace D0_Cuts{

}//D0_Cuts
namespace JPSI_Cuts{

}//JPSI_Cuts
#endif // !StMyCuts_hh
