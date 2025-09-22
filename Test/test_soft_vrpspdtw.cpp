#include <iostream>
#include <vector>
#include <numeric>
#include "../Program/InstanceCVRPLIB.h"
#include "../Program/Params.h"
#include "../Program/Split.h"
#include "../Program/Individual.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " path/to/instance.txt [veh] [cap_pen] [warp_pen]\n";
        return 1;
    }
    const std::string path = argv[1];

    // Load instance (your reader now also auto-sets isDurationConstraint)
    InstanceCVRPLIB inst(path, /*isRoundingInteger*/false);

    // Optional overrides from CLI
    int nbVeh = (inst.numVehicles > 0 ? inst.numVehicles : 999999);
    if (argc >= 3) nbVeh = std::stoi(argv[2]);
    double capPen  = (argc >= 4 ? std::stod(argv[3]) : 1000.0);
    double warpPen = (argc >= 5 ? std::stod(argv[4]) : 1000.0);

    // Algo params (set anything you need)
    AlgorithmParameters ap;
    ap.useSwapStar = 0; // coords not needed for this split path

    // Construct Params IN ONE CALL (no default ctor, no post-assignments)
    // NOTE: this matches the non-default ctor you added (including pickups).
    Params P(
        /* x_coords   */ std::vector<double>{},         // not used by this split
        /* y_coords   */ std::vector<double>{},         // not used by this split
        /* dist_mtx   */ inst.dist_mtx,                 // objective distance matrix
        /* time_mtx   */ inst.time_mtx,                 // travel time matrix (for TW/warp)
        /* service    */ inst.service_time,
        /* demands    */ inst.demands,
        /* pickups    */ inst.pickups,                  // <-- REQUIRED for SPDTW
        /* capacity   */ inst.vehicleCapacity,
        /* durLimit   */ inst.durationLimit,
        /* ready      */ inst.ready_time,
        /* due        */ inst.due_time,
        /* dispatch   */ inst.dispatchingCost,
        /* unit cost  */ inst.unitCost,
        /* nbVeh      */ nbVeh,
        /* TW flag    */ inst.isDurationConstraint,
        /* verbose    */ true,
        /* ap         */ ap
    );

    // Tune soft penalties
    P.penaltyCapacity   = capPen;
    P.penaltyMultiplier = warpPen;  // time-warp weight

    // Make a trivial chromosome 1..n
    Individual indiv(P);
    indiv.chromT.resize(P.nbClients);
    std::iota(indiv.chromT.begin(), indiv.chromT.end(), 1);

    // Run split
    Split split(P);
    split.generalSplit(indiv, P.nbVehicles);

    std::cout << "Routes used: " << indiv.eval.nbRoutes
              << " | Penalized cost: " << indiv.eval.penalizedCost << "\n";

    // Print routes
    for (int r = 0; r < P.nbVehicles; ++r) {
        if (indiv.chromR[r].empty()) continue;
        std::cout << "Route " << (r+1) << " :";
        for (int c : indiv.chromR[r]) std::cout << " " << c;
        std::cout << "\n";
    }
    return 0;
}
