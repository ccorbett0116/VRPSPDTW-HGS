#include <iostream>
#include <vector>
#include <string>
#include "InstanceCVRPLIB.h"
#include "Params.h"
#include "Individual.h"
#include "Split.h"
#include "LocalSearch.h"

static void printSolution(const Params& P, const Individual& I, const char* tag)
{
    std::cout << "---- " << tag << " ----\n";
    std::cout << "Routes used: " << I.eval.nbRoutes
              << " | Penalized cost: " << I.eval.penalizedCost << "\n";
    for (int r = 0; r < P.nbVehicles; ++r)
    {
        if (!I.chromR[r].empty())
        {
            std::cout << "Route " << (r+1) << " :";
            for (int v : I.chromR[r]) std::cout << " " << v;
            std::cout << "\n";
        }
    }
    std::cout << std::flush;
}

int main(int argc, char** argv)
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <instancePath> <nbVeh> <penCap> <penWarp>\n";
        return 1;
    }
    const std::string path = argv[1];
    const int max_nbVeh = std::stoi(argv[2]);
    const double penCap  = std::stod(argv[3]);
    const double penWarp = std::stod(argv[4]);

    // Load instance (VRPSPDTW format)
    InstanceCVRPLIB inst(path, /*isRoundingInteger=*/false);
    std::cout << "----- FLEET SIZE SPECIFIED: SET TO " << max_nbVeh << " VEHICLES\n";
    std::cout << "----- INSTANCE SUCCESSFULLY LOADED WITH "
              << inst.nbClients << " CLIENTS AND " << inst.numVehicles << " VEHICLES\n";

    // Build Params (durationLimit from depot due time)
    const double durationLimit = inst.due_time[0];
    AlgorithmParameters ap; // defaults are fine
    Params P(
        inst.x_coords, inst.y_coords,
        inst.dist_mtx, inst.time_mtx,
        inst.service_time, inst.demands, inst.pickups,
        inst.vehicleCapacity, durationLimit,
        inst.ready_time, inst.due_time,
        max_nbVeh,
        /*isDurationConstraint=*/true,
        /*verbose=*/true,
        ap
    );
    // Penalties & costs
    P.penaltyCapacity   = penCap;
    P.penaltyMultiplier = penWarp;          // soft time-window "warp" penalty
    // Use units from the instance (already read by InstanceCVRPLIB)
    P.unitCost          = inst.unitCost;
    P.dispatchingCost   = inst.dispatchingCost;

    // Build initial giant tour 1..n
    Individual indiv(P);
    indiv.chromT.clear();
    for (int i = 1; i <= P.nbClients; ++i) indiv.chromT.push_back(i);

    // Split to get an initial routed solution
    Split split(P);
    split.generalSplit(indiv, max_nbVeh);
    indiv.evaluateCompleteCost(P);
    printSolution(P, indiv, "After Split (baseline)");

    // Run Local Search (uses load -> run -> export)
    LocalSearch ls(P);                                // ctor
    ls.run(indiv, P.penaltyCapacity, P.penaltyDuration);  // calls loadInside + search  :contentReference[oaicite:0]{index=0}
    ls.exportIndividual(indiv);                       // push improved routes back

    // Re-evaluate with *exact* VRPSPDTW cost (distance + warp + cap + dispatch)
    indiv.evaluateCompleteCost(P);
    printSolution(P, indiv, "After Local Search");

    return 0;
}
