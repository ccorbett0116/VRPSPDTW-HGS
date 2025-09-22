/*MIT License
   (same header/license block as your project)
*/
#ifndef SPLIT_H
#define SPLIT_H

#include "Params.h"
#include "Individual.h"
#include <algorithm>
#include <vector>
#include <cmath>

// Per-client precomputed values used by Split
struct ClientSplit
{
    // Objective-distance legs
    double d0_x;     // objective dist: depot -> i
    double dx_0;     // objective dist: i -> depot
    double dnext;    // objective dist: i -> i+1  (i < n), else sentinel

    // Timing legs (for TW scheduling)
    double t0_x;     // travel time: depot -> i
    double tx_0;     // travel time: i -> depot
    double tnext;    // travel time: i -> i+1     (i < n), else sentinel

    // Demands / service / windows
    double demand;       // delivery (>= 0)
    double pickup;       // pickup   (>= 0), defaults to 0 if unavailable
    double serviceTime;  // service duration
    double twOpen;       // earliest start
    double twClose;      // latest start

    ClientSplit()
    : d0_x(0.), dx_0(0.), dnext(0.),
      t0_x(0.), tx_0(0.), tnext(0.),
      demand(0.), pickup(0.), serviceTime(0.),
      twOpen(0.), twClose(1.e30) {}
};

// Simple deque used by O(n) Split for pure CVRP
struct Trivial_Deque
{
    std::vector<int> myDeque;
    int indexFront;
    int indexBack;
    inline void pop_front(){ indexFront++; }
    inline void pop_back(){ indexBack--; }
    inline void push_back(int i){ indexBack++; myDeque[indexBack] = i; }
    inline int  get_front(){ return myDeque[indexFront]; }
    inline int  get_next_front(){ return myDeque[indexFront + 1]; }
    inline int  get_back(){ return myDeque[indexBack]; }
    inline int  size(){ return indexBack - indexFront + 1; }
    void reset(int firstNode) { myDeque[0] = firstNode; indexBack = 0; indexFront = 0; }

    Trivial_Deque(int nbElements, int firstNode)
    {
        myDeque.resize(nbElements);
        myDeque[0] = firstNode;
        indexBack = 0;
        indexFront = 0;
    }
};

class Split
{
private:
    // Problem params
    const Params& params;
    int maxVehicles;

    // Linear split scaffolding
    std::vector<ClientSplit>            cliSplit;
    std::vector<std::vector<double>>    potential; // potential[k][i]
    std::vector<std::vector<int>>       pred;      // pred[k][i]
    std::vector<double>                 sumDistance; // prefix of d_{k,k+1}
    std::vector<double>                 sumLoad;     // prefix of deliveries
    std::vector<double>                 sumService;  // prefix of service

    // ---- Helpers for O(n) CVRP path (now scaled by unitCost) ----
    inline double propagate(int i, int j, int k)
    {
        // Distance part must be scaled by unitCost; add dispatchingCost once per route
        return potential[k][i]
            + params.unitCost * ((sumDistance[j] - sumDistance[i + 1]) + cliSplit[i + 1].d0_x + cliSplit[j].dx_0)
            + params.penaltyCapacity * std::max<double>(sumLoad[j] - sumLoad[i] - params.vehicleCapacity, 0.)
            + params.dispatchingCost;
    }
    inline bool dominates(int i, int j, int k)
    {
        // Apply same unitCost scaling to the distance terms used in dominance tests
        return potential[k][j] + params.unitCost * cliSplit[j + 1].d0_x >
               potential[k][i] + params.unitCost * (cliSplit[i + 1].d0_x + (sumDistance[j + 1] - sumDistance[i + 1]))
             + params.penaltyCapacity * (sumLoad[j] - sumLoad[i]);
    }
    inline bool dominatesRight(int i, int j, int k)
    {
        return potential[k][j] + params.unitCost * cliSplit[j + 1].d0_x <
               potential[k][i] + params.unitCost * (cliSplit[i + 1].d0_x + (sumDistance[j + 1] - sumDistance[i + 1])) + MY_EPSILON;
    }

    // ---- Soft VRPSPDTW edge cost (Bellman) ----
    // Cost of serving segment (i+1..j):
    //   unitCost * distance + soft capacity (onboard-max) + soft time windows (time-warp) + dispatchingCost
    inline double edgeCostSoftVRPSPDTW(int i, int j) const;

    // Split variants
    int splitSimple(Individual& indiv); // unlimited fleet
    int splitLF(Individual& indiv);     // limited fleet

public:
    void generalSplit(Individual& indiv, int nbMaxVehicles);
    Split(const Params& params);
};

#endif
