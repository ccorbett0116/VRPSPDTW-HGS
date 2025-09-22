#include "Split.h"
#include <stdexcept>
#include <cmath>

void Split::generalSplit(Individual & indiv, int nbMaxVehicles)
{
    // Never fewer than trivial bin-packing lower bound
    maxVehicles = std::max<int>(nbMaxVehicles, std::ceil(params.totalDemand / params.vehicleCapacity));

    // (Re)build the DP arrays to the exact shape we will use  <-- NEW & IMPORTANT
    potential.assign(maxVehicles + 1, std::vector<double>(params.nbClients + 1, 1.e30));
    pred.assign     (maxVehicles + 1, std::vector<int   >(params.nbClients + 1, -1));

    // Reinit the per-client scaffolding (as you already do)
    for (int i = 1; i <= params.nbClients; i++)
    {
        const int cust = indiv.chromT[i - 1];

        // objective distance legs
        cliSplit[i].d0_x  = params.timeCost[0][cust];
        cliSplit[i].dx_0  = params.timeCost[cust][0];
        if (i < params.nbClients)
            cliSplit[i].dnext = params.timeCost[cust][indiv.chromT[i]];
        else
            cliSplit[i].dnext = -1.e30;

        // travel-time legs
        cliSplit[i].t0_x  = params.timeMatrix[0][cust];
        cliSplit[i].tx_0  = params.timeMatrix[cust][0];
        if (i < params.nbClients)
            cliSplit[i].tnext = params.timeMatrix[cust][indiv.chromT[i]];
        else
            cliSplit[i].tnext = -1.e30;

        // demands / pickups / service / TW
        cliSplit[i].demand      = params.cli[cust].demand;
        cliSplit[i].pickup      = params.pickups[cust];
        cliSplit[i].serviceTime = params.serviceTime[cust];
        cliSplit[i].twOpen      = params.readyTime[cust];
        cliSplit[i].twClose     = params.dueTime[cust];

        // prefix sums for O(n) CVRP path (still used in no-TW branch)
        sumLoad[i]     = sumLoad[i - 1] + cliSplit[i].demand;
        sumService[i]  = sumService[i - 1] + cliSplit[i].serviceTime;
        sumDistance[i] = sumDistance[i - 1] + cliSplit[i - 1].dnext;
    }

    // Try unlimited-fleet split; if it fails, go to limited-fleet
    if (splitSimple(indiv) == 0)
        splitLF(indiv);

    // Finish individual costs
    indiv.evaluateCompleteCost(params);
}

// ---- Soft VRPSPDTW edge cost (Bellman, patched) ----
inline double Split::edgeCostSoftVRPSPDTW(int i, int j) const
{
    double distObj = 0.0;   // objective distance (timeCost)
    double timeNow = 0.0;   // schedule time (timeMatrix)
    double warp    = 0.0;   // accumulated time-warp

    // SPDTW capacity trackers for segment (i+1..j)
    double del  = 0.0;      // sum deliveries in the segment
    double pick = 0.0;      // sum pickups in the segment
    double maxDelta = 0.0;  // max over prefixes of (pick - del)

    for (int u = i + 1; u <= j; ++u)
    {
        if (u == i + 1)
        {
            distObj += cliSplit[u].d0_x;
            timeNow  = std::max(cliSplit[u].twOpen, cliSplit[u].t0_x);
        }
        else
        {
            distObj += cliSplit[u - 1].dnext;
            timeNow  = std::max(
                timeNow + cliSplit[u - 1].serviceTime + cliSplit[u - 1].tnext,
                cliSplit[u].twOpen
            );
        }

        // Soft TW: allow late arrival but penalize
        const double latestAtU = std::min(
            cliSplit[u].twClose,
            params.durationLimit - cliSplit[u].serviceTime - cliSplit[u].tx_0
        );
        if (timeNow > latestAtU)
        {
            warp   += (timeNow - latestAtU);
            timeNow = latestAtU;  // clamp so subsequent nodes aren’t even later
        }

        // SPDTW prefixes
        del  += cliSplit[u].demand;
        pick += cliSplit[u].pickup;
        maxDelta = std::max(maxDelta, pick - del);
    }

    // Return to depot (objective distance)
    distObj += cliSplit[j].dx_0;

    // Compute final schedule duration (depot -> first .. last -> depot)
    double duration = timeNow + cliSplit[j].serviceTime + cliSplit[j].tx_0;

    // Add depot due-time warp (if we return after depot’s due time)
    if (duration > params.dueTime[0])
        warp += (duration - params.dueTime[0]);

    // Duration excess vs global limit
    const double durationExcess = std::max(0.0, duration - params.durationLimit);

    // True onboard max for SPDTW on this segment
    const double onboardMax = del + std::max(0.0, maxDelta);
    const double capOver    = std::max(0.0, onboardMax - params.vehicleCapacity);

    // Always return a cost (never 1e30). Infeasibility is encoded via penalties.
    return params.unitCost          * distObj
         + params.dispatchingCost
         + params.penaltyCapacity   * capOver
         + params.penaltyDuration   * durationExcess
         + params.penaltyMultiplier * warp;
}


int Split::splitSimple(Individual& indiv)
{
    const int n = params.nbClients;     // number of clients (excludes depot)
    const int V = params.nbVehicles;    // max number of vehicles

    // DP tables: potential[k][j] = min cost to cover clients [1..j] with k routes
    std::vector<std::vector<double>> potential(V + 1, std::vector<double>(n + 1, 1e30));
    std::vector<std::vector<int>> pred(V + 1, std::vector<int>(n + 1, -1));

    // Base case: 0 routes, only depot (j=0)
    potential[0][0] = 0.0;

    // Fill DP
    for (int k = 1; k <= V; k++) {
        for (int i = 0; i < n; i++) {
            if (potential[k-1][i] >= 1e29) continue; // unreachable

            for (int j = i + 1; j <= n; j++) {
                double cost = edgeCostSoftVRPSPDTW(i, j);
                if (cost >= 1e29) continue; // infeasible segment

                double newCost = potential[k-1][i] + cost;
                if (newCost < potential[k][j]) {
                    potential[k][j] = newCost;
                    pred[k][j] = i;
                }
            }
        }
    }

    // Find best feasible solution with ≤ V routes
    int bestK = -1;
    double bestCost = 1e30;
    for (int k = 1; k <= V; k++) {
        if (potential[k][n] < bestCost) {
            bestCost = potential[k][n];
            bestK = k;
        }
    }

    if (bestK == -1) {
        throw std::runtime_error("SplitSimple failed: no feasible solution found");
    }

    // Clear existing routes
    for (int k = 0; k < params.nbVehicles; k++) indiv.chromR[k].clear();

    // Backtrack
    int end = n;
    int nbRoutes = 0;
    int k = bestK;

    while (end > 0 && k > 0) {
        int begin = pred[k][end];
        if (begin < 0) {
            throw std::runtime_error("SplitSimple backtrack failed: no predecessor for node " + std::to_string(end));
        }

        auto &route = indiv.chromR[nbRoutes];
        route.clear();
        for (int ii = begin; ii < end; ii++)
            route.push_back(indiv.chromT[ii]);

        end = begin;
        k--;
        nbRoutes++;
    }

    if (end != 0) {
        throw std::runtime_error("SplitSimple backtrack failed: did not reach depot");
    }

    // Success
    return 1;
}



int Split::splitLF(Individual& indiv)
{
    // Initialize all layers
    for (int k = 0; k <= maxVehicles; k++)
    {
        potential[k][0] = (k == 0 ? 0.0 : 1.e30);
        for (int i = 1; i <= params.nbClients; i++) potential[k][i] = 1.e30;
    }

    if (params.isDurationConstraint)
    {
        // k-layer Bellman using soft VRPSPDTW arc costs
        for (int k = 0; k < maxVehicles; k++)
        {
            for (int i = k; i < params.nbClients && potential[k][i] < 1.e29; i++)
            {
                for (int j = i + 1; j <= params.nbClients; j++)
                {
                    const double cost = edgeCostSoftVRPSPDTW(i, j);
                    if (potential[k][i] + cost < potential[k + 1][j])
                    {
                        potential[k + 1][j] = potential[k][i] + cost;
                        pred[k + 1][j]      = i;
                    }
                }
            }
        }
    }
    else
    {
        // Original O(n) CVRP split with deque (distances scaled; fixed cost via propagate())
        Trivial_Deque queue(params.nbClients + 1, 0);
        for (int k = 0; k < maxVehicles; k++)
        {
            queue.reset(k);
            for (int i = k + 1; i <= params.nbClients && queue.size() > 0; i++)
            {
                potential[k + 1][i] = propagate(queue.get_front(), i, k);
                pred[k + 1][i]      = queue.get_front();

                if (i < params.nbClients)
                {
                    if (!dominates(queue.get_back(), i, k))
                    {
                        while (queue.size() > 0 && dominatesRight(queue.get_back(), i, k))
                            queue.pop_back();
                        queue.push_back(i);
                    }
                    while (queue.size() > 1 &&
                           propagate(queue.get_front(), i + 1, k) >
                           propagate(queue.get_next_front(), i + 1, k) - MY_EPSILON)
                        queue.pop_front();
                }
            }
        }
    }

    if (potential[maxVehicles][params.nbClients] > 1.e29)
        throw std::string("ERROR : no Split solution has been propagated until the last node");

    // Possibly fewer vehicles is cheaper
    double best = potential[maxVehicles][params.nbClients];
    int nbRoutes = maxVehicles;
    for (int k = 1; k < maxVehicles; k++)
        if (potential[k][params.nbClients] < best)
            { best = potential[k][params.nbClients]; nbRoutes = k; }

    // Build routes
    for (int k = params.nbVehicles - 1; k >= nbRoutes; k--) indiv.chromR[k].clear();

    int end = params.nbClients;
    for (int k = nbRoutes - 1; k >= 0; k--)
    {
        indiv.chromR[k].clear();
        const int begin = pred[k + 1][end];
        for (int ii = begin; ii < end; ii++)
            indiv.chromR[k].push_back(indiv.chromT[ii]);
        end = begin;
    }
    return (end == 0);
}

Split::Split(const Params& params) : params(params)
{
    // Allocate linear split structures
    cliSplit    = std::vector<ClientSplit>(params.nbClients + 1);
    sumDistance = std::vector<double>(params.nbClients + 1, 0.0);
    sumLoad     = std::vector<double>(params.nbClients + 1, 0.0);
    sumService  = std::vector<double>(params.nbClients + 1, 0.0);
    potential   = std::vector<std::vector<double>>(params.nbVehicles + 1,
                    std::vector<double>(params.nbClients + 1, 1.e30));
    pred        = std::vector<std::vector<int>>(params.nbVehicles + 1,
                    std::vector<int>(params.nbClients + 1, 0));
}
