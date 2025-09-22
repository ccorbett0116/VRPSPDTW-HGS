#include "LocalSearch.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <cassert>

#ifndef MY_EPSILON
#define MY_EPSILON 1e-9
#endif

#ifndef LS_DEBUG
#define LS_DEBUG 0
#endif

#if LS_DEBUG
#include <iostream>
static inline int cid(const Node* n){ return n? n->cour : -1; }
static inline int rid(const Route* r){ return r? r->cour : -1; }
#endif
#ifndef NDEBUG
static bool validateRouteForward(const Route* R, int nbClients) {
    if (!R || !R->depot) return false;
    const Node* n = R->depot;
    int depots = 0, steps = 0;
    do {
        if (!n) return false;
        if (n->isDepot) depots++;
        n = n->next;
        if (++steps > nbClients + 2) return false; // loop / overflow
    } while (n && depots < 2);
    return depots == 2; // exactly one full lap depot→...→depot
}

static bool validateRouteBackward(const Route* R, int nbClients) {
    if (!R || !R->depot) return false;
    const Node* n = R->depot;
    int depots = 0, steps = 0;
    do {
        if (!n) return false;
        if (n->isDepot) depots++;
        n = n->prev;
        if (++steps > nbClients + 2) return false;
    } while (n && depots < 2);
    return depots == 2;
}

static inline void assertRouteOK(const Route* R, int nbClients, const char* where) {
    if (!validateRouteForward(R, nbClients) || !validateRouteBackward(R, nbClients)) {
        std::cerr << "[LS][FATAL] Route ring corrupted at " << where
                  << " (route #" << (R ? R->cour : -1) << ")\n";
        std::abort();
    }
}
#endif

#include <iostream>
#include <sstream>
#include <iomanip>

static long long LS_move_counter = 0;

static std::string routeToStr(const Route* R) {
    if (!R || !R->depot) return "-";
    std::ostringstream oss;
    Node* n = R->depot->next;
    bool first = true;
    while (n && !n->isDepot) {
        if (!first) oss << ' ';
        oss << n->cour;
        first = false;
        n = n->next;
    }
    return oss.str();
}

static inline void LS_logMove(const char* name,
                              const Route* RU, const Route* RV,
                              const Node* U, const Node* X,
                              const Node* V, const Node* Y,
                              double oldSum, double newSum)
{
    if (LS_move_counter%1000==0) {
        std::cerr << "[LS][" << std::setw(6) << LS_move_counter << "] "
                  << name
                  << "  Δ=" << (newSum - oldSum)
                  << "  old=" << oldSum
                  << " new=" << newSum
                  << "\n    U=" << (U ? U->cour : -1)
                  << " X=" << (X ? X->cour : -1)
                  << " V=" << (V ? V->cour : -1)
                  << " Y=" << (Y ? Y->cour : -1)
                  << "\n    RouteU(" << (RU ? RU->cour : -1) << "): [" << routeToStr(RU) << "]"
                  << "\n    RouteV(" << (RV ? RV->cour : -1) << "): [" << routeToStr(RV) << "]"
                  << "\n";
    }
    LS_move_counter++;
}


// =================== Segment primitives (TW + SPDTW) ===================
// Each TWSPD_Seg summarizes a contiguous customer block in forward order.
// Also maintain per-route O(m^2) caches of forward blocks [i..j] and
// their reversed counterpart reverse([i..j]) to get true O(1) move evals
// even for reversal moves (2-opt intra).

// Exact step-by-step simulation that matches Individual evaluation perfectly
static inline TWSPD_Seg simulateExactSegment(const std::vector<int>& customers, const Params& P)
{
    if (customers.empty()) return TWSPD_Seg{};

    TWSPD_Seg S{};
    S.len = customers.size();
    S.first = customers[0];
    S.last = customers.back();

    // Simulate exactly like Individual::evaluateCompleteCost
    // Start from arrival at first customer (after depot travel)
    double currentTime = P.timeMatrix[0][S.first]; // Include depot travel time
    double totalDist = 0.0;
    double totalWarp = 0.0;
    double totalServ = 0.0;

    // SPDTW tracking
    double totalDel = 0.0, prefDel = 0.0, prefPick = 0.0, maxDelta = 0.0, minDelta = 0.0;
    for (int u : customers) {
        totalDel += P.cli[u].demand;
    }

    for (int i = 0; i < (int)customers.size(); i++) {
        int node = customers[i];

        // Travel time (internal to segment only)
        if (i > 0) {
            int prev = customers[i-1];
            totalDist += P.timeCost[prev][node];
            currentTime += P.timeMatrix[prev][node];
        }

        // Time windows (exactly like Individual evaluation)
        const double earliest = P.readyTime[node];
        const double latest = std::min(
            P.dueTime[node],
            P.durationLimit - P.serviceTime[node] - P.timeMatrix[node][0]
        );

        if (currentTime < earliest) {
            currentTime = earliest;
        } else if (currentTime > latest) {
            totalWarp += (currentTime - latest);
            currentTime = latest;
        }

        // Service time
        totalServ += P.serviceTime[node];
        currentTime += P.serviceTime[node];

        // SPDTW prefixes
        prefDel += P.cli[node].demand;
        prefPick += P.pickups[node];
        double delta = prefPick - prefDel;
        maxDelta = std::max(maxDelta, delta);
        minDelta = std::min(minDelta, delta);
    }

    // Add final depot travel time
    currentTime += P.timeMatrix[S.last][0];

    // Check final depot time window
    if (currentTime > P.dueTime[0]) {
        totalWarp += (currentTime - P.dueTime[0]);
    }

    // Fill segment fields to match route_cost_from_seg expectations
    S.dist = totalDist;
    S.Serv = totalServ;
    // dur should be the internal time (without depot legs for route_cost_from_seg)
    S.dur = currentTime - P.timeMatrix[0][S.first] - P.timeMatrix[S.last][0];
    S.TW = totalWarp;

    // Time window bounds (for composition if needed)
    S.E = P.readyTime[S.first];
    S.L = std::min(P.dueTime[S.last], P.durationLimit - P.serviceTime[S.last] - P.timeMatrix[S.last][0]);

    // SPDTW fields
    S.Del = totalDel;
    S.Pick = 0.0; for (int u : customers) S.Pick += P.pickups[u];
    S.sumDelta = S.Pick - S.Del;
    S.maxPref = maxDelta;
    S.minPref = minDelta;
    S.totalDel = totalDel;

    return S;
}


static inline double route_cost_from_seg(const Params& P, const TWSPD_Seg& S, double penaltyCapacity, double penaltyDuration, double penaltyTimeWarp)
{
    if (S.len == 0) return 0.0;  // empty route: cost 0, no dispatching

    constexpr int DEPOT = 0;

    // Distance: depot → first, internal, last → depot
    double dist = P.timeCost[DEPOT][S.first] + S.dist + P.timeCost[S.last][DEPOT];
    double currentTime = P.timeMatrix[DEPOT][S.first] + S.dur + P.timeMatrix[S.last][DEPOT];

    // Time warp
    double warp = S.TW;
    if (currentTime > P.dueTime[DEPOT])
        warp += (currentTime - P.dueTime[DEPOT]);

    // Duration excess
    double durExc = 0.0;
    if (currentTime > P.durationLimit)
        durExc = (currentTime - P.durationLimit);

    // Capacity excess (onboard max rule = Del + max(0, Pick-Del prefixes))
    double onboardMax = S.Del + std::max(0.0, S.maxPref);
    double capExc = 0.0;
    if (onboardMax > P.vehicleCapacity)
        capExc = (onboardMax - P.vehicleCapacity);

    // Penalized cost (dispatching cost once per non-empty route)
    return  P.unitCost        * dist
          + P.dispatchingCost * 1.0   // one dispatch per non-empty route
          + penaltyCapacity * capExc
          + penaltyDuration * durExc
          + penaltyTimeWarp * warp;
}



// =================== LocalSearch core ===================

LocalSearch::LocalSearch(Params& P)
    : params(P),
      searchCompleted(false),
      nbMoves(0),
      loopID(0),
      nodeU(nullptr), nodeV(nullptr), nodeX(nullptr), nodeY(nullptr),
      routeU(nullptr), routeV(nullptr),
      nodeUPrevIndex(-1), nodeUIndex(-1), nodeXIndex(-1), nodeXNextIndex(-1),
      nodeVPrevIndex(-1), nodeVIndex(-1), nodeYIndex(-1), nodeYNextIndex(-1),
      loadU(0.0), loadV(0.0), loadX(0.0), loadY(0.0),
      serviceU(0.0), serviceV(0.0), serviceX(0.0), serviceY(0.0),
      penaltyCapacityLS(P.penaltyCapacity),
      penaltyDurationLS(P.penaltyDuration),
      penaltyTimeWarpLS(P.penaltyMultiplier),
      intraRouteMove(false)
{
    // 1) Routes & depot nodes
    routes.resize(P.nbVehicles);
    for (int r = 0; r < P.nbVehicles; ++r) {
        Route& R = routes[r];
        R.cour = r;
        R.nbCustomers = 0;

        if (!R.depot) R.depot = new Node();
        R.depot->cour     = 0;
        R.depot->isDepot  = true;
        R.depot->route    = &R;
        R.depot->next     = R.depot;
        R.depot->prev     = R.depot;
        R.whenLastModified = 0;
        R.whenLastTestedSWAPStar = -1;
    }

    // 2) Clients array: one Node per customer id in [1..nbClients]
    clients.clear();
    clients.resize(P.nbClients + 1); // index 0 unused (depot handled per-route)
    for (int i = 1; i <= P.nbClients; ++i) {
        Node& n = clients[i];
        n.cour       = i;
        n.isDepot    = false;
        n.route      = nullptr;
        n.prev       = nullptr;
        n.next       = nullptr;
        n.position   = 0;
        n.whenLastTestedRI = -1;
    }

    // 3) Correlated vertex order arrays if you use them elsewhere
    orderNodes.resize(P.nbClients);
    for (int i = 0; i < P.nbClients; ++i) orderNodes[i] = i + 1;
    orderRoutes.resize(P.nbVehicles);
    for (int r = 0; r < P.nbVehicles; ++r) orderRoutes[r] = r;

    // 4) Segment caches for O(1) eval
    prefSeg.resize(P.nbVehicles);
    sufSeg.resize(P.nbVehicles);
    revSufSeg.resize(P.nbVehicles);
    rangeFwd.resize(P.nbVehicles);
    rangeRev.resize(P.nbVehicles);

    emptyRoutes.clear();
}




void LocalSearch::routeToSequence(Route* R, std::vector<int>& seq) const
{
    seq.clear();
    if (!R || !R->depot) return;

    Node* n = R->depot->next;
    while (n && !n->isDepot) {
        seq.push_back(n->cour);
        n = n->next;
    }
}

void LocalSearch::updateNodePositions(Route* R)
{
    int pos = 0;
    for (Node* n = R->depot->next; !n->isDepot; n = n->next)
    {
        n->position = ++pos;
        n->route = R;
    }
    R->nbCustomers = pos;
}

static inline size_t idx2D(int m, int i, int j) {
    assert(i >= 0 && i < m);
    assert(j >= 0 && j < m);
    return (size_t)i * (size_t)m + (size_t)j;
}

void LocalSearch::updateRouteData(Route* R)
{
    if (!R) return;

#if LS_DEBUG
    if (!validateRouteLinkedList(R)) {
        std::cerr << "ERROR: Route " << R->cour << " linked list corrupted!" << std::endl;
        printRouteLinkedList(R);
        std::abort();
    }
#endif

    // Comprehensive route integrity check
    if (!R->depot) {
        std::cerr << "ERROR: Route " << R->cour << " has null depot!" << std::endl;
        return;
    }

    updateNodePositions(R);

    const int r = R->cour;
    const int m = R->nbCustomers;

    if ((int)prefSeg.size()    <= r) prefSeg.resize(r + 1);
    if ((int)sufSeg.size()     <= r) sufSeg.resize(r + 1);
    if ((int)revSufSeg.size()  <= r) revSufSeg.resize(r + 1);
    if ((int)rangeFwd.size()   <= r) rangeFwd.resize(r + 1);
    if ((int)rangeRev.size()   <= r) rangeRev.resize(r + 1);

    auto& PREF = prefSeg[r];
    auto& SUF  = sufSeg[r];
    auto& RSUF = revSufSeg[r];

    // Only resize if needed to avoid memory reallocation
    int required_segments = std::max(1, m);
    if ((int)PREF.size() != required_segments) PREF.resize(required_segments);
    if ((int)SUF.size() != required_segments) SUF.resize(required_segments);
    if ((int)RSUF.size() != required_segments) RSUF.resize(required_segments);

    std::vector<int> seq; seq.reserve(m);
    routeToSequence(R, seq);

    // Safety check: if sequence is corrupted, don't proceed
    if (seq.size() != (size_t)m) {
        std::cerr << "ERROR: updateRouteData sequence size mismatch. Expected " << m
                  << " customers, got " << seq.size() << std::endl;
        return;
    }

    if (m == 0)
    {
        rangeFwd[r].clear();
        rangeRev[r].clear();
        R->routeCost = 0.0;
        R->whenLastModified = nbMoves;
        return;
    }

    // Check for size mismatch between m and actual sequence size
    if ((int)seq.size() != m) {
        std::cerr << "CRITICAL: Size mismatch! Route " << r << " reports m=" << m
                  << " customers but sequence has " << seq.size() << " elements" << std::endl;
        std::cerr << "Using actual sequence size for processing." << std::endl;
        // Use the actual sequence size to prevent bounds errors
        const_cast<int&>(m) = (int)seq.size();
    }

    for (size_t i = 0; i < seq.size(); i++) {
        if (seq[i] < 0 || seq[i] >= (int)params.cli.size()) {
            std::cerr << "CRITICAL: Sequence corrupted before processing! Position " << i
                      << " has invalid customer " << seq[i] << std::endl;

            // JIT DEBUG TRIGGER: Source sequence corruption detected!

            return; // Abort processing this route
        }
    }

    // Build ALL segments using exact simulation for perfect accuracy
    // This gives us O(1) lookups with zero approximation error

    // O(m^2) forward and reversed range blocks (build these first)
    // Only resize if needed to avoid repeated memory allocations
    size_t required_size = (size_t)m * (size_t)m;
    if (rangeFwd[r].size() != required_size) {
        rangeFwd[r].resize(required_size);
        rangeRev[r].resize(required_size);
    }

    // Forward segments [i..j] using exact simulation
    for (int i = 0; i < m; ++i)
    {
        // Check if seq got corrupted during processing
        if ((int)seq.size() != m) {
            std::cerr << "CRITICAL: seq vector size changed during processing! Was " << m
                      << " now " << seq.size() << ". Route=" << r << std::endl;
            return; // Abort processing
        }

        for (int j = i; j < m; ++j)
        {
            // Declare variables before any goto targets
            size_t index;
            bool skip_subseq = false;

            // Validate indices before creating subsequence
            if (i < 0 || j >= (int)seq.size() || i > j) {
                std::cerr << "ERROR: Invalid subsequence indices i=" << i << " j=" << j
                          << " for sequence size=" << seq.size() << std::endl;
                continue;
            }

            auto start_iter = seq.begin() + i;
            auto end_iter = seq.begin() + j + 1;

            // Use local vector to avoid memory corruption
            std::vector<int> subseq(start_iter, end_iter);

            // Validate subsequence before simulation
            for (size_t k = 0; k < subseq.size(); k++) {
                if (subseq[k] < 0 || subseq[k] >= (int)params.cli.size()) {
                    std::cerr << "ERROR: Corrupted subsequence at [" << i << ".." << j << "] position " << k
                              << " customer=" << subseq[k] << " (valid: 0-" << (params.cli.size()-1) << ")" << std::endl;
                    std::cerr << "Original sequence size: " << seq.size() << " subseq size: " << subseq.size() << std::endl;
                    std::cerr << "Full subsequence: [";
                    for (size_t x = 0; x < subseq.size() && x < 20; x++) {
                        std::cerr << subseq[x];
                        if (x < subseq.size() - 1 && x < 19) std::cerr << ", ";
                    }
                    if (subseq.size() > 20) std::cerr << ", ...";
                    std::cerr << "]" << std::endl;

                    // JIT DEBUG TRIGGER: Corrupted subsequence validation!

                    skip_subseq = true;
                    break; // Exit validation loop
                }
            }

            if (skip_subseq) continue; // Skip this corrupted subsequence

            // Validate 2D index before using it
            index = idx2D(m, i, j);
            if (index >= rangeFwd[r].size()) {
                std::cerr << "ERROR: 2D index out of bounds! idx2D(" << m << "," << i << "," << j
                          << ")=" << index << " >= rangeFwd[" << r << "].size()=" << rangeFwd[r].size() << std::endl;
                continue; // Skip this operation
            }

            rangeFwd[r][index] = simulateExactSegment(subseq, params);
        }
    }

    // Reversed segments reverse([i..j]) using exact simulation
    for (int i = 0; i < m; ++i)
    {
        for (int j = i; j < m; ++j)
        {
            // Declare variables before any goto targets
            size_t rev_index;
            bool skip_rev_subseq = false;

            // Validate indices before creating subsequence
            if (i < 0 || j >= (int)seq.size() || i > j) {
                std::cerr << "ERROR: Invalid reversed subsequence indices i=" << i << " j=" << j
                          << " for sequence size=" << seq.size() << std::endl;
                continue;
            }

            // Use local vector to avoid memory corruption
            std::vector<int> subseq(seq.begin() + i, seq.begin() + j + 1);
            std::reverse(subseq.begin(), subseq.end());

            // Validate subsequence before simulation
            for (size_t k = 0; k < subseq.size(); k++) {
                if (subseq[k] < 0 || subseq[k] >= (int)params.cli.size()) {
                    std::cerr << "ERROR: Corrupted reversed subsequence at [" << i << ".." << j << "] position " << k
                              << " customer=" << subseq[k] << " (valid: 0-" << (params.cli.size()-1) << ")" << std::endl;
                    std::cerr << "Original sequence size: " << seq.size() << " subseq size: " << subseq.size() << std::endl;
                    skip_rev_subseq = true;
                    break; // Exit validation loop
                }
            }

            if (skip_rev_subseq) continue; // Skip this corrupted subsequence

            // Validate 2D index before using it
            rev_index = idx2D(m, i, j);
            if (rev_index >= rangeRev[r].size()) {
                std::cerr << "ERROR: 2D index out of bounds for rangeRev! idx2D(" << m << "," << i << "," << j
                          << ")=" << rev_index << " >= rangeRev[" << r << "].size()=" << rangeRev[r].size() << std::endl;
                continue; // Skip this operation
            }

            rangeRev[r][rev_index] = simulateExactSegment(subseq, params);
        }
    }

    // Now build prefixes and suffixes from exact segments
    // Forward prefixes: PREF[i] = rangeFwd[0..i]
    for (int i = 0; i < m; ++i)
    {
        PREF[i] = rangeFwd[r][idx2D(m, 0, i)];
    }

    // Forward suffixes: SUF[i] = rangeFwd[i..m-1]
    for (int i = 0; i < m; ++i)
    {
        SUF[i] = rangeFwd[r][idx2D(m, i, m - 1)];
    }

    // Reverse suffixes: RSUF[i] = rangeRev[i..m-1]
    for (int i = 0; i < m; ++i)
    {
        RSUF[i] = rangeRev[r][idx2D(m, i, m - 1)];
    }

    R->routeCost = route_cost_from_seg(params, PREF[m-1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    R->whenLastModified = nbMoves;

    // Memory cleanup: shrink temp vector if it grew too large
    if (tempSubseq.capacity() > 100) {
        tempSubseq.clear();
        tempSubseq.shrink_to_fit();
    }
}

void LocalSearch::printRouteLinkedList(Route* R) const {
    if (!R || !R->depot) {
        std::cerr << "Invalid route or depot" << std::endl;
        return;
    }

    std::cerr << "Route " << R->cour << " (customers: " << R->nbCustomers << "): ";
    Node* n = R->depot;
    int count = 0;
    int maxCount = params.nbClients + 10;

    do {
        if (n->isDepot) {
            std::cerr << "DEPOT";
        } else {
            std::cerr << n->cour;
        }

        if (n->next && !n->next->isDepot) {
            std::cerr << " -> ";
        }

        n = n->next;
        count++;

        if (count > maxCount) {
            std::cerr << " ... (loop detected)";
            break;
        }
    } while (n != R->depot);

    std::cerr << std::endl;
}

bool LocalSearch::validateRouteLinkedList(Route* R) const {
    if (!R || !R->depot) {
        std::cerr << "ERROR: Invalid route or depot" << std::endl;
        return false;
    }

    // Check depot properties
    if (!R->depot->isDepot || R->depot->route != R) {
        std::cerr << "ERROR: Depot validation failed for route " << R->cour << std::endl;
        return false;
    }

    // ---- Forward traversal ----
    Node* n = R->depot->next;
    int fwdCount = 0;
    std::set<int> seenCustomers;

    while (n && n != R->depot && fwdCount <= params.nbClients + 5) {
        if (n->isDepot && n != R->depot) {
            std::cerr << "ERROR: Multiple depots in route " << R->cour << std::endl;
            return false;
        }

        if (!n->isDepot) {
            if (n->cour < 1 || n->cour > params.nbClients) {
                std::cerr << "ERROR: Invalid customer ID " << n->cour
                          << " in route " << R->cour << std::endl;
                return false;
            }
            if (!seenCustomers.insert(n->cour).second) {
                std::cerr << "ERROR: Duplicate customer " << n->cour
                          << " in route " << R->cour << std::endl;
                return false;
            }
            if (n->route != R) {
                std::cerr << "ERROR: Node " << n->cour
                          << " has incorrect route pointer in route " << R->cour << std::endl;
                return false;
            }
        }

        // Check consistency
        if (n->next->prev != n) {
            std::cerr << "ERROR: Pointer inconsistency at node "
                      << (n->isDepot ? "DEPOT" : std::to_string(n->cour))
                      << " in route " << R->cour << std::endl;
            return false;
        }

        n = n->next;
        fwdCount++;
    }

    if (n != R->depot) {
        std::cerr << "ERROR: Route " << R->cour << " is not a closed circular list" << std::endl;
        return false;
    }

    // ---- Backward traversal ----
    n = R->depot->prev;
    int backCount = 0;
    while (n && n != R->depot && backCount <= params.nbClients + 5) {
        if (n->next->prev != n || n->prev->next != n) {
            std::cerr << "ERROR: Broken backward link at node "
                      << (n->isDepot ? "DEPOT" : std::to_string(n->cour))
                      << " in route " << R->cour << std::endl;
            return false;
        }
        n = n->prev;
        backCount++;
    }
    if (n != R->depot) {
        std::cerr << "ERROR: Backward traversal could not return to depot in route "
                  << R->cour << std::endl;
        return false;
    }

    // ---- Count check ----
    if ((int)seenCustomers.size() != R->nbCustomers) {
        std::cerr << "ERROR: Customer count mismatch in route " << R->cour
                  << ": expected " << R->nbCustomers
                  << ", found " << seenCustomers.size() << std::endl;
        return false;
    }

    return true;
}
void LocalSearch::validateAllRoutes(const char* where) const {
    for (int r = 0; r < params.nbVehicles; ++r) {
        if (!validateRouteLinkedList(const_cast<Route*>(&routes[r]))) {
            std::cerr << "[CORRUPTION] " << where << " — route " << r << " failed\n";
            printRouteLinkedList(const_cast<Route*>(&routes[r]));
            std::abort();
        }
    }
}


void LocalSearch::run(Individual & indiv, double penaltyCapacityLS, double penaltyDurationLS, double penaltyTimeWarpLS)
{
    this->penaltyCapacityLS = penaltyCapacityLS;
    this->penaltyDurationLS = penaltyDurationLS;
    this->penaltyTimeWarpLS = penaltyTimeWarpLS;

    loadIndividual(indiv);

    for (int r = 0; r < params.nbVehicles; ++r)
        updateRouteData(&routes[r]);

    // Initialize total solution cost
    totalSolutionCost = 0.0;
    for (int r = 0; r < params.nbVehicles; ++r) {
        totalSolutionCost += routes[r].routeCost;
    }

    // DEBUG: Check initial state - DISABLED FOR STABILITY
    // std::cerr << "=== INITIAL STATE DEBUG ===" << std::endl;
    std::shuffle(orderNodes.begin(),  orderNodes.end(),  params.ran);
    std::shuffle(orderRoutes.begin(), orderRoutes.end(), params.ran);
    for (int i = 1; i <= params.nbClients; i++)
        if (params.ran() % params.ap.nbGranular == 0)
            std::shuffle(params.correlatedVertices[i].begin(), params.correlatedVertices[i].end(), params.ran);

    searchCompleted = false;
    for (loopID = 0; !searchCompleted; loopID++)
    {
        if (loopID > 1) searchCompleted = true;


        for (int posU = 0; posU < params.nbClients; posU++)
        {
            nodeU = &clients[orderNodes[posU]];
            int lastTestRINodeU = nodeU->whenLastTestedRI;
            nodeU->whenLastTestedRI = nbMoves;

            for (int posV = 0; posV < (int)params.correlatedVertices[nodeU->cour].size(); posV++)
            {
                nodeV = &clients[params.correlatedVertices[nodeU->cour][posV]];
                if (nodeU == nodeV) continue;
                if (loopID == 0 || std::max<int>(nodeU->route->whenLastModified, nodeV->route->whenLastModified) > lastTestRINodeU)
                {
                    setLocalVariablesRouteU();
                    setLocalVariablesRouteV();
                    if (move1()) continue;
                    if (move2()) continue;
                    if (move3()) continue;

                    if (nodeUIndex <= nodeVIndex && move4()) continue;
                    if (move5()) continue;
                    if (nodeUIndex <= nodeVIndex && move6()) continue;

                    if (routeU == routeV && move7()) continue;  // 2-opt (intra, reversal) //LAST BUG HERE
                    if (routeU != routeV && move8()) continue;  // 2-opt* (inter)
                    if (routeU != routeV && move9()) continue;  // cross-exchange

                    if (nodeV->prev->isDepot)
                    {
                        nodeV = nodeV->prev;
                        setLocalVariablesRouteV();
                        if (move1()) continue;
                        if (move2()) continue;
                        if (move3()) continue;
                        if (routeU == routeV && move7()) continue;
                        if (routeU != routeV && move8()) continue;
                        if (routeU != routeV && move9()) continue;
                    }
                }
            }

            if (loopID > 0 && !emptyRoutes.empty())
            {
                nodeV = routes[*emptyRoutes.begin()].depot;
                setLocalVariablesRouteU();
                setLocalVariablesRouteV();
                if (move1()) continue;
                if (move2()) continue;
                if (move3()) continue;
                if (move9()) continue;
            }
        }

        if (params.ap.useSwapStar == 1 && params.areCoordinatesProvided)
        {
            for (int rU = 0; rU < params.nbVehicles; rU++)
            {
                routeU = &routes[orderRoutes[rU]];
                int lastTestSWAPStarRouteU = routeU->whenLastTestedSWAPStar;
                routeU->whenLastTestedSWAPStar = nbMoves;
                for (int rV = 0; rV < params.nbVehicles; rV++)
                {
                    routeV = &routes[orderRoutes[rV]];
                    if (routeU->nbCustomers > 0 && routeV->nbCustomers > 0 && routeU->cour < routeV->cour
                        && (loopID == 0 || std::max<int>(routeU->whenLastModified, routeV->whenLastModified) > lastTestSWAPStarRouteU))
                        if (CircleSector::overlap(routeU->sector, routeV->sector))
                            swapStar();
                }
            }
        }
    }
    exportIndividual(indiv);
    if (loopID == 0 || nbMoves > 0) {
        std::cerr << "LS finished after " << loopID << " loops and "
                  << nbMoves << " moves. Final cost: " << totalSolutionCost << "\n";
    }
}

void LocalSearch::setLocalVariablesRouteU()
{
    routeU = nodeU->route;
    nodeX  = nodeU->next;

    // indices must be positions along the route, not customer IDs
    nodeUIndex     = nodeU->position;
    nodeUPrevIndex = nodeU->prev->isDepot ? 0 : nodeU->prev->position;

    nodeXIndex     = nodeX->isDepot ? (nodeUIndex + 1) : nodeX->position;
}

void LocalSearch::setLocalVariablesRouteV()
{
    routeV = nodeV->route;
    nodeY  = nodeV->next;

    nodeVIndex     = nodeV->position;
    nodeVPrevIndex = nodeV->prev->isDepot ? 0 : nodeV->prev->position;

    nodeYIndex     = nodeY->isDepot ? (nodeVIndex + 1) : nodeY->position;

    intraRouteMove = (routeU == routeV);
}



// =========================== MOVES ===========================

bool LocalSearch::move1()
{
    if (routeU == routeV) return false;
    if (!nodeU || !nodeV) return false;
    if (nodeU == nodeV)   return false;
    if (nodeU->isDepot || nodeV->isDepot) return false;
    if (routeU == routeV && nodeV == nodeU->prev) return false;

    Route* RU_before = routeU;
    Route* RV_before = routeV;
    const double oldTotalCost = totalSolutionCost;

    // 🔑 O(1) old costs from prefix segments
    double oldCostU = route_cost_from_seg(params, prefSeg[routeU->cour][routeU->nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    double oldCostV = route_cost_from_seg(params, prefSeg[routeV->cour][routeV->nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

    // Build sequences after move
    std::vector<int> seqU, seqV;
    routeToSequence(RU_before, seqU);
    routeToSequence(RV_before, seqV);

    // Remove nodeU from U
    std::vector<int> newSeqU = seqU;
    int posU = nodeU->position - 1;
    if (posU >= 0 && posU < (int)newSeqU.size())
        newSeqU.erase(newSeqU.begin() + posU);

    // Insert nodeU into V after nodeV
    std::vector<int> newSeqV = seqV;
    int insertPos = nodeV->position;
    if (insertPos >= 0 && insertPos <= static_cast<int>(newSeqV.size()))
        newSeqV.insert(newSeqV.begin() + insertPos, nodeU->cour);

    const TWSPD_Seg TU_after = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
    const TWSPD_Seg TV_after = simulateExactSegment(newSeqV, params);

    double newCostU = route_cost_from_seg(params, TU_after, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    double newCostV = route_cost_from_seg(params, TV_after, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

    double predictedNewCost = oldTotalCost - oldCostU - oldCostV + newCostU + newCostV;
    if (predictedNewCost >= oldTotalCost) return false;

    // Apply move
    insertNode(nodeU, nodeV);
    nbMoves++;
    searchCompleted = false;
    totalSolutionCost = predictedNewCost;

    updateRouteData(RU_before);
    if (RU_before != RV_before) updateRouteData(RV_before);

    LS_logMove("move1", RU_before, RV_before, nodeU, nullptr, nodeV, nullptr,
               oldTotalCost, predictedNewCost);

    return true;
}



bool LocalSearch::move2()
{
    if (routeU == routeV) return false;
    if (!nodeU || !nodeX || !nodeV) return false;
    if (nodeU->next != nodeX) return false;       // must be consecutive
    if (nodeX->isDepot) return false;
    if (nodeU == nodeV || nodeX == nodeV) return false;

    Route* RU_before = routeU;
    Route* RV_before = routeV;
    const double oldTotalCost = totalSolutionCost;

    // Helper to fetch current cost safely
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };
    const double oldU = curCost(routeU);
    const double oldV = curCost(routeV);

    // Convert routes to sequences
    std::vector<int> seqU, seqV;
    routeToSequence(routeU, seqU);
    routeToSequence(routeV, seqV);

    // Remove nodeU and nodeX from routeU (by position)
    std::vector<int> newSeqU = seqU;
    int posU = nodeU->position - 1;  // nodeU index
    newSeqU.erase(newSeqU.begin() + posU, newSeqU.begin() + posU + 2);

    // Insert [nodeU, nodeX] into routeV after nodeV
    std::vector<int> newSeqV = seqV;
    int insertPos = nodeV->position; // after nodeV
    newSeqV.insert(newSeqV.begin() + insertPos, {nodeU->cour, nodeX->cour});

    // Re-simulate both routes
    const TWSPD_Seg TU = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
    const TWSPD_Seg TV = simulateExactSegment(newSeqV, params);

    const double newRouteCostU = route_cost_from_seg(params, TU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newRouteCostV = route_cost_from_seg(params, TV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

    const double newTotalCost = oldTotalCost - oldU - oldV + newRouteCostU + newRouteCostV;

    if (newTotalCost >= oldTotalCost) return false;

    // Apply actual move on solution
    insertNode(nodeU, nodeV);
    insertNode(nodeX, nodeU);
    nbMoves++; searchCompleted = false;
    totalSolutionCost = newTotalCost;

    updateRouteData(RU_before);
    if (RU_before != RV_before) updateRouteData(RV_before);

    LS_logMove("move2", routeU, routeV, nodeU, nodeX, nodeV, nullptr,
               oldTotalCost, newTotalCost);

    return true;
}


bool LocalSearch::move3()
{
    if (routeU == routeV) return false;
    if (!nodeU || !nodeX || !nodeV) return false;
    if (nodeU->next != nodeX) return false;   // must be consecutive
    if (nodeX->isDepot) return false;
    if (nodeU == nodeV || nodeX == nodeV) return false;

    Route* RU_before = routeU;
    Route* RV_before = routeV;
    const double oldTotalCost = totalSolutionCost;

    // Fetch old costs from prefSeg (avoid stale routeCost)
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };
    const double oldU = curCost(routeU);
    const double oldV = curCost(routeV);

    // Convert routes to sequences
    std::vector<int> seqU, seqV;
    routeToSequence(routeU, seqU);
    routeToSequence(routeV, seqV);

    // Remove nodeU and nodeX from routeU (by index range)
    std::vector<int> newSeqU = seqU;
    int posU = nodeU->position - 1;  // index of nodeU
    newSeqU.erase(newSeqU.begin() + posU, newSeqU.begin() + posU + 2);

    // Insert [X, U] into routeV after nodeV
    std::vector<int> newSeqV = seqV;
    int insertPos = nodeV->position; // after nodeV
    newSeqV.insert(newSeqV.begin() + insertPos, {nodeX->cour, nodeU->cour});

    // Re-simulate both routes
    const TWSPD_Seg TU = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
    const TWSPD_Seg TV = simulateExactSegment(newSeqV, params);

    const double newRouteCostU = route_cost_from_seg(params, TU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newRouteCostV = route_cost_from_seg(params, TV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

    const double newTotalCost = oldTotalCost - oldU - oldV + newRouteCostU + newRouteCostV;

    if (newTotalCost >= oldTotalCost) return false;

    // Apply actual move on solution
    insertNode(nodeX, nodeV);   // put X after V
    insertNode(nodeU, nodeX);   // then put U after X
    nbMoves++; searchCompleted = false;
    totalSolutionCost = newTotalCost;

    updateRouteData(RU_before);
    if (RU_before != RV_before) updateRouteData(RV_before);

    LS_logMove("move3", routeU, routeV, nodeU, nodeX, nodeV, nullptr,
               oldTotalCost, newTotalCost);

    return true;
}

inline void LocalSearch::swapSingletonsSafe(Node* a, Node* b)
{
    if (!a || !b || a == b) return;

    // Basic validation
    if (!(a->prev && a->next && a->prev->next == a && a->next->prev == a)) std::abort();
    if (!(b->prev && b->next && b->prev->next == b && b->next->prev == b)) std::abort();

    Route* ra = a->route;
    Route* rb = b->route;

    // Check if nodes are adjacent - this requires special handling
    bool adjacent = (a->next == b) || (b->next == a);

    if (adjacent) {
        // For adjacent nodes, need a simple reversal operation
        Node* first = (a->next == b) ? a : b;
        Node* second = (a->next == b) ? b : a;

        // Save the outer neighbors
        Node* beforeFirst = first->prev;
        Node* afterSecond = second->next;

        // Reverse the pair: beforeFirst -> second -> first -> afterSecond
        beforeFirst->next = second;
        second->prev = beforeFirst;
        second->next = first;
        first->prev = second;
        first->next = afterSecond;
        afterSecond->prev = first;

        // Swap route pointers if different
        if (ra != rb) {
            a->route = rb;
            b->route = ra;
        }
    } else {
        // Non-adjacent case: original logic
        // Save neighbors
        Node* ap = a->prev;
        Node* an = a->next;
        Node* bp = b->prev;
        Node* bn = b->next;

        // --- Detach both completely ---
        ap->next = an;
        an->prev = ap;
        bp->next = bn;
        bn->prev = bp;

        // --- Insert b where a was ---
        ap->next = b;
        b->prev = ap;
        b->next = an;
        an->prev = b;

        // --- Insert a where b was ---
        bp->next = a;
        a->prev = bp;
        a->next = bn;
        bn->prev = a;

        // --- Swap route pointers if different ---
        if (ra != rb) {
            a->route = rb;
            b->route = ra;
            // nbCustomers stay balanced
        }
    }

    // Basic validation
    if (!(a->prev && a->next && a->prev->next == a && a->next->prev == a)) std::abort();
    if (!(b->prev && b->next && b->prev->next == b && b->next->prev == b)) std::abort();
}




bool LocalSearch::move4()
{
    if (!nodeU || !nodeV) return false;
    if (nodeU == nodeV) return false;
    if (nodeU->isDepot || nodeV->isDepot) return false;

    const double oldTotalCost = totalSolutionCost;

    // Build sequences
    std::vector<int> seqU, seqV;
    routeToSequence(routeU, seqU);
    routeToSequence(routeV, seqV);

    double newTotalCost;

    // Compute "old" costs from prefSeg (avoid stale route->routeCost)
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };

    const double oldU = curCost(routeU);
    const double oldV = (routeV != routeU) ? curCost(routeV) : 0.0;

    if (routeU == routeV)
    {
        // Intra-route swap: swap at positions
        std::vector<int> newSeq = seqU;
        int posU = nodeU->position - 1;
        int posV = nodeV->position - 1;
        std::swap(newSeq[posU], newSeq[posV]);

        const TWSPD_Seg T = simulateExactSegment(newSeq, params);
        const double newRouteCost = route_cost_from_seg(params, T, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

        newTotalCost = oldTotalCost - oldU + newRouteCost;
        if (newTotalCost >= oldTotalCost) return false;

        // Apply real move
        swapSingletonsSafe(nodeU, nodeV);
        nbMoves++; searchCompleted = false;
        totalSolutionCost = newTotalCost;

        updateRouteData(routeU);
        LS_logMove("move4", routeU, routeU, nodeU, nodeX, nodeV, nullptr,
                   oldTotalCost, newTotalCost);
        return true;
    }
    else
    {
        // Inter-route swap: replace by position
        std::vector<int> newSeqU = seqU;
        std::vector<int> newSeqV = seqV;

        int posU = nodeU->position - 1;
        int posV = nodeV->position - 1;

        newSeqU[posU] = nodeV->cour;
        newSeqV[posV] = nodeU->cour;

        const TWSPD_Seg TU = simulateExactSegment(newSeqU, params);
        const TWSPD_Seg TV = simulateExactSegment(newSeqV, params);

        const double newRouteCostU = route_cost_from_seg(params, TU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
        const double newRouteCostV = route_cost_from_seg(params, TV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

        newTotalCost = oldTotalCost - oldU - oldV + newRouteCostU + newRouteCostV;
        if (newTotalCost >= oldTotalCost) return false;

        // Apply real move
        swapSingletonsSafe(nodeU, nodeV);
        nbMoves++; searchCompleted = false;
        totalSolutionCost = newTotalCost;

        updateRouteData(routeU);
        updateRouteData(routeV);
        LS_logMove("move4", routeU, routeV, nodeU, nodeX, nodeV, nullptr,
                   oldTotalCost, newTotalCost);
        return true;
    }
}


bool LocalSearch::move5()
{
    if (routeU == routeV) return false;
    if (!nodeU || !nodeX || !nodeV) return false;
    if (nodeU->next != nodeX) return false;     // must be consecutive
    if (nodeX->isDepot || nodeV->isDepot) return false;
    if (nodeU == nodeV || nodeX == nodeV) return false;

    const double oldTotalCost = totalSolutionCost;

    // Safe old costs from prefSeg
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };
    const double oldU = curCost(routeU);
    const double oldV = curCost(routeV);

    // Convert routes to sequences
    std::vector<int> seqU, seqV;
    routeToSequence(routeU, seqU);
    routeToSequence(routeV, seqV);

    // --- Simulate routeU after swap ---
    std::vector<int> newSeqU = seqU;
    int posU = nodeU->position - 1; // index of U
    // erase [U,X]
    newSeqU.erase(newSeqU.begin() + posU, newSeqU.begin() + posU + 2);
    // insert V at posU
    newSeqU.insert(newSeqU.begin() + posU, nodeV->cour);

    // --- Simulate routeV after swap ---
    std::vector<int> newSeqV = seqV;
    int posV = nodeV->position - 1; // index of V
    // erase V
    newSeqV.erase(newSeqV.begin() + posV);
    // insert [U,X] at posV
    newSeqV.insert(newSeqV.begin() + posV, {nodeU->cour, nodeX->cour});

    // Re-simulate
    const TWSPD_Seg TU = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
    const TWSPD_Seg TV = simulateExactSegment(newSeqV, params);

    const double newRouteCostU = route_cost_from_seg(params, TU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newRouteCostV = route_cost_from_seg(params, TV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

    const double newTotalCost = oldTotalCost - oldU - oldV + newRouteCostU + newRouteCostV;
    if (newTotalCost >= oldTotalCost) return false;

    // Apply the actual move: swap U with V, then insert X after U
    swapNode(nodeU, nodeV);
    insertNode(nodeX, nodeU);

    nbMoves++; searchCompleted = false;
    totalSolutionCost = newTotalCost;

    updateRouteData(routeU);
    if (routeU != routeV) updateRouteData(routeV);

    LS_logMove("move5", routeU, routeV, nodeU, nodeX, nodeV, nullptr,
               oldTotalCost, newTotalCost);

    return true;
}

bool LocalSearch::move6()
{
    if (!nodeU || !nodeX || !nodeV || !nodeY) return false;
    if (routeU == routeV) return false;
    if (nodeU->next != nodeX || nodeV->next != nodeY) return false; // pairs must be consecutive
    if (nodeX->isDepot || nodeY->isDepot) return false;
    if (nodeU == nodeV || nodeX == nodeY || nodeU == nodeY || nodeX == nodeV) return false;

    Route* RU_before = routeU;
    Route* RV_before = routeV;

    const double oldTotalCost = totalSolutionCost;

    // Current route costs from segment caches (avoid stale route->routeCost)
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };
    const double oldUcost = curCost(routeU);
    const double oldVcost = curCost(routeV);

    // Build sequences
    std::vector<int> seqU, seqV;
    routeToSequence(routeU, seqU);
    routeToSequence(routeV, seqV);

    // Positions (0-based)
    const int posU = nodeU->position - 1; // [U, X] at posU, posU+1
    const int posV = nodeV->position - 1; // [V, Y] at posV, posV+1

    // --- Simulate routeU after replacing [U,X] with [V,Y] at posU ---
    std::vector<int> newSeqU = seqU;
    newSeqU.erase(newSeqU.begin() + posU, newSeqU.begin() + posU + 2);        // remove [U,X]
    newSeqU.insert(newSeqU.begin() + posU, { nodeV->cour, nodeY->cour });     // insert [V,Y]

    // --- Simulate routeV after replacing [V,Y] with [U,X] at posV ---
    std::vector<int> newSeqV = seqV;
    newSeqV.erase(newSeqV.begin() + posV, newSeqV.begin() + posV + 2);        // remove [V,Y]
    newSeqV.insert(newSeqV.begin() + posV, { nodeU->cour, nodeX->cour });     // insert [U,X]

    // Re-simulate both routes using the exact same oracle as elsewhere
    const TWSPD_Seg TU = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
    const TWSPD_Seg TV = newSeqV.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqV, params);

    const double newRouteCostU = route_cost_from_seg(params, TU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newRouteCostV = route_cost_from_seg(params, TV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

    const double predictedNewTotal =
        oldTotalCost - oldUcost - oldVcost + newRouteCostU + newRouteCostV;

    if (predictedNewTotal >= oldTotalCost) return false;

    // --- Apply the actual move on linked structure ---
    // Insert semantics assumed: insertNode(A, B) inserts A immediately AFTER B.
    // Want [V,Y] to appear where [U,X] was, and [U,X] where [V,Y] was.
    Node* Uprev = nodeU->prev; // anchor before U in routeU
    Node* Vprev = nodeV->prev; // anchor before V in routeV

    // Move V,Y pair into routeU at U's position (after Uprev)
    insertNode(nodeV, Uprev);
    insertNode(nodeY, nodeV);

    // Move U,X pair into routeV at V's position (after Vprev)
    insertNode(nodeU, Vprev);
    insertNode(nodeX, nodeU);

    nbMoves++;
    searchCompleted = false;
    totalSolutionCost = predictedNewTotal;

    // Update routes (use where the nodes ended up)
    // After the operations: V,Y now belong to former routeU; U,X to former routeV.
    updateRouteData(RU_before);
    if (RU_before != RV_before) updateRouteData(RV_before);

    LS_logMove("move6", RU_before, RV_before, nodeU, nodeX, nodeV, nodeY,
               oldTotalCost, predictedNewTotal);

    return true;
}


bool LocalSearch::move7()
{
    if (routeU != routeV) return false;
    if (!nodeU || !nodeV) return false;
    if (nodeU->isDepot || nodeV->isDepot) return false;
    if (nodeU->next == nodeV || nodeV->next == nodeU) return false;


    const int posU = nodeU->position - 1;
    const int posV = nodeV->position - 1;
    if (posU >= posV) return false;

    const double oldTotalCost = totalSolutionCost;

    // Safe old cost from prefSeg
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };
    const double oldCostU = curCost(routeU);

    // Build new sequence with [posU+1 .. posV] reversed
    std::vector<int> seq;
    routeToSequence(routeU, seq);

    std::vector<int> newSeq;
    newSeq.insert(newSeq.end(), seq.begin(), seq.begin() + posU + 1);                // prefix
    newSeq.insert(newSeq.end(), seq.rbegin() + (seq.size() - posV - 1),              // reversed middle
                  seq.rbegin() + (seq.size() - (posU + 1)));
    newSeq.insert(newSeq.end(), seq.begin() + posV + 1, seq.end());                  // suffix

    const TWSPD_Seg newSeg = simulateExactSegment(newSeq, params);
    const double newRouteCost = route_cost_from_seg(params, newSeg, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newTotalCost = oldTotalCost - oldCostU + newRouteCost;

    if (newTotalCost >= oldTotalCost) return false;

    // Apply reversal of [U.next .. V]
    Node* first = nodeU->next;
    Node* last  = nodeV;
    Node* before = first->prev;
    Node* after  = last->next;

    // Safer block collection: collect nodes from first up to and INCLUDING last
    std::vector<Node*> block;
    Node* n = first;
    while (true) {
        block.push_back(n);
        if (n == last) break;  // Stop after including the last node
        n = n->next;

        // Safety check to prevent infinite loops
        if (block.size() > 1000) {
            return false;
        }
    }

    // Cut block out
    before->next = after;
    after->prev = before;

    // CRITICAL FIX: Clear all pointers in the block to prevent corruption
    for (Node* node : block) {
        node->prev = nullptr;
        node->next = nullptr;
    }

    // Reinsert nodes in reverse order
    Node* anchor = before;
    for (int i = (int)block.size() - 1; i >= 0; --i) {
        insertNode(block[i], anchor);
        anchor = block[i];
    }

    nbMoves++;
    searchCompleted = false;
    totalSolutionCost = newTotalCost;

    updateRouteData(routeU);

    LS_logMove("move7", routeU, routeU, nodeU, nullptr, nodeV, nullptr,
               oldTotalCost, newTotalCost);
    return true;
}

bool LocalSearch::move8()
{
    if (routeU == routeV) return false;
    if (!nodeU || !nodeV) return false;
    if (nodeU->isDepot || nodeV->isDepot) return false;

    Node* X = nodeU->next; // head of U-tail
    Node* Y = nodeV->next; // head of V-tail
    if (!X || !Y || X->isDepot || Y->isDepot) return false; // forbid empty tails

    const int posU = nodeU->position - 1;
    const int posV = nodeV->position - 1;

    const double oldTotalCost = totalSolutionCost;

    // Old route costs from prefSeg (O(1), avoids stale routeCost)
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };
    const double oldCostU = curCost(routeU);
    const double oldCostV = curCost(routeV);

    // Build sequences
    std::vector<int> seqU, seqV;
    routeToSequence(routeU, seqU);
    routeToSequence(routeV, seqV);

    // New sequence for routeU: prefixU + tailV
    std::vector<int> newSeqU(seqU.begin(), seqU.begin() + posU + 1);
    newSeqU.insert(newSeqU.end(), seqV.begin() + posV + 1, seqV.end());

    // New sequence for routeV: prefixV + tailU
    std::vector<int> newSeqV(seqV.begin(), seqV.begin() + posV + 1);
    newSeqV.insert(newSeqV.end(), seqU.begin() + posU + 1, seqU.end());

    const TWSPD_Seg TU = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
    const TWSPD_Seg TV = newSeqV.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqV, params);

    const double newRouteCostU = route_cost_from_seg(params, TU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newRouteCostV = route_cost_from_seg(params, TV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

    const double newTotalCost = oldTotalCost - oldCostU - oldCostV + newRouteCostU + newRouteCostV;
    if (newTotalCost >= oldTotalCost) return false;

    // Snapshot original routes
    Route* RU0 = routeU;
    Route* RV0 = routeV;

    // Collect tail nodes BEFORE edits
    std::vector<Node*> tailNodesU;
    for (Node* n = X; !n->isDepot; n = n->next) tailNodesU.push_back(n);
    std::vector<Node*> tailNodesV;
    for (Node* n = Y; !n->isDepot; n = n->next) tailNodesV.push_back(n);

    // Apply: insert V-tail after U
    Node* anchor = nodeU;
    for (Node* n : tailNodesV) { insertNode(n, anchor); anchor = n; }

    // Then insert U-tail after V
    anchor = nodeV;
    for (Node* n : tailNodesU) { insertNode(n, anchor); anchor = n; }

    nbMoves++; searchCompleted = false;
    totalSolutionCost = newTotalCost;

    updateRouteData(RU0);
    if (RV0 != RU0) updateRouteData(RV0);

    // Sanity: routes should not be empty
    if (RU0->nbCustomers == 0 || RV0->nbCustomers == 0) {
        std::cerr << "[LS][FATAL] move8 produced empty route; aborting.\n";
        std::abort();
    }

    LS_logMove("move8", RU0, RV0, nodeU, X, nodeV, Y, oldTotalCost, newTotalCost);
    return true;
}

bool LocalSearch::move9()
{
    if (routeU == routeV) return false;
    if (!nodeU || !nodeV) return false;
    if (nodeU->isDepot || nodeV->isDepot) return false;

    Node* X = nodeU->next; // head of U-tail
    Node* Y = nodeV->next; // head of V-tail
    if (!X || !Y || X->isDepot || Y->isDepot) return false;

    const int posU = nodeU->position - 1;
    const int posV = nodeV->position - 1;

    const double oldTotalCost = totalSolutionCost;

    // Safe old costs from prefSeg
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };
    const double oldCostU = curCost(routeU);
    const double oldCostV = curCost(routeV);

    // Build new sequences
    std::vector<int> seqU, seqV;
    routeToSequence(routeU, seqU);
    routeToSequence(routeV, seqV);

    // RouteU: prefixU + reversed tailV
    std::vector<int> newSeqU(seqU.begin(), seqU.begin() + posU + 1);
    newSeqU.insert(newSeqU.end(), seqV.rbegin(), seqV.rend() - (posV + 1));

    // RouteV: prefixV + reversed tailU
    std::vector<int> newSeqV(seqV.begin(), seqV.begin() + posV + 1);
    newSeqV.insert(newSeqV.end(), seqU.rbegin(), seqU.rend() - (posU + 1));

    const TWSPD_Seg TU = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
    const TWSPD_Seg TV = newSeqV.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqV, params);

    const double newRouteCostU = route_cost_from_seg(params, TU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newRouteCostV = route_cost_from_seg(params, TV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    const double newTotalCost = oldTotalCost - oldCostU - oldCostV + newRouteCostU + newRouteCostV;

    if (newTotalCost >= oldTotalCost) return false;

    // Snapshot routes
    Route* RU0 = routeU;
    Route* RV0 = routeV;

    // Gather tails before edits
    std::vector<Node*> tailNodesU;
    for (Node* n = X; !n->isDepot; n = n->next) tailNodesU.push_back(n);
    std::vector<Node*> tailNodesV;
    for (Node* n = Y; !n->isDepot; n = n->next) tailNodesV.push_back(n);

    // Insert reversed V-tail after U
    Node* anchor = nodeU;
    for (int i = (int)tailNodesV.size() - 1; i >= 0; --i) {
        insertNode(tailNodesV[i], anchor);
        anchor = tailNodesV[i];
    }

    // Insert reversed U-tail after V
    anchor = nodeV;
    for (int i = (int)tailNodesU.size() - 1; i >= 0; --i) {
        insertNode(tailNodesU[i], anchor);
        anchor = tailNodesU[i];
    }

    nbMoves++;
    searchCompleted = false;
    totalSolutionCost = newTotalCost;

    updateRouteData(RU0);
    if (RV0 != RU0) updateRouteData(RV0);

    if (RU0->nbCustomers == 0 || RV0->nbCustomers == 0) {
        std::cerr << "[LS][FATAL] move9 produced empty route; aborting.\n";
        std::abort();
    }

    LS_logMove("move9", RU0, RV0, nodeU, nullptr, nodeV, nullptr,
               oldTotalCost, newTotalCost);
    return true;
}



// =========================== HELPERS ===========================
// Insert node U right after node V in its route.
// Detaches U from its current position if necessary.
// Insert node U right after node V (possibly same route). Safe for intra/inter.
inline void LocalSearch::insertNode(Node* U, Node* V)
{
    if (!U || !V || U == V || V->next == U) return;

    Route* oldRoute = U->route;
    Route* newRoute = V->route;

    // 1. Detach U from its current place
    Node* up = U->prev;
    Node* un = U->next;
    if (up) up->next = un;
    if (un) un->prev = up;

    // If removing from an old route, decrement its count
    if (oldRoute && oldRoute != newRoute) {
        oldRoute->nbCustomers--;
    }

    // 2. Insert U after V in V's route
    Node* vn = V->next;
    V->next = U;
    U->prev = V;
    U->next = vn;
    if (vn) vn->prev = U;

    // Update route pointer and counts
    if (oldRoute != newRoute) {
        newRoute->nbCustomers++;
    }
    U->route = newRoute;

}





// Swap two nodes U and V, possibly across routes.
inline void LocalSearch::swapNode(Node* U, Node* V)
{
    if (!U || !V || U == V) return;

    Route* rU = U->route;
    Route* rV = V->route;

    Node* Up = U->prev; Node* Un = U->next;
    Node* Vp = V->prev; Node* Vn = V->next;

    // --- Detach both ---
    Up->next = Un; Un->prev = Up;
    Vp->next = Vn; Vn->prev = Vp;

    // --- Insert V where U was ---
    Up->next = V; V->prev = Up;
    V->next = Un; Un->prev = V;

    // --- Insert U where V was ---
    Vp->next = U; U->prev = Vp;
    U->next = Vn; Vn->prev = U;

    // --- Route bookkeeping ---
    if (rU != rV) {
        U->route = rV;
        V->route = rU;
        // nbCustomers doesn’t change: each route loses one, gains one
    }
}


// Rebuild all routes from the Individual into LocalSearch data structures.
void LocalSearch::loadIndividual(const Individual& indiv)
{
    // First, completely reset all nodes and routes
    for (int i = 1; i <= params.nbClients; ++i) {
        Node& n = clients[i];
        n.route = nullptr;
        n.prev = nullptr;
        n.next = nullptr;
        n.position = 0;
        n.whenLastTestedRI = -1;
    }

    // Reset all routes to depot-only rings
    emptyRoutes.clear();
    for (int r = 0; r < params.nbVehicles; ++r) {
        Route& R = routes[r];
        R.depot->next = R.depot;
        R.depot->prev = R.depot;
        R.nbCustomers = 0;
        R.whenLastModified = nbMoves;
        if (indiv.chromR[r].empty()) emptyRoutes.insert(r);
    }

    // Link customer nodes according to chromR
    for (int r = 0; r < params.nbVehicles; ++r) {
        Route& R = routes[r];
        Node* prev = R.depot;

        const auto& routeSeq = indiv.chromR[r];
        for (int idx = 0; idx < (int)routeSeq.size(); ++idx) {
            int c = routeSeq[idx];
            Node* n = &clients[c];

            // CRITICAL FIX: Ensure node is completely detached before reusing
            if (n->prev) n->prev->next = n->next;
            if (n->next) n->next->prev = n->prev;
            n->prev = nullptr;
            n->next = nullptr;
            n->route = nullptr;

            // Insert after 'prev'
            n->route = &R;
            n->prev = prev;
            n->next = prev->next;
            prev->next->prev = n;
            prev->next = n;

            prev = n;
            R.nbCustomers++;
        }

        // Close the ring: last customer -> depot
        if (prev != R.depot) {  // Only if added customers
            prev->next = R.depot;
            R.depot->prev = prev;
        }
    }
}


// Export LocalSearch routes back into an Individual
void LocalSearch::exportIndividual(Individual& indiv)
{
    // routes → chromR
    indiv.chromR.clear();
    indiv.chromR.resize(params.nbVehicles);

    // clear adjacency (index 0 is depot)
    std::fill(indiv.successors.begin(),   indiv.successors.end(),   0);
    std::fill(indiv.predecessors.begin(), indiv.predecessors.end(), 0);

    for (int r = 0; r < params.nbVehicles; ++r)
    {
        Route& R = routes[r];
        Node* prev = R.depot;          // start at depot (0)
        for (Node* n = R.depot->next; !n->isDepot; n = n->next)
        {
            indiv.chromR[r].push_back(n->cour);

            // prev -> n arc
            indiv.successors[prev->cour]   = n->cour;
            indiv.predecessors[n->cour]    = prev->cour;

            prev = n;
        }
        // close route n_last -> depot
        if (prev != R.depot)
        {
            indiv.successors[prev->cour] = 0;
            // If you rely on depot predecessor for anything, you can set:
            // indiv.predecessors[0] = prev->cour;   // optional
        }
    }
}



// ======================= SWAP* OPERATOR =======================
// Full SWAP* (relocate + exchange) adapted from Vidal HGS-CVRP,
// using segment concatenations (rangeFwd / rangeRev) for O(1) eval.

bool LocalSearch::swapStar()
{
    auto curCost = [&](const Route* R) {
        if (!R || routes[R->cour].nbCustomers == 0) return 0.0;
        return route_cost_from_seg(params, prefSeg[R->cour][routes[R->cour].nbCustomers - 1], penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
    };

    for (int rU = 0; rU < params.nbVehicles; rU++)
    {
        Route* RU = &routes[rU];
        if (RU->nbCustomers == 0) continue;

        for (int rV = rU + 1; rV < params.nbVehicles; rV++)
        {
            Route* RV = &routes[rV];
            if (RV->nbCustomers == 0) continue;

            // pruning by geometry
            if (!CircleSector::overlap(RU->sector, RV->sector))
                continue;

            // -------------------------------
            // RELOCATE: move U from RU to RV
            // -------------------------------
            for (Node* U = RU->depot->next; !U->isDepot; U = U->next)
            {
                Node* X = U->next;
                if (X->isDepot) continue;  // nothing to relocate (U is last)

                const double oldCostRU = curCost(RU);
                const double oldCostRV = curCost(RV);

                // Build modified sequence for RU (remove U by index)
                std::vector<int> seqU;
                routeToSequence(RU, seqU);
                const int posU = U->position - 1;

                std::vector<int> newSeqU = seqU;
                newSeqU.erase(newSeqU.begin() + posU);

                const TWSPD_Seg newSegU = newSeqU.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqU, params);
                const double newRouteCostU = route_cost_from_seg(params, newSegU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

                // Try all insertion anchors V in RV (including depot to allow front insertion)
                for (Node* V = RV->depot; ; V = V->next)
                {
                    // Build modified sequence for RV: insert U **after** V
                    std::vector<int> seqV;
                    routeToSequence(RV, seqV);
                    const int posV = V->position - 1;   // depot -> -1

                    std::vector<int> newSeqV = seqV;
                    newSeqV.insert(newSeqV.begin() + (posV + 1), U->cour);

                    const TWSPD_Seg newSegV = newSeqV.empty() ? TWSPD_Seg{} : simulateExactSegment(newSeqV, params);
                    const double newRouteCostV = route_cost_from_seg(params, newSegV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

                    const double newTotalCost =
                        totalSolutionCost - oldCostRU - oldCostRV + newRouteCostU + newRouteCostV;

                    if (newTotalCost < totalSolutionCost - MY_EPSILON)
                    {
                        // Apply relocate: detach U from RU, insert after V in RV
                        // Safer to use helpers to preserve invariants
                        // Remove U from its current position:
                        U->prev->next = U->next;
                        U->next->prev = U->prev;
                        U->route = RU; // still RU until insert

                        // Insert after V
                        insertNode(U, V);

                        nbMoves++; searchCompleted = false;
                        totalSolutionCost = newTotalCost;
                        updateRouteData(RU);
                        updateRouteData(RV);
                        return true; // restart after improvement
                    }
                    if (V->next->isDepot) break;
                }
            }

            // -------------------------------------------
            // EXCHANGE: swap singleton U (in RU) with V (in RV)
            // -------------------------------------------
            for (Node* U = RU->depot->next; !U->isDepot; U = U->next)
            {
                if (U->next->isDepot) continue; // (your original guard kept this; keep behavior)
                const int posU = U->position - 1;

                const double oldCostRU = curCost(RU);
                const double oldCostRV = curCost(RV);

                for (Node* V = RV->depot->next; !V->isDepot; V = V->next)
                {
                    if (V->next->isDepot) continue; // keep original behavior
                    const int posV = V->position - 1;

                    // Build modified sequence for RU (replace U with V at posU)
                    std::vector<int> seqU, seqV;
                    routeToSequence(RU, seqU);
                    routeToSequence(RV, seqV);

                    std::vector<int> newSeqU = seqU;
                    newSeqU[posU] = V->cour;

                    // Build modified sequence for RV (replace V with U at posV)
                    std::vector<int> newSeqV = seqV;
                    newSeqV[posV] = U->cour;

                    const TWSPD_Seg newSegU = simulateExactSegment(newSeqU, params);
                    const TWSPD_Seg newSegV = simulateExactSegment(newSeqV, params);

                    const double newRouteCostU = route_cost_from_seg(params, newSegU, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);
                    const double newRouteCostV = route_cost_from_seg(params, newSegV, penaltyCapacityLS, penaltyDurationLS, penaltyTimeWarpLS);

                    const double newTotalCost =
                        totalSolutionCost - oldCostRU - oldCostRV + newRouteCostU + newRouteCostV;

                    if (newTotalCost < totalSolutionCost - MY_EPSILON)
                    {
                        // Apply swap of U and V
                        swapNode(U, V);

                        nbMoves++; searchCompleted = false;
                        totalSolutionCost = newTotalCost;
                        updateRouteData(RU);
                        updateRouteData(RV);
                        return true; // restart after improvement
                    }
                }
            }
        }
    }
    return false;
}




