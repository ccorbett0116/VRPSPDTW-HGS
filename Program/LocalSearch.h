/*MIT License

Copyright(c) 2020 Thibaut Vidal
... (keep your existing license header) ...
*/
#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include "Individual.h"

struct Node;

// Structure containing a route
struct Route
{
    int    cour;                     // Route index
    int    nbCustomers;              // Number of customers visited in the route
    int    whenLastModified;         // "When" this route has been last modified
    int    whenLastTestedSWAPStar;   // "When" the SWAP* moves for this route have been last tested
    Node * depot;                    // Pointer to the associated depot

    // Legacy fields (kept for pruning & book-keeping)
    double duration;                 // Sum of timeCost + service (legacy duration metric)
    double load;                     // Sum of deliveries only (legacy load)
    double reversalDistance;         // Difference of cost if the route is reversed
    double penalty;                  // Legacy penalty used in distance-pruning (kept)

    // Geometry (unchanged)
    double polarAngleBarycenter;
    CircleSector sector;

    // NEW: exact soft VRPSPDTW components for this route
    double distanceObj = 0.0;        // objective distance: d(0,first)+Σd(i,i+1)+d(last,0)
    double timeWarp    = 0.0;        // accumulated time-warp (lateness) on the route
    double onboardMax  = 0.0;        // SPDTW max onboard load over the route
    double routeCost   = 0.0;        // exact cost: unitCost*distanceObj + cap + warp + dispatch
};

struct Node
{
    bool   isDepot;
    int    cour;
    int    position;
    int    whenLastTestedRI;
    Node * next;
    Node * prev;
    Route* route;

    // Legacy cumulative metrics (kept for RI pruning)
    double cumulatedLoad;
    double cumulatedTime;
    double cumulatedReversalDistance;

    double deltaRemoval; // Precomputed in SWAP*
};

// Structure used in SWAP* to remember the three best insertion positions of a customer in a given route
struct ThreeBestInsert
{
    int   whenLastCalculated;
    double bestCost[3];
    Node * bestLocation[3];

    void compareAndAdd(double costInsert, Node * placeInsert)
    {
        if (costInsert >= bestCost[2]) return;
        else if (costInsert >= bestCost[1]) { bestCost[2] = costInsert; bestLocation[2] = placeInsert; }
        else if (costInsert >= bestCost[0])
        {
            bestCost[2] = bestCost[1]; bestLocation[2] = bestLocation[1];
            bestCost[1] = costInsert;  bestLocation[1] = placeInsert;
        }
        else
        {
            bestCost[2] = bestCost[1]; bestLocation[2] = bestLocation[1];
            bestCost[1] = bestCost[0]; bestLocation[1] = bestLocation[0];
            bestCost[0] = costInsert;  bestLocation[0] = placeInsert;
        }
    }

    void reset() { bestCost[0] = bestCost[1] = bestCost[2] = 1.e30; bestLocation[0] = bestLocation[1] = bestLocation[2] = NULL; }
    ThreeBestInsert() { reset(); }
};

// Structured used to keep track of the best SWAP* move
struct SwapStarElement
{
    double moveCost = 1.e30;
    Node * U = NULL;
    Node * bestPositionU = NULL;
    Node * V = NULL;
    Node * bestPositionV = NULL;
};

// Segment state for O(1) VRPSPDTW evaluation (time-warp + SPDTW capacity)
struct TWSPD_Seg {
    int    len = 0;
    int    first = 0;
    int    last = 0;

    double dist = 0.0;
    double Serv = 0.0;
    double dur  = 0.0;

    double E = 0.0, L = 1e30;
    double TW = 0.0;

    double Del = 0.0;
    double Pick = 0.0;
    double sumDelta = 0.0;
    double maxPref = 0.0;
    double minPref = 0.0;

    // NEW
    double totalDel = 0.0;  // for onboardMax = totalDel + maxDelta
};





// Main local search structure
class LocalSearch
{
private:
    double totalSolutionCost;
    Params & params;                         // Problem parameters
    bool searchCompleted;
    int  nbMoves;
    std::vector<int> orderNodes;
    std::vector<int> orderRoutes;
    std::set<int>    emptyRoutes;
    int  loopID;

    // Linked-list representation
    std::vector<Node> clients;
    std::vector<Node> depots;
    std::vector<Node> depotsEnd;
    std::vector<Route> routes;
    std::vector<std::vector<ThreeBestInsert>> bestInsertClient;
    // O(1) gate caches per route: prefix/suffix segments in tour order
    std::vector<std::vector<TWSPD_Seg>> prefSeg;
    std::vector<std::vector<TWSPD_Seg>> sufSeg;
    std::vector<std::vector<TWSPD_Seg>> revSufSeg;

    std::vector<std::vector<TWSPD_Seg>> rangeFwd;
    std::vector<std::vector<TWSPD_Seg>> rangeRev;
    // TEMP vars used in loops
    std::vector<int> tempSubseq; // Reusable temporary vector to avoid allocations
    Node * nodeU; Node * nodeX;
    Node * nodeV; Node * nodeY;
    Route * routeU; Route * routeV;
    int nodeUPrevIndex, nodeUIndex, nodeXIndex, nodeXNextIndex;
    int nodeVPrevIndex, nodeVIndex, nodeYIndex, nodeYNextIndex;
    double loadU, loadX, loadV, loadY;
    double serviceU, serviceX, serviceV, serviceY;
    double penaltyCapacityLS, penaltyDurationLS;
    bool intraRouteMove;

    // --- Legacy helpers (kept) ---
    void setLocalVariablesRouteU();
    void setLocalVariablesRouteV();
    inline double penaltyExcessDuration(double myDuration) { return std::max<double>(0., myDuration - params.durationLimit) * penaltyDurationLS; }
    inline double penaltyExcessLoad(double myLoad) { return std::max<double>(0., myLoad - params.vehicleCapacity) * penaltyCapacityLS; }

    // --- Exact VRPSPDTW helpers (NEW) ---
    void routeToSequence(Route* R, std::vector<int>& seq) const;     // builds [customers] for this route
    double costOfSequence(const std::vector<int>& seq) const;        // exact soft VRPSPDTW cost for a sequence
    bool exactImprove_move1(Node* U, Node* V) const;
    bool exactImprove_move2(Node* U, Node* V, Node* X) const;
    bool exactImprove_move3(Node* U, Node* V, Node* X) const;
    bool exactImprove_move4(Node* U, Node* V) const;
    bool exactImprove_move5(Node* U, Node* V, Node* X) const;
    bool exactImprove_move6(Node* U, Node* V, Node* X, Node* Y) const;
    bool exactImprove_move7(Node* U, Node* V) const;
    bool exactImprove_move8(Node* U, Node* V) const;
    bool exactImprove_move9(Node* U, Node* V) const;

    // Sequence editing utilities (NEW)
    static void eraseFirst(std::vector<int>& a, int val);
    static void insertAfterValue(std::vector<int>& a, int afterVal, int val);
    static void insertAfterIndex(std::vector<int>& a, int afterPos, int val);
    static int  findFirstPos(const std::vector<int>& a, int val);
    static void reverseSegment(std::vector<int>& a, int l, int r);   // inclusive l..r

    void debugSegmentVsIndividual();

    void debugSingleRoute(int routeIdx, const std::vector<int> &sequence);

    void validateMove1Debug(Node *nodeU, Node *nodeV);

    /* RELOCATE */
    bool move1();
    bool move2();
    bool move3();

    static void swapSingletonsSafe(Node *a, Node *b);

    /* SWAP */
    bool move4();
    bool move5();
    bool move6();

    /* 2-OPT, 2-OPT* */
    bool move7();
    bool move8();
    bool move9();

    /* SWAP* */
    bool swapStar();
    double getCheapestInsertSimultRemoval(Node * U, Node * V, Node *& bestPosition);
    void preprocessInsertions(Route * R1, Route * R2);

    /* STRUCTURE MUTATIONS */
    static void insertNode(Node * U, Node * V);
    static void swapNode(Node * U, Node * V);

    void updateNodePositions(Route *R);

    void updateRouteData(Route * myRoute);

    void printRouteLinkedList(Route *R) const;

    bool validateRouteLinkedList(Route *R) const;

    void validateAllRoutes(const char *where) const;

public:
    void run(Individual & indiv, double penaltyCapacityLS, double penaltyDurationLS);
    void loadIndividual(const Individual & indiv);
    void exportIndividual(Individual & indiv);
    LocalSearch(Params & params);
};

#endif
