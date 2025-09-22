#include "Individual.h"

void Individual::evaluateCompleteCost(const Params & params) {
    eval = EvalIndiv();

    for (int r = 0; r < params.nbVehicles; r++) {
        if (chromR[r].empty()) continue;

        double currentTime = 0.0;
        double routeDistance = 0.0;

        // ---- Exact SPDTW capacity tracking ----
        double totalDel = 0.0;    // total deliveries on this route
        for (int u : chromR[r]) totalDel += params.cli[u].demand;

        double prefDel = 0.0, prefPick = 0.0, maxDelta = 0.0;
        // ---------------------------------------

        for (int i = 0; i < (int)chromR[r].size(); i++) {
            int node = chromR[r][i];

            // Travel (objective vs. schedule)
            if (i == 0) {
                routeDistance += params.timeCost[0][node];
                currentTime   += params.timeMatrix[0][node];
            } else {
                int prev = chromR[r][i-1];
                routeDistance += params.timeCost[prev][node];
                currentTime   += params.timeMatrix[prev][node];
            }

            // Time windows with soft warp and depot-return feasibility at each node
            const double earliest = params.readyTime[node];
            const double latest   = std::min(
                params.dueTime[node],
                params.durationLimit - params.serviceTime[node] - params.timeMatrix[node][0]
            );
            if (currentTime < earliest) currentTime = earliest;
            else if (currentTime > latest) { eval.timeWarp += (currentTime - latest); currentTime = latest; }

            // Service time
            currentTime += params.serviceTime[node];

            // ---- SPDTW prefixes for capacity ----
            prefDel  += params.cli[node].demand;
            prefPick += params.pickups[node];       // pickups must be in Params
            maxDelta  = std::max(maxDelta, prefPick - prefDel);
            // -------------------------------------
        }

        // Return to depot (objective + schedule)
        routeDistance += params.timeCost[chromR[r].back()][0];
        currentTime   += params.timeMatrix[chromR[r].back()][0];

        // Optional depot TW check (only if you model a closing time at depot)
        if (currentTime > params.dueTime[0]) eval.timeWarp += (currentTime - params.dueTime[0]);

        // Duration excess (route-level)
        if (currentTime > params.durationLimit) eval.durationExcess += (currentTime - params.durationLimit);

        // ---- Compute SPDTW capacity excess using onboardMax ----
        const double onboardMax = totalDel + std::max(0.0, maxDelta);
        if (onboardMax > params.vehicleCapacity)
            eval.capacityExcess += (onboardMax - params.vehicleCapacity);
        // --------------------------------------------------------

        // Totals
        eval.distance += routeDistance;
        eval.nbRoutes++;
    }

	eval.penalizedCost =
	  params.unitCost        * eval.distance      // scale objective distance
	+ params.dispatchingCost * eval.nbRoutes      // fixed cost per non-empty route
	+ params.penaltyCapacity * eval.capacityExcess
	+ params.penaltyDuration * eval.durationExcess
	+ params.penaltyMultiplier * eval.timeWarp;

	eval.isFeasible =
		(eval.capacityExcess < MY_EPSILON) &&
		(!params.isDurationConstraint || eval.durationExcess < MY_EPSILON) &&
		(eval.timeWarp < MY_EPSILON) &&
		(eval.nbRoutes <= params.nbVehicles);
}


Individual::Individual(Params & params)
{
	successors = std::vector <int>(params.nbClients + 1);
	predecessors = std::vector <int>(params.nbClients + 1);
	chromR = std::vector < std::vector <int> >(params.nbVehicles);
	chromT = std::vector <int>(params.nbClients);
	for (int i = 0; i < params.nbClients; i++) chromT[i] = i + 1;
	std::shuffle(chromT.begin(), chromT.end(), params.ran);
	eval.penalizedCost = 1.e30;
}

Individual::Individual(Params & params, std::string fileName) : Individual(params)
{
	double readCost;
	chromT.clear();
	std::ifstream inputFile(fileName);
	if (inputFile.is_open())
	{
		std::string inputString;
		inputFile >> inputString;
		// Loops in the input file as long as the first line keyword is "Route"
		for (int r = 0; inputString == "Route"; r++)
		{
			inputFile >> inputString;
			getline(inputFile, inputString);
			std::stringstream ss(inputString);
			int inputCustomer;
			while (ss >> inputCustomer) // Loops as long as there is an integer to read in this route
			{
				chromT.push_back(inputCustomer);
				chromR[r].push_back(inputCustomer);
			}
			inputFile >> inputString;
		}
		if (inputString == "Cost") inputFile >> readCost;
		else throw std::string("Unexpected token in input solution");

		// Some safety checks and printouts
		evaluateCompleteCost(params);
		if ((int)chromT.size() != params.nbClients) throw std::string("Input solution does not contain the correct number of clients");
		if (!eval.isFeasible) throw std::string("Input solution is infeasible");
		if (eval.penalizedCost != readCost)throw std::string("Input solution has a different cost than announced in the file");
		if (params.verbose) std::cout << "----- INPUT SOLUTION HAS BEEN SUCCESSFULLY READ WITH COST " << eval.penalizedCost << std::endl;
	}
	else
		throw std::string("Impossible to open solution file provided in input in : " + fileName);
}
