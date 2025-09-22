#include "Genetic.h"
#include "commandline.h"
#include "LocalSearch.h"
#include "Split.h"
#include "InstanceCVRPLIB.h"
using namespace std;

int main(int argc, char *argv[])
{
	try
	{
		// Reading the arguments of the program
		CommandLine commandline(argc, argv);

		// Print all algorithm parameter values
		if (commandline.verbose) print_algorithm_parameters(commandline.ap);

		// Reading the data file and initializing some data structures
		if (commandline.verbose) std::cout << "----- READING INSTANCE: " << commandline.pathInstance << std::endl;
		InstanceCVRPLIB cvrp(commandline.pathInstance, commandline.isRoundingInteger);

		Params params(
			cvrp.x_coords,
			cvrp.y_coords,
			cvrp.dist_mtx,
			cvrp.time_mtx,       // Travel time matrix
			cvrp.service_time,
			cvrp.demands,
			cvrp.pickups,
			cvrp.vehicleCapacity,
			cvrp.durationLimit,
			cvrp.ready_time,     // Lower bounds (ready_time)
			cvrp.due_time,       // Upper bounds (due_time)
			cvrp.dispatchingCost,
			cvrp.unitCost,
			commandline.nbVeh, //cvrp.numVehicles, <-- This is the fix but dont apply until the LS is O(1 or m)
			cvrp.isDurationConstraint,
			commandline.verbose,
			commandline.ap
		);

		// Running HGS
		Genetic solver(params);
		solver.run();
		
		// Exporting the best solution
		if (solver.population.getBestFound() != NULL)
		{
			if (params.verbose) std::cout << "----- WRITING BEST SOLUTION IN : " << commandline.pathSolution << std::endl;
			solver.population.exportCVRPLibFormat(*solver.population.getBestFound(),commandline.pathSolution);
			solver.population.exportSearchProgress(commandline.pathSolution + ".PG.csv", commandline.pathInstance);
		}
	}
	catch (const string& e) { std::cout << "EXCEPTION | " << e << std::endl; }
	catch (const std::exception& e) { std::cout << "EXCEPTION | " << e.what() << std::endl; }
	return 0;
}
