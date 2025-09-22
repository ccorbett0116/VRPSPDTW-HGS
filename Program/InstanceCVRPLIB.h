//
// Created by chkwon on 3/22/22.
//

#ifndef INSTANCECVRPLIB_H
#define INSTANCECVRPLIB_H
#include <string>
#include <vector>

class InstanceCVRPLIB
{
public:
	// Distance and time matrices
	std::vector<std::vector<double>> dist_mtx;
	std::vector<std::vector<double>> time_mtx;

	// Per-customer data
	std::vector<double> service_time;
	std::vector<double> demands;
	std::vector<double> pickups;
	std::vector<double> ready_time;
	std::vector<double> due_time;

	// Coordinates (for CDP instances)
	std::vector<double> x_coords;
	std::vector<double> y_coords;

	// Global parameters
	double durationLimit = 1.e30;          // Route duration limit
	double vehicleCapacity;                // Capacity limit
	bool isDurationConstraint = false;     // Whether problem includes duration/TW constraints
	bool areCoordinatesProvided = false;   // True only for CDP instances
	int nbClients;                         // Number of clients (excluding depot)
	std::string problemName;
	std::string problemType;
	int numVehicles = 0;
	double dispatchingCost = 0.0;
	double unitCost = 0.0;

	InstanceCVRPLIB(std::string pathToInstance, bool isRoundingInteger);
};

#endif // INSTANCECVRPLIB_H
