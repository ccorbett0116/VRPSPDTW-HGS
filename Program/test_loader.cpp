#include <iostream>
#include "InstanceCVRPLIB.h"

int main() {
    try {
        InstanceCVRPLIB instance("./data/1000_1.vrpsdptw", false);

        std::cout << "Loaded instance with " << instance.nbClients << " clients.\n";
        std::cout << "Problem name     : " << instance.problemName << "\n";
        std::cout << "Problem type     : " << instance.problemType << "\n";
        std::cout << "Vehicle capacity : " << instance.vehicleCapacity << "\n";
        std::cout << "Num vehicles     : " << instance.numVehicles << "\n";
        std::cout << "Dispatching cost : " << instance.dispatchingCost << "\n";
        std::cout << "Unit cost        : " << instance.unitCost << "\n";

        std::cout << "\nFirst 5 nodes:\n";
        for (int i = 0; i < 5; ++i) {
            std::cout << "Node " << i
                << ": demand=" << instance.demands[i]
                << ", pickup=" << instance.pickups[i]
                << ", time window=[" << instance.ready_time[i]
                << "," << instance.due_time[i] << "]"
                << ", service=" << instance.service_time[i] << "\n";
        }

        std::cout << "\nSample distances (from node 0):\n";
        for (int j = 1; j <= 5; ++j) {
            std::cout << "Distance 0 to " << j << ": "
                << instance.dist_mtx[0][j]
                << " | Time: " << instance.time_mtx[0][j] << "\n";
        }
    }
    catch (const std::string& e) {
        std::cerr << "Error: " << e << "\n";
    }
}
