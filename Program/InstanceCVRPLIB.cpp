// Created by chkwon on 3/22/22.
// Extended to support CDP instances with coordinates (Euclidean distance)

#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "InstanceCVRPLIB.h"

InstanceCVRPLIB::InstanceCVRPLIB(std::string pathToInstance, bool isRoundingInteger)
{
    std::ifstream inputFile(pathToInstance);
    if (!inputFile.is_open())
        throw std::string("Cannot open instance file: " + pathToInstance);

    std::string line;
    bool isCDP = false;

    // Peek first line to detect CDP
    std::streampos startPos = inputFile.tellg();
    if (std::getline(inputFile, line)) {
        if (line.find("Cdp") == 0 || line.find("CDP") == 0) {
            isCDP = true;
        }
    }
    inputFile.clear();
    inputFile.seekg(startPos);

    if (isCDP) {
        // ---------------- CDP PARSER ----------------
        // First line: name
        std::getline(inputFile, problemName);

        // Skip header lines until capacity block
        while (std::getline(inputFile, line)) {
            if (line.find("NUMBER") != std::string::npos && line.find("CAPACITY") != std::string::npos)
                break;
        }

        // Parse nbClients, numVehicles, capacity
        int numCustomers = 0;
        inputFile >> numCustomers >> numVehicles >> vehicleCapacity;
        std::getline(inputFile, line); // flush newline
        nbClients = numCustomers-1;

        // Skip to CUSTOMER section
        while (std::getline(inputFile, line)) {
            if (line.find("CUSTOMER") != std::string::npos && line.find("CUST") != std::string::npos)
                break;
        }
        // Skip header row
        std::getline(inputFile, line);

        // Allocate storage
        demands.resize(nbClients + 1);
        pickups.resize(nbClients + 1);
        ready_time.resize(nbClients + 1);
        due_time.resize(nbClients + 1);
        service_time.resize(nbClients + 1);
        x_coords.resize(nbClients + 1);
        y_coords.resize(nbClients + 1);

        // Read customers (including depot 0)
        for (int i = 0; i <= nbClients; i++) {
            int id;
            double x, y, del, pick, ready, due, service;
            inputFile >> id >> x >> y >> del >> pick >> ready >> due >> service;
            if (id != i)
                throw std::runtime_error("Node index mismatch in CDP file at node " + std::to_string(id));

            x_coords[i] = x;
            y_coords[i] = y;
            demands[i] = del;
            pickups[i] = pick;
            ready_time[i] = ready;
            due_time[i] = due;
            service_time[i] = service;
        }

        // Build distance/time matrices (Euclidean distance = time)
        dist_mtx.resize(nbClients + 1, std::vector<double>(nbClients + 1, 0.0));
        time_mtx.resize(nbClients + 1, std::vector<double>(nbClients + 1, 0.0));
        for (int i = 0; i <= nbClients; i++) {
            for (int j = 0; j <= nbClients; j++) {
                if (i == j) {
                    dist_mtx[i][j] = 0.0;
                    time_mtx[i][j] = 0.0;
                } else {
                    double dx = x_coords[i] - x_coords[j];
                    double dy = y_coords[i] - y_coords[j];
                    double dist = std::sqrt(dx*dx + dy*dy);
                    dist_mtx[i][j] = dist;
                    time_mtx[i][j] = dist;
                }
            }
        }


        // CDP assumptions
        dispatchingCost = 2000.0;
        unitCost = 1.0;
        areCoordinatesProvided = true;
        isDurationConstraint = true; // CDP always has TWs
        return;
    }

    // ---------------- VRPSPDTW PARSER (original) ----------------
    while (std::getline(inputFile, line)) {
        if (line.empty()) continue;

        if (line.find("NAME") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                problemName = line.substr(colonPos + 1);
                problemName.erase(0, problemName.find_first_not_of(" \t"));
            }
        }
        else if (line.find("TYPE") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                problemType = line.substr(colonPos + 1);
                problemType.erase(0, problemType.find_first_not_of(" \t"));
            }
        }
        else if (line.find("DIMENSION") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                nbClients = std::stoi(line.substr(colonPos + 1)) - 1;
            }
        }
        else if (line.find("VEHICLES") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                numVehicles = std::stoi(line.substr(colonPos + 1));
            }
        }
        else if (line.find("DISPATCHINGCOST") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                dispatchingCost = std::stod(line.substr(colonPos + 1));
            }
        }
        else if (line.find("UNITCOST") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                unitCost = std::stod(line.substr(colonPos + 1));
            }
        }
        else if (line.find("CAPACITY") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                vehicleCapacity = std::stod(line.substr(colonPos + 1));
            }
        }
        else if (line.find("EDGE_WEIGHT_TYPE") == 0) {
            size_t colonPos = line.find(":");
            if (colonPos != std::string::npos) {
                std::string type = line.substr(colonPos + 1);
                type.erase(0, type.find_first_not_of(" \t"));
                if (type != "EXPLICIT") {
                    throw std::runtime_error("Expected EDGE_WEIGHT_TYPE: EXPLICIT");
                }
            }
        }
        else if (line == "NODE_SECTION") {
            break;
        }
    }

    // Allocate storage
    demands.resize(nbClients + 1);
    service_time.resize(nbClients + 1);
    pickups.resize(nbClients + 1);
    ready_time.resize(nbClients + 1);
    due_time.resize(nbClients + 1);

    // Parse NODE_SECTION lines: i,delivery,pickup,start,end,service
    for (int i = 0; i <= nbClients; ++i) {
        std::getline(inputFile, line);
        std::stringstream ss(line);
        std::string field;
        int id;
        double delivery, pickup, ready, due, service;

        std::getline(ss, field, ','); id = std::stoi(field);
        std::getline(ss, field, ','); delivery = std::stod(field);
        std::getline(ss, field, ','); pickup   = std::stod(field);
        std::getline(ss, field, ','); ready    = std::stod(field);
        std::getline(ss, field, ','); due      = std::stod(field);
        std::getline(ss, field, ','); service  = std::stod(field);

        if (id != i) throw std::string("Node index mismatch at line: " + line);

        demands[i] = delivery;
        pickups[i] = pickup;
        ready_time[i] = ready;
        due_time[i] = due;
        service_time[i] = service;
    }

    // Allocate distance/time matrices
    dist_mtx.resize(nbClients + 1, std::vector<double>(nbClients + 1, 0.0));
    time_mtx.resize(nbClients + 1, std::vector<double>(nbClients + 1, 0.0));

    // Look for DISTANCETIME_SECTION
    while (std::getline(inputFile, line)) {
        if (line == "DISTANCETIME_SECTION")
            break;
    }

    while (std::getline(inputFile, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string fromStr, toStr, distStr, timeStr;

        if (!std::getline(ss, fromStr, ',')) continue;
        if (!std::getline(ss, toStr, ',')) continue;
        if (!std::getline(ss, distStr, ',')) continue;
        if (!std::getline(ss, timeStr, ',')) continue;

        try {
            int from = std::stoi(fromStr);
            int to = std::stoi(toStr);
            double dist = std::stod(distStr);
            double time = std::stod(timeStr);
            dist_mtx[from][to] = dist;
            time_mtx[from][to] = time;
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "Invalid line in DISTANCETIME_SECTION: " << line << "\n";
            continue;
        }
    }

    // Detect duration/TW constraints
    bool hasTW = false, hasService = false, hasTime = false;

    for (int i = 0; i <= nbClients; ++i) {
        if (ready_time[i] > 0.0 || due_time[i] < 1.e29) hasTW = true;
        if (service_time[i] > 0.0)                       hasService = true;
    }

    for (int a = 0; a <= nbClients && !hasTime; ++a) {
        for (int b = 0; b <= nbClients; ++b) {
            if (time_mtx[a][b] > 0.0) { hasTime = true; break; }
        }
    }

    isDurationConstraint = (hasTW || hasService || hasTime);
    areCoordinatesProvided = false; // VRP files don't provide coords
}
