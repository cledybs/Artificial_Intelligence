#ifndef HEADER_H 
#define HEADER_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <chrono>

#define nCities 20  // Number of Cities
#define pSize 100   // Population Size
#define mRate 0.01  // Mutation Rate
#define eSize 10    // Elite Size

// City struct
struct City {
    int name;
    int x;
    int y;

    double distance(const City& city) const {
        int xDis = std::abs(x - city.x);
        int yDis = std::abs(y - city.y);
        return std::sqrt((xDis * xDis) + (yDis * yDis));
    }

    bool operator==(const City& other) const {
        return name == other.name && x == other.x && y == other.y;
    }
};

// Fitness Function
class Fitness {
public:
    Fitness(const std::vector<City>& route) : route(route), distance(0.0), fitness(0.0) {}

    double routeDistance() {
        if (distance == 0.0) {
            double pathDistance = 0.0;
            for (size_t i = 0; i < route.size(); i++) {
                const City& fromCity = route[i];
                const City& toCity = (i + 1 < route.size()) ? route[i + 1] : route[0];
                pathDistance += fromCity.distance(toCity);
            }
            distance = pathDistance;
        }
        return distance;
    }

    double routeFitness() {
        if (fitness == 0.0) {
            fitness = 1.0 / routeDistance();
        }
        return fitness;
    }

private:
    std::vector<City> route;
    double distance;
    double fitness;
};

#endif
