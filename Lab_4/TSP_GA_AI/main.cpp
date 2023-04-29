#include "header.h"

// A random solution/path
std::vector<City> createRoute(const std::vector<City>& cityList) {
    std::vector<City> route = cityList;
    std::random_shuffle(route.begin(), route.end());
    return route;
}

// Create an initial population
std::vector<std::vector<City>> initialPopulation(int popSize, const std::vector<City>& cityList) {
    std::vector<std::vector<City>> population;
    for (int i = 0; i < popSize; i++) {
        population.push_back(createRoute(cityList));
    }
    return population;
}

// Rank individuals
std::vector<std::pair<int, double>> rankRoutes(const std::vector<std::vector<City>>& population) { // CORREGIRRR
    std::vector<std::pair<int, double>> fitnessResults;
    for (int i = 0; i < population.size(); i++) {
        Fitness fitness(population[i]);
        fitnessResults.emplace_back(i, fitness.routeFitness());
    }
    std::sort(fitnessResults.begin(), fitnessResults.end(),
        [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
            return a.second > b.second;
        });
    return fitnessResults;
}

// Create a selection function
std::vector<int> selection(const std::vector<std::pair<int, double>>& popRanked, int eliteSize) {
    std::vector<int> selectionResults;
    std::vector<double> cumulativeFitness;
    double sumFitness = 0.0;
    for (auto i = popRanked.begin(); i != popRanked.end(); i++) {
        sumFitness += i->second;
        cumulativeFitness.push_back(sumFitness);
    }

    std::vector<double> cumulativePercentages;
    for (auto i = cumulativeFitness.begin(); i != cumulativeFitness.end(); i++) {
        cumulativePercentages.push_back(100*(*i)/cumulativeFitness.back());
    }

    for (int i = 0; i < eliteSize; i++) 
        selectionResults.push_back(popRanked[i].first);
    

    for (int i = 0; i < popRanked.size() - eliteSize; i++) {
        double pick = std::rand() % 100;
        for (int j = 0; j < popRanked.size(); j++) {
            if (pick <= cumulativePercentages[j]) {
                selectionResults.push_back(popRanked[j].first);
                break;
            }
        }
    }
    return selectionResults;
}

// Create a mating pool
std::vector<std::vector<City>> matingPool(const std::vector<std::vector<City>>& population, 
                                            const std::vector<int>& selectionResults) {
    std::vector<std::vector<City>> matingpool;
    for (int i = 0; i < selectionResults.size(); i++) {
        int index = selectionResults[i];
        matingpool.push_back(population[index]);
    }
    return matingpool;
}

// Create a crossover function
std::vector<City> breed(const std::vector<City>& parent1, const std::vector<City>& parent2) {
    std::vector<City> child(parent1.size());
    int geneA = std::rand() % parent1.size();
    int geneB = std::rand() % parent1.size();

    int startGene = std::min(geneA, geneB);
    int endGene = std::max(geneA, geneB);

    for (int i = startGene; i < endGene; i++) {
        child[i] = parent1[i];
    }

    int childIndex = endGene;
    for (int i = 0; i < parent2.size(); i++) {
        if (std::find(child.begin(), child.end(), parent2[i]) == child.end()) { // CORREGIRRR
            child[childIndex] = parent2[i];
            childIndex = (childIndex + 1) % child.size();
        }
    }

    return child;
}

// Create a function to run crossover over the full mating pool
std::vector<std::vector<City>> breedPopulation(std::vector<std::vector<City>>& matingpool, 
                                                int eliteSize) {
    std::vector<std::vector<City>> children;
    int length = matingpool.size() - eliteSize;
    std::vector<std::vector<City>> pool = matingpool;

    for (int i = 0; i < eliteSize; i++) 
        children.push_back(matingpool[i]);

    for (int i = 0; i < length; i++) {
        int randomIndex1 = std::rand() % pool.size();
        int randomIndex2 = std::rand() % pool.size();
        std::vector<City> child = breed(pool[randomIndex1], pool[randomIndex2]);
        children.push_back(child);
    }

    return children;
}

// Create a function to mutate a single route
std::vector<City> mutate(std::vector<City> individual, double mutationRate) {
    for (int i = 0; i < individual.size(); i++) {
        if (static_cast<double>(std::rand()) / RAND_MAX < mutationRate) {
            int swapWith = std::rand() % individual.size();
            std::swap(individual[i], individual[swapWith]);
        }
    }
    return individual;
}

// Create a function to run mutation over the entire population
std::vector<std::vector<City>> mutatePopulation(const std::vector<std::vector<City>>& population, 
                                                double mutationRate) {
    std::vector<std::vector<City>> mutatedPop;
    for (int i = 0; i < population.size(); i++) {
        std::vector<City> mutatedInd = mutate(population[i], mutationRate);
        mutatedPop.push_back(mutatedInd);
    }
    return mutatedPop;
}

// Put all steps together to create the next generation
std::vector<std::vector<City>> nextGeneration(const std::vector<std::vector<City>>& currentGen, 
                                                int eliteSize, double mutationRate) {
    std::vector<std::pair<int, double>> popRanked = rankRoutes(currentGen);
    std::vector<int> selectionResults = selection(popRanked, eliteSize);
    std::vector<std::vector<City>> matingpool = matingPool(currentGen, selectionResults);
    std::vector<std::vector<City>> children = breedPopulation(matingpool, eliteSize);
    std::vector<std::vector<City>> nextGeneration = mutatePopulation(children, mutationRate);
    return nextGeneration;
}

// Compute the Average Fitness of the Population
double averagePop(const std::vector<std::pair<int, double>>& population) {
    double sum = 0.0;
    for (auto i = population.begin(); i != population.end(); i++)
        sum += i->second;
    
    return sum / population.size();
}

// Create the genetic algorithm
std::vector<City> geneticAlgorithm(const std::vector<City>& population, int popSize, 
                                    int eliteSize, double mutationRate, int generations) {
    std::vector<std::vector<City>> pop = initialPopulation(popSize, population);
    std::vector<double> progressAv;
    std::vector<double> progress;
    progressAv.push_back(1.0 / averagePop(rankRoutes(pop)));
    progress.push_back(1.0 / rankRoutes(pop)[0].second);
    std::cout << "Generation [0] Average Fitness:\t" << progressAv[0] << std::endl;
    std::cout << "Generation [0] Best Fitness:\t" << progress[0] << std::endl;

    for (int i = 1; i <= generations; i++) {
        pop = nextGeneration(pop, eliteSize, mutationRate);
        progressAv.push_back(1.0 / averagePop(rankRoutes(pop)));
        progress.push_back(1.0 / rankRoutes(pop)[0].second);
        if (i % 20 == 0) {
            std::cout << "Generation [" << i << "] Average Fitness: " << progressAv[i] << std::endl;
            std::cout << "Generation [" << i << "] Best Fitness: " << progress[i] << std::endl;
        }
    }

    int bestRouteIndex = rankRoutes(pop)[0].first;
    std::vector<City> bestRoute = pop[bestRouteIndex];

    // Here show the Average Fitness Vs Generation 

    // Here show the Best Fitness Vs Generation

    return bestRoute;
}

// Main
int main() {
    std::vector<City> cityList;
    std::srand(std::time(nullptr));

    for (int i = 0; i < nCities; i++) {
        City city;
        city.name = i;
        city.x = std::rand() % 200;
        city.y = std::rand() % 200;
        cityList.push_back(city);
    }

    std::vector<City> bestRoute = geneticAlgorithm(cityList, nCities, eSize, mRate, pSize);

    for (const auto& city : bestRoute) {
        std::cout << "(" << city.name << ") ";
    }
    std::cout << std::endl;

    // Here show the cities and the best path manos

    return 0;
}