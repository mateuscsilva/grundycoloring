#ifndef GRUNDYVERIFIER_H
#define GRUNDYVERIFIER_H

#include <iostream>
#include <algorithm>
#include "GrundyInstance.h"

bool grundySolutionVerifier(const GrundyInstance& instance, std::vector<int> nodesOrder, 
	std::vector<int> nodesColor);

bool verifyColor(const GrundyInstance& instance, int vertex, int color, std::vector<int> vertexColor);

bool verifyNeighborhood(const GrundyInstance& instance, int vertex, std::vector<int> inducedSubgraph);

#endif