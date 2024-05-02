#include <iostream>
#include <algorithm>
#include <queue>
#include "GrundyInstance.h"
#include "GrundyVerifier.h"

bool grundySolutionVerifier(const GrundyInstance& instance, std::vector<int> nodesOrder, 
	std::vector<int> nodesColor){
	
	bool validSolution = true;
	int maxColor = -1;
	std::vector<int> vertexColor;
	vertexColor.resize(instance.getNumNodes());

	if(nodesOrder.size() != instance.getNumNodes()) { validSolution = false; }

	for(int i = 0; i < nodesOrder.size(); i++){
		validSolution = verifyColor(instance, nodesOrder[i], nodesColor[nodesOrder[i]], vertexColor);
		if(validSolution){
			vertexColor[nodesOrder[i]] = nodesColor[nodesOrder[i]];
		}else{
			break;
		}
	}
	return validSolution;
}

bool verifyColor(const GrundyInstance& instance, int vertex, int color, std::vector<int> vertexColor){
	int minColor = 0;
	bool isChoosedColor = false;
	
	while(!isChoosedColor){
		minColor++;
		isChoosedColor = true;
		for(int i=0; i < instance.nodes[vertex].size(); i++){
			if(minColor == vertexColor[instance.nodes[vertex][i]]){
				isChoosedColor = false;
				break;
			}
		}
	}

	bool correctColor = (minColor == color) ? true : false;
	return correctColor;
}

bool verifyNeighborhood(const GrundyInstance& instance, int vertex, std::vector<int> inducedSubgraph){
	bool isNeighbour = false;

	for(int i=0; i < instance.nodes[vertex].size(); i++){
		if(std::find(inducedSubgraph.begin(), inducedSubgraph.end(), instance.nodes[vertex][i]) != 
			inducedSubgraph.end()){
			isNeighbour = true;
			break;
		}
	}
	return isNeighbour;
}