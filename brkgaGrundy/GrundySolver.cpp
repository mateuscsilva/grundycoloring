#include "GrundySolver.h"
#include <iostream>
#include <cmath>
#include <set>
#include <queue>  
#include <utility>
#include <string>
#include <algorithm>

GrundySolver::GrundySolver(const GrundyInstance& instance, const std::vector< double >& chromosome,
	std::string solverTye) : color(0), nodesColor(instance.getNumNodes()) {
	
	if(solverTye.compare("grundyNP") == 0){
		color = GrundySolverNormalPriority(instance, chromosome, nodesColor);
	}
}

GrundySolver::~GrundySolver() {
}

int GrundySolver::GrundySolverNormalPriority(const GrundyInstance& instance, 
	const std::vector< double >& chromosome, std::vector<int>& nodesColor){
	
	int color = 0;
	std::pair <double, int> pairPriorityVertex;
	std::vector <int> colorNodes;
	std::priority_queue<std::pair<double, int> > notSelectedNodes;
	
	colorNodes.resize(instance.getNumNodes());
	

	for(int i = 0; i < instance.getNumNodes(); i++){ 
		std::pair <double, int> pairChromPriority;
		pairChromPriority = std::make_pair(chromosome[i], i);
		notSelectedNodes.push(pairChromPriority); 
	}
	
	while(!notSelectedNodes.empty()) {
		pairPriorityVertex = notSelectedNodes.top();
		int vertex = pairPriorityVertex.second;
		colorNodes[vertex] = chooseNodeColor(instance, colorNodes, vertex);
		notSelectedNodes.pop();	
		nodesOrder.push_back(vertex);

		if(colorNodes[vertex] > color){
			color = colorNodes[vertex];
		}
	}
	nodesColor = colorNodes;
	return color;
}

unsigned GrundySolver::getHighestColor() const { return color; }

int GrundySolver::chooseNodeColor(const GrundyInstance& instance, const std::vector< int >& colorNodes,
 	int vertex) const {
	
	int actualColor = 0;
	bool isChoosedColor = false;
	while(!isChoosedColor){
		actualColor++;
		isChoosedColor = true;
		for(int i=0; i < instance.nodes[vertex].size(); i++){
			if(actualColor == colorNodes[instance.nodes[vertex][i]]){
				isChoosedColor = false;
				break;
			}
		}
	}
	return actualColor; 
}

std::vector<int> GrundySolver::getNodesColor() const{
	return nodesColor;
}

std::vector<int> GrundySolver::getNodesOrder() const{
	return nodesOrder;
}