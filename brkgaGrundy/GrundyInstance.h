#ifndef GRUNDYINSTANCE_H
#define GRUNDYINSTANCE_H

#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

class GrundyInstance {
public:
	typedef std::runtime_error Error;

	GrundyInstance(const std::string& instanceFile) throw(Error);
	virtual ~GrundyInstance();

	std::vector< std::vector<int> > nodes;
	std::vector<int> colorNodes;
	std::vector <std::vector<double> > initialSolutions;
	std::string instanceFileName;
	std::string solverType;
	int executionId;
	int globalBestFitness = 0;
	// Getters:
	unsigned getNumNodes() const;
	void writeOutput(const std::string outputFile, const std::string instanceFile, 
		const int maxColor, const std::string solverType, bool validSolution, double executionTime,
		double pe, double pm, double rhoe, int generation, unsigned relGeneration, double relGenerationTime);
	void readInitialSolutions(const std::string& solFile) throw(Error);

	int getMaxDegree() const;
	int getNumberOfNeighbors(int vertex) const;
	void setGlobalBestFitness(int globalBestFitness);
	
private:

	std::string name;
	std::string comment;
	std::string problemType;
	int maxDegree;

	unsigned nNodes;

	void readNodes(const std::string& line, int simpleInstancePattern) throw(Error);
	bool isEOF(const std::string& line) const;
	void trim(std::string& str) const;
	std::vector<double> convertToChromossome(std::vector<int> nodesOrder);
};

#endif