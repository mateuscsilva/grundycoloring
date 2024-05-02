#ifndef GRUNDYDECODER_H
#define GRUNDYDECODER_H

#include "GrundySolver.h"
#include "GrundyInstance.h"

class GrundyDecoder {
public:
	GrundyDecoder(const GrundyInstance& instance, std::string solverType);
	virtual ~GrundyDecoder();

	// Decodes a chromosome into a solution to the TSP:
	double decode(const std::vector< double >& chromosome) const;

private:
	const GrundyInstance& instance;
	const std::string solverType;
};

#endif