#include "GrundyDecoder.h"

GrundyDecoder::GrundyDecoder(const GrundyInstance& _instance, std::string type) : 
	instance(_instance), solverType(type) {
}

GrundyDecoder::~GrundyDecoder() {
}

double GrundyDecoder::decode(const std::vector< double >& chromosome) const {
	// 1) Solve the problem (i.e., create a tour out of this chromosome):
	// Avoids race conditions by making sure we have a single TSPSolver for each thread calling
	// ::decode (as long as TSPSolver does not make use of 'static' or other gimmicks):
	GrundySolver solver(instance, chromosome, solverType);

	// 2) Extract the fitness (tour distance):
	const unsigned fitness = solver.getHighestColor();

	// 3) Return:
	return double(fitness);
}