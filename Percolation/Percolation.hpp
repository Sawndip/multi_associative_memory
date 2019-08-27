class Vertex;

namespace percolation {
	// Threshold is in the interval [0,1] and a value of t runs percolation on
	// the graph induced only by edges with presence < t.
	std::set<Vertex*> asynchronous_threshold(std::set<Vertex*> boostrap, int K,
			double threshold);
	// A function which returns a set with the active vertices after
	// percolation.
	std::set<Vertex*> asynchronous(std::set<Vertex*> boostrap, int K);

	std::set<Vertex*> asynchronous_timed_bootstrap(std::set<std::pair<double,Vertex*> > boostrap, int K);

	// Dynamic Percolation which returns which vertices are in the refractory period
	// The simulation stops after time T
	std::vector<std::pair<double, int> > dynamic_threshold(std::set<Vertex*> bootstrap, 
			int K, double T, int init_voltage, double threshold);
	std::vector<std::pair<double, int> > dynamic(std::set<Vertex*> bootstrap, int K,
			double T, int init_voltage);

	std::set<Vertex*> pattern_activation(std::set<Vertex*> bootstrap, int K_e, int K_i, int threshold, int pattern);
	std::vector<std::pair<double, int> > async_fixed_pattern_activation(std::set<Vertex*> bootstrap, 
			int K_e, int K_i, double threshold, int pattern);
	std::set<Vertex*> one_round_percolation(std::set<Vertex*> bootstrap, int K_e,
			int K_i, int threshold, int pattern);
	std::set<Vertex*> bayesian_pattern_activation(std::vector<Vertex>& sources, std::vector<Vertex>& targets, std::set<Vertex*> pattern, int patterns);
}
