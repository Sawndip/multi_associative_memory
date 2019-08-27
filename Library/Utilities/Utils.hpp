#include <vector>
#include <set>
#include <fstream>
#include <iostream>

#include "../Edges/Edge.h"
#include "../Vertices/Geo_vertex.h"
#include "../Graph_Constructs/Layered_graph.h"
#include "Connect.hpp"
#include "Patterns.hpp"
#include "Print.hpp"
#include "Reset.hpp"

class Vertex;

namespace Utils {
	// Random bootstrap with n vertices
	template<typename T>
	std::set<Vertex*> sample_bootstrap(vector<T>& vertices, int n) {
		if (n > vertices.size()) {
			clog << "You cannot have a bootstrap larger than your vertex set.\n";
		} else if (n == 0) {
			return std::set<Vertex*>();
		}
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(0, vertices.size() - 1);
		std::set<int> indices;
		if (n < 0.5 * vertices.size()) {
			while (indices.size() < n) {
				indices.insert(dis(gen));
			}
		} else {
			set<int> not_indices;
			while (not_indices.size() < vertices.size() - n) {
				not_indices.insert(dis(gen));
			}
			for (int i = 0; i < (int) vertices.size(); ++i) {
				if (not_indices.count(i) == 0) {
					indices.insert(i);
				}
			}
		}
		std::set<Vertex*> bootstrap;
		std::set<int>::iterator it;
		for (it = indices.begin(); it != indices.end(); ++it) {
			bootstrap.insert(&vertices[(*it)]);
		}
		
		return bootstrap;
	}

	template<typename T>
	std::set<Vertex*> sample_bootstrap(set<T>& vertices, int n) {
		if (n > vertices.size()) {
			clog << "You cannot have a bootstrap larger than your vertex set.\n";
		} else if (n == 0) {
			return std::set<Vertex*>();
		}
		std::vector<Vertex*> temp(vertices.begin(), vertices.end());

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(0, vertices.size() - 1);
		std::set<int> indices;
		if (n < 0.5 * vertices.size()) {
			while (indices.size() < n) {
				indices.insert(dis(gen));
			}
		} else {
			set<int> not_indices;
			while (not_indices.size() < vertices.size() - n) {
				not_indices.insert(dis(gen));
			}
			for (int i = 0; i < (int) vertices.size(); ++i) {
				if (not_indices.count(i) == 0) {
					indices.insert(i);
				}
			}
		}
		std::set<Vertex*> bootstrap;
		std::set<int>::iterator it;
		for (it = indices.begin(); it != indices.end(); ++it) {
			bootstrap.insert(temp[(*it)]);
		}
		
		return bootstrap;
	}

	template<typename T>
	int intersection_size(set<T>& A, set<T>& B) {
		int count = 0;
		typename set<T>::iterator it;
		for (it = A.begin(); it != A.end(); it++) {
			if (B.count(*it) > 0) {
				count++;
			}
		}

		return count;	
	}

	// Random bootstrap where each vertex is chosen with probability p
	std::set<Vertex*> sample_bootstrap(std::vector<Vertex>& vertices, double p);

	// Batch process
	void batch_process(std::set<Vertex*>& sources, std::set<Vertex*>& targets,
			void (*process_edge)(Edge&));

	// Counts the active vertices
	int count_active(set<Vertex*> vertices, int K);
	// Counts the active vertices outside the given pattern
	int count_active_outside(set<Vertex*> pattern, set<Vertex*> vertices, int K);

	int intersection_size(vector<Vertex>& vertices, set<Vertex*>& active);

	double squared_euclidean_distance_2D(double x1, double y1, double x2, double y2);
	double squared_euclidean_distance_3D(double x1, double y1, double z1, double x2, double y2, double z2);

	// Binary entropy function
	double entropy(double x);

	// Function to calculate information
	double calc_info(int patterns, int active, int active_pattern, int n, int k, double p);
}
