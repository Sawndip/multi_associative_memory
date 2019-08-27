#include <random>
#include <iostream>
#include <cmath>

#include "Utils.hpp"

namespace Utils {
	std::set<Vertex*> sample_bootstrap(vector<Vertex>& vertices, double p) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::binomial_distribution<> d(vertices.size(), p);

		return sample_bootstrap(vertices, d(gen));
	}

	void batch_process(set<Vertex*>& sources, set<Vertex*>& targets,
			void (*process_edge)(Edge&)) {
		random_device rd;
		mt19937 gen(rd());
		typename set<Vertex*>::iterator it1;
		for (it1 = sources.begin(); it1 != sources.end(); ++it1) {
			vector<Edge>::iterator it2;
			for (it2 = (*it1)->neighbors.begin(); it2 != (*it1)->neighbors.end();
					++it2) {
				if (targets.count((*it2).target) > 0) {
					process_edge((*it2));
				}
			}
		}
	}

	/**
	 * @param	vertices	vector of vertices from which we count
	 * @param	K			activation threshold
	 */
	int count_active(set<Vertex*> vertices, int K) {
		set<Vertex*>::iterator it;
		int num_active = 0;
		for (it = vertices.begin(); it != vertices.end(); it++) {
			if ((*it)->get_voltage() >= K && (*it)->get_type() == 1) {
				num_active++;
			}
		}

		return num_active;
	}

	int count_active_outside(set<Vertex*> pattern, set<Vertex*> vertices, int K) {
		set<Vertex*>::iterator it;
		int num_active = 0;
		for (it = vertices.begin(); it != vertices.end(); it++) {
			if (pattern.count(*it) == 0) {	
				if ((*it)->get_voltage() >= K && (*it)->get_type() == 1) {
					num_active++;
				}
			}
		}

		return num_active;
	}

	int intersection_size(vector<Vertex>& vertices, set<Vertex*>& active) {
		int count = 0;
		for (int i = 0; i < (int) vertices.size(); i++) {
			if (active.count(&vertices[i]) > 0) {
				count++;
			}
		}

		return count;
	}

	double squared_euclidean_distance_2D(double x1, double y1, double x2, double y2) {
		return pow(abs(x1 - x2), 2.0)
				+ pow(abs(y1 - y2), 2.0); // (x1 - x2)^2 + (y1 - y2)^2
	}

	double squared_euclidean_distance_3D(double x1, double y1, double z1, double x2, double y2, double z2) {
		return pow(abs(x1 - x2), 2.0)
				+ pow(abs(y1 - y2), 2.0) + pow(abs(z1 - z2), 2.0); // (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
	}

	double entropy(double x) {
		if (x == 0 || x == 1) {
			return 0;
		}

		return -x * log2(x) - (1 - x) * log2(1 - x);
	}

	double calc_info(int patterns, int active, int active_pattern, int n, int k, double p) { // p is p_aff + p_rec

		double p10 = (double)(k - active_pattern) / (double)k;
		double p01 = (double)(active - active_pattern) / (double)(n - k);
		
		double I = -(k * entropy(p10) + (n - k) * entropy(p01)) + n * entropy((double)k / (double)n * (1 - p10) + (double)(n - k) / (double)n * p01);

		return (double)(I * patterns) / (double) (pow(n, 2.0) * p);
	}
}