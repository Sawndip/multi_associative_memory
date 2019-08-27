/*
 * Connect.hpp
 *
 *  Created on: May 11, 2015
 *      Author: Trujic
 */

#include <random>

#include "Connect.hpp"
#include "Utils.hpp"

namespace Connect {
	void connect(vector<Vertex>& sources, vector<Vertex>& targets, double density) {
		random_device rd;
		mt19937 gen(rd());

		if (density < 0.5) {
			binomial_distribution<> d(targets.size() - 1, density);

			for (int i = 0; i < (int) sources.size(); ++i) {
				int num_neighbors = d(gen);
				set<int> neighbor_ids;
				uniform_int_distribution<> u_dis(0, targets.size() - 1);

				while ((int) neighbor_ids.size() < num_neighbors) {
					int rand_neighbor = u_dis(gen);
					if (&targets[rand_neighbor] != &sources[i]) {
						neighbor_ids.insert(rand_neighbor);
					}
				}

				set<int>::iterator it;
				for (it = neighbor_ids.begin(); it != neighbor_ids.end(); ++it) {
					Edge e(sources[i], targets[*it]);
					sources[i].add_neighbor(e);
				}
			}
		} else {
			binomial_distribution<> d(targets.size() - 1, 1.0 - density);

			for (int i = 0; i < (int) sources.size(); ++i) {
				int num_non_neighbors = d(gen);
				set<int> non_neighbor_ids;
				uniform_int_distribution<> u_dis(0, targets.size() - 1);

				while ((int) non_neighbor_ids.size() < num_non_neighbors) {
					int rand_non_neighbor = u_dis(gen);
					if (&targets[rand_non_neighbor] != &sources[i]) {
						non_neighbor_ids.insert(rand_non_neighbor);
					}
				}

				for (int j = 0; j < (int) targets.size(); j++) {
					if (non_neighbor_ids.count(j) == 0) {
						Edge e(sources[i], targets[j]);
						sources[i].add_neighbor(e);
					}
				}
			}
		}
	}

	void connect(vector<Vertex>& sources, vector<Vertex>& targets, double density, void (*process_edge)(Edge&)) {
		random_device rd;
		mt19937 gen(rd());

		if (density < 0.5) {
			binomial_distribution<> d(targets.size() - 1, density);

			for (int i = 0; i < (int) sources.size(); ++i) {
				int num_neighbors = d(gen);
				set<int> neighbor_ids;
				uniform_int_distribution<> u_dis(0, targets.size() - 1);

				while ((int) neighbor_ids.size() < num_neighbors) {
					int rand_neighbor = u_dis(gen);
					if (&targets[rand_neighbor] != &sources[i]) {
						neighbor_ids.insert(rand_neighbor);
					}
				}

				set<int>::iterator it;
				for (it = neighbor_ids.begin(); it != neighbor_ids.end(); ++it) {
					Edge e(sources[i], targets[*it]);
					process_edge(e);
					sources[i].add_neighbor(e);
				}
			}
		} else {
			binomial_distribution<> d(targets.size() - 1, 1.0 - density);

			for (int i = 0; i < (int) sources.size(); ++i) {
				int num_non_neighbors = d(gen);
				set<int> non_neighbor_ids;
				uniform_int_distribution<> u_dis(0, targets.size() - 1);

				while ((int) non_neighbor_ids.size() < num_non_neighbors) {
					int rand_non_neighbor = u_dis(gen);
					if (&targets[rand_non_neighbor] != &sources[i]) {
						non_neighbor_ids.insert(rand_non_neighbor);
					}
				}

				for (int j = 0; j < (int) targets.size(); j++) {
					if (non_neighbor_ids.count(j) == 0) {
						Edge e(sources[i], targets[j]);
						process_edge(e);
						sources[i].add_neighbor(e);
					}
				}
			}
		}
	}

	void connect(vector<Geo_vertex>& sources, vector<Geo_vertex>& targets, double radius, void (*process_edge)(Edge&)) {
		for (int i = 0; i < (int) sources.size(); i++) {
			for (int j = 0; j < (int) targets.size(); j++) {
				if (&sources[i] != &targets[j]) {
					if (Utils::squared_euclidean_distance_2D(sources[i].get_x(), sources[i].get_y(), 
																targets[j].get_x(), targets[j].get_y()) <= (radius * radius)) {
						Edge e(sources[i], targets[j]);
						process_edge(e);
						sources[i].add_neighbor(e);
					}
				}
			}
		}
	}

	void connect_with_in(vector<Vertex>& sources, vector<Vertex>& targets,
			double density, void (*process_edge)(Edge&)) {
		random_device rd;
		mt19937 gen(rd());

		if (density < 0.5) {
			binomial_distribution<> d(targets.size() - 1, density);

			for (int i = 0; i < (int) sources.size(); ++i) {
				int num_neighbors = d(gen);
				set<int> neighbor_ids;
				uniform_int_distribution<> u_dis(0, targets.size() - 1);

				while ((int) neighbor_ids.size() < num_neighbors) {
					int rand_neighbor = u_dis(gen);
					if (&targets[rand_neighbor] != &sources[i]) {
						neighbor_ids.insert(rand_neighbor);
					}
				}

				set<int>::iterator it;
				for (it = neighbor_ids.begin(); it != neighbor_ids.end(); ++it) {
					Edge* e = new Edge(sources[i], targets[*it]);
					process_edge(*e);
					sources[i].add_neighbor(*e);
					targets[*it].inneighbors.push_back(e);
				}
			}
		} else {
			binomial_distribution<> d(targets.size() - 1, 1.0 - density);

			for (int i = 0; i < (int) sources.size(); ++i) {
				int num_non_neighbors = d(gen);
				set<int> non_neighbor_ids;
				uniform_int_distribution<> u_dis(0, targets.size() - 1);

				while ((int) non_neighbor_ids.size() < num_non_neighbors) {
					int rand_non_neighbor = u_dis(gen);
					if (&targets[rand_non_neighbor] != &sources[i]) {
						non_neighbor_ids.insert(rand_non_neighbor);
					}
				}

				for (int j = 0; j < (int) targets.size(); j++) {
					if (non_neighbor_ids.count(j) == 0) {
						Edge* e = new Edge(sources[i], targets[j]);
						process_edge(*e);
						sources[i].add_neighbor(*e);
						targets[j].inneighbors.push_back(e);
					}
				}
			}
		}
	}
}