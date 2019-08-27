/*
 * Connect.hpp
 *
 *  Created on: May 13, 2015
 *      Author: Trujic
 */

#ifndef RESET_HPP_
#define RESET_HPP_

#include <limits>

#include "../Vertices/Vertex.h"

namespace Reset {
	template<typename T>
	void reset_vertices(std::vector<T>& vertices) {
		for (int i = 0; i < (int) vertices.size(); i++) {
			vertices[i].set_refractory_period(1.0);
			vertices[i].set_voltage(0.0);
			vertices[i].set_infinity_voltage(0.0);
			vertices[i].set_sleep_until(-1.0);
			vertices[i].spike(-1.0);
			vertices[i].set_activation_threshold(0);
			vertices[i].set_inifinity_activation_threshold(0);
			vertices[i].patterns.clear();
		}
	}

	template<typename T>
	void reset_voltages(std::vector<T>& vertices, double voltage) {
		for (int i = 0; i < (int) vertices.size(); i++) {
			vertices[i].voltage = voltage;
			vertices[i].infinity_voltage = voltage;
		}
	}

	template<typename T>
	void reset_spike_time(std::vector<T>& vertices, double spike_time) {
		for (int i = 0; i < (int) vertices.size(); i++) {
			vertices[i].spike(spike_time);
		}
	}

	template<typename T>
	void reset_edges(std::vector<T>& sources, std::vector<T>& targets, void (*process_edge)(Edge&)) {
		set<int> indices;
		for (int i = 0; i < (int) targets.size(); i++) {
			indices.insert(targets[i].get_id());
		}

		for (int i = 0; i < (int) sources.size(); i++) {
			for (int j = 0; j < (int) sources[i].neighbors.size(); j++) {
				if (indices.count(sources[i].neighbors[j].target->get_id())) {
					sources[i].neighbors[j].delivery_time = -1;
					sources[i].neighbors[j].set_pattern(std::numeric_limits<int>::max());
					process_edge(sources[i].neighbors[j]);
				}
			}
		}
	}

	template<typename T>
	void reset_patterns(std::vector<std::set<T> >& patterns_A, std::vector<std::set<T> >& patterns_B) {
		std::set<Vertex*>::iterator it;
		for (int i = 0; i < (int) patterns_A.size(); i++) {
			for (it = patterns_A[i].begin(); it != patterns_A[i].end(); it++) {
				(*it)->patterns.clear();
				for (int j = 0; j < (int) (*it)->neighbors.size(); j++) {
					if ((*it)->neighbors[j].target->get_type() == 1) {
						(*it)->neighbors[j].set_weight(0);
						(*it)->neighbors[j].set_pattern(numeric_limits<int>::max());
					}
				}
			}
		}

		for (int i = 0; i < (int) patterns_B.size(); i++) {
			for (it = patterns_B[i].begin(); it != patterns_B[i].end(); it++) {
				(*it)->patterns.clear();
				for (int j = 0; j < (int) (*it)->neighbors.size(); j++) {
					if ((*it)->neighbors[j].target->get_type() == 1) {
						(*it)->neighbors[j].set_weight(0);
						(*it)->neighbors[j].set_pattern(numeric_limits<int>::max());
					}
				}
			}
		}

		patterns_A.clear();
		patterns_B.clear();
	}

	template<typename T>
	void full_reset(std::vector<T>& sources, void (*process_edge)(Edge&)) {
		for (int i = 0; i < (int) sources.size(); i++) {
			sources[i].set_refractory_period(1.0);
			sources[i].set_voltage(0.0);
			sources[i].set_infinity_voltage(0.0);
			sources[i].set_sleep_until(-1.0);
			sources[i].spike(-1.0);
			sources[i].set_activation_threshold(0);
			sources[i].set_inifinity_activation_threshold(0);
			sources[i].patterns.clear();

			for (int j = 0; j < (int) sources[i].neighbors.size(); j++) {
				// Reset all outgoing edges
				sources[i].neighbors[j].delivery_time = -1;
				sources[i].neighbors[j].set_pattern(std::numeric_limits<int>::max());
				process_edge(sources[i].neighbors[j]);

				// Reset all neighbours of the vertex as well
				sources[i].neighbors[j].target.set_refractory_period(1.0);
				sources[i].neighbors[j].target.set_voltage(0.0);
				sources[i].neighbors[j].target.set_infinity_voltage(0.0);
				sources[i].neighbors[j].target.set_sleep_until(-1.0);
				sources[i].neighbors[j].target.spike(-1.0);
				sources[i].neighbors[j].target.set_activation_threshold(0);
				sources[i].neighbors[j].target.set_inifinity_activation_threshold(0);
				sources[i].neighbors[j].target.patterns.clear();
			}
		}
	}

	template<typename T>
	void graph_reset_vertices(T& LG) {
		for (int i = 0; i < (int) LG.layers.size(); i++) {
			// Reset excitatory vertices
			for (int j = 0; j < (int) LG.layers[i].excitatory.size(); j++) {
				LG.layers[i].excitatory[j].set_refractory_period(1.0);
				LG.layers[i].excitatory[j].set_voltage(0.0);
				LG.layers[i].excitatory[j].set_infinity_voltage(0.0);
				LG.layers[i].excitatory[j].set_sleep_until(-1.0);
				LG.layers[i].excitatory[j].spike(-1.0);
				LG.layers[i].excitatory[j].set_activation_threshold(0);
				LG.layers[i].excitatory[j].set_inifinity_activation_threshold(0);
				LG.layers[i].excitatory[j].patterns.clear();

			}

			// Reset inhibitory vertices
			for (int j = 0; j < (int) LG.layers[i].inhibitory.size(); j++) {
				LG.layers[i].inhibitory[j].set_refractory_period(1.0);
				LG.layers[i].inhibitory[j].set_voltage(0.0);
				LG.layers[i].inhibitory[j].set_infinity_voltage(0.0);
				LG.layers[i].inhibitory[j].set_sleep_until(-1.0);
				LG.layers[i].inhibitory[j].spike(-1.0);
				LG.layers[i].inhibitory[j].set_activation_threshold(0);
				LG.layers[i].inhibitory[j].set_inifinity_activation_threshold(0);
				LG.layers[i].inhibitory[j].patterns.clear();
			}
		}
	}

	template<typename T>
	void graph_reset_edges(T& LG) {
		for (int i = 0; i < (int) LG.layers.size(); i++) {
			// Reset outgoing edges of excitatory vertices
			for (int j = 0; j < (int) LG.layers[i].excitatory.size(); j++) {
				for (int k = 0; k < (int) LG.layers[i].excitatory[j].neighbors.size(); k++) {
					LG.layers[i].excitatory[j].neighbors[k].delivery_time = -1;
					LG.layers[i].excitatory[j].neighbors[k].set_pattern(std::numeric_limits<int>::max());
				}
			}

			// Reset outgoing edges of inhibitory vertices
			for (int j = 0; j < (int) LG.layers[i].inhibitory.size(); j++) {
				for (int k = 0; k < (int) LG.layers[i].inhibitory[j].neighbors.size(); k++) {
					LG.layers[i].inhibitory[j].neighbors[k].delivery_time = -1;
					LG.layers[i].inhibitory[j].neighbors[k].set_pattern(std::numeric_limits<int>::max());
				}
			}
		}
	}
}

#endif /* RESET_HPP_ */