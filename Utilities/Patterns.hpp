/*
 * Connect.hpp
 *
 *  Created on: May 13, 2015
 *      Author: Trujic
 */

#ifndef PATTERNS_HPP_
#define PATTERNS_HPP_

#include <vector>
#include <set>
#include <iostream>

#include "../Vertices/Vertex.h"

namespace Patterns {
	template<typename T>
	std::set<Vertex*> sample_pattern(std::vector<T>& vertices, int k, int pattern_num) {
		if (k > vertices.size()) {
			clog << "You cannot have a pattern larger than your vertex set.\n";
		}

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(0, vertices.size() - 1);

		std::set<int> indices;
		if (k < 0.5 * vertices.size()) {
			while (indices.size() < k) {
				indices.insert(dis(gen));
			}
		} else {
			set<int> not_indices;
			while (not_indices.size() < vertices.size() - k) {
				not_indices.insert(dis(gen));
			}
			for (int i = 0; i < (int) vertices.size(); ++i) {
				if (not_indices.count(i) == 0) {
					indices.insert(i);
				}
			}
		}
		std::set<Vertex*> pattern;
		std::set<int>::iterator it;
		for (it = indices.begin(); it != indices.end(); ++it) {
			pattern.insert(&vertices[(*it)]);
			(&vertices[(*it)])->patterns.insert(pattern_num);
		}

		return pattern;
	}

	/**
	 * @param	sources		source vertices for patterns (layer A)
	 * @param	targets		target vertices for patterns (layer B)
	 * @param	k 			size of each pattern
	 * @param	bias_A		number of patters in layer A corresponding to one pattern in layer B
	 * @param	bias_B		number of patters in layer B corresponding to one pattern in layer A
	 * @param	current_A	current number of stored patterns in layer A
	 * @param	current_B	current number of stored patterns in layer B
	 *
	 * @return	a pair of vectors of sets, first vector corresponds to patterns inserted to layer A, second
	 			vector corresponds to patterns inserted to layer B
	 */
	template<typename T>
	pair<vector<set<Vertex*> >, vector<set<Vertex*> > > reveal_patterns(
			vector<T>& sources, vector<T>& targets, int k, int bias_A, int bias_B, int current_A, int current_B) {
		vector<set<Vertex*> > sample_A;
		vector<set<Vertex*> > sample_B;
		set<Vertex*>::iterator it;

		for (int i = 0; i < bias_B; i++) {
			set<Vertex*> B = Patterns::sample_pattern(targets, k, current_B + i);
			for (int j = 0; j < bias_A; j++) {
				set<Vertex*> A = Patterns::sample_pattern(sources, k, current_A + (i * bias_A) + j);

				// All edges between Ai and Bi become strong
				for (it = A.begin(); it != A.end(); it++) {
					for (int k = 0; k < (int) (*it)->neighbors.size(); k++) {
						if (B.count((*it)->neighbors[k].target) != 0) {
							(*it)->neighbors[k].set_weight(1);
							// The edge belongs to the min of its current pattern or the new pattern
							(*it)->neighbors[k].set_pattern(min(current_B, (*it)->neighbors[k].get_pattern()));
						}
					}
				}

				// All edges inside Bi become strong
				for (it = B.begin(); it != B.end(); it++) {
					for (int k = 0; k < (int) (*it)->neighbors.size(); k++) {
						if (B.count((*it)->neighbors[k].target) != 0) {
							(*it)->neighbors[k].set_weight(1);
							// The edge belongs to the min of its current pattern or the new pattern_B
							(*it)->neighbors[k].set_pattern(min(current_B, (*it)->neighbors[k].get_pattern()));
						}
					}
				}

				sample_A.push_back(A);
			}

			sample_B.push_back(B);
		}

		return make_pair(sample_A, sample_B);
	}
}

#endif /* PATTERNS_HPP_ */