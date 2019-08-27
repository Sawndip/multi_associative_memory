/*
 * Layered_graph.cpp
 *
 *  Created on: Feb 18, 2015
 *      Author: Trujic
 */

#include "Layered_graph.h"
#include "../Utilities/Utils.hpp"
#include "../Percolation/Percolation.hpp"

#include <iostream>
#include <algorithm>

using namespace std;

Layered_graph::Layered_graph() {
}

/**
 * Constructor of the layered graph. Note that it *only* creates vertices
 * and DOES NOT connect them
 *
 * @param	excitatory	number of excitatory vertices in i-th layer
 * @param	inhibitory	number of inhibitory vertices in i-th layer
 * @param	ref_p_pos	refractory period of excitatory vertices
 * @param	ref_p_neg	refractory period of inhibitory vertices
 */
Layered_graph::Layered_graph(vector<int>& excitatory, vector<int>& inhibitory, double ref_p_pos,
		double ref_p_neg) {
	for (int i = 0; i < (int) excitatory.size(); i++) {
		Basic_graph G(excitatory[i], inhibitory[i], ref_p_pos, ref_p_neg);
		layers.push_back(G);
	}
}

/**
 * @param	sources		vector of source vertices
 * @param 	targets		vector of target vertices
 * @param	k_size		size of a pattern
 * @param	bias_A		number of patters in layer A corresponding to one pattern in layer B
 * @param	bias_B		number of patters in layer B corresponding to one pattern in layer A
 * @param	M			number of patterns to reveal
 */
void Layered_graph::reveal_patterns(vector<Vertex>& sources,
		vector<Vertex>& targets, int k_size, int bias_A, int bias_B, int M) {
	for (int i = 0; i < M; i++) {
		pair<vector<set<Vertex*> >, vector<set<Vertex*> > > patterns = Patterns::reveal_patterns(
				sources, targets, k_size, bias_A, bias_B, patterns_A.size() + 1, patterns_B.size() + 1);

		patterns_A.insert(patterns_A.end(), patterns.first.begin(), patterns.first.end());
		patterns_B.insert(patterns_B.end(), patterns.second.begin(), patterns.second.end());
	}
}

/**
 * Runs the percolation in pattern revealing manner
 *
 * @param	k_size		size of a pattern
 * @param	K_e			activation threshold for excitatory vertices
 * @param	K_i			activation threshold for inhibitory vertices
 * @param	bias_A		number of patters in layer A corresponding to one pattern in layer B
 * @param	bias_B		number of patters in layer B corresponding to one pattern in layer A
 *
 * @return	a pair of vectors of (x, y) points where x is the number of revealed patterns
 * 			and y is the number of exposed vertices
 */
pair<vector<pair<double, int> >, vector<pair<double, int> > > Layered_graph::run_pattern_percolation(set<Vertex*>& bootstrap,
		int k_size, int K_e, int K_i, int bias_A, int bias_B, double num_first, double num_second) {
	set<Vertex*> exposed;
	bootstrap.clear();

	vector<pair<double, int> > history_patterns;
	vector<pair<double, int> > history;

	int upper_bound = (num_second > 0) ? 2 : 1;

	reveal_patterns(layers[0].excitatory, layers[1].excitatory, k_size, bias_A, bias_B, upper_bound);

	if (num_second > 0) {
		set<Vertex*> sample = Utils::sample_bootstrap(patterns_A[0], k_size * num_first);
		bootstrap.insert(sample.begin(), sample.end());
		sample = Utils::sample_bootstrap(patterns_A[1], k_size * num_second);
		bootstrap.insert(sample.begin(), sample.end());
	} else {
		bootstrap.insert(patterns_A[0].begin(), patterns_A[0].end());
	}

	exposed = percolation::pattern_activation(bootstrap, K_e, K_i, 1, upper_bound);

	int pattern_active = Utils::count_active(patterns_B[0], K_e);
	int active_ex = Utils::intersection_size(layers[1].excitatory, exposed);

	Reset::reset_voltages(layers[1].excitatory, 0);
	Reset::reset_voltages(layers[1].inhibitory, 0);

	/**
	 *  IMPORTANT: Uncomment only if density of afferent edges is optimal!
	 *  Otherwise we might get stuck in an endless loop here
	 */
	while (pattern_active < k_size * 0.8) {
		Reset::reset_patterns(patterns_A, patterns_B);
		reveal_patterns(layers[0].excitatory, layers[1].excitatory, k_size, bias_A, bias_B, upper_bound);

		bootstrap.clear();
		
		if (num_second > 0) {
			set<Vertex*> sample = Utils::sample_bootstrap(patterns_A[0], k_size * num_first);
			bootstrap.insert(sample.begin(), sample.end());
			sample = Utils::sample_bootstrap(patterns_A[1], k_size * num_second);
			bootstrap.insert(sample.begin(), sample.end());
		} else {
			bootstrap.insert(patterns_A[0].begin(), patterns_A[0].end());
		}

		exposed = percolation::pattern_activation(bootstrap, K_e, K_i, 1, upper_bound);

		pattern_active = Utils::count_active(patterns_B[0], K_e);
		active_ex = Utils::intersection_size(layers[1].excitatory, exposed);

		Reset::reset_voltages(layers[1].excitatory, 0);
		Reset::reset_voltages(layers[1].inhibitory, 0);
	}

	history_patterns.push_back(make_pair(upper_bound, pattern_active));
	history.push_back(make_pair(upper_bound, active_ex));

	// Exponential Search for the upper bound on number of patterns before full percolation
	while (active_ex < 0.8 * layers[0].excitatory.size()) {
		reveal_patterns(layers[0].excitatory, layers[1].excitatory, k_size, bias_A, bias_B, upper_bound);

		upper_bound *= 2;

		exposed = percolation::pattern_activation(bootstrap, K_e, K_i, 1, upper_bound);
		pattern_active = Utils::count_active(patterns_B[0], K_e);
		active_ex = Utils::intersection_size(layers[1].excitatory, exposed);

		history_patterns.push_back(make_pair(upper_bound * bias_A, pattern_active));
		history.push_back(make_pair(upper_bound * bias_A, active_ex));

		cout << "With " << upper_bound * bias_A << " revealed patterns, activated in pattern " 
				<< pattern_active << " and in total: " << active_ex << endl;

		Reset::reset_voltages(layers[1].excitatory, 0);
		Reset::reset_voltages(layers[1].inhibitory, 0);
	}

	// Binary Search for the threshold number of revealed patterns before full percolation
	int lower_bound = upper_bound / 2;
	int pivot = (upper_bound + lower_bound + 1) / 2;
	while (lower_bound < upper_bound) {
		exposed = percolation::pattern_activation(bootstrap, K_e, K_i, 1, pivot);
		pattern_active = Utils::count_active(patterns_B[0], K_e);
		active_ex = Utils::intersection_size(layers[1].excitatory, exposed);

		history_patterns.push_back(make_pair(pivot * bias_A, pattern_active));
		history.push_back(make_pair(pivot * bias_A, active_ex));

		if (active_ex > 0.8 * layers[0].excitatory.size()) {
			upper_bound = pivot - 1;
		} else {
			lower_bound = pivot + 1;
		}

		cout << "With " << pivot * bias_A << " revealed patterns, activated in pattern "
				<< pattern_active << " and in total: " << active_ex << endl;

		Reset::reset_voltages(layers[1].excitatory, 0);
		Reset::reset_voltages(layers[1].inhibitory, 0);

		pivot = (lower_bound + upper_bound + 1) / 2;
		cout << lower_bound << " " << pivot << " " << upper_bound << endl;
	}

	return make_pair(history_patterns, history);
}

double Layered_graph::learn(int thrlow, int thrhigh, double pweak, double pstrong,
		bool supervised) {
	int num_strong = 0;
	random_device rd;
	mt19937 gen(rd());
	bernoulli_distribution dstrong(pstrong), dweak(pweak);
	for (int i = 0; i < (int) layers[0].excitatory.size(); ++i) {
		for (int j = 0; j < (int) layers[0].excitatory[i].neighbors.size();
				j++) {
			Edge *e = &layers[0].excitatory[i].neighbors[j];
			// Update memory based on last round
			if (e->target->spike_time >= 0 && e->delivery_time >= 0) {
				if (e->delivery_time > e->target->spike_time) {
					e->add_memory(0);
				} else {
					e->add_memory(1);
				}
			} else if (supervised) {
				// All the signals arrive before the vertex spikes
				e->add_memory(1);
			}
			// Update weights
			if ((e->target->spike_time >= 0 && e->delivery_time > 0)
					|| supervised) {
				if ((e->get_weight() == 1) && (e->memory_value <= thrlow)) {
					e->set_weight(1 - dweak(gen));
				} else if ((e->get_weight() == 0)
						&& (e->memory_value >= thrhigh)) {
					e->set_weight(dstrong(gen));
				}
			}
			e->delivery_time = -1; // Reset
			num_strong += e->get_weight();
		}
	}

	return num_strong / (double) layers[1].excitatory.size();
}

/**
 *	Prints the layered graph in a Mathematica-readable manner
 *
 *	@param	out		a file to which to write
 *	@param	G		a layered graph
 */
ostream& operator<<(ostream &out, Layered_graph &G) {
	for (int i = 0; i < (int) G.layers[1].excitatory.size(); i++) {
		int index = 0;

		for (int j = 0; j < (int) G.layers[1].excitatory[i].neighbors.size();
				j++) {
			while (index
					!= G.layers[1].excitatory[i].neighbors[j].target->get_id()) {
				out << "0 ";
				index++;
			}
			out << "1 ";
			index++;
		}
		while (index < (int) G.layers[1].excitatory.size()) {
			out << "0 ";
			index++;
		}
		out << "\n";
	}

	return out;
}
