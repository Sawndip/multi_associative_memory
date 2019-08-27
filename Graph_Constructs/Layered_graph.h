/*
 * Layered_graph.h
 *
 *  Created on: Feb 18, 2015
 *      Author: Trujic
 */

#ifndef LAYERED_GRAPH_H_
#define LAYERED_GRAPH_H_

#include <vector>
#include <set>
#include <fstream>

#include "Basic_graph.h"
#include "Geometric_graph.h"

using namespace std;

class Layered_graph {
public:
	vector<Basic_graph> layers;
	vector<set<Vertex*> > patterns_A;
	vector<set<Vertex*> > patterns_B;

	Layered_graph();
	Layered_graph(vector<int>& excitatory, vector<int>& inhibitory, double ref_p_pos,
			double ref_p_neg);

	void reveal_patterns(vector<Vertex>& sources, vector<Vertex>& targets, int k_size, int A, 
			int B, int M);

	pair<vector<pair<double, int> >, vector<pair<double, int> > > run_pattern_percolation(set<Vertex*>& bootstrap,
			int k_size, int K_e, int K_i, int bias_A, int bias_B, double num_first, double num_second);

	double learn(int thrlow, int thrhigh, double pweak, double pstrong, bool supervised);

	friend ostream& operator<<(ostream &out, Layered_graph &G);
};

#endif /* LAYERED_GRAPH_H_ */
