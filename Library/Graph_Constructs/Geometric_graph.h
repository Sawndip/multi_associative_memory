/*
 * Geometric_graph.h
 *
 *  Created on: Apr 14, 2015
 *      Author: Trujic
 */

#ifndef GEOMETRIC_GRAPH_H_
#define GEOMETRIC_GRAPH_H_

#include <vector>
#include "../Vertices/Geo_vertex.h"
#include "Basic_graph.h"

using namespace std;

class Geometric_graph {
public:
	vector<Geo_vertex> excitatory, inhibitory;

	Geometric_graph(int n_ex, int n_in, double spike_rate_pos, double spike_rate_neg);
	Geometric_graph(int n_ex, int n_in, double spike_rate_pos, double spike_rate_neg, int init_voltage);
	Geometric_graph(int n_ex, int n_in, double spike_rate_pos, double spike_rate_neg, double X, double Y);
};

#endif /* GEOMETRIC_GRAPH_H_ */
