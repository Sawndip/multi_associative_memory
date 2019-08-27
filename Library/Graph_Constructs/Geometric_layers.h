/*
 * Geometric_layers.h
 *
 *  Created on: Apr 15, 2015
 *      Author: Trujic
 */

#include "Layered_graph.h"
#include "Geometric_graph.h"

class Geometric_layers: public Layered_graph {
private:
	double distance;

	void connect(vector<Geo_vertex>& sources, vector<Geo_vertex>& targets, double distance, void (*process_edge)(Edge&));
public:
	vector<Geometric_graph> layers;

	Geometric_layers(vector<int>& excitatory, vector<int>& inhibitory, double distance, double ref_p_pos, double ref_p_neg);

	void connect(double radius, double distance, void (*process_edge)(Edge&));
};


