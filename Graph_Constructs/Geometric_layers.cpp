/*
 * Geometric_layers.cpp
 *
 *  Created on: Apr 15, 2015
 *      Author: Trujic
 */

#include <iostream>
 
#include "Geometric_layers.h"
#include "../Utilities/Utils.hpp"
#include "../Utilities/Connect.hpp"

using namespace std;

Geometric_layers::Geometric_layers(vector<int>& excitatory, vector<int>& inhibitory, double distance, double ref_p_pos, double ref_p_neg) {
	this->distance = distance;

	for (int i = 0; i < (int) excitatory.size(); i++) {
		Geometric_graph GG(excitatory[i], inhibitory[i], ref_p_pos, ref_p_neg);
		layers.push_back(GG);
	}
}

//TODO
//	[ ] Move the connect function to Connect utils?
void Geometric_layers::connect(double radius, double distance, void (*process_edge)(Edge&)) {
	//TODO: Add inhibition
	for (int i = 1; i < (int) layers.size(); i++) {
		clog << "Ex -> Ex\n";
		Connect::connect(layers[i].excitatory, layers[i].excitatory, radius, process_edge);
	}

	clog << "A -> B, Ex -> Ex\n";
	connect(layers[0].excitatory, layers[1].excitatory, distance, process_edge);
}

void Geometric_layers::connect(vector<Geo_vertex>& sources, vector<Geo_vertex>& targets, double distance, void (*process_edge)(Edge&)) {
	for (int i = 0; i < (int) sources.size(); i++) {
		for (int j = 0; j < (int) targets.size(); j++) {
			if (Utils::squared_euclidean_distance_3D(sources[i].get_x(), sources[i].get_y(), targets[j].get_x(), targets[j].get_y(), 
														0.0, this->distance) <= distance * distance) {
				Edge e(sources[i], targets[j]);
				process_edge(e);
				sources[i].add_neighbor(e);
			}
		}
	}
}


