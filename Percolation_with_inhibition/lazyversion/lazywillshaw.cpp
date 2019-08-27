#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <set>
#include <random>
#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <algorithm>
using namespace std;

struct Vertex { // Vertex structure
	int active;
	int voltage; // active if it reaches activation threshold. If above act_threshold then the neuron has already spiked (and will not be considered further)
	//int activation_time; // stores percolation round when vertex turned active
	//vector<int> strong_neighbours;
};

struct Pattern{ // Pattern
	set<int> verticesinpattern;
};





// single willshaw experiment
int lazywillshaw (int number_vertices, int pattern_size, int multi_association_factor, int fidelity_outside, int act_threshold) {

	// Random devices;
    std::random_device rd;
    std::mt19937 gen(rd());
	uniform_int_distribution<> uniform_vertex(0, number_vertices - 1);

	// to be searched for and returned!
	int number_patterns_rec_layer=1;

	// to store the connections matrix
	//int strong_afferent_connections[pattern_size][number_vertices];
	vector<vector <int> > strong_afferent_connections;
	int number_strong_outside_connections =0;

	// initialize afferent connections
	for (int i = 0; i < pattern_size; i++) {
		vector <int> vertex_i_afferent_connections;
		for (int j = 0; j < number_vertices; j++) {
			vertex_i_afferent_connections.push_back(0);
		}
		strong_afferent_connections.push_back(vertex_i_afferent_connections);
	}
	

	// vertex second layer set
	vector<Vertex> Secondlayervertices;

	//vertices added to second layer
	for (int i = 0; i < number_vertices; i++) {
		Vertex v;
		v.voltage = 0;
		v.active = 0;
		Secondlayervertices.push_back(v);
	}

	// defining the source pattern
	Pattern Source;
	for (int i = 0; i < pattern_size; i++) {
		Source.verticesinpattern.insert(i);
	}

	// target pattern
	Pattern Target;
	while (Target.verticesinpattern.size() < pattern_size) {
		Target.verticesinpattern.insert(uniform_vertex(gen));
	}

	// inserting the target pattern
	for(int i = 0; i < pattern_size; i++) {
		for(set<int>::iterator vertex_target_it = Target.verticesinpattern.begin(); vertex_target_it != Target.verticesinpattern.end(); vertex_target_it++) {
			strong_afferent_connections[i][*vertex_target_it] = 1;
			Secondlayervertices[*vertex_target_it].voltage = act_threshold;
			Secondlayervertices[*vertex_target_it].active = 1;
		}
	}

	// other patterns - insert until fidelity_outside is breached
	int outside_active_vertices = 0;
	int outside_active_vertices_previous_round = outside_active_vertices;
	while (outside_active_vertices < fidelity_outside - pattern_size) {
		outside_active_vertices_previous_round = outside_active_vertices;

		// try to insert one more pattern
		number_patterns_rec_layer++;


		// patterns afferent layer
		//vector<Pattern> afferent_patterns_vector;
		set<int> afferent_patterns_union; //stores the union of the s patterns added in first layer
		for (int i = 0; i < multi_association_factor; i++) {
			Pattern new_afferent_pattern;
			while (new_afferent_pattern.verticesinpattern.size() < pattern_size) {
				int vertex_index = uniform_vertex(gen);
				new_afferent_pattern.verticesinpattern.insert(vertex_index);
				afferent_patterns_union.insert(vertex_index);
			}
			//afferent_patterns_vector.push_back(new_afferent_pattern);
		}

		vector<int> intersection(pattern_size);
		//intersection vertices are those that will have additional strong connections to outside vertices!
		vector<int>::iterator it =set_intersection (afferent_patterns_union.begin(),afferent_patterns_union.end(),Source.verticesinpattern.begin(),Source.verticesinpattern.end(), intersection.begin());
		intersection.resize(it-intersection.begin());

		if(intersection.size() != 0) {
			// build corresponding recurrent pattern
			Pattern new_recurrent_pattern;
			while (new_recurrent_pattern.verticesinpattern.size() < pattern_size) {
				new_recurrent_pattern.verticesinpattern.insert(uniform_vertex(gen));
			}

			

			for (vector<int>::iterator intersection_iterator = intersection.begin(); intersection_iterator != intersection.end(); intersection_iterator++) {
				for (set<int>::iterator new_recurrent_pattern_iterator = new_recurrent_pattern.verticesinpattern.begin(); new_recurrent_pattern_iterator != new_recurrent_pattern.verticesinpattern.end(); new_recurrent_pattern_iterator++) {
					
					if(strong_afferent_connections[*intersection_iterator][*new_recurrent_pattern_iterator] == 0) {
						strong_afferent_connections[*intersection_iterator][*new_recurrent_pattern_iterator] = 1;
						number_strong_outside_connections++;

						if (Secondlayervertices[*new_recurrent_pattern_iterator].voltage < act_threshold) {
							Secondlayervertices[*new_recurrent_pattern_iterator].voltage++;
						}
						if (Secondlayervertices[*new_recurrent_pattern_iterator].voltage == act_threshold && Secondlayervertices[*new_recurrent_pattern_iterator].active == 0) {
							Secondlayervertices[*new_recurrent_pattern_iterator].active = 1;
							outside_active_vertices++;

						}
					}
				}
			}
		}
	}
	number_patterns_rec_layer--;
	return number_patterns_rec_layer;

}

int main(){
	//parameters
	int number_vertices = 400000;
	int pattern_size = 15;
	int act_threshold = pattern_size;
	int multi_association_factor = 64;
	int fidelity_outside = 2*pattern_size;
	int final_number_patterns_stored = 4*floor(number_vertices*number_vertices/pattern_size/pattern_size/multi_association_factor);

	//number tests;
	int number_trials = 30;

	double average_number_patterns = 0;
    //omp_set_num_threads(1);
	#pragma omp parallel for
	for (int i = 0; i < number_trials; i++) {
		cout << i << "\n";
		int number_patterns_rec_layer_trial = lazywillshaw(number_vertices, pattern_size, multi_association_factor, fidelity_outside, act_threshold);
		average_number_patterns += ((double)number_patterns_rec_layer_trial)/number_trials;
	}
	cout << average_number_patterns << "\n";
	return average_number_patterns;
}
