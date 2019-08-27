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
#include <tuple>
using namespace std;

struct Edge { // Edge structure
	int patternturnstrong; // if patternturnstrong < number associations used than edge is strong in the graph
	int patternthreshold; // an edge turns strong iff it appeared on at least a given number of patterns
	bool is_present; // whether edge belongs to the graph or not
};

struct Vertex { // Vertex structure
	int voltage; // active if it reaches activation threshold. If above act_threshold then the neuron has already spiked (and will not be considered further)
	int activation_time; // stores percolation round when vertex turned active
	vector<int> neighbours;
};

struct Twolayergraph { // Graph structure
	int verticesfirstlayer;
	int verticessecondlayer;
	vector<Vertex> Firstlayer;
	vector<Vertex> Secondlayer;
	vector<vector<Edge> > Afferentedges;
	vector<vector<Edge> > Recurrentedges;
	//other parameters that are necessary for percolation:
	int pattern_size;
	int multi_association_factor;
	int final_number_associations;
	int turn_strong_pattern_threshold;
	int act_threshold;
	int inhibition_threshold;
	double fidelity_outside;
	double fidelity_inside;
};

struct Pattern{ // Pattern
	set<int> verticesinpattern;
};

// This builds a two layer graph with afferent connection p_aff and recurrent connection density p_rec. It should be called to build a new graph
void graphconstructor(Twolayergraph& G, int verticesfirstlayer, int verticessecondlayer, double p_aff, double p_rec){
	// define layer sizes
	G.verticesfirstlayer = verticesfirstlayer;
	G.verticessecondlayer = verticessecondlayer;

	// Random Devices
    std::random_device rd;
    std::mt19937 generator(rd());
  	bernoulli_distribution afferentdistribution(p_aff);
  	bernoulli_distribution recurrentdistribution(p_rec);

	// Build vertices in first layer && afferent edges
	for (int i = 0; i < verticesfirstlayer; i++){
		vector<Edge> vertex_i_edges;
		Vertex v;
		v.voltage = 0;
		v.activation_time = -1;
		for (int j = 0; j < verticessecondlayer; j++) {//define (afferent) neighbours of v!
			Edge e;
			e.patternturnstrong = -1;
			e.patternthreshold = 0;
			e.is_present = false;
  			if (afferentdistribution(generator)) {
				v.neighbours.push_back(j);
				e.is_present = true;
  			}
			vertex_i_edges.push_back(e);
		}
		G.Firstlayer.push_back(v);
		G.Afferentedges.push_back(vertex_i_edges);
	}

	// Build vertices in second layer && recurrent edges
	for (int i = 0; i < verticessecondlayer; i++){
		Vertex v;
		v.voltage = 0;
		v.activation_time = -1;
		vector<Edge> vertex_i_edges;
		for (int j = 0; j < verticessecondlayer; j++) {//define recurrent neighbours of v!
			Edge e;
			e.patternturnstrong = -1;
			e.patternthreshold = 0;
			e.is_present = false;
  			if (recurrentdistribution(generator) && (i != j)) {
				v.neighbours.push_back(j);
				e.is_present = true;
  			}
			vertex_i_edges.push_back(e);
		}
		G.Secondlayer.push_back(v);
		G.Recurrentedges.push_back(vertex_i_edges);
	}
}

// This inserts patterns at a given two layer graph G, associating many to one. It should be called to insert patterns after a new graph was built.
void pattern_insertion(vector<Pattern>& Patternsfirstlayer, vector<Pattern>& Patternssecondlayer, vector<Pattern>& Patternsfirstlayernoisy, Twolayergraph& G,int number_associations, int multi_association_factor, int pattern_size, int turn_strong_pattern_threshold, double false_activation_noise){

	// Extract layer sizes;
	int verticesfirstlayer = G.verticesfirstlayer;
	int verticessecondlayer = G.verticessecondlayer;

	// Random devices;
    std::random_device rd;
    std::mt19937 gen(rd());
	uniform_int_distribution<> disfirstlayer(0, verticesfirstlayer - 1);
	uniform_int_distribution<> dissecondlayer(0, verticessecondlayer - 1);

	for (int i = 0; i < number_associations; i++){// add associations one by one to G;

		// Build patterns in first layer
		vector<Pattern> sourcepatterns;
		for (int j = 0; j < multi_association_factor; j++) {
			Pattern P;
			Pattern P_noise;
			while (P.verticesinpattern.size() < pattern_size) {
				int pattern_vertex = disfirstlayer(gen);
				P.verticesinpattern.insert(pattern_vertex);
				P_noise.verticesinpattern.insert(pattern_vertex);
			}
			while (P_noise.verticesinpattern.size() < floor(pattern_size*(1+false_activation_noise))) {
				int pattern_vertex = disfirstlayer(gen);
				P_noise.verticesinpattern.insert(pattern_vertex);
			}
			sourcepatterns.push_back(P);
			Patternsfirstlayer.push_back(P);
			Patternsfirstlayernoisy.push_back(P_noise);
		}

		// Build pattern in second layer
		Pattern targetpattern;
		while (targetpattern.verticesinpattern.size() < pattern_size) {
			targetpattern.verticesinpattern.insert(dissecondlayer(gen));
		}
		Patternssecondlayer.push_back(targetpattern);

		// Turn afferent edges strong
		for (int j = 0; j < multi_association_factor; j++) {
			for (set<int>::iterator vertexsourceit = sourcepatterns[j].verticesinpattern.begin(); vertexsourceit != sourcepatterns[j].verticesinpattern.end(); vertexsourceit++){
				for (set<int>::iterator vertextargetit = targetpattern.verticesinpattern.begin(); vertextargetit != targetpattern.verticesinpattern.end(); vertextargetit++){
					if(G.Afferentedges[*vertexsourceit][*vertextargetit].is_present && G.Afferentedges[*vertexsourceit][*vertextargetit].patternthreshold < turn_strong_pattern_threshold) {
						G.Afferentedges[*vertexsourceit][*vertextargetit].patternthreshold++;
						G.Afferentedges[*vertexsourceit][*vertextargetit].patternturnstrong = i;
					}
				}
			}
		}

		// Turn recurrent edges strong
		for (set<int>::iterator vertexsourceit = targetpattern.verticesinpattern.begin(); vertexsourceit != targetpattern.verticesinpattern.end(); vertexsourceit++){
			for (set<int>::iterator vertextargetit = targetpattern.verticesinpattern.begin(); vertextargetit != targetpattern.verticesinpattern.end(); vertextargetit++) {
				if (G.Recurrentedges[*vertexsourceit][*vertextargetit].is_present && G.Recurrentedges[*vertexsourceit][*vertextargetit].patternthreshold < turn_strong_pattern_threshold) {
					G.Recurrentedges[*vertexsourceit][*vertextargetit].patternthreshold++;
					G.Recurrentedges[*vertexsourceit][*vertextargetit].patternturnstrong = i;
				}
			}
		}
	}
}


// Given a graph where patterns were inserted runs bootstrap percolation on the strong edge set of G.
// The bootstrap is the set of activated vertices in the first layer.
// activated_vertices is the set of vertices activated in second layer.
// inhibition_threshold is the number of vertices where inhibition kills the bootstrap.
// number_associations is the number of patterns inserted in second layer.
// turn_strong_pattern_threshold is the number of times an edge should be in a pattern to be strong.
// Percolation stops in the round where inhibition threshold vertices turned active.
void percolationwithinhibition(set<int>& activated_vertices, Twolayergraph& G, set<int> bootstrap, int act_threshold, int inhibition_threshold, int number_associations, int turn_strong_pattern_threshold) {

	int activated_vertices_size = 0;
	int final_round = -1; //final round of percolation - inhibition kicks in
	queue<int> verticesinqueue;
	set<int>::iterator vertexiterator;
	for (vertexiterator = bootstrap.begin(); vertexiterator != bootstrap.end(); vertexiterator++) {
		verticesinqueue.push(*vertexiterator);
		verticesinqueue.push(1); // first round only has first layer vertices. Second round onwards has only second layer vertices. Round is stored to stop program at the correct round!
		G.Firstlayer[*vertexiterator].voltage = act_threshold + 1;
		G.Firstlayer[*vertexiterator].activation_time = 1;
	}
	while(!verticesinqueue.empty()){
		int vertex_index = verticesinqueue.front();
		verticesinqueue.pop();
		int vertex_activation_round = verticesinqueue.front();
		verticesinqueue.pop();
		if (vertex_activation_round == 1) { // first round is only in first layer
			int neighbourhood_size = G.Firstlayer[vertex_index].neighbours.size();
			for (int i = 0; i < neighbourhood_size; i++) {
				int target_vertex = G.Firstlayer[vertex_index].neighbours[i];
				if (G.Afferentedges[vertex_index][target_vertex].patternturnstrong < number_associations && G.Afferentedges[vertex_index][target_vertex].patternthreshold == turn_strong_pattern_threshold){// unnecessary to check whether edge is present
					if (G.Secondlayer[target_vertex].voltage < act_threshold) {
						G.Secondlayer[target_vertex].voltage++;
					}
					if (G.Secondlayer[target_vertex].voltage == act_threshold) {
						G.Secondlayer[target_vertex].voltage++;
						G.Secondlayer[target_vertex].activation_time = G.Firstlayer[vertex_index].activation_time + 1;
						verticesinqueue.push(target_vertex);
						verticesinqueue.push(G.Secondlayer[target_vertex].activation_time);
					}
				} 
			}
		}
		else { // other rounds are only in second layer
			if (activated_vertices_size < inhibition_threshold - 1) {
				activated_vertices.insert(vertex_index); // should I insert once it reaches right threshold potential or when it will send the spike? - Current is only when it sends spike
				activated_vertices_size++;
				int neighbourhood_size = G.Secondlayer[vertex_index].neighbours.size();
				for (int i = 0; i < neighbourhood_size; i++) {
					int target_vertex = G.Secondlayer[vertex_index].neighbours[i];
					if (G.Recurrentedges[vertex_index][target_vertex].is_present && G.Recurrentedges[vertex_index][target_vertex].patternturnstrong < number_associations && G.Recurrentedges[vertex_index][target_vertex].patternthreshold == turn_strong_pattern_threshold){// unnecessary to check whether edge is present
						if (G.Secondlayer[target_vertex].voltage < act_threshold) {
							G.Secondlayer[target_vertex].voltage++;
						}
						if (G.Secondlayer[target_vertex].voltage == act_threshold) {
							G.Secondlayer[target_vertex].voltage++;
							G.Secondlayer[target_vertex].activation_time = G.Secondlayer[vertex_index].activation_time + 1;
							verticesinqueue.push(target_vertex);
							verticesinqueue.push(G.Secondlayer[target_vertex].activation_time);
						}
					}
				}
			}
			else if (activated_vertices_size == inhibition_threshold - 1) { // Defines the final round
				activated_vertices.insert(vertex_index);
				activated_vertices_size++;
				final_round = G.Secondlayer[vertex_index].activation_time;
				// no point in sending the spikes since the vertices that get it can only be active at final_round + 1;
			}
			else {
				if (G.Secondlayer[vertex_index].activation_time <= final_round) {
					activated_vertices.insert(vertex_index);
					//activated_vertices_size++; unnecessary;
					// no point in sending the spikes.
				}
			}
		}
	}

	// Reset voltages so that percolation procedure can be run again!
	for (int i = 0; i < G.verticesfirstlayer; i++) {
		G.Firstlayer[i].voltage = 0;
		G.Firstlayer[i].activation_time = -1;
	}
	for (int i = 0; i < G.verticessecondlayer; i++) {
		G.Secondlayer[i].voltage = 0;
		G.Secondlayer[i].activation_time = -1;
	}
}


// Checks how many vertices turn active inside a pattern. Run many times to get averages
int percolation_inside_active_set_size(int pattern_size, int act_threshold, double p_aff, double p_rec) {

	// Graph construction
	Twolayergraph G;
	graphconstructor(G, pattern_size, pattern_size, p_aff, p_rec);

	// Pattern insertion
	vector<Pattern> Patternsfirstlayer;
	vector<Pattern> Patternssecondlayer;
	vector<Pattern> Patternsfirstlayernoisy;
	pattern_insertion(Patternsfirstlayer, Patternssecondlayer, Patternsfirstlayernoisy, G, 1, 1,
		pattern_size, 1, 0);

	// Percolation - Pattern recall
	set<int> activated_vertices;
	int inhibition_threshold = pattern_size;
	percolationwithinhibition(activated_vertices, G, Patternsfirstlayer[0].verticesinpattern, act_threshold, inhibition_threshold, 1, 1);
	int number_activated_vertices = activated_vertices.size();
	return number_activated_vertices;
}

// Given p_rec and fidelity_inside requirement determines p_aff up to some 10^{-precision} so that percolation inside still happens.
double percolation_inside_p_aff_search_given_p_rec(int pattern_size, int act_threshold, int number_tries, int precision, double p_rec, double fidelity_inside) {

	double p_aff = 0.0;

	double increment=1;
	for(int k=1; k< precision+1; k++){
		increment= increment/(double) 10;
		double average_active=0;
		while (average_active < fidelity_inside * pattern_size) {

			p_aff+=increment;
			//cout << " Check p_aff= "<< p_aff << " \n";
			p_aff = min(1.0,p_aff);
			average_active=0;

			for(int i=0; i < number_tries; i++){
				int activity_inside = percolation_inside_active_set_size(pattern_size, act_threshold, p_aff, p_rec);
				average_active +=  activity_inside/(double) number_tries;
			}
			//cout << "average activity is " << average_active << " \n";
		}
		p_aff -= increment;
	}
	return p_aff;
}


// Alternatively:
// Given p_aff and fidelity_inside requirement determines p_rec up to some 10^{-precision} so that percolation inside still happens.
double percolation_inside_p_rec_search_given_p_aff(int pattern_size, int act_threshold, int number_tries, int precision, double p_aff, double fidelity_inside) {

	double p_rec = 0.0;

	double increment=1;
	for(int k=1; k< precision+1; k++){
		increment= increment/(double) 10;
		double average_active=0;
		while (average_active < fidelity_inside * pattern_size) {

			p_rec+=increment;
			//cout << " Check p_aff= "<< p_aff << " \n";
			p_rec = min(1.0,p_rec);
			average_active=0;

			for(int i=0; i < number_tries; i++){
				int activity_inside = percolation_inside_active_set_size(pattern_size, act_threshold, p_aff, p_rec);
				average_active +=  activity_inside/(double) number_tries;
				//cout << p_r << ": " << exposed.size() << " " << active_ex << endl;
			}
			//cout << "average activity is " << average_active << " \n";
		}
		p_rec -= increment;
	}
	return p_rec;
}


// Binary searches the optimal number of patterns which can be stored respecting fidelity outside. Single graph trial, should be averaged by calling it many times.
tuple<int,double,double, int> maximal_number_patterns_without_breaching_fidelity_outside_single_graph_trial (vector<Pattern>& Patternsfirstlayer, vector<Pattern>& Patternssecondlayer, vector<Pattern>& Patternsfirstlayernoisy, Twolayergraph& G, int number_trials_pattern) {

		int active_vertices = 0;
		double average_active_vertices=0;
		int active_vertices_inside = 0;
		double average_active_vertices_inside = 0; // average over number_trials_binary_search
		int inhibition_threshold = G.inhibition_threshold;
		double average_active_vertices_inside_down_final = 0;
		double average_active_vertices_down_final = 0;
		// Do binary search for the maximal number of patterns we can store
		int act_threshold_down = G.act_threshold;
		int base_act_threshold = G.act_threshold;

		int num_patterns=1;

		for (int j = 0; j < number_trials_pattern; j++) {
			int rand_bootstrap_pattern= rand()% num_patterns*G.multi_association_factor;
			set<int> activated_vertices;
			percolationwithinhibition(activated_vertices, G, Patternsfirstlayernoisy[rand_bootstrap_pattern].verticesinpattern, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
			active_vertices = activated_vertices.size();
			average_active_vertices += (double) active_vertices / (double) number_trials_pattern; 
			set<int> inside_active_vertices;
			set_intersection(Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.begin(),Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.end(),activated_vertices.begin(),activated_vertices.end(), inserter(inside_active_vertices, inside_active_vertices.begin()));
			average_active_vertices_inside += inside_active_vertices.size()/(double) number_trials_pattern;	
		}
		if (average_active_vertices < G.fidelity_outside) {
		//if (average_active_vertices < G.fidelity_outside && average_active_vertices_inside >= G.fidelity_inside*G.pattern_size) {
			average_active_vertices_inside_down_final = average_active_vertices_inside;
			average_active_vertices_down_final = average_active_vertices;
		}

		//cout << average_active_vertices << "= average_active_vertices\n";
		//cout << average_active_vertices_inside << "= average_active_vertices_inside\n";
		while (average_active_vertices < G.fidelity_outside){ // 0.05 is arbitrary decision
		//while (average_active_vertices < G.fidelity_outside && average_active_vertices_inside >= (G.fidelity_inside-0.05)*G.pattern_size){ // 0.05 is arbitrary decision
			num_patterns*=2;
			G.act_threshold = base_act_threshold;
			//clog << "check number of revealed patterns is "<< num_patterns << "\n";
			average_active_vertices = 0;
			average_active_vertices_inside = 0;
			for(int j=0; j< number_trials_pattern; j++ ){
				int rand_bootstrap_pattern= rand()% num_patterns*G.multi_association_factor;
				set<int> activated_vertices;
				percolationwithinhibition(activated_vertices, G, Patternsfirstlayernoisy[rand_bootstrap_pattern].verticesinpattern, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
				active_vertices = activated_vertices.size();
				average_active_vertices += (double) active_vertices / (double) number_trials_pattern;
				set<int> inside_active_vertices;
				set_intersection(Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.begin(),Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.end(),activated_vertices.begin(),activated_vertices.end(), inserter(inside_active_vertices, inside_active_vertices.begin()));
				average_active_vertices_inside += inside_active_vertices.size()/(double) number_trials_pattern;	
			}
			bool go_on = true;
			//cout << average_active_vertices_inside << " = average_active_vertices_inside111\n";
			if (average_active_vertices_inside < G.fidelity_inside*G.pattern_size) { // how should I test the inside fidelity
				go_on = false;
			}
			while (go_on) {
				G.act_threshold++;
				//cout << G.act_threshold << "= G.act_threshold\n";
				double average_active_vertices_new_threshold = 0;
				double average_active_vertices_inside_new_threshold = 0;
				for(int j=0; j< number_trials_pattern; j++ ){
					int rand_bootstrap_pattern= rand()% num_patterns*G.multi_association_factor;
					set<int> activated_vertices;
					percolationwithinhibition(activated_vertices, G, Patternsfirstlayernoisy[rand_bootstrap_pattern].verticesinpattern, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
					active_vertices = activated_vertices.size();
					average_active_vertices_new_threshold += (double) active_vertices / (double) number_trials_pattern;
					set<int> inside_active_vertices;
					set_intersection(Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.begin(),Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.end(),activated_vertices.begin(),activated_vertices.end(), inserter(inside_active_vertices, inside_active_vertices.begin()));
					average_active_vertices_inside_new_threshold += inside_active_vertices.size()/(double) number_trials_pattern;	
				}
				//cout << average_active_vertices_inside_new_threshold << " = average_active_vertices_inside_new_threshold\n";
				if (average_active_vertices_inside_new_threshold < G.fidelity_inside*G.pattern_size) {
					go_on = false;
					G.act_threshold--;
				}
				else {
					average_active_vertices = average_active_vertices_new_threshold;
					average_active_vertices_inside = average_active_vertices_inside_new_threshold;
				}
			}
			if (average_active_vertices < G.fidelity_outside) {
				average_active_vertices_inside_down_final = average_active_vertices_inside;
				average_active_vertices_down_final = average_active_vertices;
				act_threshold_down = G.act_threshold;
			}
			//cout << average_active_vertices_inside << " = average_active_vertices_inside222\n";
		}

		int num_patterns_up=num_patterns;
		int num_patterns_down= num_patterns/2;
		int num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
		while(num_patterns_up - num_patterns_down > 1) {
			G.act_threshold = base_act_threshold;
			//clog << "check num patterns " << num_patterns_mid<< "\n";
			average_active_vertices = 0;
			average_active_vertices_inside = 0;
			num_patterns=num_patterns_mid;
			for(int j=0; j< number_trials_pattern; j++ ){
				int rand_bootstrap_pattern= rand()% num_patterns*G.multi_association_factor;
				set<int> activated_vertices;
				percolationwithinhibition(activated_vertices, G, Patternsfirstlayernoisy[rand_bootstrap_pattern].verticesinpattern, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
				active_vertices = activated_vertices.size();
				average_active_vertices += (double) active_vertices / (double) number_trials_pattern; 
				set<int> inside_active_vertices;
				set_intersection(Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.begin(),Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.end(),activated_vertices.begin(),activated_vertices.end(), inserter(inside_active_vertices, inside_active_vertices.begin()));
				average_active_vertices_inside += inside_active_vertices.size()/(double) number_trials_pattern;
			}
			//cout << average_active_vertices_inside << " = average_active_vertices_inside111\n";
			bool go_on = true;
			if (average_active_vertices_inside < G.fidelity_inside*G.pattern_size) { // how should I test the inside fidelity
				go_on = false;
			}
			while (go_on) { 
				G.act_threshold++;
				//cout << G.act_threshold << "= G.act_threshold\n";
				double average_active_vertices_new_threshold = 0;
				double average_active_vertices_inside_new_threshold = 0;
				for(int j=0; j< number_trials_pattern; j++ ){
					int rand_bootstrap_pattern= rand()% num_patterns*G.multi_association_factor;
					set<int> activated_vertices;
					percolationwithinhibition(activated_vertices, G, Patternsfirstlayernoisy[rand_bootstrap_pattern].verticesinpattern, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
					active_vertices = activated_vertices.size();
					average_active_vertices_new_threshold += (double) active_vertices / (double) number_trials_pattern;
					set<int> inside_active_vertices;
					set_intersection(Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.begin(),Patternssecondlayer[rand_bootstrap_pattern/G.multi_association_factor].verticesinpattern.end(),activated_vertices.begin(),activated_vertices.end(), inserter(inside_active_vertices, inside_active_vertices.begin()));
					average_active_vertices_inside_new_threshold += inside_active_vertices.size()/(double) number_trials_pattern;	
				}
				//cout << average_active_vertices_inside_new_threshold << " = average_active_vertices_inside_new_threshold\n";
				if (average_active_vertices_inside_new_threshold < G.fidelity_inside*G.pattern_size) {
					go_on = false;
					G.act_threshold--;
				}
				else {
					average_active_vertices = average_active_vertices_new_threshold;
					average_active_vertices_inside = average_active_vertices_inside_new_threshold;
				}
			}
			if (average_active_vertices < G.fidelity_outside) {
			//if (average_active_vertices < G.fidelity_outside && average_active_vertices_inside >= (G.fidelity_inside-0.05)*G.pattern_size) {
				num_patterns_down= num_patterns_mid;
				num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
				average_active_vertices_down_final = average_active_vertices;
				average_active_vertices_inside_down_final = average_active_vertices_inside;
				act_threshold_down = G.act_threshold;
			}
			else{ 
				num_patterns_up= num_patterns_mid;
				num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
			}
			//cout << average_active_vertices_inside << " = average_active_vertices_inside222\n";
		}

		return make_tuple(num_patterns_down, average_active_vertices_down_final, average_active_vertices_inside_down_final, act_threshold_down);
}


// finds the optimal number of patterns that can be stored respecting fidelities averaging over many graph trials.
double percolation_test_parameters_with_inhibition_p_rec_given (int number_vertices, int pattern_size, int act_threshold, double p_rec, double p_aff, double fidelity_inside, double fidelity_outside, ofstream& file_percolation, int number_trials_new_graph, int number_trials_pattern, int multi_association_factor, int final_number_associations, int turn_strong_pattern_threshold, int inhibition_threshold, double false_activation_noise) {

	
	// Find maximal number of patterns that can be stored;
	double maximal_number_patterns_stored = 0;
	file_percolation << "maximal_number_patterns_under_fidelity = [";
	for (int i = 0; i < number_trials_new_graph; i++) {
		//cout << i << "\n";
		//////////////// Graph construction
		Twolayergraph G;
		graphconstructor(G, number_vertices, number_vertices, p_aff, p_rec);
		G.pattern_size = pattern_size;
		G.act_threshold = act_threshold;
		G.multi_association_factor = multi_association_factor;
		G.final_number_associations = final_number_associations;
		G.turn_strong_pattern_threshold = turn_strong_pattern_threshold;
		G.fidelity_outside = fidelity_outside;
		G.inhibition_threshold = inhibition_threshold;
		G.fidelity_inside = fidelity_inside;

		//cout << "Graph is built\n";

		////////////////// Pattern insertion

		vector<Pattern> Patternsfirstlayer;
		vector<Pattern> Patternssecondlayer;
		vector<Pattern> Patternsfirstlayernoisy;
		pattern_insertion(Patternsfirstlayer, Patternssecondlayer, Patternsfirstlayernoisy, G, final_number_associations, multi_association_factor,
			pattern_size, turn_strong_pattern_threshold, false_activation_noise);

		//cout << "Pattern insertion is complete\n";
		
		auto tuple_num_pattern_and_num_active_and_num_active_inside = maximal_number_patterns_without_breaching_fidelity_outside_single_graph_trial (Patternsfirstlayer, Patternssecondlayer, Patternsfirstlayernoisy, G, number_trials_pattern);
		int number_patterns_in_this_graph = get<0>(tuple_num_pattern_and_num_active_and_num_active_inside);
		double average_number_of_active_vertices_in_this_graph =  get<1>(tuple_num_pattern_and_num_active_and_num_active_inside);
		double average_number_of_active_vertices_inside = get<2>(tuple_num_pattern_and_num_active_and_num_active_inside);
		int optimal_activation_threshold = get<3>(tuple_num_pattern_and_num_active_and_num_active_inside);
		maximal_number_patterns_stored += ((double)number_patterns_in_this_graph)/number_trials_new_graph;
		file_percolation << number_patterns_in_this_graph << "," << average_number_of_active_vertices_in_this_graph << "," << average_number_of_active_vertices_inside << "," << optimal_activation_threshold;
		//clog << number_patterns_in_this_graph << "," << average_number_of_active_vertices_in_this_graph << "," << number_of_active_vertices_inside<< "\n";
		if (i< number_trials_new_graph -1){
			file_percolation << ";";
		}
		else{
			file_percolation << "];\n";
		}
	}
	return maximal_number_patterns_stored;
}



int main(){

	// Parameters
	int number_vertices = 5000; // Number of excitatory  vertices
	int multi_association_factor=16;
	int turn_strong_pattern_threshold = 1;
	int number_trials_pattern = 100;
	int number_trials_new_graph = 30;
	double fidelity_inside = 0.95;
	int pattern_size = 71;
	double fidelity_outside = 2*pattern_size;
	int inhibition_threshold = 5000;
	int final_number_associations = ceil(4*number_vertices*number_vertices/(pattern_size*pattern_size*multi_association_factor));

	int number_act_threshold_values = 5;
	vector<int> act_threshold_values(number_act_threshold_values);
	act_threshold_values = {8,10,12,14,16};

	double false_activation_noise_increment = 0.2;
	int number_false_activation_noise_values = floor(1/false_activation_noise_increment) + 1;

	double p_rec_increment = 0.05;
	int number_p_rec_values = floor(1/p_rec_increment) + 1;

	for (int noise_iterator = 0; noise_iterator < number_false_activation_noise_values; noise_iterator++) {
		double false_activation_noise = false_activation_noise_increment*noise_iterator;
		for (int act_threshold_iterator = 0; act_threshold_iterator < number_act_threshold_values; act_threshold_iterator++) {
			int act_threshold = act_threshold_values[act_threshold_iterator];
			cout << "false_activation_noise = " << false_activation_noise << "; act_threshold = " << act_threshold << "\n";
			omp_set_num_threads(4);
			#pragma omp parallel for
			for (int p_rec_iterator = 0; p_rec_iterator < number_p_rec_values; p_rec_iterator++){
				double p_rec = p_rec_increment*p_rec_iterator;
				int number_tries = 4000;
				int precision = 3;
				double p_aff = percolation_inside_p_aff_search_given_p_rec(pattern_size, act_threshold, number_tries, precision, p_rec, fidelity_inside);
				cout << p_aff << " = p_aff\n";
				int count_file = p_rec_iterator + act_threshold_iterator*number_p_rec_values + noise_iterator*number_act_threshold_values*number_p_rec_values;

				ofstream file_percolation;

				file_percolation.open("percolation_" + to_string(count_file)+".m");

				file_percolation << "number_vertices=" << number_vertices << "; \n";
				file_percolation << "pattern_size=" << pattern_size << ";\n";
				file_percolation << "act_threshold=" << act_threshold << ";\n";
				file_percolation << "fidelity_inside=" << fidelity_inside << ";\n";
				file_percolation << "fidelity_outside=" << fidelity_outside << ";\n";
				file_percolation << "number_associations=" << multi_association_factor << ";\n";
				file_percolation << "number_trials_new_graph = " << number_trials_new_graph << ";\n";
				file_percolation << "number_trials_pattern = " << number_trials_pattern << ";\n";
				file_percolation << "inhibition_threshold = " << inhibition_threshold << ";\n";
				file_percolation << "number_p_rec_values = " << number_p_rec_values << ";\n";
				file_percolation << "p_rec = " << p_rec << ";\n";
				file_percolation << "p_aff = " << p_aff << ";\n";
				file_percolation << "number_false_activation_noise_values = " << number_false_activation_noise_values << ";\n";
				file_percolation << "false_activation_noise = " << false_activation_noise << ";\n";
				file_percolation << "p_rec_iterator = " << p_rec_iterator << ";\n";
				file_percolation << "noise_iterator = " << noise_iterator << ";\n";

				double maximal_number_patterns_under_fidelity = percolation_test_parameters_with_inhibition_p_rec_given (number_vertices, pattern_size, act_threshold, p_rec, p_aff, fidelity_inside, fidelity_outside, file_percolation, number_trials_new_graph, number_trials_pattern, multi_association_factor, final_number_associations, turn_strong_pattern_threshold, inhibition_threshold, false_activation_noise);

				file_percolation << "maximal_number_patterns_under_fidelity_average = " << maximal_number_patterns_under_fidelity << ";\n";

				file_percolation.close();

			}

		}
	}


	/*
	double false_activation_noise = 0.2;

	double p_rec = 0.3;
	int number_tries = 4000;
	int precision = 3;
	double p_aff = percolation_inside_p_aff_search_given_p_rec(pattern_size, act_threshold, number_tries, precision, p_rec, fidelity_inside);
	cout << p_aff << "\n";

	int count_file = 2;
	
	int final_number_associations = ceil(4*number_vertices*number_vertices/(pattern_size*pattern_size*multi_association_factor));

	ofstream file_percolation;

	file_percolation.open("percolation_" + to_string(count_file)+".m");
	
	file_percolation << "number_vertices=" << number_vertices << "; \n";
	file_percolation << "pattern_size=" << pattern_size << ";\n";
	file_percolation << "act_threshold=" << act_threshold << ";\n";
	file_percolation << "fidelity_inside=" << fidelity_inside << ";\n";
	file_percolation << "fidelity_outside=" << fidelity_outside << ";\n";
	file_percolation << "number_associations=" << multi_association_factor << ";\n";
	file_percolation << "number_trials_new_graph = " << number_trials_new_graph << ";\n";
	file_percolation << "number_trials_pattern = " << number_trials_pattern << ";\n";
	file_percolation << "inhibition_threshold = " << inhibition_threshold << ";\n";
	file_percolation << "p_rec = " << p_rec << ";\n";
	file_percolation << "p_aff = " << p_aff << ";\n";
	file_percolation << "false_activation_noise = " << false_activation_noise << ";\n";

	double maximal_number_patterns_under_fidelity = percolation_test_parameters_with_inhibition_p_rec_given (number_vertices, pattern_size, act_threshold, p_rec, p_aff, fidelity_inside, fidelity_outside, file_percolation, number_trials_new_graph, number_trials_pattern, multi_association_factor, final_number_associations, turn_strong_pattern_threshold, inhibition_threshold, false_activation_noise);

	file_percolation << "maximal_number_patterns_under_fidelity_average = " << maximal_number_patterns_under_fidelity << ";\n";

	file_percolation.close();
	*/
	return 1;
}

///////////////////////// TEST CODES ///////////////////////////////////////////////////////////////






	//////////// Test edge addition in graph ///////////////////////////////////////
	/*int number_recurrent_edges = 0;
	for (int i = 0; i < number_vertices; i++) {
		for (int j = 0; j < number_vertices; j++){
			if(G.Recurrentedges[i][j].is_present) {
				number_recurrent_edges++;
			}
		}
	}
	double fraction_recurrent_edges = ((double)number_recurrent_edges)/((double)number_vertices)/((double)number_vertices-1.0);

	cout << fraction_recurrent_edges << "\n";*/

	///////////// Test strong edges in graph (pattern insertion) ////////////////////////////////////////////
	/*int number_strong_afferent_edges = 0;
	for (int i = 0; i < number_vertices; i++) {
		for (int j = 0; j < number_vertices; j++){
			if(G.Afferentedges[i][j].is_present && G.Afferentedges[i][j].patternthreshold == 1 && G.Afferentedges[i][j].patternturnstrong != -1) {
				number_strong_afferent_edges++;
			}
		}
	}
	double fraction_strong_afferent_edges = ((double)number_strong_afferent_edges)/((double)number_vertices)/((double)number_vertices);

	cout << fraction_strong_afferent_edges << "\n";
	cout << number_strong_afferent_edges << "\n";*/








	////////////////// Test percolation with inhibition ////////////////////////////////

	/*int test_number_associations = (log(2)+0.5)*number_vertices*number_vertices/(pattern_size*pattern_size);
	int test_pattern;

	// Random devices - to pick pattern at random
    std::random_device rd;
    std::mt19937 gen(rd());
	uniform_int_distribution<> randompatternchoice(0, test_number_associations - 1);
	test_pattern = randompatternchoice(gen);
	set<int> activated_vertices;

	cout << "Percolation starting\n"; 
	percolationwithinhibition(activated_vertices, G, Patternsfirstlayer[test_pattern].verticesinpattern, act_threshold, inhibition_threshold, test_number_associations, turn_strong_pattern_threshold);

	int number_activated_vertices = activated_vertices.size();
	cout << number_activated_vertices << "\n";*/








	////////////////////// Test density search /////////////////////////////////////////////

	//p_aff = percolation_inside_p_aff_search_given_p_rec(pattern_size, act_threshold, 1000, 4, 1, 0.95);
	//cout << p_aff << "\n";
	//p_rec = percolation_inside_p_rec_search_given_p_aff(pattern_size, act_threshold, 1000, 4, 0.1, 0.95);
	//cout << p_rec << "\n";
