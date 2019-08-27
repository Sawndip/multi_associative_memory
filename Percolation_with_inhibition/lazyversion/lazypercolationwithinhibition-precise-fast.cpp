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

struct Edge { // Edge structure
	long long patternturnstrong; // if patternturnstrong < number associations used than edge is strong in the graph
	int patternthreshold; // an edge turns strong iff it appeared on at least a given number of patterns
	bool is_present; // whether edge belongs to the graph or not
	int target_vertex;
};

struct Vertex { // Vertex structure
	int voltage; // active if it reaches activation threshold. If above act_threshold then the neuron has already spiked (and will not be considered further)
	int activation_time; // stores percolation round when vertex turned active
	set<int> neighbours; // asymmetric - vertices that are not in pattern actually have every possible (that was tried to be added by a pattern) neighbour in this set. One has to check whether edges are actually present
};

struct Twolayerlazygraph { // Graph structure
	int verticesfirstlayer;
	int verticessecondlayer;
	int pattern_size;
	int inhibition_threshold;
	vector<Vertex> Firstlayer;
	vector<Vertex> Secondlayer;
	vector<vector<Edge> > Patternafferentedges;
	vector<vector<Edge> > Patternrecurrentedges; // currently asymmetric - edges that are not from the pattern to the outside are treated differently - is crucial to check whether they are actually present
	//should be kept ordered (according to target vertex) to speed program up
	//other parameters that are necessary for percolation:
	int multi_association_factor;
	long long final_number_associations;
	int turn_strong_pattern_threshold;
	int act_threshold;
	double fidelity_outside;
	double p_rec;
	double p_aff;
};


// computes probability of binomial distribution Bin(k,p) being i; 
double binomial_probability(int k, double p, int i) {
	double temp = lgamma(k + 1.0);
	temp -=  lgamma(i + 1.0) + lgamma(k-i + 1.0);
	temp += i*log(p) + (k-i)*log(1-p);
	return exp(temp);
}


// This builds a two layer graph with afferent connection p_aff and recurrent connection density p_rec. It should be called to build a new graph
// It only adds recurrent connections from the pattern considered to the outside. Also the first layer should consist only of the original pattern.
void lazygraphconstructor(Twolayerlazygraph& G, int verticesfirstlayer, int verticessecondlayer, int pattern_size, double p_aff, double p_rec) {
	// define layer sizes
	G.verticesfirstlayer = verticesfirstlayer;
	G.verticessecondlayer = verticessecondlayer;
	G.pattern_size = pattern_size;
	G.p_rec = p_rec;
	G.p_aff = p_aff;

	// Random Devices
	std::random_device rd;
	std::mt19937 generator(rd());
	bernoulli_distribution afferentdistribution(p_aff);
	bernoulli_distribution recurrentdistribution(p_rec);

	// Build vertices in first layer pattern && afferent edges
	for (int i = 0; i < verticesfirstlayer; i++) {
		vector<Edge> vertex_i_edges;
		Vertex v;
		v.voltage = 0;
		v.activation_time = -1;
		for (int j = 0; j < verticessecondlayer; j++) {//define (afferent) neighbours of v!
			Edge e;
			e.patternturnstrong = -1;
			e.patternthreshold = 0;
			e.is_present = false;
			e.target_vertex = j;
			if (afferentdistribution(generator)) {
				v.neighbours.insert(j);
				e.is_present = true;
			}
			vertex_i_edges.push_back(e);
		}
		G.Firstlayer.push_back(v);
		G.Patternafferentedges.push_back(vertex_i_edges);
	}

	// Build vertices in second layer pattern && recurrent edges
	for (int i = 0; i < verticessecondlayer; i++) {
		Vertex v;
		v.voltage = 0;
		v.activation_time = -1;
		vector<Edge> vertex_i_edges;
		if (i < pattern_size) {
			for (int j = 0; j < verticessecondlayer; j++) {//define recurrent neighbours of v!
				Edge e;
				e.patternturnstrong = -1;
				e.patternthreshold = 0;
				e.is_present = false;
				e.target_vertex = j;
				if (recurrentdistribution(generator) && (i != j)) {
					v.neighbours.insert(j);
					e.is_present = true;
				}
				vertex_i_edges.push_back(e);
			}
		}
		G.Patternrecurrentedges.push_back(vertex_i_edges);
		G.Secondlayer.push_back(v);
	}
}

// This inserts patterns at a given two layer graph G, associating many to one. It should be called to insert patterns after a new graph was built.
// program assumes first pattern stored at 1-pattern_size
// program does not sample all the patterns but instead only the ones that have non-empty intersection with the source and target patterns analyzed.
// It does so by considering that the number of patterns with intersection i is concentrated around the expected value!
void pattern_insertion(Twolayerlazygraph& G, long long number_patterns_rec_layer, int multi_association_factor, int pattern_size, int turn_strong_pattern_threshold, double p_rec){

	// Extract layer sizes;
	int verticesfirstlayer = G.verticesfirstlayer;
	int verticessecondlayer = G.verticessecondlayer;
	int number_vertices = verticessecondlayer;

	// Random devices;
	std::random_device rd;
	std::mt19937_64 gen(rd());
	uniform_int_distribution<> uniform_vertex_graph(0, number_vertices - 1);
	uniform_int_distribution<> uniform_vertex_pattern(0, pattern_size - 1);
	uniform_int_distribution<long long> uniform_pattern_number(1, number_patterns_rec_layer - 1);//pattern number will go from 1 
	bernoulli_distribution recurrentdistribution(p_rec);

	set<int> source_pattern;
	set<int> target_pattern;
	for (int i = 0; i < pattern_size; i++) {
		source_pattern.insert(i);
		target_pattern.insert(i);
	}

	// turning edges inside edges strong
	for(int i = 0; i < pattern_size; i++) {
		for(int j = 0; j < pattern_size; j++) {
			// afferent edges
			if(G.Patternafferentedges[i][j].is_present && G.Patternafferentedges[i][j].patternthreshold < turn_strong_pattern_threshold) {
						G.Patternafferentedges[i][j].patternthreshold++;
						G.Patternafferentedges[i][j].patternturnstrong = 0;
			}
			if (G.Patternrecurrentedges[i][j].is_present && G.Patternrecurrentedges[i][j].patternthreshold < turn_strong_pattern_threshold) {
				G.Patternrecurrentedges[i][j].patternthreshold++;
				G.Patternrecurrentedges[i][j].patternturnstrong = 0;
			}
		}
	}
	// calculates expected number of patterns with intersection i with the source pattern!
	long long total_number_relevant_patterns = 0;
	long long number_patterns_with_intersection_first_layer[pattern_size + 1];
	long long number_patterns_without_intersection_first_layer_with_intersection_second_layer[pattern_size + 1];
	// check if c++ is precise in this computation!
	long long number_patterns_outside_rec_layer = number_patterns_rec_layer - 1;
	double probability_vertex_is_in_intersection_first_layer = 1.0-exp(multi_association_factor*log(1.0 - ((double) pattern_size)/((double) number_vertices)));
	double probability_vertex_is_in_intersection_second_layer = ((double) pattern_size)/((double) number_vertices);
	//cout << "probability_vertex_is_in_intersection_first_layer = " << probability_vertex_is_in_intersection_first_layer << "\n";
	for (int i = 0; i <= pattern_size; i++) {
		double probability_intersection_size_i = binomial_probability(pattern_size, probability_vertex_is_in_intersection_first_layer, i);
		number_patterns_with_intersection_first_layer[i] = floor(number_patterns_outside_rec_layer*probability_intersection_size_i);
		if (i != 0 ) {
			total_number_relevant_patterns += number_patterns_with_intersection_first_layer[i];
		}
		else {
			for (int j = 1; j <= pattern_size; j++) {
				double probability_intersection_size_j = binomial_probability(pattern_size, probability_vertex_is_in_intersection_second_layer, j);
				number_patterns_without_intersection_first_layer_with_intersection_second_layer[j] = floor(number_patterns_with_intersection_first_layer[0]*probability_intersection_size_j);
				total_number_relevant_patterns += number_patterns_without_intersection_first_layer_with_intersection_second_layer[j];
				//cout << "number_patterns_rec_layer = " << number_patterns_rec_layer << "\n";
				//cout << "probability_intersection_size_j = " << probability_intersection_size_j <<  "\n";
				//cout <<  "number_patterns_with_intersection_second_layer["<< j <<"] = " << number_patterns_without_intersection_first_layer_with_intersection_second_layer[j] << "\n";
			}
		}
		//cout << "number_patterns_rec_layer = " << number_patterns_rec_layer << "\n";
		//cout << "probability_intersection_size_i = " << probability_intersection_size_i <<  "\n";
		//cout <<  "number_patterns_with_intersection_first_layer["<< i <<"] = " << number_patterns_with_intersection_first_layer[i] << "\n";
	}

	//define order which patterns were inserted -> important for binary search
	vector<long long> order_patterns_were_inserted;
	set<long long> numbers_already_used;
	long long size_numbers_already_used = 0;
	long long number_small_indices = 0;
	while (size_numbers_already_used < total_number_relevant_patterns) {
		long long new_pattern_number = uniform_pattern_number(gen);
 		//cout << "shit\n";
		long long old_size = size_numbers_already_used;
		numbers_already_used.insert(new_pattern_number);
		size_numbers_already_used = numbers_already_used.size();
		if (size_numbers_already_used > old_size) {
			order_patterns_were_inserted.push_back(new_pattern_number);
			//cout << new_pattern_number << "\n";
			if (new_pattern_number < 220000) {
				number_small_indices++;
			}
		}
	}
	//cout << number_patterns_rec_layer <<"! \n";
	//cout << total_number_relevant_patterns << "\n";
	//cout << number_small_indices << "\n";

	//turns edges outside edges strong
	long long index_relevant_patterns = 0;
	for (int i = 0; i <= pattern_size; i++) {
		//cout << "intersection size i = " << i << "\n";
		if (i != 0) {
			for (long long j = 0; j < number_patterns_with_intersection_first_layer[i]; j++) {
								//insert the next relevant pattern

				//generate the intersection in first layer
				set<int> first_layer_intersection;
				while (first_layer_intersection.size() < i) {
					first_layer_intersection.insert(uniform_vertex_pattern(gen));
				}

				//generate the corresponding pattern in second layer
				set<int> new_second_layer_pattern;
				while (new_second_layer_pattern.size() < pattern_size) {
					int vertex_index = uniform_vertex_graph(gen);
					new_second_layer_pattern.insert(vertex_index);
				}

				//turn afferent edges strong
				for (set<int>::iterator first_layer_intersection_iterator = first_layer_intersection.begin(); first_layer_intersection_iterator != first_layer_intersection.end(); first_layer_intersection_iterator++) {
					for(set<int>::iterator new_second_layer_pattern_iterator = new_second_layer_pattern.begin(); new_second_layer_pattern_iterator != new_second_layer_pattern.end(); new_second_layer_pattern_iterator++) {
						if(G.Patternafferentedges[*first_layer_intersection_iterator][*new_second_layer_pattern_iterator].is_present) {
							if (G.Patternafferentedges[*first_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternthreshold < turn_strong_pattern_threshold) {
								G.Patternafferentedges[*first_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternthreshold++;
							}
							if (G.Patternafferentedges[*first_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong == -1) {
								G.Patternafferentedges[*first_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns]; //should have patterns in "random" order
							}
							else if (G.Patternafferentedges[*first_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong > order_patterns_were_inserted[index_relevant_patterns]) {
								G.Patternafferentedges[*first_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns]; //should have patterns in "random" order ///** DOES NOT WORK IF TURN_STRONG_PATTERN > 1
							}
						}
					}
				}

				// OLD VERSION - NOT SO PRECISE - THROWS IMPORTANT EDGES AWAY
				// get intersection with second layer pattern!
				set<int> second_layer_intersection;
				//intersection vertices are those that will have additional strong connections to outside vertices!
				set_intersection(new_second_layer_pattern.begin(),new_second_layer_pattern.end(),target_pattern.begin(),target_pattern.end(), inserter(second_layer_intersection, second_layer_intersection.begin()));

				//turn recurrent edges strong
				for (set<int>::iterator second_layer_intersection_iterator = second_layer_intersection.begin(); second_layer_intersection_iterator != second_layer_intersection.end(); second_layer_intersection_iterator++) {
					for(set<int>::iterator new_second_layer_pattern_iterator = new_second_layer_pattern.begin(); new_second_layer_pattern_iterator != new_second_layer_pattern.end(); new_second_layer_pattern_iterator++) {
						if(G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].is_present) {
							if (G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternthreshold < turn_strong_pattern_threshold) {
								G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternthreshold++;
							}
							if (G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong == -1) {
								G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns]; //should have patterns in "random" order
							}
							else if (G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong > order_patterns_were_inserted[index_relevant_patterns]) {
								G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns]; //should have patterns in "random" order
							}
						}
					}
				}


				// ALTERNATIVE - ADD ALL RECURRENT EDGES IN CURRENTLY GENERATE PATTERN - MAKES PROGRAM MORE PRECISE

				set<int> second_layer_not_intersection;
				set_difference(new_second_layer_pattern.begin(), new_second_layer_pattern.end(), second_layer_intersection.begin(), second_layer_intersection.end(), inserter(second_layer_not_intersection, second_layer_not_intersection.begin()));

				for(set<int>::iterator second_layer_not_intersection_source_iterator = second_layer_not_intersection.begin(); second_layer_not_intersection_source_iterator != second_layer_not_intersection.end(); second_layer_not_intersection_source_iterator++) {
					for(set<int>::iterator second_layer_not_intersection_target_iterator = second_layer_not_intersection.begin(); second_layer_not_intersection_target_iterator != second_layer_not_intersection.end(); second_layer_not_intersection_target_iterator++) {
						if (*second_layer_not_intersection_source_iterator != *second_layer_not_intersection_target_iterator) {
							int target_index_down = 0;
							int target_index_up = G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].size() - 1;
							//cout << target_index_up << " = target_index_up\n";
							int target_index_middle = (target_index_down+target_index_up)/2;
							while(target_index_up > target_index_down) {
								if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_middle].target_vertex < *second_layer_not_intersection_target_iterator) {
									target_index_down = target_index_middle + 1;
									target_index_middle = (target_index_down+target_index_up)/2;
								}
								else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_middle].target_vertex > *second_layer_not_intersection_target_iterator) {
									target_index_up = target_index_middle - 1;
									target_index_middle = (target_index_down+target_index_up)/2;
								}
								else {
									target_index_up = target_index_middle;
									target_index_down = target_index_middle;
								}
							}
							//cout << target_index_middle << " = target_index_middle\n";

							if (target_index_up >= 0 && G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].target_vertex == *second_layer_not_intersection_target_iterator) {

								if(G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternthreshold < turn_strong_pattern_threshold) {
									G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternthreshold++;
								}
								if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong == -1) {
									G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns];
								}
								else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong > order_patterns_were_inserted[index_relevant_patterns]) {
									G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns];
								}
							}

							else {
								
								G.Secondlayer[*second_layer_not_intersection_source_iterator].neighbours.insert(*second_layer_not_intersection_target_iterator);
								Edge e;
								e.patternthreshold = 1;
								e.patternturnstrong = order_patterns_were_inserted[index_relevant_patterns];
								e.is_present = false;
								if(recurrentdistribution(gen)) {
									e.is_present = true;
								}
								e.target_vertex = *second_layer_not_intersection_target_iterator;
								vector<Edge>::iterator start_of_edges = G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].begin();
								if (target_index_up == -1) {
									G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].push_back(e);
								}
								else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].size() == 1) {
									if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].target_vertex < *second_layer_not_intersection_target_iterator) {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].push_back(e);
									}
									else {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].insert(start_of_edges+target_index_up, e);
									}
								}
								else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].target_vertex < *second_layer_not_intersection_target_iterator) {
									G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].insert(start_of_edges+target_index_up + 1, e);
								}
								else {
									G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].insert(start_of_edges+target_index_up, e);
								}
							}
						}
					}
				}

				index_relevant_patterns++;
			}
		}

		else {
			for (int s = 1; s <= pattern_size; s++) {
				//cout << "intersection size s = " << s << "\n";
				for (long long j = 1; j < number_patterns_without_intersection_first_layer_with_intersection_second_layer[s]; j++) {
					//insert the next relevant pattern

					// there are no strong afferent edges!
					//generate the intersection in second layer
					set<int> second_layer_intersection;
					while (second_layer_intersection.size() < s) {
						second_layer_intersection.insert(uniform_vertex_pattern(gen));
					}

					//generate the part of the pattern not in the intersection
					set<int> new_second_layer_pattern;
					while (new_second_layer_pattern.size() < pattern_size - s) {
						int vertex_index = uniform_vertex_graph(gen);
						new_second_layer_pattern.insert(vertex_index);
					}

					// OLD VERSION - NOT SO PRECISE - THROWS IMPORTANT EDGES AWAY
					//turn recurrent edges strong
					for (set<int>::iterator second_layer_intersection_iterator = second_layer_intersection.begin(); second_layer_intersection_iterator != second_layer_intersection.end(); second_layer_intersection_iterator++) {
						for(set<int>::iterator new_second_layer_pattern_iterator = new_second_layer_pattern.begin(); new_second_layer_pattern_iterator != new_second_layer_pattern.end(); new_second_layer_pattern_iterator++) {
							if(G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].is_present) {
								if ( G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternthreshold < turn_strong_pattern_threshold) {
									G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternthreshold++;
								}
								if (G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong == -1) {
									G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns]; //should have patterns in "random" order
								}
								else if (G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong > order_patterns_were_inserted[index_relevant_patterns]) {
									G.Patternrecurrentedges[*second_layer_intersection_iterator][*new_second_layer_pattern_iterator].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns]; //should have patterns in "random" order
								}
							}
						}
					}

					// ALTERNATIVE - ADD ALL RECURRENT EDGES IN CURRENTLY GENERATE PATTERN - MAKES PROGRAM MORE PRECISE

					set<int> second_layer_not_intersection;
					while (second_layer_not_intersection.size() < pattern_size - s) {
						int new_vertex = uniform_vertex_graph(gen);
						if (new_vertex >= pattern_size) {
							second_layer_not_intersection.insert(new_vertex);
						}
					}

					for(set<int>::iterator second_layer_not_intersection_source_iterator = second_layer_not_intersection.begin(); second_layer_not_intersection_source_iterator != second_layer_not_intersection.end(); second_layer_not_intersection_source_iterator++) {
						for(set<int>::iterator second_layer_not_intersection_target_iterator = second_layer_not_intersection.begin(); second_layer_not_intersection_target_iterator != second_layer_not_intersection.end(); second_layer_not_intersection_target_iterator++) {
							if (*second_layer_not_intersection_source_iterator != *second_layer_not_intersection_target_iterator) {
								int target_index_down = 0;
								int target_index_up = G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].size() - 1;
								//cout << target_index_up << " = target_index_up\n";
								int target_index_middle = (target_index_down+target_index_up)/2;
								while(target_index_up > target_index_down) {
									if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_middle].target_vertex < *second_layer_not_intersection_target_iterator) {
										target_index_down = target_index_middle + 1;
										target_index_middle = (target_index_down+target_index_up)/2;
									}
									else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_middle].target_vertex > *second_layer_not_intersection_target_iterator) {
										target_index_up = target_index_middle - 1;
										target_index_middle = (target_index_down+target_index_up)/2;
									}
									else {
										target_index_up = target_index_middle;
										target_index_down = target_index_middle;
									}
								}
								//cout << target_index_middle << " = target_index_middle\n";

								if (target_index_up >= 0 && G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].target_vertex == *second_layer_not_intersection_target_iterator) {

									if(G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternthreshold < turn_strong_pattern_threshold) {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternthreshold++;
									}
									if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong == -1) {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns];
									}
									else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong > order_patterns_were_inserted[index_relevant_patterns]) {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].patternturnstrong = order_patterns_were_inserted[index_relevant_patterns];
									}
								}

								else {
									
									G.Secondlayer[*second_layer_not_intersection_source_iterator].neighbours.insert(*second_layer_not_intersection_target_iterator);
									Edge e;
									e.patternthreshold = 1;
									e.patternturnstrong = order_patterns_were_inserted[index_relevant_patterns];
									e.is_present = false;
									if(recurrentdistribution(gen)) {
										e.is_present = true;
									}
									e.target_vertex = *second_layer_not_intersection_target_iterator;
									vector<Edge>::iterator start_of_edges = G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].begin();
									if (target_index_up == -1) {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].push_back(e);
									}
									else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].size() == 1) {
										if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].target_vertex < *second_layer_not_intersection_target_iterator) {
											G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].push_back(e);
										}
										else {
											G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].insert(start_of_edges+target_index_up, e);
										}
									}
									else if (G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator][target_index_up].target_vertex < *second_layer_not_intersection_target_iterator) {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].insert(start_of_edges+target_index_up + 1, e);
									}
									else {
										G.Patternrecurrentedges[*second_layer_not_intersection_source_iterator].insert(start_of_edges+target_index_up, e);
									}
								}
							}
						}
					}


					index_relevant_patterns++;
				}
			}
		}
	}
}

// Given a graph where patterns were inserted runs bootstrap percolation on the strong edge set of G.
// The bootstrap is the set of activated vertices in the first layer.
// activated_vertices is the set of vertices activated in second layer.
// inhibition_threshold is the number of vertices where inhibition kills the bootstrap.
// number_patterns_rec_layer is the number of patterns inserted in second layer.
// turn_strong_pattern_threshold is the number of times an edge should be in a pattern to be strong.
// Percolation stops in the round where inhibition threshold vertices turned active.
// As there are only recurrent edges from the pattern to the outside, the program generates then at random (with the correct probability) for every outside vertex that got active
void percolationwithinhibition(set<int>& activated_vertices, Twolayerlazygraph& G, set<int> bootstrap, int act_threshold, int inhibition_threshold, long long number_patterns_rec_layer, int turn_strong_pattern_threshold) {

	// probability that recurrent edges are strong (not stored in the graph)!
	double p_rec = G.p_rec;
	//cout << p_rec << " = p_rec\n";
	int number_vertices = G.verticessecondlayer;
	int pattern_size = G.pattern_size;
	int multi_association_factor = G.multi_association_factor;
	double probability_vertex_is_in_intersection_first_layer = 1.0-exp(multi_association_factor*log(1.0 - ((double) pattern_size)/((double) number_vertices)));
	double probability_vertex_is_in_intersection_second_layer = ((double) pattern_size)/((double) number_vertices);
	double probability_intersection_first_layer_size_0 = binomial_probability(pattern_size, probability_vertex_is_in_intersection_first_layer, 0);
	double probability_intersection_second_layer_size_0 = binomial_probability(pattern_size, probability_vertex_is_in_intersection_second_layer, 0);
	long long number_patterns_rec_layer_not_considered = number_patterns_rec_layer*probability_intersection_first_layer_size_0*probability_intersection_second_layer_size_0;

	double p_strong_rec = 1.0-exp(number_patterns_rec_layer_not_considered*log(1.0 - ((double) G.pattern_size*(G.pattern_size-1))/((double) G.verticessecondlayer*(G.verticessecondlayer-1))));
	//cout << p_strong_rec << " = p_strong_rec\n";
	// only true if pattern turn strong threshold is 1

	// Random Devices
	std::random_device rd;
	std::mt19937_64 generator(rd());
	//bernoulli_distribution recurrentstrongdistribution(p_strong_rec);
	//bernoulli_distribution recurrentdistribution(p_rec);
	uniform_int_distribution<> uniform_vertex_graph(0, number_vertices-1);
	binomial_distribution<> binomial_dist(number_vertices - 1, p_strong_rec*p_rec);

	int activated_vertices_size = 0;
	int final_round = -1; //final round of percolation - inhibition kicks in
	queue<int> verticesinqueue;
	for (set<int>::iterator vertex_iterator = bootstrap.begin(); vertex_iterator != bootstrap.end(); vertex_iterator++) {
		verticesinqueue.push(*vertex_iterator);
		verticesinqueue.push(1); // first round only has first layer vertices. Second round onwards has only second layer vertices. Round is stored to stop program at the correct round!
		G.Firstlayer[*vertex_iterator].voltage = act_threshold + 1;
		G.Firstlayer[*vertex_iterator].activation_time = 1;
	}
	while(!verticesinqueue.empty()){
		int vertex_index = verticesinqueue.front();
		verticesinqueue.pop();
		//cout << "vertex_index = " << vertex_index << "\n";
		int vertex_activation_round = verticesinqueue.front();
		verticesinqueue.pop();
		//cout << "vertex_activation_round = " << vertex_activation_round << "\n";
		if (vertex_activation_round == 1) { // first round is only in first layer
			for (set<int>::iterator i = G.Firstlayer[vertex_index].neighbours.begin(); i != G.Firstlayer[vertex_index].neighbours.end(); i++) {
				int target_vertex = *i;
				if (G.Patternafferentedges[vertex_index][target_vertex].patternturnstrong < number_patterns_rec_layer && G.Patternafferentedges[vertex_index][target_vertex].patternthreshold == turn_strong_pattern_threshold){// unnecessary to check whether edge is present
					if (G.Secondlayer[target_vertex].voltage < act_threshold) {
						G.Secondlayer[target_vertex].voltage++;
						//cout << "target_vertex = " << target_vertex << "\n";
						//cout << "G.Secondlayer[target_vertex].voltage = " << G.Secondlayer[target_vertex].voltage << "\n";
					}
					if (G.Secondlayer[target_vertex].voltage == act_threshold) {
						G.Secondlayer[target_vertex].voltage++;
						G.Secondlayer[target_vertex].activation_time = G.Firstlayer[vertex_index].activation_time + 1;
						verticesinqueue.push(target_vertex);
						verticesinqueue.push(G.Secondlayer[target_vertex].activation_time);
						//cout << "target_vertex = " << target_vertex << "\n";
						//cout << "G.Secondlayer[target_vertex].activation_time = " << G.Secondlayer[target_vertex].activation_time << "\n";
					}
				}
			}
		}
		else { // other rounds are only in second layer
			if (activated_vertices_size < inhibition_threshold - 1) {
				activated_vertices.insert(vertex_index); // should I insert once it reaches right threshold potential or when it will send the spike? - Current is only when it sends spike
				activated_vertices_size++;
				if (vertex_index >= pattern_size) {
					int neighbourhood_size = binomial_dist(generator);
					int old_size = 0;
					set<int> neighbourhood;
					while (neighbourhood.size() < neighbourhood_size) {
						int target_vertex = uniform_vertex_graph(generator);
						//set<int>::iterator target_vertex_index = known_neighbourhood.find(target_vertex);
						if(target_vertex != vertex_index) {
							old_size = neighbourhood.size();
							neighbourhood.insert(target_vertex);
							if(old_size < neighbourhood.size()){ //  to sample without replacement
								if (G.Secondlayer[target_vertex].voltage < act_threshold) {
									G.Secondlayer[target_vertex].voltage++;
									//cout << "target_vertex = " << target_vertex << "\n";
									//cout << "G.Secondlayer[target_vertex].voltage = " << G.Secondlayer[target_vertex].voltage << "\n";
								}
								if (G.Secondlayer[target_vertex].voltage == act_threshold) {
									G.Secondlayer[target_vertex].voltage++;
									G.Secondlayer[target_vertex].activation_time = G.Secondlayer[vertex_index].activation_time + 1;
									verticesinqueue.push(target_vertex);
									verticesinqueue.push(G.Secondlayer[target_vertex].activation_time);
									//cout << "target_vertex = " << target_vertex << "\n";
									//cout << "G.Secondlayer[target_vertex].activation_time = " << G.Secondlayer[target_vertex].activation_time << "\n";
								}
							}
						}
					}
					//set<int> known_neighbourhood;
					for (vector<Edge>::iterator vertex_neighbourhood_iterator = G.Patternrecurrentedges[vertex_index].begin(); vertex_neighbourhood_iterator != G.Patternrecurrentedges[vertex_index].end(); vertex_neighbourhood_iterator++) {
						if((*vertex_neighbourhood_iterator).is_present) {
							int target_vertex = (*vertex_neighbourhood_iterator).target_vertex;
							//known_neighbourhood.insert(target_vertex);
							set<int>::iterator target_vertex_index = neighbourhood.find(target_vertex);
							if (target_vertex_index == neighbourhood.end()) { // edges only send spike once!
								if ((*vertex_neighbourhood_iterator).patternturnstrong < number_patterns_rec_layer && (*vertex_neighbourhood_iterator).patternthreshold == turn_strong_pattern_threshold) {
									if (G.Secondlayer[target_vertex].voltage < act_threshold) {
										G.Secondlayer[target_vertex].voltage++;
									}
									if (G.Secondlayer[target_vertex].voltage == act_threshold) {
										G.Secondlayer[target_vertex].voltage++;
										G.Secondlayer[target_vertex].activation_time = G.Secondlayer[vertex_index].activation_time + 1;
										verticesinqueue.push(target_vertex);
										//cout << target_vertex << "heheh\n";
										verticesinqueue.push(G.Secondlayer[target_vertex].activation_time);
									}
								}
							}
						}
					}
				}
				else {
					//int neighbourhood_size = G.Secondlayer[vertex_index].neighbours.size();
					for (set<int>::iterator i = G.Secondlayer[vertex_index].neighbours.begin(); i != G.Secondlayer[vertex_index].neighbours.end(); i++) {
						int target_vertex = *i;
						if (G.Patternrecurrentedges[vertex_index][target_vertex].is_present && G.Patternrecurrentedges[vertex_index][target_vertex].patternturnstrong < number_patterns_rec_layer && G.Patternrecurrentedges[vertex_index][target_vertex].patternthreshold == turn_strong_pattern_threshold) {// unnecessary to check whether edge is present
							if (G.Secondlayer[target_vertex].voltage < act_threshold) {
								G.Secondlayer[target_vertex].voltage++;
								//cout << "target_vertex = " << target_vertex << "\n";
								//cout << "G.Secondlayer[target_vertex].voltage = " << G.Secondlayer[target_vertex].voltage << "\n";
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
int percolation_inside_active_set_size(int pattern_size, set<int>& bootstrap, int act_threshold, double p_aff, double p_rec) {

	// Graph construction
	Twolayerlazygraph G;
	lazygraphconstructor(G, pattern_size, pattern_size, pattern_size, p_aff, p_rec);

	//cout << "Graph is built\n" ;
	// Pattern insertion
	pattern_insertion(G, 1, 1, pattern_size, 1, p_rec);

	//cout << "Patterns are built\n" ;
	// Percolation - Pattern recall
	set<int> activated_vertices;
	int inhibition_threshold = pattern_size;

	percolationwithinhibition(activated_vertices, G, bootstrap, act_threshold, inhibition_threshold, 1, 1);
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
			set<int> bootstrap;
			for (int i = 0; i < pattern_size; i++) {
				bootstrap.insert(i);
			}
			for(int i=0; i < number_tries; i++){
				int activity_inside = percolation_inside_active_set_size(pattern_size,bootstrap, act_threshold, p_aff, p_rec);
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
			set<int> bootstrap;
			for (int i = 0; i < pattern_size; i++) {
				bootstrap.insert(i);
			}

			for(int i=0; i < number_tries; i++){
				int activity_inside = percolation_inside_active_set_size(pattern_size, bootstrap, act_threshold, p_aff, p_rec);
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
long long maximal_number_patterns_without_breaching_fidelity_outside_single_graph_trial (Twolayerlazygraph& G, int number_trials_binary_search) {

		double active_vertices_average = 0;
		int inhibition_threshold = G.inhibition_threshold;

		// Do binary search for the maximal number of patterns we can store
		set<int> bootstrap;
		for(int i = 0; i < G.pattern_size; i++) {
			bootstrap.insert(i);
		}

		long long num_patterns=1;
		int i = 0;
		while (active_vertices_average < G.fidelity_outside){
			num_patterns *= 2;
			//clog << "check number of revealed patterns is "<< num_patterns << "\n";
			int active_vertices_total = 0;
			for (int iterator = 0; iterator < number_trials_binary_search; iterator++) {
				set<int> activated_vertices;
				percolationwithinhibition(activated_vertices, G, bootstrap, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
				active_vertices_total += activated_vertices.size();
				//cout << activated_vertices.size() << " \n";;
			}
			active_vertices_average = ((double) active_vertices_total)/((double)number_trials_binary_search);
			cout << active_vertices_average << " \n";
			i++;
		}

		long long num_patterns_up=num_patterns;
		long long num_patterns_down= num_patterns/2;
		long long num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
		while(num_patterns_up - num_patterns_down > 1) {
			//clog << "check num patterns " << num_patterns_mid<< "\n";
			num_patterns=num_patterns_mid;
			int active_vertices_total = 0;
			for (int iterator = 0; iterator < number_trials_binary_search; iterator++) {
				set<int> activated_vertices;
				percolationwithinhibition(activated_vertices, G, bootstrap, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
				active_vertices_total += activated_vertices.size();
			}
			active_vertices_average = active_vertices_total/number_trials_binary_search;
			cout << active_vertices_average << "\n";
			if (active_vertices_average < G.fidelity_outside) {
				num_patterns_down= num_patterns_mid;
				num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
			}
			else{
				num_patterns_up= num_patterns_mid;
				num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
			}
		}


		return num_patterns_down;
		// Stores whole activity vector 
		/*for(int j=0; j< num_patterns_down*number_associations; j++ ){
			num_patterns = num_patterns_down;
			int bootstrap_pattern= j;
			//clog << "Shit\n";
			exposed = percolation::pattern_activation(LG.patterns_A[bootstrap_pattern], act_threshold, z_i, 1, num_patterns);
			//clog << "Shit2\n";

			int pattern_active = Utils::count_active(LG.patterns_B[bootstrap_pattern/number_associations], act_threshold);
			//inside_error_probability_aux += ((double)pattern_size - (double)pattern_active)/(double)pattern_size/(double)number_trials_pattern;
			vector_inside_activity[i][j] = pattern_active;
			active_ex = Utils::intersection_size(LG.layers[1].excitatory, exposed);
			//average_active_ex+= (double) active_ex / (double) number_trials_pattern; 
			//outside_error_probability_aux += ((double)active_ex - (double)pattern_active)/((double)(number_vertices-pattern_size))/(double)number_trials_pattern;
			vector_outside_activity[i][j] = active_ex - pattern_active;


			Reset::reset_voltages(LG.layers[1].excitatory, 0);
			Reset::reset_voltages(LG.layers[1].inhibitory, 0);
		}

		max_number_patterns_under_fidelity[i]= num_patterns_down*number_associations;
		//inside_error_probability[i] = inside_error_probability_aux;
		//outside_error_probability[i] = outside_error_probability_aux;

		clog << "Max number of patterns inserted before fidelity_outside is reached: " << num_patterns_down << "\n";

		
		//sort(history.first.begin(), history.first.end());
		//sort(history.second.begin(), history.second.end());
	   */

}

// finds the optimal number of patterns that can be stored respecting fidelities averaging over many graph trials.
// returns the average number of patterns over the trials and stores the actual sequence in the file in matlab format.
double percolation_test_parameters_with_inhibition_p_rec_given_binary_search_many_tries (int number_vertices, int pattern_size, int act_threshold, double p_rec, double p_aff, double fidelity_inside, double fidelity_outside, ofstream& file_percolation, int number_trials_new_graph, int number_trials_binary_search, int multi_association_factor, long long final_number_associations, int turn_strong_pattern_threshold, int inhibition_threshold) {

	// Find maximal number of patterns that can be stored;
	double maximal_number_patterns_stored = 0;
	file_percolation << "maximal_number_patterns_under_fidelity = [";
	for (int j = 0; j < number_trials_new_graph; j++) {
		cout << j << "\n";
		//////////////// Graph construction
		double active_vertices_average = 0;

		// Do binary search for the maximal number of patterns we can store
		set<int> bootstrap;
		for(int i = 0; i < pattern_size; i++) {
			bootstrap.insert(i);
		}

		long long num_patterns=1;
		int i = 0;
		while (active_vertices_average < fidelity_outside){
			num_patterns *= 2;
			clog << "check number of revealed patterns is "<< num_patterns << "\n";
			int active_vertices_total = 0;
			for (int iterator = 0; iterator < number_trials_binary_search; iterator++) {
				Twolayerlazygraph G;
				lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
				G.pattern_size = pattern_size;
				G.act_threshold = act_threshold;
				G.multi_association_factor = multi_association_factor;
				G.final_number_associations = final_number_associations;
				G.turn_strong_pattern_threshold = turn_strong_pattern_threshold;
				G.fidelity_outside = fidelity_outside;
				G.inhibition_threshold = inhibition_threshold;
				//cout << "Graph is built\n";

				////////////////// Pattern insertion

				pattern_insertion(G, num_patterns, multi_association_factor, pattern_size, turn_strong_pattern_threshold, p_rec);

				//cout << "Pattern insertion is complete\n";

				set<int> activated_vertices;
				percolationwithinhibition(activated_vertices, G, bootstrap, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
				active_vertices_total += activated_vertices.size();
				//cout << activated_vertices.size() << " \n";;
			}
			active_vertices_average = active_vertices_total/number_trials_binary_search;
			//cout << active_vertices_average << " \n";;
			i++;
		}

		long long num_patterns_up=num_patterns;
		long long num_patterns_down= num_patterns/2;
		long long num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
		while(num_patterns_up - num_patterns_down > 1) {
			clog << "check num patterns " << num_patterns_mid<< "\n";
			num_patterns=num_patterns_mid;
			int active_vertices_total = 0;
			for (int iterator = 0; iterator < number_trials_binary_search; iterator++) {
				Twolayerlazygraph G;
				lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
				G.pattern_size = pattern_size;
				G.act_threshold = act_threshold;
				G.multi_association_factor = multi_association_factor;
				G.final_number_associations = final_number_associations;
				G.turn_strong_pattern_threshold = turn_strong_pattern_threshold;
				G.fidelity_outside = fidelity_outside;
				G.inhibition_threshold = inhibition_threshold;
				//cout << "Graph is built\n";

				////////////////// Pattern insertion

				pattern_insertion(G, num_patterns, multi_association_factor, pattern_size, turn_strong_pattern_threshold, p_rec);

				//cout << "Pattern insertion is complete\n";
				
				set<int> activated_vertices;
				percolationwithinhibition(activated_vertices, G, bootstrap, G.act_threshold, inhibition_threshold, num_patterns, G.turn_strong_pattern_threshold);
				active_vertices_total += activated_vertices.size();
			}
			active_vertices_average = active_vertices_total/number_trials_binary_search;
			//cout << active_vertices_average << "\n";
			if (active_vertices_average < fidelity_outside) {
				num_patterns_down= num_patterns_mid;
				num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
			}
			else{
				num_patterns_up= num_patterns_mid;
				num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
			}
		}

		maximal_number_patterns_stored += ((double) num_patterns_down)/number_trials_new_graph;
		file_percolation << num_patterns_down << ",";
	}
	file_percolation << "]\n";
	return maximal_number_patterns_stored;
}

// finds the optimal number of patterns that can be stored respecting fidelities averaging over many graph trials.
// returns the average number of patterns over the trials and stores the actual sequence in the file in matlab format.
double percolation_test_parameters_with_inhibition_p_rec_given (int number_vertices, int pattern_size, int act_threshold, double p_rec, double p_aff, double fidelity_inside, double fidelity_outside, ofstream& file_percolation, int number_trials_new_graph, int number_trials_binary_search, int multi_association_factor, long long final_number_associations, int turn_strong_pattern_threshold, int inhibition_threshold) {

	// Find maximal number of patterns that can be stored;
	double maximal_number_patterns_stored = 0;
	file_percolation << "maximal_number_patterns_under_fidelity = [";
	for (int i = 0; i < number_trials_new_graph; i++) {
		cout << i << "\n";
		//////////////// Graph construction
		Twolayerlazygraph G;
		lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
		G.pattern_size = pattern_size;
		G.act_threshold = act_threshold;
		G.multi_association_factor = multi_association_factor;
		G.final_number_associations = final_number_associations;
		G.turn_strong_pattern_threshold = turn_strong_pattern_threshold;
		G.fidelity_outside = fidelity_outside;
		G.inhibition_threshold = inhibition_threshold;

		cout << "Graph is built\n";

		////////////////// Pattern insertion

		pattern_insertion(G, final_number_associations, multi_association_factor, pattern_size, turn_strong_pattern_threshold, p_rec);

		cout << "Pattern insertion is complete\n";

		long long number_patterns_in_this_graph = maximal_number_patterns_without_breaching_fidelity_outside_single_graph_trial(G, number_trials_binary_search);
		maximal_number_patterns_stored += ((double) number_patterns_in_this_graph)/number_trials_new_graph;
		file_percolation << number_patterns_in_this_graph << ",";
	}
	file_percolation << "];\n";
	return maximal_number_patterns_stored;
}

// finds the optimal number of patterns that can be stored respecting fidelities over many graph trials.
// stores the sequence of number of patterns obtained
void percolation_test_parameters_with_inhibition_p_rec_given_optimizer (vector<int>& maximal_number_patterns_under_fidelity, int number_vertices, int pattern_size, int act_threshold, double p_rec, double p_aff, double fidelity_inside, double fidelity_outside, ofstream& file_percolation, int number_trials_new_graph, int number_trials_binary_search, int multi_association_factor, long long final_number_associations, int turn_strong_pattern_threshold, int inhibition_threshold) {


	// Find maximal number of patterns that can be stored;
	//double maximal_number_patterns_stored = 0;
	for (int i = 0; i < number_trials_new_graph; i++) {
		cout << i << "\n";
		//////////////// Graph construction
		Twolayerlazygraph G;
		lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
		G.pattern_size = pattern_size;
		G.act_threshold = act_threshold;
		G.multi_association_factor = multi_association_factor;
		G.final_number_associations = final_number_associations;
		G.turn_strong_pattern_threshold = turn_strong_pattern_threshold;
		G.fidelity_outside = fidelity_outside;
		G.inhibition_threshold = inhibition_threshold;

		cout << "Graph is built\n";

		////////////////// Pattern insertion

		pattern_insertion(G, final_number_associations, multi_association_factor, pattern_size, turn_strong_pattern_threshold, p_rec);

		cout << "Pattern insertion is complete\n";

		long long number_patterns_in_this_graph = maximal_number_patterns_without_breaching_fidelity_outside_single_graph_trial(G, number_trials_binary_search);
		//maximal_number_patterns_stored += ((double)number_patterns_in_this_graph)/number_trials_new_graph;
		maximal_number_patterns_under_fidelity.push_back(number_patterns_in_this_graph);
	}
	//return maximal_number_patterns_stored;
}

int main(){
	// Parameters
	/*int count_file = 1;


	int num_parallel = 8;
    omp_set_num_threads(4);
	#pragma omp parallel for
	for (int i = 0; i < num_parallel; i++) {
		long long number_vertices = 5000; // Number of excitatory  vertices
		int pattern_size = 18; // Size of a pattern
		int multi_association_factor=1;
		int turn_strong_pattern_threshold = 1;
		int number_trials_new_graph = 30;
		int number_trials_binary_search = 1; // this only makes sense for the percolation_test_binary_search_many_tries
		double fidelity_inside = 0.95;

		int act_threshold = 9;
		int inhibition_threshold = pattern_size;
		double fidelity_outside = 2*pattern_size; // The program assumes that fidelity outside is the bound on the number of vertices active outside not the fraction.

		long long final_number_associations = ceil(2*number_vertices*number_vertices/(pattern_size*pattern_size*multi_association_factor));
		cout << multi_association_factor*final_number_associations << "!\n";

		ofstream file_percolation;
		count_file = i+20;
		file_percolation.open("percolation_" + to_string(count_file)+".m");
		file_percolation << "number_vertices=" << number_vertices << "; \n";
		file_percolation << "pattern_size=" << pattern_size << ";\n";
		file_percolation << "act_threshold=" << act_threshold << ";\n";
		//file_percolation << "fraction_rec_aff=" << fraction_rec_aff << ";\n";
		//file_percolation << "num_p_r_values=" << num_parallel << " ;\n";
		file_percolation << "fidelity_inside=" << fidelity_inside << ";\n";
		file_percolation << "fidelity_outside=" << fidelity_outside << ";\n";
		file_percolation << "number_associations=" << multi_association_factor << ";\n";
		//file_percolation << "num_num_association_values=" << num_num_association_values << "; \n";
		//file_percolation << "s_iterator=" << s_iterator+1 << "; \n";
		file_percolation << "number_trials_new_graph = " << number_trials_new_graph << ";\n";
		file_percolation << "number_trials_binary_search = " << number_trials_binary_search << ";\n";
		file_percolation << "inhibition_threshold = " << inhibition_threshold << ";\n";

		// optimal p_rec search
		// precision 0.1 <-- improve for better precision

		double p_rec = 0.1*i;
		int number_tries = 10000;
		int precision = 4;
		double p_aff = percolation_inside_p_aff_search_given_p_rec(pattern_size, act_threshold, number_tries, precision, p_rec, fidelity_inside);
		//double p_aff = 0.8915;
		cout << p_aff << "\n";

		file_percolation << "p_rec = " << p_rec << ";\n";
		file_percolation << "p_aff = " << p_aff << ";\n";

		double maximal_number_patterns_under_fidelity = percolation_test_parameters_with_inhibition_p_rec_given(number_vertices, pattern_size, act_threshold, p_rec, p_aff, fidelity_inside, fidelity_outside, file_percolation, number_trials_new_graph, number_trials_binary_search, multi_association_factor, final_number_associations, turn_strong_pattern_threshold, inhibition_threshold);
	    // ROUNDING DOWN TOO MUCH ON AVERAGE ^ FIX IT

		//number_trials_binary_search = 100;
		//double maximal_number_patterns_under_fidelity = percolation_test_parameters_with_inhibition_p_rec_given_binary_search_many_tries(number_vertices, pattern_size, act_threshold, p_rec, p_aff, fidelity_inside, fidelity_outside, file_percolation, number_trials_new_graph, number_trials_binary_search, multi_association_factor, final_number_associations, turn_strong_pattern_threshold, inhibition_threshold);
		
		file_percolation << "maximal_number_patterns_under_fidelity_average = " << maximal_number_patterns_under_fidelity << "\n";
		//return maximal_number_patterns_under_fidelity;	
	}*/
	

	int count_file = 10;

	long long number_vertices = 100000; // Number of excitatory  vertices
	int pattern_size = 13; // Size of a pattern
	int multi_association_factor=64;
	int turn_strong_pattern_threshold = 1;
	int number_trials_new_graph = 30;
	int number_trials_binary_search = 1; // this only makes sense for the percolation_test_binary_search_many_tries
	double fidelity_inside = 0.95;

	int act_threshold = 13;
	int inhibition_threshold = pattern_size;
	double fidelity_outside = 2*pattern_size; // The program assumes that fidelity outside is the bound on the number of vertices active outside not the fraction.

	long long final_number_associations = ceil(4*number_vertices*number_vertices/(pattern_size*pattern_size*multi_association_factor));
	cout << multi_association_factor*final_number_associations << "!\n";

	ofstream file_percolation;

	file_percolation.open("percolation_" + to_string(count_file)+".m");
	file_percolation << "number_vertices=" << number_vertices << "; \n";
	file_percolation << "pattern_size=" << pattern_size << ";\n";
	file_percolation << "act_threshold=" << act_threshold << ";\n";
	//file_percolation << "fraction_rec_aff=" << fraction_rec_aff << ";\n";
	//file_percolation << "num_p_r_values=" << num_parallel << " ;\n";
	file_percolation << "fidelity_inside=" << fidelity_inside << ";\n";
	file_percolation << "fidelity_outside=" << fidelity_outside << ";\n";
	file_percolation << "number_associations=" << multi_association_factor << ";\n";
	//file_percolation << "num_num_association_values=" << num_num_association_values << "; \n";
	//file_percolation << "s_iterator=" << s_iterator+1 << "; \n";
	file_percolation << "number_trials_new_graph = " << number_trials_new_graph << ";\n";
	file_percolation << "number_trials_binary_search = " << number_trials_binary_search << ";\n";
	file_percolation << "inhibition_threshold = " << inhibition_threshold << ";\n";

	// optimal p_rec search
	// precision 0.1 <-- improve for better precision

	double p_rec = 0.5;
	int number_tries = 10000;
	int precision = 4;
	double p_aff = percolation_inside_p_aff_search_given_p_rec(pattern_size, act_threshold, number_tries, precision, p_rec, fidelity_inside);
	//double p_aff = 1;
	cout << p_aff << "\n";

	file_percolation << "p_rec = " << p_rec << ";\n";
	file_percolation << "p_aff = " << p_aff << ";\n";

	double maximal_number_patterns_under_fidelity = percolation_test_parameters_with_inhibition_p_rec_given(number_vertices, pattern_size, act_threshold, p_rec, p_aff, fidelity_inside, fidelity_outside, file_percolation, number_trials_new_graph, number_trials_binary_search, multi_association_factor, final_number_associations, turn_strong_pattern_threshold, inhibition_threshold);
    // ROUNDING DOWN TOO MUCH ON AVERAGE ^ FIX IT

	//number_trials_binary_search = 100;
	//double maximal_number_patterns_under_fidelity = percolation_test_parameters_with_inhibition_p_rec_given_binary_search_many_tries(number_vertices, pattern_size, act_threshold, p_rec, p_aff, fidelity_inside, fidelity_outside, file_percolation, number_trials_new_graph, number_trials_binary_search, multi_association_factor, final_number_associations, turn_strong_pattern_threshold, inhibition_threshold);
	
	file_percolation << "maximal_number_patterns_under_fidelity_average = " << maximal_number_patterns_under_fidelity << ";\n";
	return maximal_number_patterns_under_fidelity;

	
	/*Twolayerlazygraph G;
	double p_aff = 1;
	double p_rec = 0.1;
	//final_number_associations = 46;

	lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
	pattern_insertion(G, final_number_associations, multi_association_factor, pattern_size, turn_strong_pattern_threshold, p_rec);
	int number_strong_afferent_edges = 0;
	for (int i = 0; i < pattern_size; i++) {
		for (int j = 0; j < number_vertices; j++){
			if(G.Patternafferentedges[i][j].is_present && G.Patternafferentedges[i][j].patternthreshold == 1 && G.Patternafferentedges[i][j].patternturnstrong != -1) {
				number_strong_afferent_edges++;
			}
		}
	}
	double fraction_strong_afferent_edges = ((double)number_strong_afferent_edges)/((double)number_vertices)/((double)pattern_size);

	cout << fraction_strong_afferent_edges << "\n";
	cout << number_strong_afferent_edges << "\n";

	int number_strong_recurrent_edges = 0;
	int number_extra_strong_recurrent_edges = 0;
	int everything_ok = 1;
	for (int i = 0; i < number_vertices; i++) {
		if (i  < pattern_size) {
			for (int j = 0; j < number_vertices; j++){
				if(G.Patternrecurrentedges[i][j].is_present && G.Patternrecurrentedges[i][j].patternthreshold == 1 && G.Patternrecurrentedges[i][j].patternturnstrong != -1) {
					number_strong_recurrent_edges++;
				}
			}			
		}
		else {
			int last_target_vertex = -1;
			for (vector<Edge>::iterator i_neighbourhood = G.Patternrecurrentedges[i].begin(); i_neighbourhood != G.Patternrecurrentedges[i].end(); i_neighbourhood++) {
				if((*i_neighbourhood).is_present && (*i_neighbourhood).patternthreshold == 1 && (*i_neighbourhood).patternturnstrong != -1 && (*i_neighbourhood).patternturnstrong < final_number_associations) {
					number_extra_strong_recurrent_edges++;
					cout << "souce vertex = " << i << "\n";
					cout << "target vertex = " << (*i_neighbourhood).target_vertex << "\n";
					if(last_target_vertex == -1) {
						last_target_vertex = (*i_neighbourhood).target_vertex;
					}
					else {
						if(last_target_vertex > (*i_neighbourhood).target_vertex) {
							everything_ok = 0;
							cout << everything_ok << "\n";
							cout << "souce vertex = " << i << "\n";
							cout << "target vertex = " << (*i_neighbourhood).target_vertex << "\n";
							cout << "last_target vertex = " << last_target_vertex << "\n";
							return 1;
						}
						last_target_vertex = (*i_neighbourhood).target_vertex;
					}
				}
			}
		}
	}
	double fraction_strong_recurrent_edges = ((double)number_strong_recurrent_edges)/((double)number_vertices-1)/((double)pattern_size);
	double fraction_extra_strong_recurrent_edges = ((double)number_extra_strong_recurrent_edges)/((double)number_vertices-1)/((double)number_vertices-pattern_size);
	
	cout << everything_ok << "\n";
	cout << fraction_strong_recurrent_edges << "\n";
	cout << fraction_extra_strong_recurrent_edges << "\n";
	cout << number_strong_recurrent_edges << "\n";
	cout << number_extra_strong_recurrent_edges << "\n";*/

}

///////////////////////// TEST CODES ///////////////////////////////////////////////////////////////


	//////////// Test edge addition in graph ///////////////////////////////////////
	/*Twolayerlazygraph G;
	double p_aff = 0.8;
	double p_rec = 0.1;
	lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
	int number_recurrent_edges = 0;
	for (int i = 0; i < pattern_size; i++) {
		for (int j = 0; j < number_vertices; j++){
			if(G.Patternrecurrentedges[i][j].is_present) {
				number_recurrent_edges++;
			}
		}
	}
	double fraction_recurrent_edges = ((double)number_recurrent_edges)/((double)pattern_size)/((double)number_vertices-1);

	cout << fraction_recurrent_edges << "\n";

	int number_afferent_edges = 0;
	for (int i = 0; i < pattern_size; i++) {
		for (int j = 0; j < number_vertices; j++){
			if(G.Patternafferentedges[i][j].is_present) {
				number_afferent_edges++;
			}
		}
	}
	double fraction_afferent_edges = ((double)number_afferent_edges)/((double)number_vertices)/((double)pattern_size);

	cout << fraction_afferent_edges << "\n";*/



	///////////// Test strong edges in graph (pattern insertion) ////////////////////////////////////////////
	/*
	Twolayerlazygraph G;
	double p_aff = 1;
	double p_rec = 0.1;
	//final_number_associations = 46;

	lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
	pattern_insertion(G, final_number_associations, multi_association_factor, pattern_size, turn_strong_pattern_threshold, p_rec);
	int number_strong_afferent_edges = 0;
	for (int i = 0; i < pattern_size; i++) {
		for (int j = 0; j < number_vertices; j++){
			if(G.Patternafferentedges[i][j].is_present && G.Patternafferentedges[i][j].patternthreshold == 1 && G.Patternafferentedges[i][j].patternturnstrong != -1) {
				number_strong_afferent_edges++;
			}
		}
	}
	double fraction_strong_afferent_edges = ((double)number_strong_afferent_edges)/((double)number_vertices)/((double)pattern_size);

	cout << fraction_strong_afferent_edges << "\n";
	cout << number_strong_afferent_edges << "\n";

	int number_strong_recurrent_edges = 0;
	int number_extra_strong_recurrent_edges = 0;
	for (int i = 0; i < number_vertices; i++) {
		if (i  < pattern_size) {
			for (int j = 0; j < number_vertices; j++){
				if(G.Patternrecurrentedges[i][j].is_present && G.Patternrecurrentedges[i][j].patternthreshold == 1 && G.Patternrecurrentedges[i][j].patternturnstrong != -1) {
					number_strong_recurrent_edges++;
				}
			}			
		}
		else {
			for (vector<Edge>::iterator i_neighbourhood = G.Patternrecurrentedges[i].begin(); i_neighbourhood != G.Patternrecurrentedges[i].end(); i_neighbourhood++) {
				if((*i_neighbourhood).is_present && (*i_neighbourhood).patternthreshold == 1 && (*i_neighbourhood).patternturnstrong != -1 && (*i_neighbourhood).patternturnstrong < final_number_associations) {
					number_extra_strong_recurrent_edges++;
					//cout << "souce vertex = " << i << "\n";
					//cout << "target vertex = " << (*i_neighbourhood).target_vertex << "\n";
				}
			}
		}
	}
	double fraction_strong_recurrent_edges = ((double)number_strong_recurrent_edges)/((double)number_vertices-1)/((double)pattern_size);
	double fraction_extra_strong_recurrent_edges = ((double)number_extra_strong_recurrent_edges)/((double)number_vertices-1)/((double)number_vertices-pattern_size);

	cout << fraction_strong_recurrent_edges << "\n";
	cout << fraction_extra_strong_recurrent_edges << "\n";
	cout << number_strong_recurrent_edges << "\n";
	cout << number_extra_strong_recurrent_edges << "\n";
	*/


	////////////////// Test percolation with inhibition ////////////////////////////////
	/*Twolayerlazygraph G;
	double p_aff = 1;
	double p_rec = 0;

	lazygraphconstructor(G, pattern_size, number_vertices, pattern_size, p_aff, p_rec);
	pattern_insertion(G, final_number_associations, multi_association_factor, pattern_size, turn_strong_pattern_threshold, p_rec);

	int test_number_patterns_inserted = (log(2)-0.6)*number_vertices*number_vertices/(pattern_size*pattern_size);
	test_number_patterns_inserted = 220000;
	// Random devices - to pick pattern at random
	set<int> bootstrap;
	for(int i = 0; i < G.pattern_size; i++) {
		bootstrap.insert(i);
	}
	set<int> activated_vertices;

	cout << "Percolation starting\n"; 
	percolationwithinhibition(activated_vertices, G, bootstrap, act_threshold, inhibition_threshold, test_number_patterns_inserted, turn_strong_pattern_threshold);
	int number_activated_vertices = activated_vertices.size();
	cout << number_activated_vertices << "\n";*/
