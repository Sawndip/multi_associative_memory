/*
 * florian_marcelo_inhibition_async.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: Trujic
 */
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <cmath>
#include <algorithm>
#include <stdlib.h> 
#include <omp.h>

#include "../../../../../Library/Percolation/Percolation.hpp"
#include "../../../../../Library/Graph_Constructs/Layered_graph.h"
#include "../../../../../Library/Utilities/Utils.hpp"

#define _USE_MATH_DEFINES

using namespace std;

random_device rd;
mt19937 gen(rd());
exponential_distribution<> d(1.0);

void process_edge(Edge& e) {
	if (e.source->get_type() == -1 || e.target->get_type() == -1) {
		e.set_weight(1);
	} else {
		e.set_weight(0);
	}
	e.get_delay = [] {
		return d(gen);
	};
}
// fidelity inside is fraction of verticies to turn active in average
void percolation_test_parameters( 
	int number_vertices, int pattern_size, int act_threshold,  double fidelity_inside, double fidelity_outside,
	ofstream& file_percolation, int number_trials_new_graph, int number_trials_pattern, int number_associations, int final_number_patterns_stored) {

	
// Binary search for best recurrent density
	double p_rec=0;
	double p_aff=0;
	double p_rec_down = 0.0; 
	double p_rec_up=1.0;
	int z_i=0; // wont be used, just to call the function
	int p_rec_precision=3;
	double p_rec_precision_inc_inc=0.1;
	double p_rec_precision_inc=1;
	double best_p_rec=0;
	double best_p_aff=0;
	double max_number_patterns=0;
	for(int a=0; a< p_rec_precision; a++){
		p_rec_precision_inc *= p_rec_precision_inc_inc;
		clog << "new precision" << p_rec_precision_inc << " \n ";
		int number_p_rec_values= floor( (p_rec_up-p_rec_down)/p_rec_precision_inc+1);
		clog << "number_p_rec_values is "<< number_p_rec_values << " \n";
		clog << p_rec_up << "  " << p_rec_down << " " << p_rec_precision_inc  <<" \n";
		vector<double> p_rec_values(floor( (p_rec_up-p_rec_down)/p_rec_precision_inc+1), 0);
		vector<double> maximal_number_pattern_values(floor( (p_rec_up-p_rec_down)/p_rec_precision_inc+1), 0);
		//vector<vector<int>> p_rec_values(number_trials_new_graph, vector<int>(number_p_rec_values));
		//vector<vector<int>> maximal_number_pattern_values(number_p_rec_values, vector<int>(number_trials_new_graph));

		best_p_rec=0;
		best_p_aff=0;
		max_number_patterns=0;
		for(int b=0; b< number_p_rec_values ; b++){
			p_rec= p_rec_down + b* p_rec_precision_inc;
			clog << "test rec density "<< p_rec <<" \n ";



			clog << "\n --- OPTIMAL DENSITY SEARCH --- \n";
			set<Vertex*> exposed;

			//Search for right density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			vector<int> pattern_ex({pattern_size, pattern_size});
			vector<int> pattern_in({0, 0});


				int num_tries = 10000;


				int precision=4; // determine the correct affarent density up to 10^-precision
				double increment=1;
				p_aff=0;
				for(int k=1; k< precision+1; k++){
					increment= increment/(double) 10;
					double average_active=0;
					while (average_active < fidelity_inside * pattern_size) {  
					
						p_aff+=increment;
						//cout << " Check p_aff= "<< p_aff << " \n";
						
						p_rec = min(1.0, p_rec);
						average_active=0;


						for(int i=0; i < num_tries; i++){
							Layered_graph LG(pattern_ex, pattern_in, 1.0, 1.0);
							Connect::connect(LG.layers[1].excitatory, LG.layers[1].excitatory, p_rec, process_edge);
							// Connect::connect(LG.layers[1].excitatory, LG.layers[1].inhibitory, p_r, process_edge);
							// Connect::connect(LG.layers[1].inhibitory, LG.layers[1].excitatory, 1.0, process_edge);
							Connect::connect(LG.layers[0].excitatory, LG.layers[1].excitatory, p_aff, process_edge);

							LG.reveal_patterns(LG.layers[0].excitatory, LG.layers[1].excitatory, pattern_size, 1, 1, 1);

							exposed = percolation::pattern_activation(LG.patterns_A[0], act_threshold, z_i, 1, 1);

							average_active +=  Utils::count_active(exposed, act_threshold)/(double) num_tries;
							//cout << p_r << ": " << exposed.size() << " " << active_ex << endl;


						}
						
						
					}
					p_aff-=increment;

				}
				
				
				clog << "Optimal p_a for which almost all vertices in Bi are exposed: "
				<< p_aff << "\n";
		
			

			//End Search for right density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


			vector<int> excitatory({number_vertices, number_vertices});
			vector<int> inhibitory({0, 10});

			vector<int> max_patterns_trials(number_trials_new_graph);
			vector<int> max_number_patterns_under_fidelity(number_trials_new_graph);
			double active_ex=0;
			for (int i = 0; i < number_trials_new_graph; i++) {
				active_ex = 0;
				clog << "graph trial " << i<< " \n";



				Layered_graph LG(excitatory, inhibitory, 1.0, 1.0);
				

				
				Connect::connect(LG.layers[1].excitatory, LG.layers[1].excitatory, p_rec, process_edge);
				
				Connect::connect(LG.layers[1].excitatory, LG.layers[1].inhibitory, 0, process_edge);
				
				Connect::connect(LG.layers[1].inhibitory, LG.layers[1].excitatory, 0, process_edge);
				
				Connect::connect(LG.layers[0].excitatory, LG.layers[1].excitatory, p_aff, process_edge);
				

				//clog << " --- PATTERN REVEALING NO. " << i + 1 << " --- \n";

				
				
				LG.reveal_patterns(LG.layers[0].excitatory, LG.layers[1].excitatory, pattern_size, number_associations, 1, final_number_patterns_stored);
				
				clog << "Search for max number patterns under fidelity \n";
				pair<vector<pair<double, int> >, vector<pair<double, int> > > data; 
				
				

				double average_active_ex=0;


				// Do binary search for the maximal number of patterns we can store such that
				
				int num_patterns=1;
				while (average_active_ex < fidelity_outside){
					num_patterns*=2;
					//clog << "check number of revealed patterns is "<< num_patterns << "\n";
					average_active_ex=0;
						for(int j=0; j< number_trials_pattern; j++ ){
							int rand_bootstrap_pattern= rand()% num_patterns*number_associations;
							//clog << "Shit\n";
							exposed = percolation::pattern_activation(LG.patterns_A[rand_bootstrap_pattern], act_threshold, z_i, 1, num_patterns);
							//clog << "Shit2\n";

							int pattern_active = Utils::count_active(LG.patterns_B[rand_bootstrap_pattern/number_associations], act_threshold);
							active_ex = Utils::intersection_size(LG.layers[1].excitatory, exposed);
							average_active_ex+= (double) active_ex / (double) number_trials_pattern; 
						

							data.first.push_back(make_pair(num_patterns * number_associations, pattern_active));
							data.second.push_back(make_pair(num_patterns * number_associations, active_ex));

							Reset::reset_voltages(LG.layers[1].excitatory, 0);
							Reset::reset_voltages(LG.layers[1].inhibitory, 0);
						}

				}

				int num_patterns_up=num_patterns;
				int num_patterns_down= num_patterns/2;
				int num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
				while( num_patterns_up - num_patterns_down > 1) {
					//clog << "check num patterns " << num_patterns_mid;
					average_active_ex=0;
					num_patterns=num_patterns_mid;
						for(int j=0; j< number_trials_pattern; j++ ){
							int rand_bootstrap_pattern= rand()% num_patterns*number_associations;
							//clog << "Shit\n";
							exposed = percolation::pattern_activation(LG.patterns_A[rand_bootstrap_pattern], act_threshold, z_i, 1, num_patterns);
							//clog << "Shit2\n";

							int pattern_active = Utils::count_active(LG.patterns_B[rand_bootstrap_pattern/number_associations], act_threshold);
							active_ex = Utils::intersection_size(LG.layers[1].excitatory, exposed);
							average_active_ex+= (double) active_ex / (double) number_trials_pattern; 
						

							data.first.push_back(make_pair(num_patterns * number_associations, pattern_active));
							data.second.push_back(make_pair(num_patterns * number_associations, active_ex));

							Reset::reset_voltages(LG.layers[1].excitatory, 0);
							Reset::reset_voltages(LG.layers[1].inhibitory, 0);
						}
					if (average_active_ex < fidelity_outside) {
						num_patterns_down= num_patterns_mid;
						num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
					}	
					else{ 
						num_patterns_up= num_patterns_mid;
						num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
					}
				}
				maximal_number_pattern_values[b]+= (double) num_patterns_down  / (double) number_trials_new_graph  * number_associations;

				clog << "Max number of patterns inserted before fidelity_outside is reached: " << maximal_number_pattern_values[b] << "\n";

				
				//sort(history.first.begin(), history.first.end());
				//sort(history.second.begin(), history.second.end());
			   
				
				
			}
			clog << "Average number Mmax "<< maximal_number_pattern_values[b] << "\n";
			clog << " maximum so far " << max_number_patterns << "\n";
			if (max_number_patterns < maximal_number_pattern_values[b]){
				max_number_patterns=maximal_number_pattern_values[b];
				best_p_rec= p_rec;
				best_p_aff= p_aff;

			}
			if (max_number_patterns > 1.05* maximal_number_pattern_values[b]){// to speed this up
				b=number_p_rec_values;
			}
		
			
		}
		p_rec_down=max(0.0, best_p_rec - p_rec_precision_inc);
		p_rec_up= min(1.0, best_p_rec + p_rec_precision_inc);	
		clog << "p_rec_up" << p_rec_up << " \n";
		clog << "p_rec_down" << p_rec_down << " \n";
	}
	file_percolation << "mat_s_prec_paff_mmax = [ "<< number_associations << " , " << best_p_rec << " , " << best_p_aff << 
		" , " << max_number_patterns << " ];";
	clog << "mat_s_prec_paff_mmax = [ "<< number_associations << " , " << best_p_rec << " , " << best_p_aff << 
		" , " << max_number_patterns << " ]; \n";
}

int main(){
	int num_parallel=12;
	//int count_file=1;
	#pragma omp parallel for
	for (int i=0; i < num_parallel;i++){
		int number_vertices = 5000; // Number of excitatory  vertices
		int act_threshold=1;
		int pattern_size=1;
		double alpha = 0.5;
		if (i<6) {
			pattern_size = 50*(i+1); // Size of a pattern
			act_threshold = 10; // Activation threshold for excitatory vertices			
		}
		else{
			pattern_size = 50*(i-5); // Size of a pattern
			act_threshold = 20; // Activation threshold for excitatory vertices						
		}
		//int number_associations_additive_increment=1;
		int number_associations=8;
		int num_num_association_values=1;
		int count_file=i+1;
			//for(int s_iterator=0; s_iterator< num_num_association_values; s_iterator++){
			//	if(s_iterator==0){
			//		number_associations = 1; // Bias towards number of patterns in A
			//	}
			//	else { 
			//		number_associations += number_associations_additive_increment;
			//	}
				double fidelity_inside=0.95; // between 0 and 1
				int number_trials_pattern=100;
				int number_trials_new_graph=30;
				int final_number_patterns_stored = (number_vertices/pattern_size)*(number_vertices/pattern_size)/number_associations;
				double fidelity_outside=2*pattern_size;



				ofstream file_percolation;


				file_percolation.open("percolation_" + to_string(count_file)+".m");
				




				file_percolation << "number_vertices=" << number_vertices << "; \n";
				file_percolation << "pattern_size=" << pattern_size << ";\n";
				file_percolation << "alpha=" << alpha << ";\n";
				file_percolation << "act_threshold=" << act_threshold << ";\n";
				//file_percolation << "fraction_rec_aff=" << fraction_rec_aff << ";\n";
				
				file_percolation << "fidelity_inside=" << fidelity_inside << ";\n";
				file_percolation << "fidelity_outside=" << fidelity_outside << ";\n";
				file_percolation << "number_associations=" << number_associations << ";\n";
				file_percolation << "num_num_association_values=" << num_num_association_values << "; \n";

				//file_percolation << "s_iterator=" << s_iterator+1 << "; \n";
				file_percolation << "number_trials_new_graph=" << number_trials_new_graph << ";\n";
				file_percolation << "number_trials_pattern=" << number_trials_pattern << ";\n";
				//file_percolation << " mat_s_prec_paff_mmax= []; \n";


				clog << "number_vertices=" << number_vertices << "; \n";
				clog << "pattern_size=" << pattern_size << ";\n";
				clog << "alpha=" << alpha << ";\n";
				clog << "act_threshold=" << act_threshold << ";\n";
				//clog << "fraction_rec_aff=" << fraction_rec_aff << ";\n";
				clog << "fidelity_inside=" << fidelity_inside << ";\n";
				clog << "number_associations=" << number_associations << ";\n";

				percolation_test_parameters( 
				number_vertices, pattern_size, act_threshold,  
				 fidelity_inside, fidelity_outside, file_percolation,  
				 number_trials_new_graph,  number_trials_pattern, number_associations, final_number_patterns_stored);


				file_percolation.close();
				//count_file++; // This should be at the end of the innermost forloop
			//}
	}
}