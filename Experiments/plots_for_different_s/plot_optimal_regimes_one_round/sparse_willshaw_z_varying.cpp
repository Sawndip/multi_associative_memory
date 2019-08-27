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
	int number_vertices, int pattern_size, int act_threshold, double p_rec, double fidelity_inside, double fidelity_outside,
	ofstream& file_percolation, int number_trials_new_graph, int number_trials_pattern, int number_associations, int final_number_patterns_stored) {

	
	//double p_rec = 0.0; // Density of recurrent edges
	double p_aff = 0; // will have to search for density of afferent edges
	int z_i=0; // wont be used, just to call the function

	clog << "\n --- OPTIMAL DENSITY SEARCH --- \n";
	set<Vertex*> exposed;
	int active_ex = 0;


//Search for right density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	vector<int> pattern_ex({pattern_size, pattern_size});
	vector<int> pattern_in({0, 0});

	int num_tries = 10000;
	
	int precision=4; // determine the correct density up to 10^-precision
	double increment=1;
	for(int k=1; k< precision+1; k++){
		increment= increment/(double) 10;
		double average_active=0;
		while (average_active < fidelity_inside * pattern_size) {  
		
			p_aff+=increment;
			cout << " Check p_aff= "<< p_aff << " \n";
			p_aff = min(1.0,p_aff);
			//p_rec = min(1.0, p_rec);
			average_active=0;


			for(int i=0; i < num_tries; i++){
				Layered_graph LG(pattern_ex, pattern_in, 1.0, 1.0);
				Connect::connect(LG.layers[1].excitatory, LG.layers[1].excitatory, p_rec, process_edge);
				// Connect::connect(LG.layers[1].excitatory, LG.layers[1].inhibitory, p_r, process_edge);
				// Connect::connect(LG.layers[1].inhibitory, LG.layers[1].excitatory, 1.0, process_edge);
				Connect::connect(LG.layers[0].excitatory, LG.layers[1].excitatory, p_aff, process_edge);

				LG.reveal_patterns(LG.layers[0].excitatory, LG.layers[1].excitatory, pattern_size, 1, 1, 1);

				exposed = percolation::one_round_percolation(LG.patterns_A[0], act_threshold, z_i, 1, 1);

				//cout << "p_a " << p_aff <<"  p_r " << p_rec << " " << Utils::count_active(exposed, act_threshold) << " \n" ;
				
				average_active +=  Utils::count_active(exposed, act_threshold)/(double) num_tries;
				//cout << p_r << ": " << exposed.size() << " " << active_ex << endl;


			}
			cout << "average activity is " << average_active << " \n";
			
		}
		p_aff -= increment;

	}


	clog << "Optimal p_r for which almost all vertices in Bi are exposed: "
			<< p_rec << "\n";
	clog << "Optimal p_a for which almost all vertices in Bi are exposed: "
			<< p_aff << "\n";
	file_percolation << "p_rec=" << p_rec << ";\n";
	file_percolation << "p_aff=" << p_aff << ";\n";

//End Search for right density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	vector<int> excitatory({number_vertices, number_vertices});
	vector<int> inhibitory({0, 10});

	vector<int> max_patterns_trials(number_trials_new_graph);
	vector<int> max_number_patterns_under_fidelity(number_trials_new_graph);
	//vector<double> inside_error_probability(number_trials_new_graph);
	//vector<double> outside_error_probability(number_trials_new_graph);
	vector<vector<int>> vector_inside_activity(number_trials_new_graph, vector<int>(final_number_patterns_stored*number_associations));
	vector<vector<int>> vector_outside_activity(number_trials_new_graph, vector<int>(final_number_patterns_stored*number_associations));
	for (int i = 0; i < number_trials_new_graph; i++) {
		for (int j = 0; j < final_number_patterns_stored*number_associations; j++) {
			vector_inside_activity[i][j] = 0;
			vector_outside_activity[i][j] = 0;
		}
	}
	for (int i = 0; i < number_trials_new_graph; i++) {
		active_ex = 0;



		Layered_graph LG(excitatory, inhibitory, 1.0, 1.0);



		Connect::connect(LG.layers[1].excitatory, LG.layers[1].excitatory, p_rec, process_edge);

		Connect::connect(LG.layers[1].excitatory, LG.layers[1].inhibitory, 0, process_edge);

		Connect::connect(LG.layers[1].inhibitory, LG.layers[1].excitatory, 0, process_edge);

		Connect::connect(LG.layers[0].excitatory, LG.layers[1].excitatory, p_aff, process_edge);


		//clog << " --- PATTERN REVEALING NO. " << i + 1 << " --- \n";


		LG.reveal_patterns(LG.layers[0].excitatory, LG.layers[1].excitatory, pattern_size, number_associations, 1, final_number_patterns_stored);
		
		clog << "Search for max number patterns under fidelity \n";
		
		

		double average_active_ex=0;


		// Do binary search for the maximal number of patterns we can store such that
		
		int num_patterns=1;
		while (average_active_ex < fidelity_outside){
			num_patterns*=2;
			clog << "check number of revealed patterns is "<< num_patterns << "\n";
			average_active_ex=0;
				for(int j=0; j< number_trials_pattern; j++ ){
					int rand_bootstrap_pattern= rand()% num_patterns*number_associations;
					//clog << "Shit\n";
					exposed = percolation::one_round_percolation(LG.patterns_A[rand_bootstrap_pattern], act_threshold, z_i, 1, num_patterns);
					//clog << "Shit2\n";

					//int pattern_active = Utils::count_active(LG.patterns_B[rand_bootstrap_pattern/number_associations], act_threshold);
					active_ex = Utils::intersection_size(LG.layers[1].excitatory, exposed);
					average_active_ex+= (double) active_ex / (double) number_trials_pattern; 


					Reset::reset_voltages(LG.layers[1].excitatory, 0);
					Reset::reset_voltages(LG.layers[1].inhibitory, 0);
				}

		}

		int num_patterns_up=num_patterns;
		int num_patterns_down= num_patterns/2;
		int num_patterns_mid= (num_patterns_up+num_patterns_down)/2;
		while( num_patterns_up - num_patterns_down > 1) {
			clog << "check num patterns " << num_patterns_mid<< "\n";
			average_active_ex=0;
			num_patterns=num_patterns_mid;
				for(int j=0; j< number_trials_pattern; j++ ){
					int rand_bootstrap_pattern= rand()% num_patterns*number_associations;
					//clog << "Shit\n";
					exposed = percolation::one_round_percolation(LG.patterns_A[rand_bootstrap_pattern], act_threshold, z_i, 1, num_patterns);
					//clog << "Shit2\n";

					//int pattern_active = Utils::count_active(LG.patterns_B[rand_bootstrap_pattern/number_associations], act_threshold);
					active_ex = Utils::intersection_size(LG.layers[1].excitatory, exposed);
					average_active_ex+= (double) active_ex / (double) number_trials_pattern; 
				


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


		// Computes activity vectors for exact information
		/*for(int j=0; j< num_patterns_down*number_associations; j++ ){
			num_patterns = num_patterns_down;
			int bootstrap_pattern= j;
			//clog << "Shit\n";
			exposed = percolation::one_round_percolation(LG.patterns_A[bootstrap_pattern], act_threshold, z_i, 1, num_patterns);
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
		}*/

		max_number_patterns_under_fidelity[i]= num_patterns_down*number_associations;
		//inside_error_probability[i] = inside_error_probability_aux;
		//outside_error_probability[i] = outside_error_probability_aux;

		clog << "Max number of patterns inserted before fidelity_outside is reached: " << num_patterns_down << "\n";

		
		//sort(history.first.begin(), history.first.end());
		//sort(history.second.begin(), history.second.end());
	   
		
		
	}

	file_percolation << "maximal_number_patterns_under_fidelity=[";
	for( int i=0; i< number_trials_new_graph-1; i++){
		file_percolation << max_number_patterns_under_fidelity[i] << ","; 
	}
	file_percolation << max_number_patterns_under_fidelity[number_trials_new_graph-1] << ";\n"; 
	

	/*file_percolation << "inside_activity=[";
	for( int i=0; i< number_trials_new_graph; i++){
		for (int j = 0; j < final_number_patterns_stored*number_associations-1; j++) {
			file_percolation << vector_inside_activity[i][j] << ",";
		}
		if(i<number_trials_new_graph-1){
			file_percolation << vector_inside_activity[i][final_number_patterns_stored*number_associations-1] << ";";
		}
		else {
			file_percolation << vector_inside_activity[i][final_number_patterns_stored*number_associations-1] << "]";			
		}
	}

	file_percolation << "outside_activity=[";

	for( int i=0; i< number_trials_new_graph; i++){
		for (int j = 0; j < final_number_patterns_stored*number_associations-1; j++) {
			file_percolation << vector_outside_activity[i][j] << ",";
		}
		if(i<number_trials_new_graph-1){
			file_percolation << vector_outside_activity[i][final_number_patterns_stored*number_associations-1] << ";";
		}
		else {
			file_percolation << vector_outside_activity[i][final_number_patterns_stored*number_associations-1] << "]";			
		}
	}*/
}

int main(){

	int num_parallel=5;
	//double p_r_increment = 0.01;
	//int count_file=1;
    omp_set_num_threads(1);
    //omp_set_num_threads(5);
	#pragma omp parallel for
	for (int i=0; i < num_parallel;i++){
		//int p_r_iterator = i;
		//int num_p_r_values = 50;
		//int num_fraction_rec_aff_values=50;
		//double fract_rec_aff_multiplicative_increment= pow(2, 0.5);
		//double fraction_rec_aff= 0.25/ fract_rec_aff_multiplicative_increment;
		//double fraction_rec_aff_precision = 0.025;
		//for(int p_r_iterator=0; p_r_iterator< num_p_r_values; p_r_iterator++){
			//fraction_rec_aff*= fract_rec_aff_multiplicative_increment;
			//if (fraction_rec_aff<1){
			//	fraction_rec_aff += fraction_rec_aff_precision;
			//}
			//else{
			//	fraction_rec_aff += 10*fraction_rec_aff_precision;
			//}
			//double p_r = p_r_iterator*p_r_increment;
		int number_vertices = 5000; // Number of excitatory  vertices
		int act_threshold=12;
		int pattern_size=12;
		double p_r = 0.25*i;

			//int number_associations_multiplicative_increment=2;
			int number_associations=128;
			//int num_num_association_values=7;
			//for(int s_iterator=0; s_iterator< num_num_association_values; s_iterator++){
			//	if(s_iterator==0){
			//		number_associations = 1; // Bias towards number of patterns in A
			//	}
			//	else {
			//		number_associations*= number_associations_multiplicative_increment;
			//	}
			
			int count_file = i+1;
			double fidelity_inside=0.95; // between 0 and 1
			int number_trials_pattern=100;
			int number_trials_new_graph=10;
			int final_number_patterns_stored = 4*(number_vertices/pattern_size)*(number_vertices/pattern_size)/number_associations;
			double fidelity_outside=2*pattern_size;



			ofstream file_percolation;

			file_percolation.open("percolation_" + to_string(count_file)+".m");
					




			file_percolation << "number_vertices=" << number_vertices << "; \n";
			file_percolation << "pattern_size=" << pattern_size << ";\n";
			file_percolation << "act_threshold=" << act_threshold << ";\n";
			//file_percolation << "fraction_rec_aff=" << fraction_rec_aff << ";\n";
			//file_percolation << "num_p_r_values=" << num_parallel << " ;\n";
			file_percolation << "fidelity_inside=" << fidelity_inside << ";\n";
			file_percolation << "fidelity_outside=" << fidelity_outside << ";\n";
			file_percolation << "number_associations=" << number_associations << ";\n";
			//file_percolation << "num_num_association_values=" << num_num_association_values << "; \n";
			//file_percolation << "s_iterator=" << s_iterator+1 << "; \n";
			file_percolation << "number_trials_new_graph=" << number_trials_new_graph << ";\n";
			file_percolation << "number_trials_pattern=" << number_trials_pattern << ";\n";
				


			clog << "number_vertices=" << number_vertices << "; \n";
			clog << "pattern_size=" << pattern_size << ";\n";
			clog << "act_threshold=" << act_threshold << ";\n";
			clog << "p_r=" << p_r << ";\n";
			clog << "fidelity_inside=" << fidelity_inside << ";\n";
			clog << "number_associations=" << number_associations << ";\n";

			percolation_test_parameters(
				number_vertices, pattern_size, act_threshold,  p_r, 
				 fidelity_inside, fidelity_outside, file_percolation,  number_trials_new_graph,  number_trials_pattern, number_associations, final_number_patterns_stored);

			file_percolation.close();
				//count_file++; // This should be at the end of the innermost forloop
			//}
		//}
		//}
	}}