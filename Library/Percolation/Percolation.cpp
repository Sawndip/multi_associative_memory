#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <random>

#include "Percolation.hpp" //Includes Layers
#include "../Vertices/Vertex.h"
#include "../Edges/Edge.h"
#include "../Signals/Transmitting_Edge.h"
#include "../Utilities/Utils.hpp"

// Method for (-Inf, Inf) problem in Bayesian activation model
double F(double x) {
	if (x > 0) {
		return x;
	} else {
		return 1;
	}
}

// Method for (-Inf, Inf) problem in Bayesian activation model
double G(double x) {
	if (x > 0) {
		return 0;
	} else {
		return 1;
	}
}

namespace percolation {
	std::set<Vertex*> asynchronous_threshold(std::set<Vertex*> bootstrap, int K,
			double threshold) {
		std::set<Vertex*> exposed;
		priority_queue<Transmitting_Edge, vector<Transmitting_Edge> > Q;
		std::set<Vertex*>::iterator it;
		for(it = bootstrap.begin(); it != bootstrap.end(); ++it) {
			if((*it)->spike_time < 0) {
				(*it)->spike(0);
			}
			(*it)->expose((*it)->spike_time, Q, threshold);
			exposed.insert(*it);
		}
		
		// Run percolation
		while (!Q.empty()) {
			Transmitting_Edge e = Q.top();
			Q.pop();
			e.update(Q, exposed, K, threshold);
		}

		return exposed;
	}

	std::set<Vertex*> asynchronous(std::set<Vertex*> bootstrap, int K) {
		return asynchronous_threshold(bootstrap, K, 1);
	}

	std::set<Vertex*> asynchronous_timed_bootstrap(
			std::set<pair<double, Vertex*> > bootstrap, int K) {
		std::set<Vertex*> exposed;
		priority_queue<Transmitting_Edge, vector<Transmitting_Edge> > Q;
		std::set<pair<double, Vertex*> >::iterator it;
		for (it = bootstrap.begin(); it != bootstrap.end(); ++it) {
			(*it).second->spike((*it).first);
			(*it).second->expose((*it).first, Q, 1);
			exposed.insert((*it).second);
		}

		// Run percolation
		while (!Q.empty()) {
			Transmitting_Edge e = Q.top();
			Q.pop();
			e.update(Q, exposed, K, 1);
		}
		return exposed;
	}

	std::vector<std::pair<double, int> > dynamic_threshold(
			std::set<Vertex*> bootstrap, int K, double T, int init_voltage,
			double threshold) {
		double t = 0.0;
		vector<pair<double, int> > history;
		priority_queue<Transmitting_Edge, vector<Transmitting_Edge> > Q;
		std::set<Vertex*>::iterator it;

		for (it = bootstrap.begin(); it != bootstrap.end(); ++it) {
			(*it)->expose(0, Q, threshold);
			double rp = (*it)->sample_refractory_period();
			(*it)->set_sleep_until(t + rp);
			(*it)->voltage = init_voltage;
			history.push_back(make_pair(0, 1));
			history.push_back(make_pair(rp, -1));
		}

		// Run percolation
		while (!Q.empty() && t < T) {
			Transmitting_Edge e = Q.top();
			Q.pop();
			t = e.time;
			Vertex *target = e.edge->target, *source = e.edge->source;
			if (target->awake(t)) {
				target->voltage = target->voltage
						+ (source->type * e.edge->get_weight());
				if (target->voltage >= K) {
					clog << t << "\n";
					target->expose(e.time, Q, threshold);
					double rp = target->sample_refractory_period();
					// TODO: Wrap the next two lines into a "reset" function?
					target->set_sleep_until(t + rp);
					target->voltage = init_voltage;
					history.push_back(make_pair(t, 1));
					history.push_back(make_pair(t + rp, -1));
				} else if (target->voltage < 0) {
					target->voltage = 0;
				}
			}
		}
		
		return history;
	}

	std::vector<std::pair<double, int> > dynamic(std::set<Vertex*> bootstrap, int K,
			double T, int init_voltage) {
		return dynamic_threshold(bootstrap, K, T, init_voltage, 1);
	}

	/**
	 * @param	bootstrap		bootstrap set for the percolation
	 * @param	K_e				activation threshold for excitatory vertices
	 * @param	K_i				activation threshold for inhibitory vertices
	 * @param	threshold		threshold for edge working or not
	 * @param	pattern			number of revealed patterns
	 *
	 * @return	set<Vertex*> 	set of exposed vertices
	 */
	std::set<Vertex*> pattern_activation(std::set<Vertex*> bootstrap, int K_e,
			int K_i, int threshold, int pattern) {
		std::set<Vertex*> exposed;
		priority_queue<Transmitting_Edge, vector<Transmitting_Edge> > Q;
		std::set<Vertex*>::iterator it;
		for (it = bootstrap.begin(); it != bootstrap.end(); ++it) {
			(*it)->spike(0);
			(*it)->expose(0, 0, Q, threshold, pattern);
		}

		while (!Q.empty()) {
			Transmitting_Edge e = Q.top();
			Q.pop();
			e.update(Q, exposed, K_e, K_i, threshold, pattern);
		}

		return exposed;
	}

/* ONE ROUND PERCOLATION */
	std::set<Vertex*> one_round_percolation(std::set<Vertex*> bootstrap, int K_e,
			int K_i, int threshold, int pattern) {
		std::set<Vertex*> exposed;
		priority_queue<Transmitting_Edge, vector<Transmitting_Edge> > Q;
		std::set<Vertex*>::iterator it;
		for (it = bootstrap.begin(); it != bootstrap.end(); ++it) {
			(*it)->spike(0);
			(*it)->expose(0, 0, Q, threshold, pattern);
		}
		priority_queue<Transmitting_Edge, vector<Transmitting_Edge> > Qtwo;
		while (!Q.empty()) {
			Transmitting_Edge e = Q.top();
			Q.pop();
			e.update(Qtwo, exposed, K_e, K_i, threshold, pattern);
		}
		while (!Qtwo.empty()) {
			Transmitting_Edge e = Qtwo.top();
			Qtwo.pop();
			e.updateoneround(Qtwo, exposed, K_e, K_i, threshold, pattern);
		}

		return exposed;
	}
	/**
	 * For a fixed number of revealed patterns stores the timestamp for each vertex activation
	 *
	 * @param	bootstrap		bootstrap set for the percolation
	 * @param	K_e				activation threshold for excitatory vertices
	 * @param	K_i				activation threshold for inhibitory vertices
	 * @param	threshold		threshold for edge working or not
	 * @param	pattern			number of revealed patterns
	 *
	 * @return	set<Vertex*> 	set of exposed vertices
	 */
	vector<pair<double, int> > async_fixed_pattern_activation(
			std::set<Vertex*> bootstrap, int K_e, int K_i, double threshold,
			int pattern) {
		vector<pair<double, int> > history;
		set<Vertex*> exposed;
		priority_queue<Transmitting_Edge, vector<Transmitting_Edge> > Q;

		std::set<Vertex*>::iterator it;
		for (it = bootstrap.begin(); it != bootstrap.end(); ++it) {
			(*it)->expose(0, 0, Q, threshold, pattern);
		}

		int exposed_ex = 0;
		double t = 0.0;
		while (!Q.empty()) {
			Transmitting_Edge e = Q.top();
			Q.pop();

			t = e.time;
			Vertex *target = e.edge->target, *source = e.edge->source;
			if (exposed.count(target) == 0 && (e.edge->get_pattern() <= pattern)) {
				target->voltage += source->type * e.edge->get_weight();
				if ((target->type == -1 && target->voltage >= K_i)
						|| (target->type == 1 && target->voltage >= K_e)) {
					target->expose(t, e.round + 1, Q, threshold, pattern);
					exposed.insert(target);
					if (target->get_type() == 1) {
						exposed_ex++;
						history.push_back(make_pair(t, exposed_ex));
					}
				}
			}
		}

		return history;
	}

	std::set<Vertex*> bayesian_pattern_activation(std::vector<Vertex>& sources, std::vector<Vertex>& targets, std::set<Vertex*> pattern, int patterns) {
		std::set<Vertex*> exposed;

		for (int j = 0; j < (int) targets.size(); j++) {
			double edge_contribution = 0;
			double infinity_edge_contribution = 0;

			for (int i = 0; i < (int) sources.size(); i++) {
				double M11 = Utils::intersection_size(sources[i].patterns, targets[j].patterns); // Number of patterns i and j are in
				double M01 = targets[j].patterns.size() - M11; // Number of patterns in which i is not in and j is in
				double M00 = (patterns - sources[i].patterns.size()) - M01; // Number of patterns i and j are not in
				double M10 = (patterns - targets[j].patterns.size()) - M00; // Number of patterns in which i is in and j is not in

				double weight = log((F(M11) * F(M00)) / (F(M10) * F(M01)));
				double infinity_weight = G(M10) + G(M01) - G(M11) - G(M00);

				if (infinity_weight > 2 || infinity_weight < -2) {
					clog << "ERROR!\n";
				}

				edge_contribution += log(F(M01) / F(M00));
				infinity_edge_contribution += (G(M00) - G(M01));

				targets[j].voltage += (pattern.count(&sources[i]) * weight);
				targets[j].infinity_voltage += (pattern.count(&sources[i]) * infinity_weight);
			}

			double M0 = patterns - targets[j].patterns.size();
			double M1 = targets[j].patterns.size();

			// IMPORTANT: Use sources.size() instead of targets[j].inneighbors.size() when working with a complete graph
			// 			  for which no inneighbours are stored!
			targets[j].activation_threshold = -((targets[j].inneighbors.size() - 1) * log(F(M0) / F(M1)) + edge_contribution);
			targets[j].infinity_activation_threshold = -((targets[j].inneighbors.size() - 1) * (G(M1) - G(M0)) + infinity_edge_contribution);

			targets[j].infinity_voltage = (targets[j].infinity_voltage == -0) ? 0 : targets[j].infinity_voltage;
			targets[j].infinity_activation_threshold = (targets[j].infinity_activation_threshold == -0) ? 0 : targets[j].infinity_activation_threshold;

			if (targets[j].infinity_voltage > targets[j].infinity_activation_threshold) {
				exposed.insert(&targets[j]);
			} else if (targets[j].infinity_voltage == targets[j].infinity_activation_threshold &&
				targets[j].voltage >= targets[j].activation_threshold) {
				exposed.insert(&targets[j]);
			}
		}

		return exposed;
	}
}