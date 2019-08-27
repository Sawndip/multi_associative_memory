#include <random>
#include <algorithm>
#include <iostream>
#include <limits>

#include "../Vertices/Vertex.h"
#include "Edge.h"

Edge::Edge(Vertex& source, Vertex& target) {
	this->source = &source;
	this->target = &target;
	std::exponential_distribution<> d(1.0 / source.get_transmission_delay());
	if (source.get_type() == -1 || target.get_type() == -1) {
		pattern = 0;
	} else {
		pattern = numeric_limits<int>::max();
	}
	delivery_time = -1;
}

// Note that this is fixed, we might want to change this!
/*double Edge::get_delay() {
	random_device rd;
	mt19937 gen(rd());
	//std::exponential_distribution<> d(1.0);
	uniform_real_distribution<double> d(0.0, 1.0);
	return d(gen);
	//return double_properties[DELAY];
}*/

/*void Edge::set_delay(double delay) {
	double_properties[DELAY] = delay;
}*/

bool Edge::is_strong() {
	return double_properties[WEIGHT] > 0.0;
}

void Edge::set_weight(double w) {
	double_properties[WEIGHT] = w;
}

double Edge::get_weight() {
	return double_properties[WEIGHT];
}

bool Edge::is_reliable() {
	random_device rd;
	mt19937 gen(rd());
	if (d(gen)) {
		return true;
	} else {
		delivery_time = -1;
		return false;
	}
}

void Edge::set_reliability(double p) {
	global_reliability = p;
	d = bernoulli_distribution(p);
}

bool Edge::is_reliable_local() {
	if (double_properties.count(LOCAL_RELIABILITY) > 0) {
		random_device rd;
		mt19937 gen(rd());
		//d_local = bernoulli_distribution(double_properties[LOCAL_RELIABILITY]);
		return d_local(gen);
	} else {
		return true;
	}
}

void Edge::set_reliability_local(double p) {
	double_properties[LOCAL_RELIABILITY] = p;
	d_local = bernoulli_distribution(p);
}

bool Edge::is_present(double threshold) {
	return double_properties[PRESENCE] < threshold;
}

void Edge::set_presence(double presence) {
	double_properties[PRESENCE] = presence;
}

bool Edge::works(double threshold) {
	if (double_properties.count(WEIGHT) == 0) {
		return is_present(threshold) && is_reliable_local();
	} else {
		return is_present(threshold) && is_reliable_local() && is_strong();
	}
}

bool Edge::works() {
	if (double_properties.count(WEIGHT) == 0) {
		return is_reliable_local();
	} else {
		return is_reliable_local() && is_strong();
	}
}

void Edge::set_pattern(int k) {
	this->pattern = k;
}

int Edge::get_pattern() {
	return pattern;
}

void Edge::add_memory(int val) {
	memory_value = memory_value - memory[memory_idx] + val;
	memory[memory_idx] = val;
	memory_idx = (memory_idx + 1) % memory.size();
}

void Edge::add_substance(int val) {
	substance = substance * exp(-substance_decay) + val;
}

void Edge::set_substance_decay(double val) {
	substance_decay = val;
}

