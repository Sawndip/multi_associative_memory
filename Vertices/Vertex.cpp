#include "Vertex.h"
#include "../Edges/Edge.h"
#include "../Signals/Transmitting_Edge.h"
#include <iostream>
#include <algorithm>

using namespace std;

Vertex::Vertex(int id, int type, double refractory_period) {
	this->id = id;
	this->voltage = 0;
	this->infinity_voltage = 0;
	this->type = type;
	this->refractory_period = refractory_period;
	this->neighbors = vector<Edge>();
	this->inneighbors = vector<Edge*>();
	this->sleep_until = -1.0;
	this->spike_time = -1.0;
	this->patterns = set<int>();
	this->activation_threshold = 0;
	this->infinity_activation_threshold = 0;
}

Vertex::Vertex(int type, double refractory_period) {
	this->voltage = 0;
	this->infinity_voltage = 0;
	this->type = type;
	this->refractory_period = refractory_period;
	this->neighbors = vector<Edge>();
	this->inneighbors = vector<Edge*>();
	this->sleep_until = -1.0;
	this->spike_time = -1.0;
	this->patterns = set<int>();
	this->activation_threshold = 0;
	this->infinity_activation_threshold = 0;
}

Vertex::Vertex(int id, int type, double refractory_period, int voltage) {
	this->id = id;
	this->type = type;
	this->refractory_period = refractory_period;
	this->voltage = voltage;
	this->infinity_voltage = voltage;
	this->neighbors = vector<Edge>();
	this->inneighbors = vector<Edge*>();
	this->sleep_until = -1.0;
	this->spike_time = -1.0;
	this->patterns = set<int>();
	this->activation_threshold = 0;
	this->infinity_activation_threshold = 0;
}

int Vertex::get_type() {
	return type;
}

void Vertex::set_type(int type) {
	this->type = type;
}

double Vertex::get_voltage() {
	return voltage;
}

void Vertex::set_voltage(double voltage) {
	this->voltage = voltage;
}

double Vertex::get_infinity_voltage() {
	return infinity_voltage;
}

void Vertex::set_infinity_voltage(double voltage) {
	this->infinity_voltage = voltage;
}

int Vertex::get_id() {
	return id;
}

void Vertex::set_id(int id) {
	this->id = id;
}

double Vertex::get_refractory_period() {
	return refractory_period;
}

void Vertex::set_refractory_period(double sr) {
	refractory_period = sr;
}

double Vertex::sample_refractory_period() {
	random_device rd;
	mt19937 gen(rd());
	exponential_distribution<double> d(refractory_period);
	return d(gen);
}

double Vertex::get_transmission_delay() {
	return transmission_delay;
}

void Vertex::set_transmission_delay(double td) {
	transmission_delay = td;
}

double Vertex::sample_transmission_delay() {
	random_device rd;
	mt19937 gen(rd());
	exponential_distribution<double> d(transmission_delay);
	return d(gen);
}

double Vertex::get_sleep_until() {
	return sleep_until;
}

void Vertex::set_sleep_until(double t) {
	sleep_until = t;
}

bool Vertex::awake(double t) {
	return t > sleep_until;
}

vector<double>* Vertex::get_history() {
	return &spike_times;
}

void Vertex::add_neighbor(Edge& e) {
	this->neighbors.push_back(e);
}

void Vertex::spike(double t) {
	spike_time = t;
}

double Vertex::get_activation_threshold() {
	return activation_threshold;
}

void Vertex::set_activation_threshold(double threshold) {
	this->activation_threshold = threshold;
}

double Vertex::get_infinity_activation_threshold() {
	return infinity_activation_threshold;
}

void Vertex::set_infinity_activation_threshold(double threshold) {
	this->infinity_activation_threshold = threshold;
}

void Vertex::expose(double time, priority_queue<Transmitting_Edge>& Q,
		double threshold) {
	this->spike(time);
	for (int i = 0; i < (int) neighbors.size(); ++i) {
		double delivery_time = time + this->neighbors[i].get_delay();
		this->neighbors[i].delivery_time = delivery_time;
		if (this->neighbors[i].works(threshold)) {
			Transmitting_Edge tmp = Transmitting_Edge(delivery_time, this->neighbors[i]);
			Q.push(tmp);
		}
	}
}

/****ONE ROUND PERCOLATION*/
void Vertex::exposeoneround(double time, int round, priority_queue<Transmitting_Edge>& Q,
		double threshold, int pattern) {
	this->spike(time);
	for (int i = 0; i < (int) neighbors.size(); ++i) {
		double delivery_time = time + this->neighbors[i].get_delay();
		this->neighbors[i].delivery_time = delivery_time;
		if (this->neighbors[i].works(threshold)) {
			Transmitting_Edge tmp = Transmitting_Edge(delivery_time, this->neighbors[i]);
			//Q.push(tmp);
		}
	}
}
/**
 * @param	round		round in synchronous percolation
 * @param	Q			priority queue for edges
 * @param	threshold	threshold for edge working or not
 * @param	pattern		the number of revealed patterns
 */
void Vertex::expose(double time, int round, priority_queue<Transmitting_Edge>& Q,
		double threshold, int pattern) {
	for (int i = 0; i < (int) neighbors.size(); ++i) {
		double delivery_time = time + this->neighbors[i].get_delay();
		this->neighbors[i].delivery_time = delivery_time;
		if (this->neighbors[i].works(threshold) && this->neighbors[i].get_pattern() <= pattern) {
			Transmitting_Edge tmp = Transmitting_Edge(delivery_time, round, this->neighbors[i]);
			Q.push(tmp);
		}
	}
}
