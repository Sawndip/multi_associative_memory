#ifndef VERTEX_H
#define VERTEX_H
#include <vector>
#include <set>
#include <random>
#include <queue>
#include "../Edges/Edge.h"

using namespace std;

class Transmitting_Edge;

class Vertex {
	friend class Edge;
	friend class Transmitting_Edge;
private:
	int id;
	double refractory_period;
	double transmission_delay;
	double sleep_until;
	vector<double> spike_times;
public:
	double voltage; // This is just more convenient, we might change it later
	double infinity_voltage;
	int type;
	double activation_threshold;
	double infinity_activation_threshold;
	double spike_time;
	vector<Edge> neighbors;
	vector<Edge*> inneighbors;
	set<int> patterns;

	Vertex(int type, double refractory_period);
	Vertex(int id, int type, double refractory_period);
	Vertex(int id, int type, double refractory_period, int voltage);

	int get_type();
	void set_type(int type);
	double get_voltage();
	void set_voltage(double voltage);
	double get_infinity_voltage();
	void set_infinity_voltage(double voltage);
	int get_id();
	void set_id(int id);
	double get_refractory_period();
	void set_refractory_period(double sr);
	double sample_refractory_period();
	double get_transmission_delay();
	void set_transmission_delay(double td);
	double sample_transmission_delay();
	double get_sleep_until();
	void set_sleep_until(double t);
	bool awake(double t);
	void spike(double t);
	double get_activation_threshold();
	void set_activation_threshold(double threshold);
	double get_infinity_activation_threshold();
	void set_infinity_activation_threshold(double threshold);


	vector<double>* get_history();

	void add_neighbor(Edge& e);

	template<typename Distribution>
	void add_neighbor(Edge& e, Distribution& presence) {
		random_device rd;
		mt19937 gen(rd);
		add_neighbor(e, presence(gen));
	}

	void expose(double time, priority_queue<Transmitting_Edge>& Q,
			double threshold);
	void expose(double time, int round, priority_queue<Transmitting_Edge>& Q,
			double threshold, int k);
	/****ONE ROUND PERCOLATION*/
	void exposeoneround(double time, int round, priority_queue<Transmitting_Edge>& Q,
			double threshold, int k);

};

#endif
