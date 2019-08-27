#ifndef EDGE_H
#define EDGE_H
#include<vector>
#include<random>
#include<utility>
#include<unordered_map>
#include<functional>

using namespace std;

class Vertex;

enum Property {
	WEIGHT, DELAY, PRESENCE, LOCAL_RELIABILITY
};

class Edge {
	friend class Transmitting_Edge;
private:
	int pattern;
protected:
	unordered_map<Property, double, hash<int> > double_properties;
	double global_reliability = 1.0;
public:
	bernoulli_distribution d, d_local;
	Vertex *source, *target;
	double delivery_time;

	vector<int> memory;
	function<double()> get_delay = [] () {return 1.0;};

	int memory_idx = 0, memory_value = 0;
    double substance = 0;
    double substance_decay = 0;

	Edge(Vertex& source, Vertex& target);

	bool is_strong();

	void set_weight(double w);
	double get_weight();

	bool is_reliable();
	void set_reliability(double p);

	bool is_reliable_local();
	void set_reliability_local(double p);

	bool is_present(double threshold);
	void set_presence(double p);

	bool works(double threshold);
	bool works(); // Ignoring presence

	void set_pattern(int k);
	int get_pattern();

	void add_memory(int val);

	void set_substance_decay(double val);
    void add_substance(int val);
};

#endif
