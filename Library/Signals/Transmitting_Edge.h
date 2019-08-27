#ifndef TRANSMITTING_EDGE_H
#define TRANSMITTING_EDGE_H

#include<set>
#include<queue>
#include "../Edges/Edge.h"

class Transmitting_Edge {
public:
	Edge *edge;
	double time;
	int round;

	Transmitting_Edge(double time, Edge& edge);
	Transmitting_Edge(double time, int round, Edge& edge);

	void update(priority_queue<Transmitting_Edge>& Q, set<Vertex*>& Exposed,
			int K, double threshold);

	void update(priority_queue<Transmitting_Edge>& Q, set<Vertex*>& Exposed,
			int K_e, int K_i, double threshold, int pattern);

	friend bool operator<(const Transmitting_Edge& lhs,
			const Transmitting_Edge& rhs);
	void updateoneround(priority_queue<Transmitting_Edge>& Q,
		set<Vertex*>& Exposed, int K_e, int K_i, double threshold,
		int pattern);

};

#endif
