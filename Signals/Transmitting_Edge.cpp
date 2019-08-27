#include "Transmitting_Edge.h"
#include "../Vertices/Vertex.h"

Transmitting_Edge::Transmitting_Edge(double time, Edge& edge) {
	this->time = time;
	this->edge = &edge;
	this->round = 0;
}

Transmitting_Edge::Transmitting_Edge(double time, int round, Edge& edge) {
	this->time = time;
	this->edge = &edge;
	this->round = round;
}

void Transmitting_Edge::update(priority_queue<Transmitting_Edge>& Q,
		set<Vertex*>& Exposed, int K, double threshold) {
	Vertex *target = this->edge->target, *source = this->edge->source;
	if (Exposed.count(target) == 0) {
		target->voltage += source->type * this->edge->get_weight();
		if (target->voltage >= K) {
            target->spike(this->time);
			target->expose(this->time, Q, threshold);
			Exposed.insert(target);
		}
	}
}

/**
 * @param	Q			priority queue for edges
 * @param	exposed		thus far exposed vertices
 * @param	K_e			activation threshold for excitatory vertices
 * @param	K_i			activation threshold for inhibitory vertices
 * @param	threshold	threshold for edge working or not
 * @param	pattern		the number of revealed patterns
 */
void Transmitting_Edge::update(priority_queue<Transmitting_Edge>& Q,
		set<Vertex*>& Exposed, int K_e, int K_i, double threshold,
		int pattern) {
	Vertex *target = this->edge->target, *source = this->edge->source;
	if (Exposed.count(target) == 0 && this->edge->get_pattern() <= pattern) {
		target->voltage += source->type * this->edge->get_weight();
		if ((target->type == -1 && target->voltage >= K_i)
				|| (target->type == 1 && target->voltage >= K_e)) {
			target->spike(this->time);
			target->expose(this->time, this->round + 1, Q, threshold, pattern);
			Exposed.insert(target);
		}
	}
}

/**
/**** ONE ROUND UPDATE
 * @param	Q			priority queue for edges
 * @param	exposed		thus far exposed vertices
 * @param	K_e			activation threshold for excitatory vertices
 * @param	K_i			activation threshold for inhibitory vertices
 * @param	threshold	threshold for edge working or not
 * @param	pattern		the number of revealed patterns
 */
void Transmitting_Edge::updateoneround(priority_queue<Transmitting_Edge>& Q,
		set<Vertex*>& Exposed, int K_e, int K_i, double threshold,
		int pattern) {
	Vertex *target = this->edge->target, *source = this->edge->source;
	if (Exposed.count(target) == 0 && this->edge->get_pattern() <= pattern) {
		target->voltage += source->type * this->edge->get_weight();
		if ((target->type == -1 && target->voltage >= K_i)
				|| (target->type == 1 && target->voltage >= K_e)) {
			target->spike(this->time);
			target->exposeoneround(this->time, this->round + 1, Q, threshold, pattern);
			Exposed.insert(target);
		}
	}
}

// This function is just for the priority queue to work without needing to
// specify a comparison function.
bool operator<(const Transmitting_Edge& lhs, const Transmitting_Edge& rhs) {
	if (lhs.time == rhs.time) {
		if (lhs.round == rhs.round) {
			return lhs.edge->target->get_type() > rhs.edge->target->get_type();
		} else {
			return lhs.round > rhs.round;
		}
	} else {
		return lhs.time > rhs.time;
	}
}

