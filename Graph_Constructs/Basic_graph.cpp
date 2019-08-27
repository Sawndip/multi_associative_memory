#include <random>
#include <vector>

#include "Basic_graph.h"

using namespace std;

Basic_graph::Basic_graph(int n_ex, int n_in, double ref_p_pos,
		double ref_p_neg) {
	for (int j = 0; j < n_ex; ++j) {
		Vertex tmp = Vertex(j, 1, ref_p_pos);
		this->excitatory.push_back(tmp);
	}
	for (int i = 0; i < n_in; ++i) {
		Vertex tmp = Vertex(n_ex + i, -1, ref_p_neg);
		this->inhibitory.push_back(tmp);
	}
}

Basic_graph::Basic_graph(int n_ex, int n_in, double ref_p_pos, double ref_p_neg,
		int init_voltage) {
	for (int j = 0; j < n_ex; ++j) {
		Vertex tmp = Vertex(j, 1, ref_p_pos);
		tmp.voltage = init_voltage;
		this->excitatory.push_back(tmp);
	}
	for (int i = 0; i < n_in; ++i) {
		Vertex tmp = Vertex(n_ex + i, -1, ref_p_neg);
		tmp.voltage = init_voltage;
		this->inhibitory.push_back(tmp);
	}
}