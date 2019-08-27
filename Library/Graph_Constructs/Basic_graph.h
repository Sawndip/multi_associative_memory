#ifndef BASIC_GRAPH_H
#define BASIC_GRAPH_H

#include <vector>
#include <random>
#include <queue>
#include <set>

#include "../Vertices/Vertex.h"
#include "../Edges/Edge.h"
#include "../Signals/Transmitting_Edge.h"

using namespace std;

class Basic_graph{
public:
    vector<Vertex> excitatory, inhibitory;

    Basic_graph(int n_ex, int n_in, double spike_rate_pos, double spike_rate_neg);
    Basic_graph(int n_ex, int n_in, double spike_rate_pos, double spike_rate_neg, int init_voltage);
};

#endif