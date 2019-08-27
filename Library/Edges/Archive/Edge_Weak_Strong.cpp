/*
#include "Edge_Weak_Strong.h"
#include<random>
#include "../Vertices/Vertex.h"

Edge_Weak_Strong::Edge_Weak_Strong(Vertex& source, Vertex& target) : Edge(source,target) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::exponential_distribution<> d(source.get_spike_rate());
  this->delay = d(gen);
  this->weight = 0.0;
}

double Edge_Weak_Strong::get_delay() {
    return delay;
}

bool Edge_Weak_Strong::is_strong() {
    return weight > 0.0;
}

void Edge_Weak_Strong::set_weight(double w) {
    weight = w;
}

double Edge_Weak_Strong::get_weight() {
    return weight;
}

bool Edge_Weak_Strong::is_reliable() {
    return true;
}

bool Edge_Weak_Strong::is_present(double threshold) {
    return true;
}

*/
