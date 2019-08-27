/*
#ifndef EDGE_WEAK_STRONG_H
#define EDGE_WEAK_STRONG_H
#include<vector>
#include<random>
#include "Edge.h"
using namespace std;

class Vertex;

class Edge_Weak_Strong : protected Edge {
  protected:
    double weight;
    double delay;
  public:
    Edge_Weak_Strong(Vertex& source, Vertex& target) : Edge(source,target);

    virtual Edge_Weak_Strong* clone() const {
      return Edge_Weak_Strong(*this);
    }
    double get_delay();
    bool is_strong();
    void set_weight(double w);
    double get_weight();
    bool is_reliable();
    bool is_present(double threshold);
};

#endif
*/
