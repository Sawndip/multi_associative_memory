/*
 * Geo_vertex.h
 *
 *  Created on: Apr 14, 2015
 *      Author: Trujic
 */

#ifndef GEO_VERTEX_H_
#define GEO_VERTEX_H_

#include <vector>
#include "Vertex.h"

using namespace std;

class Geo_vertex: public Vertex {
private:
	double x, y; // Position of the vertex in the plane
public:
	Geo_vertex(int type, double refractory_peroid, double x, double y);
	Geo_vertex(int id, int type, double refractory_period, double x, double y);
	Geo_vertex(int id, int type, double refractory_peroid, int voltage, double x, double y);

	double get_x();
	void set_x(double x);
	double get_y();
	void set_y(double y);
};

#endif /* GEO_VERTEX_H_ */
