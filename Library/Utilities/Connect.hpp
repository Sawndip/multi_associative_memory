/*
 * Connect.hpp
 *
 *  Created on: May 11, 2015
 *      Author: Trujic
 */

#ifndef CONNECT_HPP_
#define CONNECT_HPP_

#include <vector>
#include <map>

class Edge;
class Vertex;
class Geo_vertex;

namespace Connect {
	void connect(std::vector<Vertex>& sources, std::vector<Vertex>& targets, double density);
	void connect(std::vector<Vertex>& sources, std::vector<Vertex>& targets, double density, void (*process_edge)(Edge&));
	void connect(std::vector<Geo_vertex>& sources, std::vector<Geo_vertex>& targets, double radius, void (*process_edge)(Edge&));
	void connect_with_in(std::vector<Vertex>& sources, std::vector<Vertex>& targets, double density, void (*process_edge)(Edge&));
}

#endif /* CONNECT_HPP_ */