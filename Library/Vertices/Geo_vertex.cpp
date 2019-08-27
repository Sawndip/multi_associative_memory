/*
 * Geo_vertex.cpp
 *
 *  Created on: Apr 14, 2015
 *      Author: Trujic
 */

#include "Geo_vertex.h"

#include <vector>
#include <random>

#include <iostream>

using namespace std;

Geo_vertex::Geo_vertex(int type, double refractory_period, double x, double y) : Vertex(type, refractory_period) {
	this->x = x;
	this->y = y;
}

Geo_vertex::Geo_vertex(int id, int type, double refractory_period, double x, double y) : Vertex(id, type, refractory_period) {
	this->x = x;
	this->y = y;
}

Geo_vertex::Geo_vertex(int id, int type, double refractory_period, int voltage, double x, double y) : Vertex(id, type, refractory_period, voltage) {
	this->x = x;
	this->y = y;
}

double Geo_vertex::get_x() {
	return this->x;
}

void Geo_vertex::set_x(double x) {
	this->x = x;
}

double Geo_vertex::get_y() {
	return this->y;
}

void Geo_vertex::set_y(double y) {
	this->y = y;
}
