/*
 * Geometric_graph.cpp
 *
 *  Created on: Apr 14, 2015
 *      Author: Trujic
 */

#include <vector>
#include <random>
#include <iostream>

#include "Geometric_graph.h"
#include "../Utilities/Utils.hpp"

using namespace std;

/**
 * Geometric graph on a unit square
 */
Geometric_graph::Geometric_graph(int n_ex, int n_in, double ref_p_pos, double ref_p_neg) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> d(0.0, 1.0);

	for (int i = 0; i < n_ex; i++) {
		double x = d(gen);
		double y = d(gen);
		Geo_vertex v = Geo_vertex(i, 1, ref_p_pos, x, y);
		this->excitatory.push_back(v);
	}

	for (int i = 0; i < n_ex; i++) {
		double x = d(gen);
		double y = d(gen);
		Geo_vertex v = Geo_vertex(n_ex + i, -1, ref_p_pos, x, y);
		this->inhibitory.push_back(v);
	}
}

/**
 * Geometric graph on a unit square with initial voltage of vertices set
 */
Geometric_graph::Geometric_graph(int n_ex, int n_in, double ref_p_pos, double ref_p_neg, int init_voltage) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> d(0.0, 1.0);

	for (int i = 0; i < n_ex; i++) {
		double x = d(gen);
		double y = d(gen);
		Geo_vertex v = Geo_vertex(i, 1, ref_p_pos, x, y);
		this->excitatory.push_back(v);
	}

	for (int i = 0; i < n_ex; i++) {
		double x = d(gen);
		double y = d(gen);
		Geo_vertex v = Geo_vertex(n_ex + i, -1, ref_p_pos, x, y);
		this->inhibitory.push_back(v);
	}
}

/**
 * Geometric graph on a X by Y square
 */
Geometric_graph::Geometric_graph(int n_ex, int n_in, double ref_p_pos, double ref_p_neg, double X, double Y) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> d_x(0.0, X);
	uniform_real_distribution<> d_y(0.0, Y);

	for (int i = 0; i < n_ex; i++) {
		double x = d_x(gen);
		double y = d_y(gen);
		Geo_vertex v = Geo_vertex(i, 1, ref_p_pos, x, y);
		this->excitatory.push_back(v);
	}

	for (int i = 0; i < n_ex; i++) {
		double x = d_x(gen);
		double y = d_y(gen);
		Geo_vertex v = Geo_vertex(n_ex + i, -1, ref_p_pos, x, y);
		this->inhibitory.push_back(v);
	}
}
