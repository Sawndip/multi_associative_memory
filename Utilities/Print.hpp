/*
 * Print.hpp
 *
 *  Created on: May 5, 2015
 *      Author: Trujic
 */

#ifndef PRINT_HPP_
#define PRINT_HPP_

#include <vector>
#include <map>
#include <fstream>
#include <iostream>

namespace Print {
	// Print a vector
	template<typename T>
	void print(std::vector<T>& result) {
		std::cout << "{" << result[0];
		for (int i = 1; i < result.size(); ++i) {
			std::cout << ", " << result[i];
		}
		std::cout << "}\n";
	}

	// Print a vector with name
	template<typename T>
	void print(std::string s, std::vector<T>& result) {
		std::cout << s << " = {" << result[0];
		for (int i = 1; i < result.size(); ++i) {
			std::cout << ", " << result[i];
		}
		std::cout << "};\n";
	}

	// Print a vector with name to a file
	template<typename T>
	void print(std::string s, std::vector<T>& result, std::ofstream& file) {
		file << s << " = {" << result[0];
		for (int i = 1; i < result.size(); ++i) {
			file << ", " << result[i];
		}
		file << "};\n";
	}

	// Print a vector of pairs
	template<typename T, typename K>
	void print(std::vector<std::pair<T, K> >& result) {
		std::cout << "{{" << result[0].first << ", " << result[0].second << "}";
		for (int i = 1; i < result.size(); ++i) {
			std::cout << ", {" << result[i].first << ", " << result[i].second << "}";
		}
		std::cout << "};\n";
	}

	// Print a vector of pairs to a file
	template<typename T, typename K>
	void print(std::vector<std::pair<T, K> >& result, std::ofstream& file) {
		file << "{{" << result[0].first << ", " << result[0].second << "}";
		for (int i = 1; i < result.size(); ++i) {
			file << ", {" << result[i].first << ", " << result[i].second << "}";
		}
		file << "};\n";
	}

	// Print a vector of pairs with name
	template<typename T, typename K>
	void print(std::string s, std::vector<std::pair<T, K> >& result) {
		std::cout << s << " = {{" << result[0].first << ", " << result[0].second << "}";
		for (int i = 1; i < result.size(); ++i) {
			std::cout << ", {" << result[i].first << ", " << result[i].second << "}";
		}
		std::cout << "};\n";
	}

	// Print a vector of vectors with name
	template<typename T>
	void print(std::string s, std::vector<std::vector<T> >& result) {
		std::cout << s << "{";
		Print::print(result[0]);
		std::cout << "}";
		for (int i = 1; i < result.size(); ++i) {
			std::cout << ", ";
			Print::print(result[i]);
		}
		std::cout << "};\n";
	}

	// Print a vector of pairs of vectors
	template<typename T, typename K>
	void print(std::vector<std::pair<T, std::vector<K> > >& result) {
		std::cout << "data = {{" << result[0].first << ", {" << result[0].second[0];
		for (int i = 1; i < result[0].second.size(); ++i) {
			std::cout << "," << result[0].second[i];
		}
		std::cout << "}}";

		for (int i = 1; i < result.size(); ++i) {
			std::cout << ", {" << result[i].first << ", {" << result[i].second[0];
			for (int j = 0; j < result[i].second.size(); ++j) {
				std::cout << ", " << result[i].second[j];
			}
			std::cout << "}}";
		}
		std::cout << "};\n";
	}

	// Print a vector of pairs of vectors to a file
	template<typename T, typename K>
	void print(std::vector<std::pair<T, std::vector<K> > >& result, std::ofstream& file) {
		file << "data = {{" << result[0].first << ", {" << result[0].second[0];
		for (int i = 1; i < result[0].second.size(); ++i) {
			file << "," << result[0].second[i];
		}
		file << "}}";

		for (int i = 1; i < result.size(); ++i) {
			file << ", {" << result[i].first << ", {" << result[i].second[0];
			for (int j = 0; j < result[i].second.size(); ++j) {
				file << ", " << result[i].second[j];
			}
			file << "}}";
		}
		file << "};\n";
	}

	// Print a map of vectors with name
	template<typename K, typename T>
	void print(std::string s, std::map<K, std::vector<T> >& result) {
		auto it = result.begin();

		std::cout << s << " = {{" << (*it).first << ", {" << (*it).second[0];
		for (int i = 1; i < (*it).second.size(); ++i) {
			std::cout << ", " << (*it).second[i];
		}
		std::cout << "}}";
		it++;
		for (; it != result.end(); ++it) {
			std::cout << ", {" << (*it).first << ", {" << (*it).second[0];
			for (int j = 1; j < (*it).second.size(); ++j) {
				std::cout << ", " << (*it).second[j];
			}
			std::cout << "}}";
		}
		std::cout << "};\n";
	}

	//TODO: Find a way to remove this method?
	// Special print method
	template<typename T>
	void print(int idx, std::vector<T>& result, std::ofstream& file) {
		file << "data" << idx << " = {" << result[0].first;
		for (int i = 1; i < result.size(); ++i) {
			file << ", " << result[i].first;
		}
		file << "};\n";
		file << "data" << idx << "den = {" << result[0].second.first;
		for (int i = 1; i < result.size(); ++i) {
			file << ", " << result[i].second.first;
		}
		file << "};\n";
		file << "data" << idx << "chg = {" << result[0].second.second;
		for (int i = 1; i < result.size(); ++i) {
			file << ", " << result[i].second.second;
		}
		file << "};\n";
	}
}

#endif /* PRINT_HPP_ */