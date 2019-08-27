# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CC = g++
CFLAGS  = -Wall -O3 -std=c++11 -fopenmp
#-stdlib=libc++

# the build target executable:
TARGET = find_optimal_p_for_each_z_and_k_s_fixed

all: $(TARGET)

$(TARGET): $(TARGET).cpp Vertex.o Geo_vertex.o Edge.o Transmitting_Edge.o Layered_graph.o Basic_graph.o Percolation.o Utils.o Connect.o
	$(CC) $(CFLAGS)  -o $(TARGET) $(TARGET).cpp Layered_graph.o Basic_graph.o Percolation.o Transmitting_Edge.o Edge.o Vertex.o Geo_vertex.o Utils.o Connect.o

Layered_graph.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Graph_Constructs/Layered_graph.cpp

Basic_graph.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Graph_Constructs/Basic_graph.cpp

Percolation.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Percolation/Percolation.cpp

Transmitting_Edge.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Signals/Transmitting_Edge.cpp

Vertex.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Vertices/Vertex.cpp

Geo_vertex.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Vertices/Geo_vertex.cpp

Edge.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Edges/Edge.cpp

Connect.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Utilities/Connect.cpp

Utils.o:
	$(CC) $(CFLAGS)  -c ../../../../../Library/Utilities/Utils.cpp

clean:
	$(RM) $(TARGET)
