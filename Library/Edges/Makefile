# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -c -g -Wall -O3 -std=c++11 #-stdlib=libc++

# the build target executable:
TARGET = Edge

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET).o $(TARGET).cpp

clean:
	$(RM) $(TARGET)
