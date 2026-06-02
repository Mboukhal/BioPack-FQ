CXX      := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O2 -march=native -Iinclude
LDFLAGS  := -pthread
# Uncomment when implementing compression/hashing:
# LDFLAGS  := -lzstd -lxxhash -lcrypto -pthread

SRCS     := $(wildcard src/*.cpp)
OBJS     := $(SRCS:.cpp=.o)
TARGET   := biopack

.PHONY: all clean test run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

test: $(TARGET)
	$(CXX) $(CXXFLAGS) -o test/test_runner test/test_basic.cpp src/encoder.cpp
	./test/test_runner

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) test/test_runner
