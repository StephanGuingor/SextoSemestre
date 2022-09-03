CXX = g++
SRC_DIRS = .
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
EXEC = main

run: build
	./$(EXEC)

build: $(SRCS)
	$(CXX) $(SRCS) -o $(EXEC) --std=c++20

.PHONY: clean
clean:
	rm $(EXEC)