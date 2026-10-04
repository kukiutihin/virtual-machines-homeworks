CXX = clang++
CXXFLAGS = -O2 -std=c++20

ifeq ($(shell uname), Darwin)
	CXXFLAGS += -isysroot $$(xcrun --show-sdk-path)
endif

cache: cache-characteristics/main.cpp 
	$(CXX) $(CXXFLAGS) cache-characteristics/main.cpp -o cache 

cache-debug: cache-characteristics/main.cpp 
	$(CXX) $(CXXFLAGS) cache-characteristics/main.cpp -g -o cache 

.PHONY : clean
clean :
	@rm -rf cache cache.dSYM/
