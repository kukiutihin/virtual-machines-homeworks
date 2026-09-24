CXX = clang++
CXXFLAGS = -O2 

ifeq ($(shell uname), Darwin)
	CXXFLAGS += -stdlib=libc++ -isysroot $$(xcrun --show-sdk-path)
endif

cache: cache-characteristics/main.cpp 
	$(CXX) $(CXXFLAGS) cache-characteristics/main.cpp -o cache 

.PHONY : clean
clean :
	-rm cache 
