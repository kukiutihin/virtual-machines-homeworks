CXX = clang++
CXXFLAGS = -O2 -stdlib=libc++

cache: cache-characteristics/main.cpp 
	$(CXX) $(CXXFLAGS) cache-characteristics/main.cpp -o cache 

.PHONY : clean
clean :
	-rm cache 
