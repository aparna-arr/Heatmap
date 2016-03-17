CXX := g++
CXXFLAGS := -Wall -std=c++11 -g
OBJECTS := Main.o Heatmap.o Row.o

heatmap2: $(OBJECTS) 
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@	

Main.o: Main.cpp Heatmap.h Row.h
	$(CXX) $(CXXFLAGS) -c Main.cpp

Heatmap.o: Heatmap.cpp Heatmap.h Row.h
	$(CXX) $(CXXFLAGS) -c Heatmap.cpp

Row.o: Row.cpp Row.h
	$(CXX) $(CXXFLAGS) -c Row.cpp
	
clean:
	rm -f $(OBJECTS) heatmap2
