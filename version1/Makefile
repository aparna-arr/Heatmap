heatmap: CalculateHeatmap.o Heatmap.o
	g++ -Wall -o $@ CalculateHeatmap.o Heatmap.o
CalculateHeatmap.o: CalculateHeatmap.cpp Heatmap.h
	g++ -Wall -c CalculateHeatmap.cpp
Heatmap.o: Heatmap.cpp Heatmap.h
	g++ -Wall -c Heatmap.cpp
clean: 
	rm -f Heatmap.o CalculateHeatmap.o heatmap
