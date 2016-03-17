#ifndef HEATMAP_H
#define HEATMAP_H

#include<iostream>
#include<string>
#include<vector>
#include<iterator>
#include<fstream>
#include<unordered_map>
#include<sstream>
#include"Row.h"

typedef struct {
	int start;
	int end;
	double value;
} WigPeak;


void convert(char * arg, std::string& str);
void convert(char * arg, int& i);

void getInput(int argc, char * argv[], std::string &wigfile, std::string &peakfile, int &binsize, int &windowsize, int &shift, int &heatwidth, int &sortwidth, int &sortplace);

void readInPeaks(std::vector<Row*>&heatmap, std::string peakfile, int binsize, int windowsize, int shift, int heatwidth);

void readInWig(std::unordered_map<std::string, std::vector<WigPeak>>& wigpeaks, std::string wigfile);

int bsearch(std::vector<WigPeak>& array, int start);
#endif
