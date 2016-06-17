#ifndef ROW_H
#define ROW_H

#include<iostream>
#include<string>
#include<fstream>

double round(double num);

class Row 
{
	public:
	Row();
	Row(int center, int bin_size, int size, std::string chrom, int window, int shift);	
	bool map(int start, int bins, double value); // return True to go on to next wig peak
	double smooth(int sortplace, int sortwidth); // returns sort number
	std::string getChr(void);	
	int getStartPos(void);
	void printToFP(std::ofstream& file);

	private:
	int binsize;
	int start_pos;
	int num_bins;
	int smooth_win;
	int smooth_shift;
	std::string chr;

	double * array;
};


#endif
