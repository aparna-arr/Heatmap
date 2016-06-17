#ifndef HEATMAP_H
#define HEATMAP_H

#include<iostream>
#include<string>
#include<stack>
#include<vector>
#include<iterator>
#include<exception>
#include<sstream>
#include<fstream>

class Peak
{
	public:
	Peak();
	Peak(int length, int start_pos);
	Peak(int length, int start_pos, std::string chrom);

	void addSignal(int offset, int length, double value);

	int getStart(void);
	int getEnd(void);
	int getLen(void);	
	std::string getChr(void);

	void print(std::ofstream &out);

	private:
	int start;
	int end;
	int len;
	std::string chr;

	double * points; // each bp of peak in window
};

typedef struct WigPeak 
{
	int start;
	int end;
	double value;
	std::string chr;
} WigPeak;

class UserOpts
{
	public:
	UserOpts() {};	
	UserOpts(int argc, char * argv[]);

	void printUsage(void);
	std::vector<std::string> handleOpts(int argc, char * argv[]);
	std::string getWigFilename(void);
	std::string getBedFilename(void);

	private:
	std::string wigFile;
	std::string bedFile;
};

void calc(std::vector<WigPeak> * wigPeaks, std::vector<Peak> &bedPeaks);
std::stack< std::vector<WigPeak> * > readInWig(std::string filename);
std::vector< std::vector<Peak> > readInBed(std::string filename);
#endif
