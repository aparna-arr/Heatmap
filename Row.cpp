#include"Row.h"
using namespace std;

Row::Row() {};

Row::Row(int center, int bin_size, int size, std::string chrom, int window, int shift) 
{
	start_pos = center - (size/2)*bin_size;

//	cerr << "debug: Row:Row() center is " << center << " startpos is " << start_pos << " size is " << size << endl;
	
	num_bins = size; 
	smooth_win = window;
	smooth_shift = shift;
	chr = chrom;

	binsize = bin_size;
	
	array = new double[size];

//	cerr << "debug(): Row::Row(): binsize is " << binsize << endl;
//	cerr << "size is " << size << endl;
	for (int i = 0; i < size; i++)
		array[i] = 0.0;		
}
	
bool Row::map(int start, int bins, double value)
{
//	cerr << "debug: map(): begin" << endl;

//	cerr << "binsize is " << binsize << " bins is " << bins << endl;

	if (start - start_pos > num_bins * binsize)
	{
//		cerr << "debug: map(): FALSE end" << endl;
		return false;
	}

	if (start + bins * binsize <= start_pos)
	{
//		cerr << "debug: map(): TRUE end" << endl;
		return true;	
	}

	int bin_start = (start - start_pos) / binsize;

	if (bin_start < 0)
		bin_start = 0;

	int bin_end = num_bins;

	if (bins + bin_start < num_bins - 1)
		bin_end = bin_start + bins;

//	cerr << "in map(): binstart is " << bin_start << " binend is " << bin_end << " value is " << value << endl;
	
	for (int i = bin_start; i < bin_end; i++)
		array[i] += value;	

//	cerr << "debug: map(): end" << endl;
	return true;
}

double Row::smooth(int sortplace, int sortwidth)
{
//	cerr << "start smooth" << endl;
//	cerr << "num_bins is " << num_bins << " smooth_shift is " << smooth_shift << endl;
//	cerr << "sortwidth is " << sortwidth << " sortplace is " << sortplace << endl;


//	cerr << "START DEBUG smooth()" << endl;
//	for (int d = 0; d < num_bins; d++)
//		cerr << array[d] << " ";

//	cerr << endl;
//	cerr << "END debug loop" << endl;

	int sortstart, sortend;

	if (sortplace == 0)
	{
		sortstart = num_bins / 2 - sortwidth / 2;
		sortend = num_bins / 2 + sortwidth / 2;
	}
// fast! O(n)!!!
	double ring_buf[smooth_win];
	
	// init
	for (int i = 0; i < smooth_win; i++)
		ring_buf[i] = 0.0;

	double curr_sum = 0.0;

	for (int i = 0; i < num_bins; i+= smooth_shift)
	{
		double tmp = ring_buf[i % smooth_win];
		curr_sum = curr_sum - tmp + array[i];
		
		ring_buf[i % smooth_win] = array[i];
		
//		cerr << "smooth(): curr_sum is " << curr_sum << " iter i is " << i << endl;
		
		array[i] = (double)curr_sum / smooth_win;
//		cerr << "smooth(): stored " << array[i] << endl;
	}

	for (int j = num_bins - num_bins % smooth_shift; j < num_bins; j++)
	{
		double tmp = ring_buf[j % smooth_win];
		curr_sum = curr_sum - tmp + array[j];
		
		ring_buf[j % smooth_win] = 0.0;
	
//		cerr << "smooth(): curr_sum is " << curr_sum << " iter j is " << j << endl;

		array[j] = (double)curr_sum / smooth_win;
//		cerr << "smooth(): stored " << array[j] << endl;
	}

	double sort_num = 0.0;

//	cerr << "sortstart is " << sortstart << " sortend is " << sortend << " num_bins is " << num_bins << endl;
	for (int s = sortstart; s < sortend; s++)
		sort_num += array[s];

	
/*
// from here on, EXTREMELY slow
// Attempt 1
	double runningAvg = 0.0;
	// init
	for (int i = 0; i < smooth_win; i++)
	{
		runningAvg += array[i];	
	}

	int i;
	for (i = 0; i < num_bins - smooth_shift; i+= smooth_shift)
	{
		for (int j = i; j < smooth_shift; j++)
		{	
			double tmp = array[i]; 
			array[i] = runningAvg / smooth_win;
			runningAvg -= tmp;
		}
		
		for (int k = i + smooth_win; k < i + smooth_win + smooth_shift ; i++)
		{
			runningAvg += array[k];
		}
	}
	
	// imperfect solution
	for (;i < num_bins; i++)
		array[i] = runningAvg;
*/	
//	cerr << "end smooth" << endl;
	return sort_num;
}	

std::string Row::getChr(void)
{
	return chr;
}

int Row::getStartPos(void)
{
	return start_pos;
}

void Row::printToFP(std::ofstream& file)
{
	for (int i = 0; i < num_bins; i++)
		file << array[i] << " ";

	file << endl;
}
