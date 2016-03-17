#include"Row.h"
using namespace std;

Row::Row() {};

Row::Row(int center, int bin_size, int size, std::string chrom, int window, int shift) 
{
	start_pos = center - (size/2)*bin_size;
	num_bins = size; 
	smooth_win = window;
	smooth_shift = shift;
	chr = chrom;
	binsize = bin_size;
	array = new double[size];

	for (int i = 0; i < size; i++)
		array[i] = 0.0;		
}
	
bool Row::map(int start, int bins, double value)
{
	if (start - start_pos > num_bins * binsize)
		return false;

	if (start + bins * binsize <= start_pos)
		return true;	

	int bin_start = (start - start_pos) / binsize;

	if (bin_start < 0)
		bin_start = 0;

	int bin_end = num_bins;

	if (bins + bin_start < num_bins - 1)
		bin_end = bin_start + bins;

	for (int i = bin_start; i < bin_end; i++)
		array[i] += value;	

	return true;
}

double Row::smooth(int sortplace, int sortwidth)
{
	int sortstart, sortend;

	if (sortplace == 0)
	{
		sortstart = num_bins / 2 - sortwidth / 2;
		sortend = num_bins / 2 + sortwidth / 2;
	}

// fast! O(n)!!!
	double ring_buf[smooth_win];
	
	for (int i = 0; i < smooth_win; i++)
		ring_buf[i] = 0.0;

	double curr_sum = 0.0;

	for (int i = 0; i < num_bins; i+= smooth_shift)
	{
		double tmp = ring_buf[i % smooth_win];
		curr_sum = curr_sum - tmp + array[i];
		ring_buf[i % smooth_win] = array[i];
		array[i] = (double)curr_sum / smooth_win;
	}

	for (int j = num_bins - num_bins % smooth_shift; j < num_bins; j++)
	{
		double tmp = ring_buf[j % smooth_win];
		curr_sum = curr_sum - tmp + array[j];
		ring_buf[j % smooth_win] = 0.0;
		array[j] = (double)curr_sum / smooth_win;
	}

	double sort_num = 0.0;

	for (int s = sortstart; s < sortend; s++)
		sort_num += array[s];

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
