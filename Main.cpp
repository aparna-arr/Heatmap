#include"Row.h"
#include"Heatmap.h"

using namespace std;

int main(int argc, char * argv[])
{
	vector<Row*> heatmap;
	unordered_map<string, vector<WigPeak>> wigpeaks;
	
	string wigfile, peakfile;
	int binsize, windowsize, shift, heatwidth, sortwidth, sortplace;

	try{
		getInput(argc, argv, wigfile, peakfile, binsize, windowsize, shift, heatwidth, sortwidth, sortplace);
	
		readInPeaks(heatmap, peakfile, binsize, windowsize, shift, heatwidth);
	
		readInWig(wigpeaks, wigfile);
	
		double sortNums[heatmap.size()];
		Row ** ptr_ar[heatmap.size()];
	
		string prev_chr = "INIT";
		for (vector<Row*>::iterator iter = heatmap.begin(); iter != heatmap.end(); iter++)
		{
			ptr_ar[iter - heatmap.begin()] = &(*iter);		
	
			string chr = (*iter)->getChr();
	
			if (prev_chr.compare("INIT") != 0 && chr.compare(prev_chr) != 0)
				wigpeaks.erase(prev_chr);
			
			int index = bsearch(wigpeaks[chr], (*iter)->getStartPos());
		
			while(index < (signed) wigpeaks[chr].size() && (*iter)->map(wigpeaks[chr][index].start, (wigpeaks[chr][index].end - wigpeaks[chr][index].start + 1) / binsize, wigpeaks[chr][index].value))
				index++;
	
			sortNums[iter - heatmap.begin()] = (*iter)->smooth(sortplace, sortwidth/binsize);	
			prev_chr = chr;
		}	
		
		wigpeaks.clear();
		
		Row ** tmpElement;
		double tmpNum;
	
		for (unsigned int i = 0; i < heatmap.size(); i++)
		{
			int j = i;
			while(j > 0 && sortNums[j-1] > sortNums[j])
			{
				tmpNum = sortNums[j];
				tmpElement = ptr_ar[j];
	
				sortNums[j] = sortNums[j - 1];
				ptr_ar[j] = ptr_ar[j - 1];
	
				sortNums[j-1] = tmpNum;
				ptr_ar[j-1] = tmpElement;
	
				j--;
			}
		}
	
		ofstream file;
		file.open("heatmap2.txt");
	
		for (int i = (signed)heatmap.size() - 1; i >= 0; i--)
			(*ptr_ar[i])->printToFP(file);
	
		file.close();
	}
	catch(int e)
	{
		if (e == 1)
		{
			cout << "usage: heatmap2 wigfile bedfile bin_size window_size shift heatmap_width sort_place sort_width" << endl;
			cout << "EVERYTHING must be even and in multiples of binsize! sort_place can only equal 0 (center) right now." << endl;
			return 1;
		}
	}
}	
