#include"Row.h"
#include"Heatmap.h"

using namespace std;

int main(int argc, char * argv[])
{
	// read in rows 
	vector<Row*> heatmap;
	unordered_map<string, vector<WigPeak>> wigpeaks;
	
	string wigfile, peakfile;
	int binsize, windowsize, shift, heatwidth, sortwidth, sortplace;

	try{
		getInput(argc, argv, wigfile, peakfile, binsize, windowsize, shift, heatwidth, sortwidth, sortplace);
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

//	cerr << "debug: Main(): before readInPeaks" << endl;
//	cerr << "debug: Main(): binsize is " << binsize << endl;
	readInPeaks(heatmap, peakfile, binsize, windowsize, shift, heatwidth);

//	cerr << "debug: Main(): after readInPeaks" << endl;

//	cerr << "debug: Main(): before readInWig" << endl;

	readInWig(wigpeaks, wigfile);

//	cerr << "debug: Main(): after readInWig" << endl;
	double sortNums[heatmap.size()];
	Row ** ptr_ar[heatmap.size()];

	string prev_chr = "INIT";
	for (vector<Row*>::iterator iter = heatmap.begin(); iter != heatmap.end(); iter++)
	{

		// initialize ptr array
		ptr_ar[iter - heatmap.begin()] = &(*iter);		

		string chr = (*iter)->getChr();

//		cerr << "debug: Main(): before looping through wig peaks" << endl;

//		for (vector<WigPeak>::iterator d = wigpeaks[chr].begin(); d != wigpeaks[chr].end(); d++)
//			cerr << "chr " << chr << " start " << (*d).start << " value " << (*d).value << endl;
//		cerr << "debug: Main(): after looping through wig peaks" << endl;

		if (prev_chr.compare("INIT") != 0 && chr.compare(prev_chr) != 0)
			wigpeaks.erase(prev_chr);
		
//		cerr << "Main(): chr is " << chr << endl;
//		cerr << "debug: Main(): before index search" << endl;
		int index = bsearch(wigpeaks[chr], (*iter)->getStartPos());
//		cerr << "debug: Main(): after index search startPos is " << (*iter)->getStartPos() << " wigpeak start is " << wigpeaks[chr][index].start << endl;
		
//		cerr << "index is " << index << endl;
//		cerr << "chr is " << chr << endl;
//		cerr << "wigpeaks size is " << wigpeaks[chr].size() << endl;
		//while loop here
		//return false only if peak is bigger than 
		while(index < (signed) wigpeaks[chr].size() && (*iter)->map(wigpeaks[chr][index].start, (wigpeaks[chr][index].end - wigpeaks[chr][index].start + 1) / binsize, wigpeaks[chr][index].value))
			index++;

//		cerr << "debug: after while mapping loop" << endl;
//		cerr << "Main(): sortplace is " << sortplace << " sortwidth is " << sortwidth << endl;
		sortNums[iter - heatmap.begin()] = (*iter)->smooth(sortplace, sortwidth/binsize);	
		prev_chr = chr;
	}	
	
	wigpeaks.clear();
	
	// ptr sort
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

	for (unsigned int i = 0; i < heatmap.size(); i++)
	{
		(*ptr_ar[i])->printToFP(file);
	}

	file.close();
	// for each peak in heatmap(in order)
	//	if new chromosome
	//		free prev wig chrom
	//
	// 	hash wig chr
	// 	bsearch starting wigpeak
	//	map wig regions
	//	smooth

	// free whole wig map
	// sort vector heatmap by ptr sort
	// print result(ptrsort-array, &heatmap)
}	
