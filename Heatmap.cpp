#include"Heatmap.h"
using namespace std;

void getInput(int argc, char * argv[], std::string &wigfile, std::string &peakfile, int &binsize, int &windowsize, int &shift, int &heatwidth, int &sortwidth, int &sortplace)
{
	if (argc < 9)
		throw(1);

	convert(argv[1], wigfile);
	convert(argv[2], peakfile);
	convert(argv[3], binsize);
	convert(argv[4], windowsize);
	convert(argv[5], shift);
	convert(argv[6], heatwidth);
	convert(argv[7], sortplace);
	convert(argv[8], sortwidth);
}

void convert(char * arg, string& str) 
{
	stringstream ss(arg);

	if (!(ss >> str))
		throw(1);
}

void convert(char * arg, int& i) 
{
	stringstream ss(arg);

	if (!(ss >> i))
		throw(1);
}

void readInPeaks(std::vector<Row*>&heatmap, std::string peakfile, int binsize, int windowsize, int shift, int heatwidth)
{
	ifstream file;
	file.open(peakfile);

	string line;

	while(getline(file, line) && !file.bad())
	{
		string chr;
		int start, end;

		stringstream ss(line);
		ss >> chr >> start >> end;

		int center = (start + end) / 2;

		heatmap.push_back(new Row(center, binsize, heatwidth/binsize, chr, windowsize/binsize, shift/binsize));		
	}	

	file.close();
}

void readInWig(std::unordered_map<std::string, vector<WigPeak>>& wigpeaks, std::string wigfile)
{
	ifstream file;	
	file.open(wigfile);

	string line;
	string prev_chr = "INIT";
	string chr;
	int span;

	while(getline(file, line) && !file.bad())
	{
		stringstream ss(line);

		if (line.find("variableStep") != string::npos)
		{
			int pos_start = line.find("=");	
			int pos_end = line.find(" ", pos_start+1);
			
			chr = line.substr(pos_start + 1, pos_end - (pos_start + 1));

			pos_start = line.find("=", pos_end);

			string spanstr = line.substr(pos_start + 1);
			stringstream conv(spanstr);
			conv >> span;

			if (chr.compare(prev_chr) != 0)
			{
				vector<WigPeak> tmp;
				wigpeaks.insert({{chr, tmp}});
				prev_chr = chr;
			}
		}
		else
		{		
			int pos;
			double value;
	
			ss >> pos >> value;
	
			WigPeak tmp;
			tmp.start = pos;
			tmp.end = pos + span;
			tmp.value = round(value);
			wigpeaks[chr].push_back(tmp);	
		}
	}
	
	file.close();
}

int bsearch(vector<WigPeak>& array, int start)
{
	int min_index = 0;
	int max_index = array.size() - 1;
	int med_index = 0;

	while(min_index <= max_index)
	{
		med_index = (max_index + min_index) / 2;
		
		if (array[med_index].start < start)
			min_index = med_index + 1;
		else if (array[med_index].start > start)
			max_index = med_index - 1;
		else
			return med_index;
	}

	if (med_index > 0 && array[med_index - 1].end > start)
		med_index--;
	else if (med_index < (signed) array.size() - 1 && array[med_index + 1].end < start)
		med_index++;

	return med_index;
}
