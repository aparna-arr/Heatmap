#include"Heatmap.h"
using namespace std;

Peak::Peak(void)
{
	
}

Peak::Peak(int length, int start_pos)
{
	len = length;
	start = start_pos;
	end = start_pos + len;

	points = new double[len];

	for (int i = 0; i < len; i++)
		points[i] = 0.0;
}

Peak::Peak(int length, int start_pos, string chrom)
{
	len = length;
	start = start_pos;
	end = start_pos + len;
	chr = chrom;

	points = new double[len];

	for (int i = 0; i < len; i++)
		points[i] = 0.0;
}

void Peak::addSignal(int offset, int length, double value)
{	
//	cerr << "DEBUG in addSignal offset is " << offset << " length is " << length << " len is " << len << endl;

//	cerr << "offset is " << offset << " length is " << length << " value is " << value << endl;
	for (int i = offset; i < length + offset; i++)
	{
//		cerr << "\tpoints[i] is " << points[i] << endl;
		points[i] += (double)value;
//		cerr << "\tpoints[i] is " << points[i] << endl;
	}
}

int Peak::getStart(void)
{
	return start;
}

int Peak::getEnd(void)
{
	return end;
}

int Peak::getLen(void)
{
	return len;
}

string Peak::getChr(void)
{
	return chr;
}

void Peak::print(ofstream &out)
{
	for (int i = 0; i < len; i++)
		out << points[i] << "\t";

	out << "\n";
}

/* expects wigPeaks and bedPeaks to be of the SAME CHROMOSOME and SORTED!*/

void calc(vector<WigPeak> * wigPeaks, vector<Peak> &bedPeaks)
{
//	cerr << "DEBUG calc start" << endl;
	vector<Peak>::iterator bedIter = bedPeaks.begin();
	vector<WigPeak>::iterator wigIter = wigPeaks->begin();
//	cerr << "DEBUG got iters" << endl;
	while(bedIter != bedPeaks.end())
	{
//		cerr << "\tDEBUG in first while" << endl;
		while(wigIter != wigPeaks->end() && wigIter->end < bedIter->getStart())
			wigIter++;	
//		cerr << "\tDEBUG after first inner while" << endl;
	
		while(wigIter != wigPeaks->end() && wigIter->start < bedIter->getEnd())
		{
			int offset = 0;
			int len = 0;

			if (wigIter->start > bedIter->getStart())
			{
//				cerr << "\t\tDEBUG in if" << endl;
				offset = wigIter->start - bedIter->getStart();

				if (wigIter->end >= bedIter->getEnd())
					len = bedIter->getEnd() - wigIter->start;
				else
					len = wigIter->end - wigIter->start;

//				cerr << "\t\tDEBUG end of if" << endl;
			}
			else
			{
//				cerr << "\t\tDEBUG in else" << endl;
				offset = 0;
				
				if (wigIter->end > bedIter->getEnd())
					len = bedIter->getLen();
				else
					len = wigIter->end - bedIter->getStart();
//				cerr << "\t\tDEBUG end of else" << endl;
			}

			bedIter->addSignal(offset, len, wigIter->value);
			wigIter++;
		}
//		cerr << "\tDEBUG after second inner while" << endl;
		bedIter++;
	}
}

UserOpts::UserOpts(int argc, char * argv[])
{
	wigFile = "";
	bedFile = "";

	if (argc < 3)
	{	
		printUsage();
		throw 1;
	}
	
	vector<string> args = handleOpts(argc, argv);

	for (vector<string>::iterator iter = args.begin(); iter != args.end(); iter++)
	{
		if (iter == args.begin()) // on wigfile
			wigFile = *iter;
		else
		{
			// only accepting 1 bedfile right now
			bedFile = *iter;
		}
	}	
}

vector<string> UserOpts::handleOpts(int argc, char * argv[])
{
	vector<string> args;

	for (int i = 1; i < argc; i++)
	{
		// no opts handled right now
		args.push_back(string(argv[i]));
	}

	return args;
}

void UserOpts::printUsage(void)
{
	cout << "usage: heatmap [opts] <wigfile> <processed bedfile of peaks>" << endl;
	cout << "opts are: " << endl;
}

string UserOpts::getWigFilename(void)
{
	return wigFile;
}

string UserOpts::getBedFilename(void)
{
	return bedFile;
}

stack< vector<WigPeak> * > readInWig(string filename)
{
	ifstream file(filename.c_str());
	string line;

	int currSpan = 0;
	string currChr = "INIT";

	vector<WigPeak> * currVector = new vector<WigPeak>;
	stack< vector<WigPeak> * > allChroms;

	while(getline(file, line))
	{
		stringstream linestream(line);
		stringstream test(line);
		int pos;	
		double val;

		if (!(test >> pos))
		{
			string substring;	
			linestream >> substring; // variableStep
			linestream >> substring; // chrom=
			
			string chr = substring.substr(substring.find("=") + 1, substring.length());

			linestream >> substring; // span=

			string span = substring.substr(substring.find("=") + 1, substring.length());

			if (chr != currChr)
			{
				if (currChr != "INIT")
				{
					allChroms.push(currVector);
					currVector = new vector<WigPeak>;
				}	
				currChr = chr;
			} // if
			stringstream iToS(span);	
			iToS >> currSpan;
		} // if		
		else
		{
			linestream >> pos >> val;
			WigPeak tmp;
			tmp.start = pos;
			tmp.end = pos + currSpan;
			tmp.chr = currChr;
			tmp.value = val;

			currVector->push_back(tmp);
		}
	}
	file.close();
	allChroms.push(currVector);
	return allChroms;
}

vector< vector<Peak> > readInBed(string filename)
{
	ifstream file(filename.c_str());
	string line;

	string currChr = "INIT";
	vector<Peak> * currVector = new vector<Peak>;
	vector< vector<Peak> > allPeaks;

	while(getline(file, line))
	{
		int start, end;
		string chr;
		stringstream linestream(line);

		linestream >> chr >> start >> end;
	
		if (chr != currChr)
		{
			if (currChr != "INIT")
			{
				allPeaks.push_back(*currVector);	
				currVector = new vector<Peak>;
			}
			currChr = chr;
		}	
		currVector->push_back(Peak(end - start, start, chr));
	}	
	allPeaks.push_back(*currVector);	

	file.close();	
	return allPeaks;
}
