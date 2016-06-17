#include"Heatmap.h"
using namespace std;

int main(int argc, char * argv[])
{	
	try {
		UserOpts input(argc, argv);

		cout << "Reading in wig" << endl;	
	
		stack< vector<WigPeak> * > wig = readInWig(input.getWigFilename());

		cout << "Reading in bed" << endl;
		
		vector< vector<Peak> > bed = readInBed(input.getBedFilename());

		cout << "done" << endl;

		// great, now have to find common chrs ...

		// O(n*n) horrible search for chroms
		vector< vector<Peak> > results;

		while (!wig.empty())
		{
			vector< vector<Peak> >::iterator iter;

			for (iter = bed.begin(); iter != bed.end(); iter++)
			{
				if (iter->begin()->getChr() == ((wig.top())->begin())->chr)
				{
					break;
				}
			}
			
			if (iter == bed.end())
			{
				wig.pop();
				continue;
			}

			// found a common chr
			// do stuff
		
			cout << "Found a chr!" << endl;

			cout << "starting calc()" << endl;

			calc(wig.top(), *iter);
	
			cout << "done" << endl;
			
			results.push_back(*iter);

			wig.pop();
			bed.erase(iter);
		}

		cout << "printing results" << endl;

		ofstream outfile("heatmap_outfile.txt");
			
		outfile << "peaks\t";

		for (int i = 0; i < results.begin()->begin()->getLen(); i++)
			outfile << i << "\t";

		outfile << "\n"; 

		for (vector< vector<Peak> >::iterator it = results.begin(); it != results.end(); it++)
		{
			for (vector<Peak>::iterator iter = it->begin(); iter != it->end(); iter++)
			{
				outfile << iter->getChr() << "_" << iter->getStart() << "\t";
				iter->print(outfile);		
			}
		}

		outfile.close();

		cout << "done" << endl;
	}
	catch (int e) {
		if (e == 1)
		{
			return 1;
		}
	}			
}
