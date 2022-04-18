/*
 * 190221_RNA_trimming_and_combining.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: anna
 *
 *      Inputs: fastq R1 & R2
 *      argv[1]: R1: 18 bases, where first 8 are UMI and next 10 are a unique index
 *      argv[2]: R2: 52 bases sequence
 *      argv[3]: An outfile;
 *      argv[4]: PCR_well;
 *
 *      Output: a single fastq that includes the 18 base sequence in the name of the 52 base one
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <ostream>
#include <stdio.h>
#include <string.h>

using namespace std;

int main(int argc, char* argv[]){
	string PCR_well = argv[4];

	ifstream(fastq_R1);
	fastq_R1.open(argv[1]);

	ifstream(fastq_R2);
	fastq_R2.open(argv[2]);

	ofstream(combined_fastq);
	combined_fastq.open(argv[3]);

	string line_R2;
	string partial_line_R2;
	string line_R1;
	int four_line_counter = 4;

	string R2_firstline;
	string R2_secondline;
	string R2_thirdline;
	string R2_fourthline;
	string R1_secondline;

	bool polyA = false;
	bool starts_with_polyA = false;
	int pos_of_polyA;

	while(getline(fastq_R2, line_R2))
	{
		getline(fastq_R1, line_R1);
		if (line_R2[0] == '@'  && four_line_counter == 4)
		{
			four_line_counter = 1;
			//partial_line_R2 = line_R2.substr(0, line_R2.find(" "));
		}
		else if(four_line_counter == 1)
		{
			R2_firstline = "@" + PCR_well + "_" + line_R1.substr(0,8) + "_" + line_R1.substr(8,10);
			R2_secondline = line_R2.substr(0,52);
			if (R2_secondline.find("AAAAAAAAAAAA") != string::npos)
			{
				if (R2_secondline.find("AAAAAAAAAAAA") != 0)
				{
					pos_of_polyA = R2_secondline.find("AAAAAAAAAAAA");
					R2_secondline = R2_secondline.substr(0, pos_of_polyA);
					polyA = true;
				}
				else
				{
					starts_with_polyA = true;
				}
			}
			four_line_counter++;
		}
		else if(four_line_counter == 2)
		{
			R2_thirdline = line_R2;
			four_line_counter++;
		}
		else if(four_line_counter == 3)
		{
			if (polyA == false && starts_with_polyA == false)
			{
				R2_fourthline = line_R2.substr(0,52);
				combined_fastq << R2_firstline << endl;
				combined_fastq << R2_secondline << endl;
				combined_fastq << R2_thirdline << endl;
				combined_fastq << R2_fourthline << endl;
				four_line_counter++;
			}
			else if (polyA == true)
			{
				R2_fourthline = line_R2.substr(0, pos_of_polyA);
				combined_fastq << R2_firstline << endl;
				combined_fastq << R2_secondline << endl;
				combined_fastq << R2_thirdline << endl;
				combined_fastq << R2_fourthline << endl;
				four_line_counter++;
				polyA = false;
			}
			else if (starts_with_polyA == true)
			{
				four_line_counter++; //don't put these lines in the new file at all.
				starts_with_polyA = false;
			}
		}
	}

	return 0;
}


