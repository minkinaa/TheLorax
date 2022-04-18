/*
 * 190223_sciRNA_make_real_RT_index_list.cpp
 *
 *  Created on: Feb 23, 2019
 *      Author: anna
 *
 *      This program takes in the metrics file from 190223_sciRNA_remove_duplicates and makes a new file
 *      of RT indexes that correspond to real cells
 *
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

int main(int argc, char* argv[])
{
	ifstream(metrics_file);
	metrics_file.open(argv[1]);
	int UMI_cutoff = stoi(argv[2]);
	string line;
	stringstream line_ss;
	string elem_in_line;
	int unique_umi_count;
	string index;
	int counter;

	while(getline(metrics_file, line))
	{
		index = line.substr(0,10);
		line_ss << line;
		counter = 0;
		while(getline(line_ss, elem_in_line, '\t'))
		{
			if (counter == 1)
			{
				unique_umi_count = stoi(elem_in_line);
				if (unique_umi_count > UMI_cutoff)
				{
					cout << index << endl;
				}
			}
			counter++;
		}
		line_ss.clear();
	}

	return 0;
}


