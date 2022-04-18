/*
 * 190223_sciRNA_remove_noncells_from_sam.cpp
 *
 *  Created on: Feb 23, 2019
 *      Author: anna
 *
 *      This program removes all cells below a certain UMI count based on a list of "real cell"
 *      RT indexes that are provided as input.
 *
 *      argv[1]: list of real indices
 *      argv[2]: sam file collapsed by umi
 *      argv[3]: sam outfile
 *      argv[4]: line index outfile
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

struct cell{
	string RT_index;
	vector <string> sam_lines;
};

int main(int argc, char* argv[])
{
	vector<cell> vector_of_cells;
	map<string,int> RT_index_to_pos_in_vect_map;
	cell temp_cell;

	vector<int>pos_of_each_cell_start_vect;

	ifstream(real_cell_RT_index_list);
	real_cell_RT_index_list.open(argv[1]);
	string index_line;

	while(getline(real_cell_RT_index_list,index_line))
	{
		temp_cell.RT_index = index_line;
		vector_of_cells.push_back(temp_cell);
		RT_index_to_pos_in_vect_map[index_line] = vector_of_cells.size() - 1;
	}

	ifstream(unique_reads_sam);
	unique_reads_sam.open(argv[2]);
	string sam_line;
	stringstream sam_line_ss;
	string curr_index;

	while(getline(unique_reads_sam, sam_line))
	{
		curr_index = sam_line.substr(sam_line.find('_') + 10, 10);
		if (RT_index_to_pos_in_vect_map.find(curr_index) != RT_index_to_pos_in_vect_map.end())
		{
			vector_of_cells[RT_index_to_pos_in_vect_map[curr_index]].sam_lines.push_back(sam_line);
		}
	}

	ofstream(sam_outfile);
	sam_outfile.open(argv[3]);

	ofstream(index_outfile);
	index_outfile.open(argv[4]);
	string temp_line;

	int counter = 1;
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		index_outfile << counter << endl;
		//cout << counter << ":" << temp_line << endl;
		for (int j = 0; j < vector_of_cells[i].sam_lines.size(); j++)
		{
			sam_outfile << vector_of_cells[i].sam_lines[j] << endl;
			temp_line = vector_of_cells[i].sam_lines[j];
			counter++;
		}
	}



	return 0;
}


