/*
 * 191226_CROPt_process_atac_bedfile.cpp
 *
 *  Created on: Dec 26, 2019
 *      Author: anna
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
	string cell_name;
	map<string,int> chr_and_interval_map;
	int total_read_count;
	int dedup_read_count;
	int dedup_no_MT;
};

int main (int argc, char* argv[])
{
	ifstream(bed_file);
	bed_file.open(argv[1]);

	ofstream(main_out);
	main_out.open(argv[2]);

	ofstream(map_to_diff_chr_output);
	map_to_diff_chr_output.open(argv[3]);

	ofstream(cell_umi_count);
	cell_umi_count.open(argv[4]);

	int min_mapq = stoi(argv[5]);
	int max_len = stoi(argv[6]);

	map<string,int> map_of_cells;
	vector<cell> vector_of_cells;
	cell temp_cell_struct;
	int temp_pos_in_vect;

	string line;
	string elem;
	stringstream ss;
	int counter;

	string temp_first_chr;
	string temp_second_chr;
	string temp_start;
	string temp_end;
	string temp_cell;
	string temp_chr_start_end;
	int temp_mapq;
	int temp_length;

	while(getline(bed_file, line))
	{
		ss << line;
		counter = 0;
		while(getline(ss, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_first_chr = elem;
				counter++;
			}
			else if (counter == 1)
			{
				temp_start = elem;
				counter++;
			}
			else if (counter == 3)
			{
				temp_second_chr = elem;
				counter++;
			}
			else if (counter == 5)
			{
				temp_end = elem;
				temp_length = stoi(temp_end) - stoi(temp_start);
				counter++;
			}
			else if (counter == 6)
			{
				temp_cell = elem.substr(0,elem.find("."));
				counter++;
			}
			else if (counter == 7)
			{
				temp_mapq = stoi(elem);
				counter++;
			}
			else
			{
				counter++;
			}
		}
		if (temp_mapq <= min_mapq || temp_length > max_len ||  temp_first_chr != temp_second_chr)
		{
			map_to_diff_chr_output << line << endl;
		}
		else
		{
			if (map_of_cells.find(temp_cell) == map_of_cells.end())
			{
				temp_cell_struct.cell_name = temp_cell;
				temp_chr_start_end = temp_first_chr + ":" + temp_start + "-" + temp_end;
				temp_cell_struct.chr_and_interval_map[temp_chr_start_end] = 1;
				temp_cell_struct.dedup_read_count = 1;
				temp_cell_struct.total_read_count = 1;
				if (temp_first_chr != "MT"){
					temp_cell_struct.dedup_no_MT = 1;
				} else {
					temp_cell_struct.dedup_no_MT = 0;
				}

				vector_of_cells.push_back(temp_cell_struct);
				map_of_cells[temp_cell] = vector_of_cells.size() - 1;
				temp_cell_struct = cell();
			}
			else
			{
				temp_pos_in_vect = map_of_cells[temp_cell];
				temp_chr_start_end = temp_first_chr + ":" + temp_start + "-" + temp_end;
				if (vector_of_cells[temp_pos_in_vect].chr_and_interval_map.find(temp_chr_start_end) != vector_of_cells[temp_pos_in_vect].chr_and_interval_map.end())
				{
					vector_of_cells[temp_pos_in_vect].chr_and_interval_map[temp_chr_start_end]++;
					vector_of_cells[temp_pos_in_vect].total_read_count++;
				}
				else
				{
					vector_of_cells[temp_pos_in_vect].chr_and_interval_map[temp_chr_start_end] = 1;
					vector_of_cells[temp_pos_in_vect].total_read_count++;
					vector_of_cells[temp_pos_in_vect].dedup_read_count++;
					if (temp_first_chr != "MT"){
						vector_of_cells[temp_pos_in_vect].dedup_no_MT++;
					}
				}
			}
		}
		ss.clear();
	}

	int temp_col_pos;
	int temp_dash;
	string temp_chr;


	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		for (auto it = vector_of_cells[i].chr_and_interval_map.begin(); it != vector_of_cells[i].chr_and_interval_map.end(); it++)
		{
			temp_col_pos = it -> first.find(":");
			temp_dash = it -> first.find("-");
			temp_chr = it -> first.substr(0,temp_col_pos);
			temp_start = it -> first.substr(temp_col_pos+1, temp_dash - temp_col_pos - 1);
			temp_end = it -> first.substr(temp_dash + 1);
			main_out << vector_of_cells[i].cell_name << "\t" << temp_chr << "\t" << temp_start << "\t" << temp_end << endl;
		}
		cell_umi_count << vector_of_cells[i].cell_name << "\t" << vector_of_cells[i].dedup_read_count << "\t" << vector_of_cells[i].total_read_count << "\t" << vector_of_cells[i].dedup_no_MT << endl;
	}

	return 0;
}



