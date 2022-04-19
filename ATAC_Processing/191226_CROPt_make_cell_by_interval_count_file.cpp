/*
 * 191226_CROPt_make_cell_by_interval_count_file.cpp
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
#include <math.h>

using namespace std;

struct cell{
	string cell_name;
	map<pair<string, int>, int> interval_in_map_num_to_count;
	int total_count;
	int new_cell_num;
};

int main(int argc, char* argv[])
{
	int interval_size = stoi(argv[1]);

	ifstream(hg38_chr_info_file);
	hg38_chr_info_file.open(argv[2]);

	ifstream(cell_chr_start_end_file);
	cell_chr_start_end_file.open(argv[3]);

	//ifstream(real_cell_file);
	//real_cell_file.open(argv[4]);

	int real_cell_cutoff = stoi(argv[4]);

	int min_count_per_interval = stoi(argv[5]);

	ofstream(matrix_outfile);
	matrix_outfile.open(argv[6]);

	ofstream(intervals_outfile);
	intervals_outfile.open(argv[7]);

	ofstream(barcodes_outfile);
	barcodes_outfile.open(argv[8]);

	ofstream(cell_interval_count_outfile);
	cell_interval_count_outfile.open(argv[9]);


	//generate interval_vector;
	map<pair<string,int>,int> interval_to_pos_in_vect;
	map<pair<string,int>,int> total_count_per_interval_map; //we'll only want to include those intervals which aren't 0;

	string line;
	string elem;
	stringstream ss;
	int counter;

	string temp_chr;
	int temp_length;

	int num_intervals_for_chromosome;

	getline(hg38_chr_info_file, line); //this is top line;

	cout << "making it here?" << endl;

	while(getline(hg38_chr_info_file, line))
	{
		ss << line;
		counter = 0;
		while(getline(ss, elem, ' '))
		{
			if (counter == 0)
			{
				temp_chr = elem;
				counter++;
			}
			else if (counter == 1)
			{
				cout << temp_length << endl;
				temp_length = stoi(elem);
				counter++;
			}
			else
			{
				counter++;
			}
		}
		num_intervals_for_chromosome = ceil(double(temp_length)/double(interval_size));
		cout << "chr" << temp_chr << ":" << num_intervals_for_chromosome << " bins" << endl;
		int temp_start = 1;
		for (int i = 0; i < num_intervals_for_chromosome + 1; i++)
		{
			interval_to_pos_in_vect[make_pair(temp_chr,temp_start)] = i; //so this map goes from chr & start pos to position in vect;
			temp_start+=interval_size;
		}
		ss.clear();
	}

	string temp_cell_name;
	map<string,int> cell_to_pos_in_vect;
	cell temp_cell_struct;
	vector<cell> vector_of_cells;
	int temp_start;
	int temp_interval;
	int temp_pos_in_vect;
	int line_counter = 0;

	while(getline(cell_chr_start_end_file, line))
	{
		if (line_counter%500000 == 0)
		{
			cout << line_counter << " : " << total_count_per_interval_map.size() << endl;
		}
		line_counter++;

		ss << line;
		counter = 0;
		while(getline(ss, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_cell_name = elem;
				//cout << temp_cell_name << endl;
				counter++;
			}
			else if (counter == 1)
			{
				temp_chr = elem;
				counter++;
			}
			else if (counter == 2)
			{
				temp_start = stoi(elem);
				//cout << temp_cell_name << ":" << temp_start << "div by 5000  " << endl;
				counter++;
			}
			else
			{
				counter++;
			}
		}
		if (cell_to_pos_in_vect.find(temp_cell_name) == cell_to_pos_in_vect.end())
		{
			temp_cell_struct.cell_name = temp_cell_name;
			temp_interval = floor(double(temp_start)/double(interval_size))*interval_size + 1;
			//cout << temp_interval << endl;
			temp_cell_struct.total_count = 1;
			temp_cell_struct.interval_in_map_num_to_count[make_pair(temp_chr,temp_interval)] = 1;
			vector_of_cells.push_back(temp_cell_struct);
			cell_to_pos_in_vect[temp_cell_name] = vector_of_cells.size() - 1;
			temp_cell_struct = cell();
			if (total_count_per_interval_map.find(make_pair(temp_chr,temp_interval)) == total_count_per_interval_map.end())
			{
				total_count_per_interval_map[make_pair(temp_chr,temp_interval)] = 1;
				//cout << temp_chr << ":" << temp_interval << endl;
			}
			else
			{
				total_count_per_interval_map[make_pair(temp_chr,temp_interval)]++;
			}
		}
		else //if cell in map
		{
			temp_pos_in_vect = cell_to_pos_in_vect[temp_cell_name];
			temp_interval = floor(double(temp_start)/double(interval_size))*interval_size + 1;
			//cout << temp_interval << endl;
			if (vector_of_cells[temp_pos_in_vect].interval_in_map_num_to_count.find(make_pair(temp_chr,temp_interval)) != vector_of_cells[temp_pos_in_vect].interval_in_map_num_to_count.end()) //if interval already in map
			{
				vector_of_cells[temp_pos_in_vect].interval_in_map_num_to_count[make_pair(temp_chr,temp_interval)]++;
				vector_of_cells[temp_pos_in_vect].total_count++;
				total_count_per_interval_map[make_pair(temp_chr,temp_interval)]++;
				//cout << temp_interval << " : " << vector_of_cells[temp_pos_in_vect].interval_in_map_num_to_count[make_pair(temp_chr,temp_interval)] << endl;
			}
			else
			{
				vector_of_cells[temp_pos_in_vect].interval_in_map_num_to_count[make_pair(temp_chr,temp_interval)] = 1;
				vector_of_cells[temp_pos_in_vect].total_count++;
				if (total_count_per_interval_map.find(make_pair(temp_chr,temp_interval)) == total_count_per_interval_map.end())
				{
					total_count_per_interval_map[make_pair(temp_chr,temp_interval)] = 1;
				}
				else
				{
					total_count_per_interval_map[make_pair(temp_chr,temp_interval)]++;
				}
			}
		}
		ss.clear();
	}

	//cout total number of intervals w/ count > some cutoff & in above map, total counts in those intervals, &
	cout << "Total num intervals w/ count > 1= " << total_count_per_interval_map.size() << endl;
	int total_number_over_10 = 0;
	int total_number_over_100 = 0;

	map<int, pair<string,int> > new_map_of_numbers_to_real_intervals;
	map<pair<string,int>, int> new_map_of_intervals_to_numbers;
	int real_interval_counter = 1;
	int total_count_for_real_intervals = 0;

	for(auto it = total_count_per_interval_map.begin(); it != total_count_per_interval_map.end(); it++)
	{
		if (interval_to_pos_in_vect.find(it -> first) != interval_to_pos_in_vect.end() && it -> second > min_count_per_interval)
		{
			new_map_of_numbers_to_real_intervals[real_interval_counter] = it -> first;
			new_map_of_intervals_to_numbers[it -> first] = real_interval_counter;
			intervals_outfile << it -> first.first << "." << it -> first.second << endl;
			real_interval_counter++;
			total_count_for_real_intervals+= it -> second;
		}
	}

	int total_good_cell_counter = 1;
	map<int, string> cell_num_to_real_cell;

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		if (vector_of_cells[i].total_count > real_cell_cutoff)
		{
			cell_num_to_real_cell[total_good_cell_counter] = vector_of_cells[i].cell_name;
			vector_of_cells[i].new_cell_num = total_good_cell_counter;
			barcodes_outfile << vector_of_cells[i].cell_name << endl;
			total_good_cell_counter++;
		}
	}

	//make matrix outfile
	matrix_outfile << "%%MatrixMarket matrix coordinate real general" << endl;
	matrix_outfile << "%" << endl;
	matrix_outfile << real_interval_counter - 1 << " " << total_good_cell_counter - 1 << " " << total_count_for_real_intervals << endl;

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		if (vector_of_cells[i].total_count > real_cell_cutoff)
		{
			for (auto it = vector_of_cells[i].interval_in_map_num_to_count.begin(); it != vector_of_cells[i].interval_in_map_num_to_count.end(); it++)
			{
				if(new_map_of_intervals_to_numbers.find(it -> first) != new_map_of_intervals_to_numbers.end())
				{
					matrix_outfile << new_map_of_intervals_to_numbers[it -> first] << " " << vector_of_cells[i].new_cell_num << " " << it -> second << endl;
					cell_interval_count_outfile << vector_of_cells[i].cell_name << "\t" << it -> first.first << "\t" << it -> first.second << "\t" << it -> second << endl;
				}
			}
		}
	}

	//cout << "Total num intervals w/ count > 10 = " << total_number_over_10 << endl;
	//cout << "Total num intervals w/ count > 100 = " << total_number_over_100 << endl;

	return 0;
}



