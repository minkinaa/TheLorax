/*
 * 191017_two_group_distances_updated.cpp
 *
 *  Created on: Oct 17, 2019
 *      Author: anna
 */

/*
 *	argv[1]: cigar_file
 *  argv[2] weight_same_edit="5"
	argv[3]weight_both_70M="1"
	argv[4]weight_one_X="0"
	argv[5]weight_at_least_one_AMB="0"
	argv[6]weight_diff_edits="-5"
	argv[7]weight_one_70_one_edited="-2"

	New:
	argv[8]weight if both are missing
	argv[9]weight if both are ambiguous

	argv10 --> lineage group file, old cells;
	raw_file.open(argv[11])

	argv[12]: cigar file from cells to map.
	#num_nearest_neighbors = stoi(argv[13]);
	argv[13]: all_cell_lineage_group_outfile;
	argv[14]: unique_editing_patterns_outfile;



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
	string cell_id;
	vector<string> editing_patt_vect;
	vector<int> vector_of_distances;
	vector<string> names_of_nearest_neighbors;
	vector<string> lin_groups_vector; //these are numbers but can probably be treated as strings
	int max_distance;
};

void print_cell(cell temp_cell)
{
	cout << temp_cell.cell_id << endl;
	cout << temp_cell.editing_patt_vect.size() << endl;
	for (int i = 0; i < temp_cell.editing_patt_vect.size(); i++)
	{
		cout << temp_cell.editing_patt_vect[i] << "\t";
	}
}

void print_distances(cell temp_cell)
{
	cout << temp_cell.cell_id << endl;
	cout << temp_cell.vector_of_distances.size() << endl;
	for (int i = 0; i < temp_cell.vector_of_distances.size(); i++)
	{
		cout << temp_cell.vector_of_distances[i] << "\t";
	}
}

int calculate_distance(string unedited_pattern, vector<string> cell_1_vect, vector<string> cell_2_vect, int weight_same_edit, int weight_both_70M, int weight_one_X, int weight_at_least_one_AMB, int weight_diff_edits, int weight_one_70_one_edited, int weight_both_X, int weight_both_AMB){
	int dist = 0;
	string edit_1;
	string edit_2;
	for (int i = 0; i < cell_1_vect.size(); i++)
	{
		edit_1 = cell_1_vect[i];
		edit_2 = cell_2_vect[i];
		if (edit_1 == unedited_pattern && edit_2 == unedited_pattern)
		{
			dist = dist + weight_both_70M;
		}
		else if ((edit_1 == unedited_pattern || edit_2 == unedited_pattern) && edit_1 != "X" && edit_2 != "X" && edit_1 != "AMB" && edit_2 != "AMB")
		{
			dist = dist + weight_one_70_one_edited;
		}
		else if (edit_1 == "X" && edit_2 == "X")
		{
			dist = dist + weight_both_X;
		}
		else if (edit_1 == "X" || edit_2 == "X")
		{
			dist = dist + weight_one_X;
		}
		else if (edit_1 == "AMB" && edit_2 == "AMB")
		{
			dist = dist + weight_both_AMB;
		}
		else if (edit_1 == "AMB" || edit_2 == "AMB")
		{
			dist = dist + weight_at_least_one_AMB;
		}
		else if (edit_1 == edit_2)
		{
			dist = dist + weight_same_edit;
		}
		else if (edit_1 != edit_2)
		{
			dist = dist + weight_diff_edits;
		}
		else
		{
			cout << "Missed a case!!" << endl << "Edit_1 =" << edit_1 << endl << "Edit_2 =" << edit_2 << endl;
		}
	}

	return dist;
}

int temp_max;
string temp_max_LG;
string find_LG_w_max_value_in_map(map<string, int> LG_map)
{
	temp_max = 0;
	for(auto it = LG_map.begin(); it != LG_map.end(); it++)
	{
		//cout << "inside map function " << it -> second << endl;
		if (it -> second > temp_max)
		{
			temp_max_LG = it -> first;
			temp_max = it -> second;
		}
	}
	return (temp_max_LG);
}



int main(int argc, char* argv[])
{
	string unedited = "70M:";
	int weight_same_edit = stoi(argv[2]);
	int weight_both_70M = stoi(argv[3]);
	int weight_one_X = stoi(argv[4]);
	int weight_at_least_one_AMB = stoi(argv[5]);
	int weight_diff_edits = stoi(argv[6]);
	int weight_one_70_one_edited = stoi(argv[7]);
	int weight_both_X = stoi(argv[8]);
	int weight_both_AMB = stoi(argv[9]);

	//int num_nearest_neighbors = stoi(argv[13]);

	vector<cell> vector_of_cells_to_map_to;
	map<string,int> pos_in_old_cell_vect_map;
	cell temp_cell;

	vector<cell> vector_of_new_cells;

	map<string,string> all_editing_patterns_map;

	ifstream(cigar_file_to_map_to);
	cigar_file_to_map_to.open(argv[1]);

	ifstream(cigar_file_of_new_cells);
	cigar_file_of_new_cells.open(argv[12]);

	string line;
	stringstream temp;
	string elem;
	int counter = 0;

	getline(cigar_file_to_map_to, line);
	while(getline(cigar_file_to_map_to, line))
	{
		temp << line;
		counter = 0;
		while(getline(temp, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_cell.cell_id = elem;
				counter ++;
			}
			else
			{
				temp_cell.editing_patt_vect.push_back(elem);
				if (all_editing_patterns_map.find(elem) == all_editing_patterns_map.end())
				{
					all_editing_patterns_map[elem] = elem;
				}
			}
		}
		vector_of_cells_to_map_to.push_back(temp_cell);
		pos_in_old_cell_vect_map[temp_cell.cell_id] = vector_of_cells_to_map_to.size() - 1;
		temp.clear();
		temp_cell = cell();
	}
	cout << "number_or_cells_to_map_to = " << vector_of_cells_to_map_to.size() << endl;

	//make vector of new cells
	getline(cigar_file_of_new_cells, line);
	while(getline(cigar_file_of_new_cells, line))
	{
		temp << line;
		counter = 0;
		while(getline(temp, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_cell.cell_id = elem;
				counter ++;
			}
			else
			{
				temp_cell.editing_patt_vect.push_back(elem);
				if (all_editing_patterns_map.find(elem) == all_editing_patterns_map.end())
				{
					all_editing_patterns_map[elem] = elem;
				}
			}
		}
		vector_of_new_cells.push_back(temp_cell);
		temp.clear();
		temp_cell = cell();
	}
	cout << "number_or_cells_to_map = " << vector_of_new_cells.size() << endl;


	//CALCULATE DISTACE MATRIX
	cout << "Calculating Distances ..." << endl;
	int distance;
	int temp_max_distance;
	for (int i = 0; i <  vector_of_new_cells.size(); i++)
	{
		if (i%1000 == 0)
		{
			cout << i << endl;
		}
		temp_max_distance = 0;
		for (int j = 0; j < vector_of_cells_to_map_to.size(); j++)
		{
			distance = calculate_distance(unedited, vector_of_new_cells[i].editing_patt_vect, vector_of_cells_to_map_to[j].editing_patt_vect, weight_same_edit, weight_both_70M, weight_one_X, weight_at_least_one_AMB, weight_diff_edits, weight_one_70_one_edited, weight_both_X, weight_both_AMB);
			vector_of_new_cells[i].vector_of_distances.push_back(distance);
			if (distance > temp_max_distance)
			{
				temp_max_distance = distance;
			}
		}
		vector_of_new_cells[i].max_distance = temp_max_distance;
	}

	cout << "Determining which cells are closest by distance" << endl;
	//Make list of cells which are temp_max_distance away
	for (int i = 0; i < vector_of_new_cells.size(); i++)
	{
		temp_max_distance = vector_of_new_cells[i].max_distance;
		for (int j = 0; j < vector_of_cells_to_map_to.size(); j++)
		{
			//cout << "i=" << i << " j=" << j << " temp_max_distance=" <<  temp_max_distance << endl;
			//cout << vector_of_new_cells[i].vector_of_distances.size() << endl;
			if (vector_of_new_cells[i].vector_of_distances[j] == temp_max_distance)
			{
				vector_of_new_cells[i].names_of_nearest_neighbors.push_back(vector_of_cells_to_map_to[j].cell_id);
			}
		}
		//cout << "Num nearest neighbors: " << vector_of_new_cells[i].names_of_nearest_neighbors.size() << endl;
	}

	//Make vector for old cells with lineage group info
	ifstream(old_cell_lineage_groups);
	old_cell_lineage_groups.open(argv[10]);

	getline(old_cell_lineage_groups, line);
	string lineage_group_file_header = line;

	int temp_pos_in_vect;
	while(getline(old_cell_lineage_groups, line))
	{
		counter = 0;
		temp << line;
		while(getline(temp, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_pos_in_vect = pos_in_old_cell_vect_map[elem];
				counter++;
			}
			else
			{
				vector_of_cells_to_map_to[temp_pos_in_vect].lin_groups_vector.push_back(elem);
			}
		}
		temp.clear();
	}
	cout << "LG vect size=" << vector_of_cells_to_map_to[1].lin_groups_vector.size() << endl;

	map<string,int> temp_lineage_group_map; //lineage group to count

	string temp_lineage_group;
	string max_lineage_group;
	//Determine which lineage groups new cells belong to:
	cout << "Determining which lineage groups new cells belong to" << endl;
	for (int i = 0; i < vector_of_new_cells.size(); i++)
	{
		if (i%500 == 0)
		{
			cout << i << endl;
		}
		for (int j = 0; j < vector_of_cells_to_map_to[1].lin_groups_vector.size(); j++) //iterate over 30 lineage groups
		{
			//cout << "j= " << j << endl;
			temp_lineage_group_map.clear();
			for (int k = 0; k < vector_of_new_cells[i].names_of_nearest_neighbors.size(); k++)
			{
				temp_pos_in_vect = pos_in_old_cell_vect_map[vector_of_new_cells[i].names_of_nearest_neighbors[k]];
				temp_lineage_group = vector_of_cells_to_map_to[temp_pos_in_vect].lin_groups_vector[j];
				if (temp_lineage_group_map.find(temp_lineage_group) == temp_lineage_group_map.end())
				{
					temp_lineage_group_map[temp_lineage_group] = 1;
				}
				else
				{
					temp_lineage_group_map[temp_lineage_group]++;
				}
			}
			max_lineage_group = find_LG_w_max_value_in_map(temp_lineage_group_map);
			//cout << max_lineage_group << endl;
			vector_of_new_cells[i].lin_groups_vector.push_back(max_lineage_group);
		}
	}


	ofstream(raw_file);
	raw_file.open(argv[11]);

	for (int i = 0; i < vector_of_new_cells.size(); i++)
	{
		raw_file << vector_of_new_cells[i].cell_id;
		for (int j = 0; j < vector_of_new_cells[i].vector_of_distances.size(); j++)
		{
			raw_file << "\t" << vector_of_new_cells[i].vector_of_distances[j];
		}
		raw_file << endl;
	}
	raw_file.close();

	ofstream(all_cell_lineage_groups_outfile);
	all_cell_lineage_groups_outfile.open(argv[13]);

	all_cell_lineage_groups_outfile << lineage_group_file_header << endl;
	//Add old cells
	for (int i = 0; i < vector_of_cells_to_map_to.size(); i++)
	{
		all_cell_lineage_groups_outfile << vector_of_cells_to_map_to[i].cell_id << "\t";
		for (int j = 0; j < vector_of_cells_to_map_to[i].lin_groups_vector.size(); j++)
		{
			all_cell_lineage_groups_outfile << vector_of_cells_to_map_to[i].lin_groups_vector[j] << "\t";
		}
		all_cell_lineage_groups_outfile << endl;
	}

	//Add new cells
	for (int i = 0; i < vector_of_new_cells.size(); i++)
	{
		all_cell_lineage_groups_outfile << vector_of_new_cells[i].cell_id << "\t";
		for (int j = 0; j < vector_of_new_cells[i].lin_groups_vector.size(); j++)
		{
			all_cell_lineage_groups_outfile << vector_of_new_cells[i].lin_groups_vector[j] << "\t";
		}
		all_cell_lineage_groups_outfile << endl;
	}

	//Print all unique editing patters
	ofstream(all_editing_patterns);
	all_editing_patterns.open(argv[14]);

	for (auto it = all_editing_patterns_map.begin(); it != all_editing_patterns_map.end(); it++)
	{
		all_editing_patterns << it -> first << endl;
	}

	return 0;
}









/*
 * 190927_calc_distances_between_two_groups.cpp
 *
 *  Created on: Sep 27, 2019
 *      Author: anna
 */







