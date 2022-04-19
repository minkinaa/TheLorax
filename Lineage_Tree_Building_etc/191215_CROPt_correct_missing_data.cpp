/*
 * 191215_CROPt_correct_missing_data.cpp
 *
 *  Created on: Dec 15, 2019
 *      Author: anna
 *
 *  argv[1]: cigar_file
 *  argv[2] weight_same_edit="5"
	argv[3]weight_both_70M="1"
	argv[4]weight_one_X="0"
	argv[5]weight_at_least_one_AMB="0"
	argv[6]weight_diff_edits="-5"
	argv[7]weight_one_70_one_edited="-2"

	New:
	argv[8]weight if both are missing
	argv[9]weight if both are ambiguous
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
	vector<string> corrected_editing_patt_vect;
	vector<int> vector_of_distances;
	vector<string> names_of_nearest_neighbors;
	//vector<string> lin_groups_vector; //these are numbers but can probably be treated as strings
	int max_distance;
};

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

int target_pos;
map<string, int> map_of_edit_patt_to_counts;
int max_count;
string temp_edit;
int temp_pos_of_nn;
int counter_for_not_fixed = 0;
int counter_for_corrected = 0;
string correct_edit(int target_num, int cell_num, vector<cell> &vect_of_cells, map<string,int>& map_to_pos_in_vect)
{
	bool temp_print = false;
	map_of_edit_patt_to_counts.clear();
	max_count = 0;
	for (int i = 0; i < vect_of_cells[cell_num].names_of_nearest_neighbors.size(); i++)
	{
		temp_pos_of_nn = map_to_pos_in_vect[vect_of_cells[cell_num].names_of_nearest_neighbors[i]];
		temp_edit = vect_of_cells[temp_pos_of_nn].editing_patt_vect[target_num];
		if (map_of_edit_patt_to_counts.find(temp_edit) == map_of_edit_patt_to_counts.end())
		{
			map_of_edit_patt_to_counts[temp_edit] = 1;
			if (max_count < map_of_edit_patt_to_counts[temp_edit])
			{
				max_count = map_of_edit_patt_to_counts[temp_edit];
			}
		}
		else
		{
			map_of_edit_patt_to_counts[temp_edit]++;
			if (max_count < map_of_edit_patt_to_counts[temp_edit])
			{
				max_count = map_of_edit_patt_to_counts[temp_edit];
			}
		}
	}

	//just for printing stuff for now:
	int temp_count = 0;
	for (auto it = map_of_edit_patt_to_counts.begin(); it != map_of_edit_patt_to_counts.end(); it++)
	{
		if (temp_print == true)
		{
			cout << it -> first << " : " << it -> second << endl;
		}
		if (it -> second == max_count)
		{
			temp_count++;
		}
	}
	/*
	if (temp_count > 1)
	{
		cout << vect_of_cells[cell_num].cell_id << " : target:" << target_num << endl;
		for (auto it = map_of_edit_patt_to_counts.begin(); it != map_of_edit_patt_to_counts.end(); it++)
		{
			if (it -> second == max_count)
			{
				cout << it -> first << " : " << it -> second << endl;
			}
		}
	}
	*/
	//end of printing stuff

	for (auto it = map_of_edit_patt_to_counts.begin(); it != map_of_edit_patt_to_counts.end(); it++)
	{
		if (it -> second == max_count && temp_count == 1)
		{
			counter_for_corrected++;
			return it -> first;
		}
		else if (temp_count > 1)
		{
			counter_for_not_fixed++;
			return vect_of_cells[cell_num].editing_patt_vect[target_num];
		}
	}
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

	vector<cell>vector_of_cells;
	map<string, int> pos_in_vect_of_cells_map;
	cell temp_cell;

	map<string,string> all_editing_patterns_map;

	ifstream(cigar_file);
	cigar_file.open(argv[1]);

	string line;
	stringstream temp;
	string elem;
	int counter = 0;

	getline(cigar_file, line);
	string targets_in_order = line;

	//Make
	while(getline(cigar_file, line))
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
		temp_cell.max_distance = 0;
		vector_of_cells.push_back(temp_cell);
		pos_in_vect_of_cells_map[temp_cell.cell_id] = vector_of_cells.size() - 1;
		temp.clear();
		temp_cell = cell();
	}
	cout << "number_or_cells_to_map_to = " << vector_of_cells.size() << endl;

	int distance;
	int self_dist;
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		for (int j = 0; j < vector_of_cells.size(); j++)
		{
			distance = calculate_distance(unedited, vector_of_cells[i].editing_patt_vect, vector_of_cells[j].editing_patt_vect, weight_same_edit, weight_both_70M, weight_one_X, weight_at_least_one_AMB, weight_diff_edits, weight_one_70_one_edited, weight_both_X, weight_both_AMB);
			vector_of_cells[i].vector_of_distances.push_back(distance);
			if (distance > vector_of_cells[i].max_distance && i != j)
			{
				vector_of_cells[i].max_distance = distance;
			}
		}
	}



	cout << "Determining which cells are closest by distance" << endl;
	//Make list of cells which are temp_max_distance away
	int temp_max_distance;
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		temp_max_distance = vector_of_cells[i].max_distance;
		for (int j = 0; j < vector_of_cells.size(); j++)
		{

			if (vector_of_cells[i].vector_of_distances[j] == temp_max_distance && i != j)
			{
				vector_of_cells[i].names_of_nearest_neighbors.push_back(vector_of_cells[j].cell_id);
			}
		}
		/*
		cout << vector_of_cells[i].names_of_nearest_neighbors.size() << endl;
		if (vector_of_cells[i].names_of_nearest_neighbors.size() < 4)
		{
			cout << vector_of_cells[i].cell_id << " : " << vector_of_cells[i].names_of_nearest_neighbors[0] << " : " << vector_of_cells[i].max_distance << endl;
		}
		*/
	}

	//ok, now we have a list of NN cells;
	string temp_corrected_edit;
	int total_Xs = 0;
	//int counter_for_corrected;
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		for (int j = 0; j < vector_of_cells[i].editing_patt_vect.size(); j++)
		{
			if (vector_of_cells[i].editing_patt_vect[j] != "X" && vector_of_cells[i].editing_patt_vect[j] != "AMB")
			{
				vector_of_cells[i].corrected_editing_patt_vect.push_back(vector_of_cells[i].editing_patt_vect[j]);
			}
			else if (vector_of_cells[i].editing_patt_vect[j] == "X")
			{
				total_Xs++;
				temp_corrected_edit = correct_edit(j, i, vector_of_cells, pos_in_vect_of_cells_map);
				if (vector_of_cells[i].cell_id == "Group_85")
				{
					cout << "wft....." << j << " : " << temp_corrected_edit << endl;
				}
				vector_of_cells[i].corrected_editing_patt_vect.push_back(temp_corrected_edit);
			}
			else if (vector_of_cells[i].editing_patt_vect[j] == "AMB")
			{
				temp_corrected_edit = correct_edit(j, i, vector_of_cells, pos_in_vect_of_cells_map);
				vector_of_cells[i].corrected_editing_patt_vect.push_back(temp_corrected_edit);
			}
		}
		if (vector_of_cells[i].corrected_editing_patt_vect.size() != vector_of_cells[i].editing_patt_vect.size())
		{
			cout << "Yo, try again!!" << endl;
		}
	}

	cout << "Total Xs: " << total_Xs << endl;
	cout << "unfixed: " << counter_for_not_fixed << endl;
	cout << "fixed: " << counter_for_corrected << endl;

	ofstream(new_cigar_file);
	new_cigar_file.open(argv[10]);

	new_cigar_file << targets_in_order << endl;

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		new_cigar_file << vector_of_cells[i].cell_id << "\t";
		for (int j = 0; j < vector_of_cells[i].corrected_editing_patt_vect.size(); j++)
		{
			if (j < vector_of_cells[i].corrected_editing_patt_vect.size() - 1)
			{
				new_cigar_file << vector_of_cells[i].corrected_editing_patt_vect[j] << "\t";
			}
			else
			{
				new_cigar_file << vector_of_cells[i].corrected_editing_patt_vect[j];
			}
		}
		new_cigar_file << endl;
	}








	return 0;
}




