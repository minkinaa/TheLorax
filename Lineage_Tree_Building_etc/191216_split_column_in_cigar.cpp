/*
 * 191216_split_column_in_cigar.cpp
 *
 *  Created on: Dec 16, 2019
 *      Author: anna
 *
 *
 *      argv[1]: targ_1 editing patterns split by ,
 *      argv[2]: targ_2 editing patterns split by ,
 *      argv[3]: cigar file
 *      argv[4]: modified main output file (inlcuding just target of interest): cell \t edit \t umi_count ****also subset to include real cells!! [easy in R]
 *      argv[5]: target of interest
 *      argv[6]: new_cigar_file
 *
 *
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
#include <algorithm>

using namespace std;

struct cell{
	string cell_id;
	vector<string> editing_patt_vect;
	map<string,int> edit_to_count_map;
};

vector<string> temp_vector;
vector<string> parse_edits_input(string edits_w_commas)
{
	temp_vector.clear();
	int num_commas = count(edits_w_commas.begin(), edits_w_commas.end(), ',');

	string temp_edit;
	string left_over_targets = edits_w_commas;
	for (int i = 0; i < num_commas + 2; i++)
	{
		temp_edit = left_over_targets.substr(0,left_over_targets.find(","));
		temp_vector.push_back(temp_edit);

		if (i < num_commas + 1)
		{
			left_over_targets = left_over_targets.substr(left_over_targets.find(",") + 1);
		}
	}
	return temp_vector;
}

vector<string> temp_two_editing_patterns;
map<int,int> umi_count_to_number_of_occurrences_map;
vector<string> find_two_best_editing_patterns(map<string,int>& map_of_editing_patterns)
{
	umi_count_to_number_of_occurrences_map.clear();
	temp_two_editing_patterns.clear();
	for (auto it = map_of_editing_patterns.begin(); it != map_of_editing_patterns.end(); it++) //this makes a map: UMI_count -> # of occurrences
	{
		if (umi_count_to_number_of_occurrences_map.find(it -> second) == umi_count_to_number_of_occurrences_map.end())
		{
			umi_count_to_number_of_occurrences_map[it -> second] = 1;
		}
		else
		{
			umi_count_to_number_of_occurrences_map[it -> second]++;
		}
	}
	//print for test
	/*
	for (auto it = umi_count_to_number_of_occurrences_map.rbegin(); it != umi_count_to_number_of_occurrences_map.rend(); it++)
	{
		cout << it -> first << "->" << it->second <<",";
	}
	cout << endl;
	*/
	auto it2 = umi_count_to_number_of_occurrences_map.rbegin(); //start from end!!
	if (it2 -> second == 2) //if top count occurs twice, the take top two patterns
	{
		for (auto it = map_of_editing_patterns.begin(); it != map_of_editing_patterns.end(); it++)
		{
			if (it -> second == it2 -> first)
			{
				temp_two_editing_patterns.push_back(it -> first);
			}
		}
		return temp_two_editing_patterns;
	}
	else if (it2 -> second > 2)
	{
		temp_two_editing_patterns.push_back("AMB");
		temp_two_editing_patterns.push_back("AMB");
		return temp_two_editing_patterns;
	}
	else if (it2 -> second == 1) //top UMI occurs only once
	{
		for (auto it = map_of_editing_patterns.begin(); it != map_of_editing_patterns.end(); it++)
		{
			if (it -> second == it2 -> first) //if umi count matches top umi count
			{
				temp_two_editing_patterns.push_back(it -> first);
				//it2++;  //COME BACK HERE -- this is hella confusing!!!
			}
		} //this will take care of first;
		it2++;
		if (it2 -> second == 1)
		{
			for (auto it = map_of_editing_patterns.begin(); it != map_of_editing_patterns.end(); it++)
			{
				if (it -> second == it2 -> first) //if umi count matches top umi count
				{
					temp_two_editing_patterns.push_back(it -> first);
					return temp_two_editing_patterns;
					//it2++;  //COME BACK HERE -- this is hella confusing!!!
				}
			}
		}
		else
		{
			temp_two_editing_patterns.push_back("AMB");
			return temp_two_editing_patterns;
		}
	}
	//return temp_two_editing_patterns;
}



int main(int argc, char* argv[])
{
	string target_of_interest = argv[5];

	map<string,string> edits_in_position_1;
	map<string,string> edits_in_position_2;

	string edits_in_1 = argv[1];
	string edits_in_2 = argv[2];

	vector<string> vector_of_edits_1 = parse_edits_input(edits_in_1);
	vector<string> vector_of_edits_2 = parse_edits_input(edits_in_2);

	for (int i = 0; i < vector_of_edits_1.size(); i++)
	{
		edits_in_position_1[vector_of_edits_1[i]] = vector_of_edits_1[i];
	}

	for (int i = 0; i < vector_of_edits_2.size(); i++)
	{
		edits_in_position_2[vector_of_edits_2[i]] = vector_of_edits_2[i];
	}

	ifstream(cigar_file);
	cigar_file.open(argv[3]);

	string line;
	stringstream temp;
	string elem;
	int counter = 0;

	vector <string> vector_of_targets_in_order;
	map<string,int> map_pos_of_all_targets;
	int pos_in_vect_of_AMB_target;

	getline(cigar_file, line);
	string targets_in_order = line;
	temp << line;
	counter = 0;
	while(getline(temp, elem, '\t'))
	{
		vector_of_targets_in_order.push_back(elem);
		map_pos_of_all_targets[elem] = vector_of_targets_in_order.size() - 1;
		if (elem == target_of_interest)
		{
			pos_in_vect_of_AMB_target = map_pos_of_all_targets[elem];
		}
	}
	temp.clear();

	vector<cell>vector_of_cells;
	map<string, int> pos_in_vect_of_cells_map;
	cell temp_cell;
	map<string,string> all_editing_patterns_map;

	while(getline(cigar_file, line))
	{
		//cout << line << endl;
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
		vector_of_cells.push_back(temp_cell);
		pos_in_vect_of_cells_map[temp_cell.cell_id] = vector_of_cells.size() - 1;
		temp.clear();
		temp_cell = cell();
	}

	//run through

	string temp_edit;
	int temp_pos_of_target_in_vect;

	ifstream(mod_main_out);
	mod_main_out.open(argv[4]);
	int temp_pos_in_vect;
	string temp_cell_name;
	string temp_editing_pattern;
	int temp_umi_count;
	counter = 0;

	//cout << "making it here??" << endl;

	while(getline(mod_main_out, line))
	{
		//cout << line << endl;
		temp << line;
		counter = 0;
		while(getline(temp,elem, '\t'))
		{
			if (counter == 0)
			{
				temp_cell_name = elem;
				//cout << temp_cell_name << endl;
				temp_pos_in_vect = pos_in_vect_of_cells_map[temp_cell_name];
				//cout << temp_pos_in_vect << endl;
				counter++;
			}
			else if (counter == 1)
			{
				temp_edit = elem;
				counter++;
			}
			else if (counter == 2)
			{
				temp_umi_count = stoi(elem);
				//cout << temp_umi_count << endl;
				if(vector_of_cells[temp_pos_in_vect].edit_to_count_map.find(temp_edit) == vector_of_cells[temp_pos_in_vect].edit_to_count_map.end())
				{
					vector_of_cells[temp_pos_in_vect].edit_to_count_map[temp_edit] = temp_umi_count;
				}
				else
				{
					vector_of_cells[temp_pos_in_vect].edit_to_count_map[temp_edit]+=temp_umi_count;
				}
				counter++;
			}
			//cout << vector_of_cells[temp_pos_in_vect].edit_to_count_map.size() << endl;
		}
		temp.clear();
	}

	//deal w/ maps which are larger than 2
	vector<string> temp_two_editing_patterns;
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		//cout << "ever make it here??" << endl;
		if (vector_of_cells[i].editing_patt_vect[pos_in_vect_of_AMB_target] != "X")
		{
			if (vector_of_cells[i].edit_to_count_map.size() > 2)
			{
				temp_two_editing_patterns = find_two_best_editing_patterns(vector_of_cells[i].edit_to_count_map);
				if (temp_two_editing_patterns.size() != 2)
				{
					//cout << vector_of_cells[i].cell_id << endl;
					//cout << temp_two_editing_patterns[0] << ", " << temp_two_editing_patterns[1] << ", " << temp_two_editing_patterns[2]<< endl;
					cout << "yeah, no way this code runs right:" << temp_two_editing_patterns.size() << endl;
				}
				vector_of_cells[i].edit_to_count_map.clear();
				vector_of_cells[i].edit_to_count_map[temp_two_editing_patterns[0]] = 0;
				vector_of_cells[i].edit_to_count_map[temp_two_editing_patterns[1]] = 0;
			}
		}
	}

	//correct cigar string;
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		//cout << vector_of_cells[i].edit_to_count_map.size() << endl;
		if (vector_of_cells[i].editing_patt_vect[pos_in_vect_of_AMB_target] != "X")
		{
			if (vector_of_cells[i].edit_to_count_map.size() == 1)
			{
				//cout << "size 1" << endl;
				auto it = vector_of_cells[i].edit_to_count_map.begin();
				if(edits_in_position_2.find(it -> first) != edits_in_position_2.end())
				{
					vector_of_cells[i].editing_patt_vect[pos_in_vect_of_AMB_target] = "X";
					vector_of_cells[i].editing_patt_vect.push_back(it -> first);
				}
				else
				{
					vector_of_cells[i].editing_patt_vect.push_back("X");
				}
			}
			else if (vector_of_cells[i].edit_to_count_map.size() == 2)
			{
				//cout << "size 2" << endl;
				auto it = vector_of_cells[i].edit_to_count_map.begin();
				if (edits_in_position_1.find(it -> first) != edits_in_position_1.end())
				{
					vector_of_cells[i].editing_patt_vect[pos_in_vect_of_AMB_target] = it -> first;
					it++;
					vector_of_cells[i].editing_patt_vect.push_back(it -> first); //put second element in the back!! (this is confusing!!)
					//edits_in_position_2[it -> first] = it -> first;
				}
				else if (edits_in_position_2.find(it -> first) != edits_in_position_2.end())
				{
					vector_of_cells[i].editing_patt_vect.push_back(it -> first);
					it++;
					vector_of_cells[i].editing_patt_vect[pos_in_vect_of_AMB_target] = it -> first;
					//edits_in_position_1[it -> first] = it -> first;
				}
				else
				{
					if (it -> first == "70M:")
					{
						vector_of_cells[i].editing_patt_vect.push_back(it -> first);
						it++;
						vector_of_cells[i].editing_patt_vect[pos_in_vect_of_AMB_target] = it -> first;
					}
					else
					{
						vector_of_cells[i].editing_patt_vect[pos_in_vect_of_AMB_target] = it -> first;
						it++;
						vector_of_cells[i].editing_patt_vect.push_back(it -> first); //put second element in the back!! (this is confusing!!)
					}
				}

			}
		}
		else
		{
			vector_of_cells[i].editing_patt_vect.push_back("X");
		}
	}

	ofstream(new_cigar_file);
	new_cigar_file.open(argv[6]);

	new_cigar_file << targets_in_order << "\t" << target_of_interest + "_2" << endl;

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		//cout << vector_of_cells[i].cell_id << endl;
		cout << vector_of_cells[i].editing_patt_vect.size() << endl;
		new_cigar_file << vector_of_cells[i].cell_id << "\t";
		for (int j = 0; j < vector_of_cells[i].editing_patt_vect.size(); j++)
		{
			if (j < vector_of_cells[i].editing_patt_vect.size() - 1)
			{
				new_cigar_file << vector_of_cells[i].editing_patt_vect[j] << "\t";
			}
			else
			{
				new_cigar_file << vector_of_cells[i].editing_patt_vect[j];
			}
		}
		new_cigar_file << endl;
	}

	return 0;
}
