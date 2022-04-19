/*
 * 201109_building_a_tree_3_record_all_changes.cpp
 *
 *  Created on: Nov 9, 2020
 *      Author: anna

 *      pseudocell_cigar_file.open(argv[1]);
 *      edit_num_string_outfile.open(argv[2]);
		target_edits_over_time_outfile.open(argv[3]);
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

struct pseudocell{
	string cell_name;
	string dynamic_concat_edit_string;
	vector<string> edit_string;
	vector<int> tree_dynamic_changes_num_vect; // output1
	vector<string> edits_associated_w_each_change; //output2
	vector<string> targets_associated_w_each_edit; //output2
	vector<int>targets_POS_associated_w_each_edit; //output2
};

struct target_and_edit{
	int target_pos;
	string edit;

	//adding this in order to keep track of any other edits which are as abundant (in either group) as the recorded one.
	//vector<int> any_other_shared_edit_target_pos;
	//vector<string> any_other_shared_edit_target_associated_edits;

	//vector<int> second_group_any_shared_edit_pos;
	//vector<string> second_group_any_shared_edit_pattern;
};

string temp_edit;
map<string,int> temp_map_of_edits_to_counts_per_target;
int temp_highest_count;
int global_highest_count;
string temp_highest_count_edit;
vector<int> highest_counts_per_target_vect;
vector<string> high_count_EDIT_per_target_vect;
target_and_edit find_most_abundant_edit(vector<pseudocell> &vect_of_cells){

	int num_cells = vect_of_cells.size();
	target_and_edit temp_target_edit_pair = target_and_edit();
	highest_counts_per_target_vect.clear();
	high_count_EDIT_per_target_vect.clear();
	global_highest_count = 0;

	map<string,string> already_recorded_edits; //here we'll make a list of edits which have already been recorded for this group; can do this from any cell in group
	string temp_POS_EDIT_string;
	for (int i = 0; i < vect_of_cells[0].edits_associated_w_each_change.size(); i++){
		temp_POS_EDIT_string = vect_of_cells[0].edits_associated_w_each_change[i] + "_" + to_string(vect_of_cells[0].targets_POS_associated_w_each_edit[i]);
		already_recorded_edits[temp_POS_EDIT_string] = temp_POS_EDIT_string;
	}

	for (int i = 0; i < vect_of_cells[0].edit_string.size(); i++) //for each target
	{
		temp_map_of_edits_to_counts_per_target.clear();
		for (int j = 0; j < vect_of_cells.size(); j++) // for each cell
		{
			temp_edit = vect_of_cells[j].edit_string[i];
			if (temp_map_of_edits_to_counts_per_target.find(temp_edit) != temp_map_of_edits_to_counts_per_target.end()) //if in map
			{
				temp_map_of_edits_to_counts_per_target[temp_edit]++;
			}
			else
			{
				temp_map_of_edits_to_counts_per_target[temp_edit] = 1;
			}
		}
		//Now temp_map_of_edits_to_counts_per_target contains all edits for that target and their counts
		//Let's choose the highest one and save it
		temp_highest_count = 0;
		temp_highest_count_edit = "NONE";
		for (auto it = temp_map_of_edits_to_counts_per_target.begin(); it != temp_map_of_edits_to_counts_per_target.end(); it++)
		{
			temp_POS_EDIT_string = string(it -> first) + "_" + to_string(i);
			if (it -> first != "70M:" && it -> first != "X" && it -> first != "AMB")
			{
				//if (it -> second > temp_highest_count && it -> second < num_cells)
				if (it -> second > temp_highest_count && already_recorded_edits.find(temp_POS_EDIT_string) == already_recorded_edits.end())
				{
					temp_highest_count = it -> second;
					temp_highest_count_edit = it -> first;
				}
			}
		}

		//now have a highest count and edit for that target saved; let's store those...just in vectors, why not
		highest_counts_per_target_vect.push_back(temp_highest_count);
		high_count_EDIT_per_target_vect.push_back(temp_highest_count_edit);
		if (temp_highest_count > global_highest_count)
		{
			global_highest_count = temp_highest_count;
		}
	}

	bool recorded_global = false;

	for(int i = 0; i < highest_counts_per_target_vect.size(); i++)
	{
		if (highest_counts_per_target_vect[i] == global_highest_count)
		{
			temp_target_edit_pair.target_pos = i;
			temp_target_edit_pair.edit = high_count_EDIT_per_target_vect[i];
			recorded_global = true;
			return (temp_target_edit_pair);
		}
		/*
		else if (highest_counts_per_target_vect[i] == global_highest_count && recorded_global == true){
			temp_target_edit_pair.any_other_shared_edit_target_pos.push_back(i);
			temp_target_edit_pair.any_other_shared_edit_target_associated_edits.push_back(high_count_EDIT_per_target_vect[i]);
		}
		*/
	}
	//return (temp_target_edit_pair);
	cout << "HEY, WE SHOULD NEVER MAKE IT THIS FAR!!" << 	endl;
}

void print_map(map<int,string>map_to_print){
	for (auto it = map_to_print.begin(); it != map_to_print.end(); it++)
	{
		cout << it -> second << ",";
	}
	cout << endl << endl;
}

int main(int argc, char* argv[]){

	ifstream(pseudocell_cigar_file);
	pseudocell_cigar_file.open(argv[1]);
	stringstream ss;
	string elem;
	int elem_counter;
	string line;

	pseudocell temp_pseudocell;
	vector<pseudocell> vect_of_pseudocells;

	//get target names
	vector<string> vect_of_target_names;
	string header;
	getline(pseudocell_cigar_file, header);
	ss << header;
	while(getline(ss,elem,'\t')){
		vect_of_target_names.push_back(elem);
	}
	ss.clear();

	//add in edit string information to each pseudocell struct
	int cell_counter = 0;
	while(getline(pseudocell_cigar_file, line)){
		elem_counter = 0;
		ss << line;
		while(getline(ss, elem, '\t'))
		{
			if (elem_counter == 0)
			{
				temp_pseudocell.cell_name = elem;
				temp_pseudocell.tree_dynamic_changes_num_vect.push_back(1); //so every cell
				temp_pseudocell.dynamic_concat_edit_string = to_string(1);
				elem_counter++;

			}
			else
			{
				temp_pseudocell.edit_string.push_back(elem);
			}
		}
		vect_of_pseudocells.push_back(temp_pseudocell);
		temp_pseudocell = pseudocell();
		cell_counter++;
		ss.clear();
	}


	// OK, THIS IS WHERE THE FUN BEGINS!!!!!  //
	// take #1: //
	map<string, map<int, string> > concat_string_to_map_of_pos_and_cell_names; //dynamic map of edit strings to list of positions
	//may want to change this to string -> int for second map?? We'll see...

	//fill map w/ initial edit strings
	string temp_cell_name;

	string temp_concat_string;
	for (int i = 0; i < vect_of_pseudocells.size(); i++)
	{
		temp_cell_name = vect_of_pseudocells[i].cell_name;
		temp_concat_string = vect_of_pseudocells[i].dynamic_concat_edit_string;
		concat_string_to_map_of_pos_and_cell_names[temp_concat_string][i] = temp_cell_name;
	}
	//by now you should have a map of dyn_group_strings pointing to a list of cells sharing those strings
	//currently, main map should have one just one entry.

	//do iteratively until every cell is associated w/ unique path
	string temp_dyn_num_string;
	map<int,string> temp_map_of_posit_assoc_w_dyn_num_string;
	vector<pseudocell> partial_vect_of_pseudocells;
	target_and_edit temp_best_target_and_edit;
	int new_map_size = concat_string_to_map_of_pos_and_cell_names.size();
	int old_map_size = 0;
	int counter_of_number_of_times_maps_are_equal = 0;
	//while(new_map_size > old_map_size){ //do all this until every cell is associated w/ unique path
	while(counter_of_number_of_times_maps_are_equal < vect_of_target_names.size()){
		//cout << concat_string_to_map_of_pos_and_cell_names.size() <<  " : " << vect_of_pseudocells.size() << endl;
		for (auto it = concat_string_to_map_of_pos_and_cell_names.begin(); it != concat_string_to_map_of_pos_and_cell_names.end(); it++) //for each unique dyn_num_string (or single group)
		{
			temp_best_target_and_edit = target_and_edit();
			temp_dyn_num_string = it -> first;
			temp_map_of_posit_assoc_w_dyn_num_string = it -> second;
			for (auto it2 = temp_map_of_posit_assoc_w_dyn_num_string.begin(); it2 != temp_map_of_posit_assoc_w_dyn_num_string.end(); it2++) //make new vector of associated cells
			{
				partial_vect_of_pseudocells.push_back(vect_of_pseudocells[it2 -> first]);
			}
			//Find the most abundant target
			temp_best_target_and_edit = find_most_abundant_edit(partial_vect_of_pseudocells);
			int temp_target_pos = temp_best_target_and_edit.target_pos;
			string temp_target_edit = temp_best_target_and_edit.edit;

			//vector<int> temp_additional_target_pos = temp_best_target_and_edit.any_other_shared_edit_target_pos;
			//vector<string> temp_additional_target_edit_patterns = temp_best_target_and_edit.any_other_shared_edit_target_associated_edits;

			//split on that target; use the list of CELL positions from temp_map_of_posit_assoc_w_dyn_num_string
			for (auto it2 = temp_map_of_posit_assoc_w_dyn_num_string.begin(); it2 != temp_map_of_posit_assoc_w_dyn_num_string.end(); it2++) //run through only cells associated w/ this group
			{
				if(vect_of_pseudocells[it2 -> first].edit_string[temp_target_pos] == temp_target_edit)
				{
					int len_of_numbers_vect = vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect.size();
					int next_val = vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect[len_of_numbers_vect - 1] + 1;
					vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect.push_back(next_val);
					vect_of_pseudocells[it2 -> first].dynamic_concat_edit_string = vect_of_pseudocells[it2 -> first].dynamic_concat_edit_string + "_" + to_string(next_val);
					vect_of_pseudocells[it2 -> first].edits_associated_w_each_change.push_back(temp_target_edit);
					vect_of_pseudocells[it2 -> first].targets_associated_w_each_edit.push_back(vect_of_target_names[temp_target_pos]);
					vect_of_pseudocells[it2 -> first].targets_POS_associated_w_each_edit.push_back(temp_target_pos);

				}
				else{
					int len_of_numbers_vect = vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect.size();
					int next_val = vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect[len_of_numbers_vect - 1]; //keep number the same
					vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect.push_back(next_val);
					vect_of_pseudocells[it2 -> first].dynamic_concat_edit_string = vect_of_pseudocells[it2 -> first].dynamic_concat_edit_string + "_" + to_string(next_val);
				}

				//cout << vect_of_pseudocells[it2 -> first].dynamic_concat_edit_string << " : " << vect_of_pseudocells[it2 -> first].cell_name << endl;
				//cout << vect_of_pseudocells[it2 -> first].cell_name << "\t" << vect_of_pseudocells[it2 -> first].dynamic_concat_edit_string << endl;
			}

			partial_vect_of_pseudocells.clear(); //so ready for next iteration...
		}


		//(5) you should clear and remake the map w/ the new strings
		old_map_size = concat_string_to_map_of_pos_and_cell_names.size();
		concat_string_to_map_of_pos_and_cell_names.clear();
		for (int i = 0; i < vect_of_pseudocells.size(); i++) //remake map
		{
			temp_cell_name = vect_of_pseudocells[i].cell_name;
			temp_concat_string = vect_of_pseudocells[i].dynamic_concat_edit_string;
			concat_string_to_map_of_pos_and_cell_names[temp_concat_string][i] = temp_cell_name;
		}
		new_map_size = concat_string_to_map_of_pos_and_cell_names.size();
		if (new_map_size == old_map_size)
		{
			counter_of_number_of_times_maps_are_equal++;
		}
	}
	//print stuff
	ofstream(edit_num_string_outfile);
	edit_num_string_outfile.open(argv[2]);

	int num_columns = vect_of_pseudocells[0].tree_dynamic_changes_num_vect.size();
	for (int i = 0; i < (num_columns -1); i++){
		edit_num_string_outfile << "X" << to_string(i+1) << "\t";
	}
	edit_num_string_outfile << "X" << to_string(num_columns) << endl;


	ofstream(target_edits_over_time_outfile);
	target_edits_over_time_outfile.open(argv[3]);
	target_edits_over_time_outfile << "Group\tChangeNum\tTarget\tTargetPos\tEdit" << endl;

	for (auto it = concat_string_to_map_of_pos_and_cell_names.begin(); it != concat_string_to_map_of_pos_and_cell_names.end(); it++) //for each unique dyn_num_string (or single group)
	{
		temp_map_of_posit_assoc_w_dyn_num_string = it -> second;
		for (auto it2 = temp_map_of_posit_assoc_w_dyn_num_string.begin(); it2 != temp_map_of_posit_assoc_w_dyn_num_string.end(); it2++) //run through only cells associated w/ this group
		{
			//THIS IS GOOD, Bring it back:
			edit_num_string_outfile << vect_of_pseudocells[it2 -> first].cell_name << "\t";
			for (int i = 0; i < vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect.size()-1; i++)
			{
				edit_num_string_outfile << vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect[i] << "\t";
			}
			edit_num_string_outfile << vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect[vect_of_pseudocells[it2 -> first].tree_dynamic_changes_num_vect.size()-1] << endl;


			for(int i = 0; i < vect_of_pseudocells[it2 -> first].targets_POS_associated_w_each_edit.size(); i++)
			{
				target_edits_over_time_outfile << vect_of_pseudocells[it2 -> first].cell_name << "\t";
				target_edits_over_time_outfile << i + 1 << "\t";
				target_edits_over_time_outfile << vect_of_pseudocells[it2 -> first].targets_associated_w_each_edit[i] << "\t";
				target_edits_over_time_outfile << vect_of_pseudocells[it2 -> first].targets_POS_associated_w_each_edit[i] << "\t";
				target_edits_over_time_outfile << vect_of_pseudocells[it2 -> first].edits_associated_w_each_change[i] << endl;
			}

			//cout << vect_of_pseudocells[it2 -> first].cell_name << "\t" << vect_of_pseudocells[it2 -> first].dynamic_concat_edit_string << endl;
		}
	}

	return 0;
}





