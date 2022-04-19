/*
 * 200806_make_coordinates_for_tree_plot.cpp
 *
 *  Created on: Aug 6, 2020
 *      Author: anna
 *
 *      argv[1]: table
 *      argv[2]: outfile
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
	vector<int> lin_group_vect;

};

int main(int argc, char* argv[]){

	//vector<int> x_coord;
	//vector<int> y_coord;
	//vector<int> last_x_coord;
	//vector<int> last_y_coord;

	//map<string, pair<double, double>> map_from_Col_LG_to_x_y_coordinates;
	map<string, vector<double> > map_from_Col_LG_to_x_y_coordinates;

	cell temp_cell;
	vector<cell> vector_of_cells;

	ifstream(cutdf_file);
	cutdf_file.open(argv[1]);
	string line;
	stringstream ss;
	string elem;
	int counter = 0;

	getline(cutdf_file, line);

	while(getline(cutdf_file, line)){
		ss << line;
		counter = 0;
		while(getline(ss, elem, '\t')){
			if (counter == 0){
				temp_cell.cell_name = elem;
				counter++;
			}
			else {
				temp_cell.lin_group_vect.push_back(stoi(elem));
			}
		}
		vector_of_cells.push_back(temp_cell);
		temp_cell = cell();
		ss.clear();
	}

	map_from_Col_LG_to_x_y_coordinates["0_1"].push_back(0);
	map_from_Col_LG_to_x_y_coordinates["0_1"].push_back(double(vector_of_cells.size())/2.0);
	map_from_Col_LG_to_x_y_coordinates["0_1"].push_back(0);
	map_from_Col_LG_to_x_y_coordinates["0_1"].push_back(double(vector_of_cells.size())/2.0);

	vector<int> temp_vect_of_groups_in_column;
	vector<int> temp_vect_of_groups_in_PREVIOUS_column;
	map<int, int> map_of_freq_of_each_group_in_column;
	int temp_lg_num;
	int temp_previous_lg_num;
	double number_of_groups_left;
	double temp_x;
	double temp_y;
	double temp_previous_x;
	double temp_previous_y;
	double temp_freq;
	string previous_point;
	int temp_previous_lg;

	int temp_x_pos;

	double stored_x;
	double stored_y;
	double stored_previous_x;
	double stored_previous_y;
	bool found_last_point;
	bool last_y_repeats;

	int counter_num_points_erased = 0;

	map<string, string> map_of_points_to_cell_names;

	//cout << "we made it!!" << endl;

	for (int i = 1; i < vector_of_cells[0].lin_group_vect.size(); i++){ //ignoring first column
		//cout << i << endl;
		temp_vect_of_groups_in_column.clear();
		temp_vect_of_groups_in_PREVIOUS_column.clear();
		map_of_freq_of_each_group_in_column.clear();
		number_of_groups_left = double(vector_of_cells.size());

		//make a vector of groups in that column, and a map of the number of times each appears
		for(int j = 0; j < vector_of_cells.size(); j++){
			//cout << "j=" << j << endl;
			temp_lg_num = vector_of_cells[j].lin_group_vect[i];
			if (map_of_freq_of_each_group_in_column.find(temp_lg_num) != map_of_freq_of_each_group_in_column.end()){
				//if in map;
				map_of_freq_of_each_group_in_column[temp_lg_num]++;
			}
			else {
				temp_vect_of_groups_in_column.push_back(temp_lg_num);
				map_of_freq_of_each_group_in_column[temp_lg_num] = 1;
				temp_previous_lg_num = vector_of_cells[j].lin_group_vect[i-1];
				temp_vect_of_groups_in_PREVIOUS_column.push_back(temp_previous_lg_num);
			}
		}

		//now, for each lg in each column, determine coordinates and previous coordinates.
		for (int j = 0; j < temp_vect_of_groups_in_column.size(); j++)
		{
			temp_x = double(i);
			temp_freq = double(map_of_freq_of_each_group_in_column[temp_vect_of_groups_in_column[j]]);
			temp_y = double(number_of_groups_left) - (temp_freq*.5);

			temp_previous_lg = temp_vect_of_groups_in_PREVIOUS_column[j];

			//find_previous_point_in_map
			found_last_point = false;
			temp_x_pos = i;

			while(found_last_point == false){
				previous_point = to_string(temp_x_pos-1)+"_"+ to_string(temp_previous_lg);
				cout << previous_point << endl;
				if (map_from_Col_LG_to_x_y_coordinates.find(previous_point) == map_from_Col_LG_to_x_y_coordinates.end()){
					temp_x_pos--;
				}
				else
				{
					//cout << "found it??" << endl;
					found_last_point = true;
				}
			}

			temp_previous_y = map_from_Col_LG_to_x_y_coordinates[previous_point][1];
			temp_previous_x = map_from_Col_LG_to_x_y_coordinates[previous_point][0];

			last_y_repeats = false;
			if(map_from_Col_LG_to_x_y_coordinates[previous_point][1] == map_from_Col_LG_to_x_y_coordinates[previous_point][3]){
				last_y_repeats = true;
			}

			if (temp_previous_y == temp_y && last_y_repeats == true){
				temp_previous_x = map_from_Col_LG_to_x_y_coordinates[previous_point][2];
				temp_previous_y = map_from_Col_LG_to_x_y_coordinates[previous_point][3];
				map_from_Col_LG_to_x_y_coordinates.erase(previous_point);
				counter_num_points_erased++;
				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_x);
				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_y);
				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_previous_x);
				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_previous_y);

			}
			else if (temp_previous_y != temp_y || last_y_repeats == false)
			{

				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_x);
				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_y);
				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_previous_x);
				map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_previous_y);

			}

			//map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_x);
			//map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_y);
			//map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_previous_x);
			//map_from_Col_LG_to_x_y_coordinates[to_string(i)+"_"+to_string(temp_vect_of_groups_in_column[j])].push_back(temp_previous_y);

			number_of_groups_left = number_of_groups_left - temp_freq;
		}
	}

	cout << "ERASED: " << counter_num_points_erased << endl;

	//make cell name map:
	int last_pos_in_vect = vector_of_cells[0].lin_group_vect.size()-1;
	for (int i = 0; i < vector_of_cells.size(); i++){
		map_of_points_to_cell_names[to_string(last_pos_in_vect)+"_"+to_string(vector_of_cells[i].lin_group_vect[last_pos_in_vect])] = vector_of_cells[i].cell_name;
	}


	ofstream(outfile);
	outfile.open(argv[2]);

	outfile << "x1\ty1\tx2\ty2\tcell_name" << endl;

	string temp_cell_name;

	for (auto it =  map_from_Col_LG_to_x_y_coordinates.begin(); it != map_from_Col_LG_to_x_y_coordinates.end(); it++){
		if (map_of_points_to_cell_names.find(it -> first) != map_of_points_to_cell_names.end()){
			temp_cell_name = map_of_points_to_cell_names[it -> first];
		}
		else {
			temp_cell_name = "NA";
		}

		outfile << it -> second[2] << "\t" << it -> second[3] << "\t" << it -> second[0] << "\t" << it -> second[1] << "\t" << temp_cell_name << endl;
	}

	return 0;
}


