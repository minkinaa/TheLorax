/*
 * 190219_make_cell_by_barcode_matrix.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: anna
 *
 *      argv[1] output from 190216_count_UMIs_per_barcode
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
	vector<int> barcode_unique_umi_counts_vect;
	//map<string, int> cell_barcodes_to_counts_map;
};

int main(int argc, char* argv[])
{
	vector<cell> vector_of_cells;
	map<string,int> cell_to_pos_in_vect_map;
	cell temp_cell_struct;

	ifstream(umis_per_barcode_input);
	umis_per_barcode_input.open(argv[1]);
	string line;

	vector<string> barcode_vector;
	map<string,int> barcode_to_pos_in_vect_map;

	string gt_barcode;
	string gt_barcode_umis;

	//filling barcode vector;
	while(getline(umis_per_barcode_input, line)){
		gt_barcode = line.substr(line.find("_") + 12, 10);
		//gt_barcode_umis = line.substr(line.find("_") + 23, line.find("|") - (line.find("_") + 23));
		if (barcode_to_pos_in_vect_map.find(gt_barcode) ==  barcode_to_pos_in_vect_map.end())
		{
			barcode_vector.push_back(gt_barcode);
			barcode_to_pos_in_vect_map[gt_barcode] = barcode_vector.size() - 1;
		}
	}
	umis_per_barcode_input.close();

	//filling vector for structs;

	string temp_cell_id;

	umis_per_barcode_input.open(argv[1]);
	while (getline(umis_per_barcode_input, line)){
		temp_cell_id = line.substr(0, line.find("_") + 11);
		gt_barcode = line.substr(line.find("_") + 12, 10);
		gt_barcode_umis = line.substr(line.find("_") + 23, line.find("|") - (line.find("_") + 23));
		if (cell_to_pos_in_vect_map.find(temp_cell_id) != cell_to_pos_in_vect_map.end())
		{
			vector_of_cells[cell_to_pos_in_vect_map[temp_cell_id]].barcode_unique_umi_counts_vect[barcode_to_pos_in_vect_map[gt_barcode]] = stoi(gt_barcode_umis);
		}
		else
		{
			temp_cell_struct.cell_id = temp_cell_id;
			for (int i = 0; i < barcode_vector.size(); i++){
				temp_cell_struct.barcode_unique_umi_counts_vect.push_back(0);
			}
			temp_cell_struct.barcode_unique_umi_counts_vect[barcode_to_pos_in_vect_map[gt_barcode]] = stoi(gt_barcode_umis);
			vector_of_cells.push_back(temp_cell_struct);
			cell_to_pos_in_vect_map[temp_cell_id] = vector_of_cells.size() -1;
			temp_cell_struct = cell();
		}
	}

	cout << "Barcodes:";

	for (int i = 0; i < barcode_vector.size(); i++)
	{
		cout << "\t" << barcode_vector[i];
	}
	cout << endl;

	//cout << barcode_vector.size() << endl;
	//cout << vector_of_cells[1].barcode_unique_umi_counts_vect.size() << endl;


	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		cout << vector_of_cells[i].cell_id;
		for (int j = 0; j < barcode_vector.size(); j++)
		{
			cout << "\t" << vector_of_cells[i].barcode_unique_umi_counts_vect[j];
		}
		cout << endl;
	}

	return 0;
}




