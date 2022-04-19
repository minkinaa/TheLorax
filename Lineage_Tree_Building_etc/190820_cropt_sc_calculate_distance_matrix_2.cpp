/*
 * 190820_cropt_sc_calculate_distance_matrix_2.cpp
 *
 *  Created on: Aug 20, 2019
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
	normalized_file.open(argv[8])

	New:
	argv[8]weight if both are missing
	argv[9]weight if both are ambiguous

	normalized_file.open(argv[10])
	raw_file.open(argv[11])

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
	vector<float> normalized_vector_of_distances;
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

void print_normalized_distance(cell temp_cell)
{
	cout << temp_cell.cell_id << endl;
	cout << temp_cell.normalized_vector_of_distances.size() << endl;
	for (int i = 0; i < temp_cell.normalized_vector_of_distances.size(); i++)
	{
		cout << temp_cell.normalized_vector_of_distances[i] << "\t";
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

	vector<cell> vector_of_cells;
	cell temp_cell;

	ifstream(cigar_file);
	cigar_file.open(argv[1]);

	string line;
	stringstream temp;
	string elem;
	int counter = 0;

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
			}
		}
		vector_of_cells.push_back(temp_cell);
		temp.clear();
		temp_cell = cell();
	}

	//print_cell(vector_of_cells[1]);
	//cout << vector_of_cells.size();

	//CALCULATE DISTACE MATRIX
	int distance;
	int self_dist;
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		self_dist = calculate_distance(unedited, vector_of_cells[i].editing_patt_vect, vector_of_cells[i].editing_patt_vect, weight_same_edit, weight_both_70M, weight_one_X, weight_at_least_one_AMB, weight_diff_edits, weight_one_70_one_edited, weight_both_X, weight_both_AMB);
		for (int j = 0; j < vector_of_cells.size(); j++)
		{
			distance = calculate_distance(unedited, vector_of_cells[i].editing_patt_vect, vector_of_cells[j].editing_patt_vect, weight_same_edit, weight_both_70M, weight_one_X, weight_at_least_one_AMB, weight_diff_edits, weight_one_70_one_edited, weight_both_X, weight_both_AMB);
			vector_of_cells[i].vector_of_distances.push_back(distance);
			vector_of_cells[i].normalized_vector_of_distances.push_back(float(distance)/float(self_dist));
		}
	}

	ofstream(normalized_file);
	normalized_file.open(argv[10]);

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		normalized_file << vector_of_cells[i].cell_id;
		for (int j = 0; j < vector_of_cells[i].normalized_vector_of_distances.size(); j++)
		{
			normalized_file << "\t" << vector_of_cells[i].normalized_vector_of_distances[j];
		}
		normalized_file << endl;
	}
	normalized_file.close();

	ofstream(raw_file);
	raw_file.open(argv[11]);

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		raw_file << vector_of_cells[i].cell_id;
		for (int j = 0; j < vector_of_cells[i].vector_of_distances.size(); j++)
		{
			raw_file << "\t" << vector_of_cells[i].vector_of_distances[j];
		}
		raw_file << endl;
	}
	raw_file.close();



	 //print_normalized_distance(vector_of_cells[2]);




	return 0;
}







