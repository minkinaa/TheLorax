/*
 * 200114_ASEs_for_lineage_groups.cpp
 *  Created on: Oct 16, 2019
 *      Author: anna
 *
 *      argv1 = chromosome
 *      argv2 = length
 *      argv3 = modified sam file for particular chromosome [cell_id, position, cigar, sequence, tab delimited]
 *      argv4 = outfile;
 *      argv5 = min_cutoff;
 *      argv6 = input_file_of_cells,  we'll change this to a list of cells! 1,2,5:8
 *      argv7 = chr_start_pos
 *      argv8 = chr_end_pos
 *		argv9 = cells_and_lineage_groups_file, let's assume a 2 column file, no header
 *		argv10 = list of heterozygous positions
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

struct seq_and_pos_struct{
	string sequence;
	vector<int> positions_vect;
	//bool to_print;
};

string temp_cig_letters;
string return_letters_in_cigar_string(string cigar_string, string possible_letters)
{
	temp_cig_letters = "";
	for (int i = 0; i < cigar_string.length(); i++)
	{
		if (possible_letters.find(cigar_string.substr(i,1)) != possible_letters.npos)
		{
			temp_cig_letters = temp_cig_letters + cigar_string.substr(i,1);
		}
	}
	return (temp_cig_letters);
}

void print_vector(vector<int>vect)
{
	for (int i = 0; i < vect.size(); i++)
	{
		cout << vect[i] << "\t";
	}
	cout << endl;
}


seq_and_pos_struct temp_struct_of_seq_and_pos;
seq_and_pos_struct generate_vector_from_cigar_string(string cigar_string, string seq, string possible_letters, int starting_position){
	temp_struct_of_seq_and_pos = seq_and_pos_struct();
	//temp_struct_of_seq_and_pos.to_print = false;
	string letters_in_cigar = return_letters_in_cigar_string(cigar_string, possible_letters);

	int num_matching;
	string temp_num_matching;

	string new_cigar_string;

	int curr_pos_genome = starting_position;

	int curr_pos_in_sequence = 0;
	string new_sequence = "";


	if (letters_in_cigar == "M")
	{
		num_matching = stoi(cigar_string.substr(0, cigar_string.size() - 1));
		for (int i = curr_pos_genome; i < curr_pos_genome + num_matching; i++)
		{
			temp_struct_of_seq_and_pos.positions_vect.push_back(i);
		}
		temp_struct_of_seq_and_pos.sequence = seq;
		return (temp_struct_of_seq_and_pos);
	}
	else //more than just matching characters
	{
		for (int i = 0; i < cigar_string.size(); i++) //find next letter & save number
		{
			if (possible_letters.find(cigar_string.substr(i,1)) == possible_letters.npos) //if NOT a letter
			{
				temp_num_matching = temp_num_matching + cigar_string.substr(i,1);
			}
			else // if letter
			{
				num_matching = stoi(temp_num_matching);
				if (cigar_string.substr(i,1) == "M" || cigar_string.substr(i,1) == "X" || cigar_string.substr(i,1) == "=")
				{
					for (int j = 0; j < num_matching; j++)
					{
						temp_struct_of_seq_and_pos.positions_vect.push_back(curr_pos_genome);
						new_sequence = new_sequence + seq.substr(curr_pos_in_sequence, 1);
						curr_pos_genome = curr_pos_genome + 1;
						curr_pos_in_sequence = curr_pos_in_sequence + 1;
					}
				}
				else if (cigar_string.substr(i,1) == "D" || cigar_string.substr(i,1) == "N")
				{
					curr_pos_genome = curr_pos_genome + num_matching;
					//cout << curr_pos_genome << endl;
				}
				else if (cigar_string.substr(i,1) == "I" || cigar_string.substr(i,1) == "S")
				{
					curr_pos_in_sequence = curr_pos_in_sequence + num_matching;
				}
				else if (cigar_string.substr(i,1) == "H" || cigar_string.substr(i,1) == "P")
				{
					cout << "H or P here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				}
				else
				{
					cout << "THIS SHOULDNT HAPPEN" << endl;
				}
				temp_num_matching = "";
			}

		}
	}

	temp_struct_of_seq_and_pos.sequence = new_sequence;

	if(temp_struct_of_seq_and_pos.positions_vect.size() != temp_struct_of_seq_and_pos.sequence.size())
	{
		cout << "SIZES DONT MATCH!!" << endl;
		cout << "start=" << starting_position << endl;
		cout << cigar_string << endl;
		cout << "new cigar: " << temp_struct_of_seq_and_pos.sequence << endl;
		print_vector(temp_struct_of_seq_and_pos.positions_vect);
	}

	return (temp_struct_of_seq_and_pos);
}

string temp_char;
int temp_pos;
void add_base_counts_to_vectors(string sequence, vector<int>vect_of_positions, vector<int>&vect_A, vector<int>&vect_C, vector<int>&vect_G, vector<int>&vect_T, vector<int>& total_count_vect)
{
	for (int i = 0; i < sequence.size(); i++)
	{
		temp_char = sequence.substr(i,1);
		temp_pos = vect_of_positions[i];
		if (temp_char == "A")
		{
			vect_A[temp_pos]++;
		}
		else if (temp_char == "C")
		{
			vect_C[temp_pos]++;
		}
		else if (temp_char == "G")
		{
			vect_G[temp_pos]++;
		}
		else if (temp_char == "T")
		{
			vect_T[temp_pos]++;
		}
		total_count_vect[temp_pos]++;
	}
}

string make_cell_id(string cell_id_w_umi)
{
	string temp_cell_id = "";
	temp_cell_id = temp_cell_id + cell_id_w_umi.substr(0, cell_id_w_umi.find("_")) + cell_id_w_umi.substr(cell_id_w_umi.find("_") + 9);

	return (temp_cell_id);
}

int num_commas;
int num_col;
int comma_pos;
vector<int> temp_output_vector;
vector<int> parse_string(string input_group_list)
{
	temp_output_vector.clear();
	num_commas = count(input_group_list.begin(), input_group_list.end(), ',');
	//num_col = count(input_group_list.begin(), input_group_list.end(), ':');
	string temp_grp = "";
	for (int i = 0; i < num_commas; i++)
	{
		comma_pos = input_group_list.find(",");
		temp_grp = input_group_list.substr(0, comma_pos);
		if (temp_grp.find(":") == temp_grp.npos)
		{
			temp_output_vector.push_back(stoi(temp_grp));
		}
		else
		{
			int temp_start_pos = stoi(temp_grp.substr(0, temp_grp.find(":")));
			int temp_end_pos = stoi(temp_grp.substr(temp_grp.find(":") + 1));
			for (int i = temp_start_pos; i < temp_end_pos + 1; i++)
			{
				temp_output_vector.push_back(i);
			}
		}
		input_group_list = input_group_list.substr(comma_pos + 1);
	}
	temp_grp = input_group_list;
	if (temp_grp.find(":") == temp_grp.npos)
	{
		temp_output_vector.push_back(stoi(temp_grp));
	}
	else
	{
		int temp_start_pos = stoi(temp_grp.substr(0, temp_grp.find(":")));
		int temp_end_pos = stoi(temp_grp.substr(temp_grp.find(":") + 1));
		for (int i = temp_start_pos; i < temp_end_pos + 1; i++)
		{
			temp_output_vector.push_back(i);
		}
	}
	cout << "Groups to include: ";
	for (int i = 0; i < temp_output_vector.size(); i++)
	{
		cout << temp_output_vector[i] << " ";
	}
	cout << endl;
	return (temp_output_vector);
	//will need to add a do stuff to last value thing;
}


int main (int argc, char* argv[])
{
	int min_pos = stoi(argv[7]);
	int max_pos = stoi(argv[8]);
	//NEW: read in cells;
	//vector<string> vector_of_cells_to_include();
	map<string,string> map_of_cell_ids;

/*
	ifstream(cells_to_include);
	cells_to_include.open(argv[6]);

	string cell_id;

	while(getline(cells_to_include, cell_id))
	{
		if(map_of_cell_ids.find(cell_id) == map_of_cell_ids.end())
		{
			map_of_cell_ids[cell_id] = cell_id;
			cout << cell_id <<endl;
		}
	}
*/
	string lin_groups_to_include = argv[6];
	vector<int> lin_groups_to_include_vector = parse_string(lin_groups_to_include);
	map<int,int> groups_to_include_map;
	for (int i = 0; i < lin_groups_to_include_vector.size(); i++)
	{
		cout << "adding to map:" << lin_groups_to_include_vector[i] << endl;
		groups_to_include_map[lin_groups_to_include_vector[i]] = lin_groups_to_include_vector[i];
	}

	ifstream(cells_to_include);
	cells_to_include.open(argv[9]);

	string cell_id;
	int temp_lineage_group;
	string line;
	stringstream temp;
	string elem;

	while(getline(cells_to_include, line))
	{
		cell_id = line.substr(0,line.find("\t"));
		temp_lineage_group = stoi(line.substr(line.find("\t") + 1));
		if (groups_to_include_map.find(temp_lineage_group) != groups_to_include_map.end())
		{
			map_of_cell_ids[cell_id] = cell_id;
			//cout << cell_id << endl;
			//cout << temp_lineage_group << endl;
		}
	}
	cout << "Number of cells included: " << map_of_cell_ids.size() << endl;

	int min_cutoff = stoi(argv[5]);
	string pos_letters_in_cigar = "DHIMNPSX=";

	string chromosome = argv[1];
	int chr_length = stoi(argv[2]);

	vector<int> vector_A(chr_length, 0);
	vector<int> vector_T(chr_length, 0);
	vector<int> vector_C(chr_length, 0);
	vector<int> vector_G(chr_length, 0);
	vector<int> vector_total_count(chr_length, 0);

	ifstream(mod_sam_file);
	mod_sam_file.open(argv[3]);

	int counter = 0;

	int temp_start;
	string temp_cigar;
	string temp_sequence;
	string new_sequence;
	vector<int> temp_vector_positions_in_cigar;

	int sam_file_counter = 0;

	seq_and_pos_struct temp_seq_pos_struct;

	while(getline(mod_sam_file, line))
	{
		if(sam_file_counter%100000 == 0)
		{
			cout << sam_file_counter << endl;
		}
		sam_file_counter++;
		//cout << "here!" << endl;
		counter = 0;
		temp << line;
		while(getline(temp, elem, '\t'))
		{
			if (counter == 0)
			{
				cell_id = make_cell_id(elem);
				counter++;
			}
			else if (counter == 1)
			{
				temp_start = stoi(elem);
				counter++;
			}
			else if (counter == 2)
			{
				temp_cigar = elem;
				//cout <<"FROM MAIN CIGAR:" << temp_cigar << endl;
				counter++;
			}
			else if (counter == 3)
			{
				temp_sequence = elem;
				counter ++;
			}
		}
		temp.clear();
		temp_seq_pos_struct = generate_vector_from_cigar_string(temp_cigar, temp_sequence, pos_letters_in_cigar, temp_start);
		temp_vector_positions_in_cigar = temp_seq_pos_struct.positions_vect;
		new_sequence = temp_seq_pos_struct.sequence;
		if (map_of_cell_ids.find(cell_id) != map_of_cell_ids.end())
		{
			add_base_counts_to_vectors(new_sequence, temp_vector_positions_in_cigar, vector_A, vector_C, vector_G, vector_T, vector_total_count);
		}
	}

	ifstream(positions_to_print);
	positions_to_print.open(argv[10]);

	map<int,int> positions_to_print_map;
	while(getline(positions_to_print, line))
	{
		positions_to_print_map[stoi(line)] = stoi(line);
	}

	ofstream(outfile);
	outfile.open(argv[4]);
	outfile << "chr" << "\t" << "pos" << "\t" << "A" << "\t" << "C" << "\t" << "G" << "\t" << "T" << endl;

	for (int i = min_pos; i < max_pos + 1; i++)
	{
		if (vector_total_count[i] > min_cutoff && positions_to_print_map.find(i) != positions_to_print_map.end())
		{
			outfile << chromosome  << "\t"<< i << "\t" << vector_A[i] << "\t" << vector_C[i] << "\t" << vector_G[i] << "\t" << vector_T[i] << endl;
		}
	}

	return 0;
}

