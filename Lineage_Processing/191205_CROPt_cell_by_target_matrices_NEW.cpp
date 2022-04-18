/*
 * 191205_CROPt_cell_by_target_matrices_NEW.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: anna
 */

/*
 * 190323_cropt_cell_by_target_matrices.cpp
 *
 *  Created on: Mar 23, 2019
 *      Author: anna
 *
 *      argv[1]: target list;
 *      argv[2]: output from last program (i.e. "cl36_Alltet_MAIN_OUTPUT_maxdist4")
 *      argv[3]: E_or_U_outfile;
 *      argv[4]: Most_abundant_umi_count_outfile;
 *      argv[5]: Total_umi_count_outfile;
 *      argv[6]: cigar_outfile;
 *      argv[7]: conflicts outfile;
 *      argv[8]: correlation outfile;
 *      argv[9]: low count cutoff (used if collision, but one super low count UMI)
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

using namespace std;

struct cell{
	string cell_id;
	vector<int> barcode_unique_umi_counts_vect;
	vector<int> barcode_most_abundant_umi_counts_vect;
	vector<string> barcode_edited_or_unedited_vect;
	vector<string> barcode_cigars;
	vector<string> edited_seq_conflicts;
	vector<int> correlation_vector;
	vector<int> total_reads_per_umi_vect; //new
};

int same_base_counter;
bool check_hamming_distance(string str1, string str2, int max_hamming_dist)
{
	same_base_counter = 0;
	for (int i = 0; i < str1.size(); i++)
	{
		if(str1.at(i) == str2.at(i))
		{
			same_base_counter++;
		}
	}

	if(str1.length() - same_base_counter <= max_hamming_dist)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool found_matching_seq;
string find_correct_target(vector<string> target_list, string target)
{
	found_matching_seq = false;
	for(int i = 0; i < target_list.size(); i++)
	{
		found_matching_seq = check_hamming_distance(target_list[i], target, 1);
		if (found_matching_seq == true)
		{
			return target_list[i];
		}
	}
	if (found_matching_seq == false)
	{
		return "None";
	}
}

int correlation;
int determine_cell_cell_correlation(vector<string> cell1_cigar_vect,vector<string> cell2_cigar_vect)
{
	correlation = 0;
	for (int i = 0; i < cell1_cigar_vect.size(); i++)
	{
		if(cell1_cigar_vect[i] == "X" || cell2_cigar_vect[i] == "X") //if either not captured, don't add anything
		{
			correlation+=0;
		}
		else if (cell1_cigar_vect[i] == "70M:" && cell2_cigar_vect[i] == "70M:")
		{
			correlation+=1;
		}
		else if (cell1_cigar_vect[i] == cell2_cigar_vect[i] && cell2_cigar_vect[i] != "70M:")
		{
			correlation+=10;
		}
		else if(cell1_cigar_vect[i] != cell2_cigar_vect[i])
		{
			correlation-=10;
		}
		else
		{
			cout << cell1_cigar_vect[i] << " : " << cell2_cigar_vect[i] << endl;
		}
	}
	return correlation;
}


int main(int argc, char* argv[])
{
	int low_count_cutoff = stoi(argv[9]);

	vector<cell> vector_of_cells;
	map<string,int> cell_to_pos_in_vect_map;
	cell temp_cell_struct;

	vector<string> target_vector;
	map<string,int> target_to_pos_in_vect_map;

	ifstream(target_seq_input);
	target_seq_input.open(argv[1]);
	string line;
	string temp_target;

	//cout << "made it??" << endl;
	while(getline(target_seq_input, line))
	{
		temp_target = line.substr(0, 10);
		//cout << temp_target << endl;
		target_vector.push_back(temp_target);
		target_to_pos_in_vect_map[temp_target] = target_vector.size() - 1;
	}

	for (int i = 0; i < target_vector.size(); i++)
	{
		temp_cell_struct.barcode_edited_or_unedited_vect.push_back("X");
		temp_cell_struct.barcode_most_abundant_umi_counts_vect.push_back(0);
		temp_cell_struct.barcode_unique_umi_counts_vect.push_back(0);
		temp_cell_struct.total_reads_per_umi_vect.push_back(0);
		temp_cell_struct.barcode_cigars.push_back("X");
	}

	cell zero_filled_cell_struct = temp_cell_struct;


	ifstream(edited_unedited_input);
	edited_unedited_input.open(argv[2]);
	//string line;

	int temp_umi_count;
	string temp_cell_id;
	string temp_target_id;
	string temp_E_or_U;
	string temp_cigar_string;
	int temp_pos_fisrt_underscore;
	int temp_pos_of_pipe;
	string temp_cell_target_id;
	int temp_target_pos;
	int temp_cell_pos;
	int temp_total_read_count;
	bool real_target = false;

	int counter = 0;

	while(getline(edited_unedited_input, line))
	{
		counter++;
		temp_pos_fisrt_underscore = line.find("_");
		temp_pos_of_pipe = line.find("|");
		//cout << "counter: " << counter << endl;
		//cout << line.substr(0,temp_pos_fisrt_underscore) << endl;
		//temp_umi_count = stoi(line.substr(0,temp_pos_fisrt_underscore));
		temp_umi_count = stoi(line.substr(line.find(")") + 1,temp_pos_fisrt_underscore - (line.find(")") + 1)));
		temp_cell_target_id = line.substr(temp_pos_fisrt_underscore + 1, temp_pos_of_pipe-temp_pos_fisrt_underscore-1);
		temp_E_or_U = line.substr(temp_pos_of_pipe + 1, 1);
		temp_cigar_string = line.substr(line.find("-")+1);
		temp_cell_id = temp_cell_target_id.substr(0,temp_cell_target_id.find_last_of("_"));
		temp_target_id = temp_cell_target_id.substr(temp_cell_target_id.find_last_of("_") + 1);
		temp_total_read_count = stoi(line.substr(1,line.find(")")-1));
		real_target = false;

		if(cell_to_pos_in_vect_map.find(temp_cell_id) != cell_to_pos_in_vect_map.end())
		{
			temp_cell_pos = cell_to_pos_in_vect_map[temp_cell_id];
			if (target_to_pos_in_vect_map.find(temp_target_id) != target_to_pos_in_vect_map.end())
			{
				real_target = true;
				temp_target_pos = target_to_pos_in_vect_map[temp_target_id];
			}
			else
			{
				temp_target_id = find_correct_target(target_vector, temp_target_id);
				if (target_to_pos_in_vect_map.find(temp_target_id) != target_to_pos_in_vect_map.end())
				{
					real_target = true;
					temp_target_pos = target_to_pos_in_vect_map[temp_target_id];
				}
			}
			if(real_target == true && vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] == 0)
			{
				vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos] = temp_cigar_string;
				vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] = temp_E_or_U;
				vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] = temp_umi_count;
				vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos] = temp_umi_count;
				vector_of_cells[temp_cell_pos].total_reads_per_umi_vect[temp_target_pos] = temp_total_read_count;
			}
			else if(real_target == true) //if target already recorded in cell
			{
				if (vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] == "U" && temp_E_or_U == "U") //if both unedited
				{
					vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos]+=temp_umi_count;
					vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+=temp_umi_count;
				}
				else if (vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] == "E" || temp_E_or_U == "E")
				{
					if (temp_umi_count == 1 && vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] > 1)
					{
						vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
					}
					else if(vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] == 1 && temp_umi_count > 1)
					{
						vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos] = temp_cigar_string;
						vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] = temp_E_or_U;
						vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] = temp_umi_count;
						vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
						vector_of_cells[temp_cell_pos].total_reads_per_umi_vect[temp_target_pos] = temp_total_read_count;
					}
					else if(temp_umi_count > vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos])
					{
						vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(line);
						vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(to_string(vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos]) + '_' + temp_cell_target_id + "|" + vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] + "-" + vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos]);
						vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] = "A";
						vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] = temp_umi_count;
						vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
						vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos] = temp_cigar_string;
						vector_of_cells[temp_cell_pos].total_reads_per_umi_vect[temp_target_pos] = temp_total_read_count;
					}
					else if(temp_umi_count < vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos])
					{
						vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(line);
						vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(to_string(vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos]) + '_' + temp_cell_target_id + "|" + vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] + "-" + vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos]);
						vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] = "A";
						vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
					}
					else if (temp_umi_count == 1 && vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] == 1) //this is the new one!!
					{
						//cout << temp_cell_id + "_" + temp_target_id << endl;
						if (temp_cigar_string == vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos])
						{
							vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos]+=temp_umi_count;
							vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
						}
						else if (temp_cigar_string != vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos])
						{
							if (vector_of_cells[temp_cell_pos].total_reads_per_umi_vect[temp_target_pos] > low_count_cutoff && temp_total_read_count > low_count_cutoff)
							{
								//cout << "yep, here" << endl;
								vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(line);
								vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(to_string(vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos]) + '_' + temp_cell_target_id + "|" + vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] + "-" + vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos]);
								vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] = "A";
								vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
							}
							else if (vector_of_cells[temp_cell_pos].total_reads_per_umi_vect[temp_target_pos] > low_count_cutoff)
							{
								vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
							}
							else if (temp_total_read_count > low_count_cutoff)
							{
								vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos] = temp_cigar_string;
								vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] = temp_E_or_U;
								vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos] = temp_umi_count;
								vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
								vector_of_cells[temp_cell_pos].total_reads_per_umi_vect[temp_target_pos] = temp_total_read_count;
							}
						}
					}
					else if (temp_umi_count == vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos])
					{
						vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(line);
						vector_of_cells[temp_cell_pos].edited_seq_conflicts.push_back(to_string(vector_of_cells[temp_cell_pos].barcode_most_abundant_umi_counts_vect[temp_target_pos]) + '_' + temp_cell_target_id + "|" + vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] + "-" + vector_of_cells[temp_cell_pos].barcode_cigars[temp_target_pos]);
						vector_of_cells[temp_cell_pos].barcode_edited_or_unedited_vect[temp_target_pos] = "A";
						vector_of_cells[temp_cell_pos].barcode_unique_umi_counts_vect[temp_target_pos]+= temp_umi_count;
					}
					else
					{
						cout << "this shouldn't happen" << endl;
					}
				}
			}
		}
		else
		{
			temp_cell_struct.cell_id = temp_cell_id;
			if (target_to_pos_in_vect_map.find(temp_target_id) != target_to_pos_in_vect_map.end())
			{
				real_target = true;
				temp_target_pos = target_to_pos_in_vect_map[temp_target_id];
				temp_cell_struct.barcode_edited_or_unedited_vect[temp_target_pos] = temp_E_or_U;
				temp_cell_struct.barcode_cigars[temp_target_pos] = temp_cigar_string;
				temp_cell_struct.barcode_most_abundant_umi_counts_vect[temp_target_pos] = temp_umi_count;
				temp_cell_struct.barcode_unique_umi_counts_vect[temp_target_pos] = temp_umi_count;
				temp_cell_struct.total_reads_per_umi_vect[temp_target_pos] = temp_total_read_count;
				vector_of_cells.push_back(temp_cell_struct);
				cell_to_pos_in_vect_map[temp_cell_id] = vector_of_cells.size() - 1;
				temp_cell_struct = zero_filled_cell_struct;
			}
		}
	}

	for(int i = 0; i < vector_of_cells.size(); i++)
	{
		for(int j = 0; j <  vector_of_cells.size(); j++)
		{
			vector_of_cells[i].correlation_vector.push_back(determine_cell_cell_correlation(vector_of_cells[i].barcode_cigars, vector_of_cells[j].barcode_cigars));
		}
	}



	ofstream E_or_U_outfile;
	E_or_U_outfile.open(argv[3]);

	ofstream Most_abundant_umi_count_outfile;
	Most_abundant_umi_count_outfile.open(argv[4]);

	ofstream Total_umi_count_outfile;
	Total_umi_count_outfile.open(argv[5]);

	ofstream cigar_outfile;
	cigar_outfile.open(argv[6]);

	ofstream conflicts_outfile;
	conflicts_outfile.open(argv[7]);

	ofstream correlation_outfile;
	correlation_outfile.open(argv[8]);

	E_or_U_outfile << "Targets:";
	Most_abundant_umi_count_outfile << "Targets:";
	Total_umi_count_outfile << "Targets:";
	cigar_outfile << "Targets:";
	//correlation_outfile << "Targets:";

	for (int i = 0; i < target_vector.size(); i++)
	{
		E_or_U_outfile << "\t" << target_vector[i];
		Most_abundant_umi_count_outfile << "\t" << target_vector[i];
		Total_umi_count_outfile << "\t" << target_vector[i];
		cigar_outfile << "\t" << target_vector[i];
		//correlation_outfile << "\t" << target_vector[i];
	}

	E_or_U_outfile << endl;
	Most_abundant_umi_count_outfile << endl;
	Total_umi_count_outfile << endl;
	cigar_outfile << endl;
	//correlation_outfile << endl;

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		E_or_U_outfile << vector_of_cells[i].cell_id;
		Most_abundant_umi_count_outfile << vector_of_cells[i].cell_id;
		Total_umi_count_outfile << vector_of_cells[i].cell_id;
		cigar_outfile << vector_of_cells[i].cell_id;
		//correlation_outfile << vector_of_cells[i].cell_id;

		for (int j = 0; j < target_vector.size(); j++)
		{
			E_or_U_outfile << "\t" <<  vector_of_cells[i].barcode_edited_or_unedited_vect[j];
			Most_abundant_umi_count_outfile << "\t" << vector_of_cells[i].barcode_most_abundant_umi_counts_vect[j];
			Total_umi_count_outfile << "\t" << vector_of_cells[i].barcode_unique_umi_counts_vect[j];
			cigar_outfile << "\t" << vector_of_cells[i].barcode_cigars[j];
			//correlation_outfile << "\t" << vector_of_cells[i].correlation_vector[j];
		}
		E_or_U_outfile << endl;
		Most_abundant_umi_count_outfile << endl;
		Total_umi_count_outfile << endl;
		cigar_outfile << endl;
		//correlation_outfile << endl;

		for (int k = 0; k < vector_of_cells[i].edited_seq_conflicts.size(); k++)
		{
			conflicts_outfile << vector_of_cells[i].edited_seq_conflicts[k] << endl;
		}
	}

	correlation_outfile << "Cells:";
	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		correlation_outfile << "\t" << vector_of_cells[i].cell_id;
	}
	correlation_outfile << endl;

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		correlation_outfile << vector_of_cells[i].cell_id;
		for (int j = 0; j < vector_of_cells.size(); j++)
		{
			correlation_outfile << "\t" << vector_of_cells[i].correlation_vector[j];
		}
		correlation_outfile << endl;
	}






/*
	 *      argv[3]: E_or_U_outfile;
	 *      argv[4]: Most_abundant_umi_count_outfile;
	 *      argv[5]: Total_umi_count_outfile;
	 *      argv[6]: cigar_outfile;
	 *      argv[7]: conflicts outfile;
*/



	return 0;
}





