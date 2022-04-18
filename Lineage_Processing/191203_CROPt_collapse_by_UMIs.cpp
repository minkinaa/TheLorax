/*
 * 191203_CROPt_collapse_by_UMIs.cpp
 *
 *  Created on: Dec 3, 2019
 *      Author: anna
 *
 *      argv[1]: ${PCR_well}_UMI_RT_BC_seq_output file
 *      argv[2]: file of barcodes
 *      argv[3]: file or RT indices
 *      argv[4]: outfile
 *      argv[5]: num umis per seq outfile
 *      argv[6]: well name (i.e. A10)
 *      argv[7]: plate_num
 *      argv[8]: fastq output
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

using namespace std;

struct unique_UMI_RT_BC_combo{
	string UMI_RT_BC;
	map<string,int> seq_to_count_map;
	int total_count;
	string consensus_seq;
	int consensus_count;
	bool real_seq;
};

struct seq_count_pair{
	string seq;
	int count;
};

struct num_and_count_umis_pair{
	int number_of_umis;
	int total_read_count;
};

struct RT_BC_combo_struct{
	string RT_BC;
	map<string, seq_count_pair> umi_to_count_map;
	bool has_low_count;
	map<string, num_and_count_umis_pair> sequence_to_umi_count_map;
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

bool found_matching;
string correct_sequence (string orig_seq, map<string,string> &map_of_real_seqs)
{
	found_matching = false;
	for(auto it = map_of_real_seqs.begin(); it != map_of_real_seqs.end(); it++)
	{
		//cout << "comparing to" << it-> first << endl;
		found_matching = check_hamming_distance(orig_seq, it -> first, 1);
		if (found_matching == true)
		{
			return it -> first;
			it = map_of_real_seqs.end();
		}
	}

	if (found_matching == false)
	{
		return ("NONE");
	}
	else
	{
		return ("SHOULD-NOT-HAPPEN!!");
	}
}

int main(int argc, char* argv[])
{
	string well_name = argv[6];
	string plate_num = argv[7];
	plate_num = "p" + plate_num;
	int low_count_cutoff = 2;
	int low_counts_to_correct_cutoff = 10;
	map <string, string> real_barcode_map;
	map <string, string> real_RT_index_map;

	ifstream real_barcode_file;
	real_barcode_file.open(argv[2]);

	ifstream real_RT_file;
	real_RT_file.open(argv[3]);

	string line;

	while(getline(real_barcode_file, line))
	{
		real_barcode_map[line.substr(0,10)] = line.substr(0,10);
		//cout << line << endl;
	}

	while(getline(real_RT_file, line))
	{
		real_RT_index_map[line.substr(0,10)] = line.substr(0,10);
		//cout << line.substr(0,10) << endl;
		//cout << line.substr(0,10) << "wtf" << endl;
	}

	ifstream UMI_RT_BC_file;
	UMI_RT_BC_file.open(argv[1]);
	string temp_UMI_RT_BC;
	string temp_sequence;
	string corrected_temp_UMI_RT_BC;
	string temp_UMI;
	string temp_RT;
	string temp_BC;
	bool real_combo;
	//int pos_of_first_underscore;

	vector<unique_UMI_RT_BC_combo> unique_UMI_RT_BC_combo_vect;
	map<string, int> map_BC_to_vect_pos;
	unique_UMI_RT_BC_combo temp_struct;
	int pos_in_vect;


	while(getline(UMI_RT_BC_file, line))
	{
		temp_UMI_RT_BC = line.substr(0,line.find_last_of("_"));
		//cout << temp_UMI_RT_BC << endl;
		temp_sequence = line.substr(31);
		temp_UMI = temp_UMI_RT_BC.substr(0,8);
		temp_RT = temp_UMI_RT_BC.substr(9,10);
		temp_BC = temp_UMI_RT_BC.substr(20,10);
		real_combo = false;

		if (real_RT_index_map.find(temp_RT) == real_RT_index_map.end()) //if not in map;
		{
			//cout << "here??" << endl;
			temp_RT = correct_sequence(temp_RT, real_RT_index_map);
		}
		if (real_barcode_map.find(temp_BC) == real_barcode_map.end()) //if not in map;
		{
			//cout << "BC:" << temp_BC << endl;
			temp_BC = correct_sequence(temp_BC, real_barcode_map);
			//cout << "Corrected BC:" << temp_BC << endl;
		}

		if (temp_RT != "NONE" && temp_BC != "NONE")
		{
			corrected_temp_UMI_RT_BC = temp_UMI + "_" + temp_RT + "_" + temp_BC;
			real_combo = true;
		}

		if (real_combo == true && map_BC_to_vect_pos.find(corrected_temp_UMI_RT_BC) == map_BC_to_vect_pos.end()) // if not in map;
		{
			temp_struct.UMI_RT_BC = corrected_temp_UMI_RT_BC;
			temp_struct.seq_to_count_map[temp_sequence] = 1;
			temp_struct.total_count = 1;
			//temp_struct.seq_to_count_map[temp_UMI] = 1;
			unique_UMI_RT_BC_combo_vect.push_back(temp_struct);
			map_BC_to_vect_pos[temp_struct.UMI_RT_BC] = unique_UMI_RT_BC_combo_vect.size() - 1;
			temp_struct = unique_UMI_RT_BC_combo();
		}
		else if (real_combo == true && map_BC_to_vect_pos.find(corrected_temp_UMI_RT_BC) != map_BC_to_vect_pos.end())
		{
			pos_in_vect = map_BC_to_vect_pos[corrected_temp_UMI_RT_BC];
			if (unique_UMI_RT_BC_combo_vect[pos_in_vect].seq_to_count_map.find(temp_sequence) != unique_UMI_RT_BC_combo_vect[pos_in_vect].seq_to_count_map.end())
			{
				unique_UMI_RT_BC_combo_vect[pos_in_vect].seq_to_count_map[temp_sequence]++;
				temp_struct.total_count++;
			}
			else
			{
				unique_UMI_RT_BC_combo_vect[pos_in_vect].seq_to_count_map[temp_sequence] = 1;
				temp_struct.total_count++;
			}
		}
	}

	int temp_max;
	string temp_best_seq;
	// determine consensus sequence for each UMI
	for (int i = 0; i < unique_UMI_RT_BC_combo_vect.size(); i++)
	{
		temp_max = 0;
		for (auto it = unique_UMI_RT_BC_combo_vect[i].seq_to_count_map.begin(); it != unique_UMI_RT_BC_combo_vect[i].seq_to_count_map.end(); it++)
		{
			if (it -> second > temp_max)
			{
				temp_max = it -> second;
				temp_best_seq = it -> first;
			}
		}
		unique_UMI_RT_BC_combo_vect[i].consensus_seq = temp_best_seq;
		unique_UMI_RT_BC_combo_vect[i].consensus_count = temp_max;
	}

	//run through vector, make a new vector of RT_BC -> UMI + seq/count
	RT_BC_combo_struct temp_RT_BC_struct;
	vector<RT_BC_combo_struct>RT_BC_combo_struct_vect;
	map<string, int> RT_BC_to_pos_in_vect_map;
	string temp_RT_BC;
	string temp_consensus;
	int temp_consensus_count;

	for (int i = 0; i < unique_UMI_RT_BC_combo_vect.size(); i++)
	{
		temp_RT_BC = unique_UMI_RT_BC_combo_vect[i].UMI_RT_BC.substr(9,21);
		temp_consensus = unique_UMI_RT_BC_combo_vect[i].consensus_seq;
		temp_consensus_count = unique_UMI_RT_BC_combo_vect[i].consensus_count;
		temp_UMI = unique_UMI_RT_BC_combo_vect[i].UMI_RT_BC.substr(0,8);
		if (RT_BC_to_pos_in_vect_map.find(temp_RT_BC) == RT_BC_to_pos_in_vect_map.end() && temp_consensus_count > low_count_cutoff)
		{
			temp_RT_BC_struct.RT_BC = temp_RT_BC;
			temp_RT_BC_struct.umi_to_count_map[temp_UMI].count =temp_consensus_count;
			temp_RT_BC_struct.umi_to_count_map[temp_UMI].seq = temp_consensus;
			temp_RT_BC_struct.has_low_count = false;
			if (temp_consensus_count < low_counts_to_correct_cutoff)
			{
				temp_RT_BC_struct.has_low_count = true;
			}
			RT_BC_combo_struct_vect.push_back(temp_RT_BC_struct);
			RT_BC_to_pos_in_vect_map[temp_RT_BC] = RT_BC_combo_struct_vect.size() - 1;
			temp_RT_BC_struct = RT_BC_combo_struct();
		}
		else if (RT_BC_to_pos_in_vect_map.find(temp_RT_BC) != RT_BC_to_pos_in_vect_map.end() && temp_consensus_count > low_count_cutoff)
		{
			pos_in_vect = RT_BC_to_pos_in_vect_map[temp_RT_BC];
			RT_BC_combo_struct_vect[pos_in_vect].umi_to_count_map[temp_UMI].count = temp_consensus_count;
			RT_BC_combo_struct_vect[pos_in_vect].umi_to_count_map[temp_UMI].seq = temp_consensus;
			if (temp_consensus_count < low_counts_to_correct_cutoff)
			{
				RT_BC_combo_struct_vect[pos_in_vect].has_low_count = true;
			}
		}
	}

	//UMI counts per RT_BC pair - maybe not yet??
	string low_count_UMI;
	string best_match_UMI;
	bool found_a_match;
	for (int i = 0; i < RT_BC_combo_struct_vect.size(); i++)
	{
		if (RT_BC_combo_struct_vect[i].has_low_count == true)
		{
			//cout << "Found_low_count" << endl;
			for (auto it = RT_BC_combo_struct_vect[i].umi_to_count_map.begin(); it != RT_BC_combo_struct_vect[i].umi_to_count_map.end(); it++)
			{
				if (it -> second.count < low_counts_to_correct_cutoff)
				{
					//cout << it -> second.count << endl;
					//cout << it -> first << endl;
					found_a_match = false;
					auto it2 = RT_BC_combo_struct_vect[i].umi_to_count_map.begin();
					//cout << "map size: " << RT_BC_combo_struct_vect[i].umi_to_count_map.size() << endl;
					while (found_a_match == false && it2 != RT_BC_combo_struct_vect[i].umi_to_count_map.end())
					{
						//cout << it2 -> first << endl;
						if (it -> first != it2 -> first)
						{
							//cout << "comparing to " << it2 -> first << endl;
							found_a_match = check_hamming_distance(it -> first, it2 -> first, 1);
							if (found_a_match == true & it2 -> second.count < low_counts_to_correct_cutoff)
							{
								found_a_match = false;
							}
						}
						if (found_a_match == false)
						{
							it2++;
						}
					}
					if (found_a_match == true)
					{
						//cout << "found match!!" << endl;
						//cout << it -> first << ":" << it2 -> first << endl;
						//cout << it2 -> second.count << endl;
						it2 -> second.count = it2 -> second.count + it -> second.count;
						//cout << it2 -> second.count << endl;
						it -> second.count = 0;
						//cout << it -> second.count << endl;
					}
				}
			}
		}
	}

	//print collapsed by UMIs
	ofstream(output_file);
	output_file.open(argv[4]);
	string temp_consensus_seq;

	bool temp_bool = false;
	for (int i = 0; i < RT_BC_combo_struct_vect.size(); i++)
	{
		temp_bool = false;
		if (RT_BC_combo_struct_vect[i].RT_BC == "TGCGAAGATC_GCTTATCCGC")
		{
			//cout << "OK HERE WE ARE" << endl;
			temp_bool = true;
		}
		for(auto it = RT_BC_combo_struct_vect[i].umi_to_count_map.begin(); it != RT_BC_combo_struct_vect[i].umi_to_count_map.end(); it++)
		{
			if (temp_bool == true)
			{
				//cout << it -> first << " : " << it -> second.count << endl;
			}
			if (it ->second.count > low_count_cutoff)
			{
				output_file << it -> first << "_" << RT_BC_combo_struct_vect[i].RT_BC << "_" << it -> second.count << "|" << it -> second.seq << endl;
			}
		}
	}

	//Now count UMIs
	string temp_seq;
	for (int i = 0; i < RT_BC_combo_struct_vect.size(); i++)
	{
		for (auto it = RT_BC_combo_struct_vect[i].umi_to_count_map.begin(); it != RT_BC_combo_struct_vect[i].umi_to_count_map.end(); it++)
		{
			temp_seq = it ->second.seq;
			if (RT_BC_combo_struct_vect[i].sequence_to_umi_count_map.find(temp_seq) == RT_BC_combo_struct_vect[i].sequence_to_umi_count_map.end())
			{
				if (it -> second.count > low_count_cutoff)
				{
					RT_BC_combo_struct_vect[i].sequence_to_umi_count_map[temp_seq].number_of_umis = 1;
					RT_BC_combo_struct_vect[i].sequence_to_umi_count_map[temp_seq].total_read_count = it -> second.count;
				}
			}
			else if (RT_BC_combo_struct_vect[i].sequence_to_umi_count_map.find(temp_seq) != RT_BC_combo_struct_vect[i].sequence_to_umi_count_map.end())
			{
				if (it -> second.count > low_count_cutoff)
				{
					RT_BC_combo_struct_vect[i].sequence_to_umi_count_map[temp_seq].number_of_umis++;
					RT_BC_combo_struct_vect[i].sequence_to_umi_count_map[temp_seq].total_read_count += it -> second.count;
				}
			}
		}
	}

	ofstream(second_outfile);
	second_outfile.open(argv[5]);

	ofstream(fastq_output);
	fastq_output.open(argv[8]);

	for (int i = 0; i < RT_BC_combo_struct_vect.size(); i++)
	{
		for (auto it = RT_BC_combo_struct_vect[i].sequence_to_umi_count_map.begin(); it != RT_BC_combo_struct_vect[i].sequence_to_umi_count_map.end(); it++)
		{
			second_outfile << plate_num + well_name << "_" << RT_BC_combo_struct_vect[i].RT_BC << "_" << it -> second.number_of_umis << "-" << it -> second.total_read_count << "|" << it -> first << endl;
			fastq_output << ">" << it -> second.total_read_count << "_" << it -> second.number_of_umis << "_" << plate_num + well_name << "_" << RT_BC_combo_struct_vect[i].RT_BC << endl;
			fastq_output << it -> first << endl;
		}
	}
	return 0;
}





