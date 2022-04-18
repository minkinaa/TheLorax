/*
 * 191205_CROPt_fasta_to_MAIN_OUTPUT.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: anna
 */

/*
 *
 *      argv[1] = fasta input file;
 *      argv[2] = outfile for matrix plot
 *      argv[3] = edited matrix outfile
 *      argv[4] = unedited matrix outfile
 *      argv[5] = max dist from cut site !!!!!!
 *      argv[6] = outfile for next program
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

struct cell_target_combo{
	int total_umi_count;
	string PCR_RT_Targ;
	string ref_seq;
	string target_seq;
	string long_alignment;
	string long_aligment_short_insertions;
	string cigar_string;
	int umis;
	int num_sequences;
	bool edited;
	int num_matches_at_start;
	int num_matches_at_end;
	int num_Is;
	string IsNs;
	string cigar_plus_IsNs;

	bool at_least_one_edited;
	int num_unedited;
	int num_edited; //????????
	bool correct_size;
	bool split_edit;
	map<string,int> edited_CIGAR_to_count;

	bool print_this = true;
};

string make_alignment_string(string ref, string target)
{
	string temp_alignment;
	for (int i = 0; i < ref.size(); i++)
	{
		if (ref.substr(i,1) == target.substr(i,1)){
			temp_alignment = temp_alignment + "M";
		}
		else if (ref.substr(i,1) != "-" && target.substr(i,1) == "-")
		{
			temp_alignment = temp_alignment + "D";
		}
		else if (ref.substr(i,1) == "-" && target.substr(i,1) != "-")
		{
			temp_alignment = temp_alignment + "I";
		}
		else if (ref.substr(i,1) != target.substr(i,1))
		{
			temp_alignment = temp_alignment + "N";
		}
	}
	return temp_alignment;
}
string generate_cigar_string(string long_alignment)
{
	int same_base_counter = 1;
	string temp_cigar = "";
	string curr_base;
	string last_base;
	last_base = long_alignment.substr(0,1);
	for (int i = 1; i < long_alignment.size(); i++)
	{
		curr_base = long_alignment.substr(i,1);
		if (curr_base == last_base)
		{
			same_base_counter++;
		}
		else
		{
			temp_cigar = temp_cigar + to_string(same_base_counter) + last_base;
			same_base_counter = 1;
			last_base = curr_base;
		}
	}
	temp_cigar = temp_cigar + to_string(same_base_counter) + last_base;
	return temp_cigar;
}

string matches_at_beginning = "";
int determine_num_matches_at_start(string cigar_str)
{
	int i = 0;
	matches_at_beginning = "";
	while (cigar_str.substr(i,1) != "M" && cigar_str.substr(i,1) != "N" && cigar_str.substr(i,1) != "D" && cigar_str.substr(i,1) != "I")
	{
		matches_at_beginning = matches_at_beginning + cigar_str.substr(i,1);
		i++;
	}
	if (cigar_str.substr(i,1) == "M")
	{
		return (stoi(matches_at_beginning));
		//return (matches_at_beginning);
	}
	else
	{
		return (0);
	}
}

string matches_at_end = "";
int determine_num_matches_at_end(string cigar_str)
{
	int i = cigar_str.size() - 1;
	matches_at_end = "";
	if (cigar_str.substr(i,1) == "M")
	{
		i--;
		while (cigar_str.substr(i,1) != "M" && cigar_str.substr(i,1) != "N" && cigar_str.substr(i,1) != "D" && cigar_str.substr(i,1) != "I")
		{
			matches_at_end = cigar_str.substr(i,1) + matches_at_end;
			i--;
		}
		return (stoi(matches_at_end));
	}
	else if (cigar_str.substr(i,1) != "M")
	{
		return (0);
	}
}

string short_alignment;
int leading_insertions_counter;
string convert_alignment_to_short_insertions(string long_align)
{
	short_alignment = "";
	leading_insertions_counter = 0;
	if (long_align.find("I") == long_align.npos)
	{
		return long_align;
	}
	else
	{
		for (int i = 0; i < long_align.size(); i++)
		{
			if (long_align.substr(i,1) != "I")
			{
				short_alignment = short_alignment + long_align.substr(i,1);
			}
			else if (long_align.substr(i,1) == "I" && short_alignment.size() > 0)
			{
				short_alignment.replace(short_alignment.size()-1, 1, "I");
			}
			else if (long_align.substr(i,1) == "I" && short_alignment.size() == 0)
			{
				leading_insertions_counter++;
			}
		}
		if (leading_insertions_counter > 0)
		{
			//cout << long_align << endl;
			short_alignment.replace(0,1, "I");
		}
	}
	return short_alignment;
}

string bases_near_cutsite;
bool determine_if_edited_from_short_alignment(string short_alignment, int last_base_before_cut_site, int max_dist_from_cutsite){
	bases_near_cutsite = short_alignment.substr(last_base_before_cut_site - max_dist_from_cutsite + 1, max_dist_from_cutsite*2);
	if(bases_near_cutsite.find("N") == bases_near_cutsite.npos && bases_near_cutsite.find("I") == bases_near_cutsite.npos && bases_near_cutsite.find("D") == bases_near_cutsite.npos)
	{
		return false;
	}
	else
	{
		return true;
	}
}

int curr_start_pos;
int curr_pos;
int matching_counter;
bool at_least_five_matches;
int matches_in_front;
int matches_in_back;
string short_seq;
string del_ins_cigar;
string new_long_alignment;

string correct_long_alignment_string(string short_alignment, string long_alignment, int last_base_before_cut_site, int max_dist_from_cutsite, int max_matches){
	curr_start_pos = last_base_before_cut_site;
	curr_pos = last_base_before_cut_site;
	at_least_five_matches = false;
	matching_counter = 0;
	new_long_alignment = "";
	matches_in_front = 0;
	matches_in_back = 0;

	while(at_least_five_matches == false)
	{
		if(short_alignment.substr(curr_pos,1)  == "M")
		{
			matching_counter++;
		}
		else if(short_alignment.substr(curr_pos,1)  != "M")
		{
			matching_counter = 0;
			curr_start_pos = curr_pos - 1;
		}

		if (matching_counter == max_matches)
		{
			matches_in_front = curr_start_pos;
			at_least_five_matches = true;
		}
		if (curr_pos == 0 && at_least_five_matches == false)
		{
			at_least_five_matches = true;
			matches_in_front = 0;
		}

		curr_pos --;
	}

	curr_start_pos = last_base_before_cut_site + 1;
	curr_pos = last_base_before_cut_site + 1;
	at_least_five_matches = false;
	matching_counter = 0;

	while(at_least_five_matches == false)
	{
		if(short_alignment.substr(curr_pos,1)  == "M")
		{
			matching_counter++;
		}
		else if(short_alignment.substr(curr_pos,1)  != "M")
		{
			matching_counter = 0;
			curr_start_pos = curr_pos + 1;
		}

		if (matching_counter == max_matches)
		{
			matches_in_back = short_alignment.size() - curr_start_pos;
			at_least_five_matches = true;
		}
		if (curr_pos == short_alignment.size() - 1)
		{
			matches_in_back = 0;
			at_least_five_matches = true;
		}
		curr_pos ++;
	}
	//cout << "matches_front: " << matches_in_front << endl;

	for (int i = 0; i < matches_in_front; i++)
	{
		new_long_alignment = new_long_alignment + "M";
	}
	new_long_alignment = new_long_alignment + long_alignment.substr(matches_in_front, long_alignment.size() - matches_in_front - matches_in_back);
	for (int i = 0; i < matches_in_back; i++)
	{
		new_long_alignment = new_long_alignment + "M";
	}
/*
	if (new_long_alignment != long_alignment && long_alignment != "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDMMMDMMMMMMMMMMMMMMMMMMMMMMMMMMM")
	{
		cout << "corrected! " << long_alignment << " : " << new_long_alignment << endl;
	}
	if (new_long_alignment.size() != long_alignment.size())
	{
		cout << "shit, off by something" << endl;
	}


*/
	if (matches_in_front == 0)
	{
		cout << new_long_alignment << endl;
	}

	return new_long_alignment;


	//short_seq = short_alignment.substr(matches_in_front, short_alignment.size() - matches_in_front - matches_in_back);
	//del_ins_cigar = generate_cigar_string(string long_alignment)

}

int non_match_counter;
int count_non_matches_in_cigar_string(string cigar){
	non_match_counter = 0;
	for (int i = 0; i < cigar.size(); i++)
	{
		if (cigar.substr(i,1) == "I" || cigar.substr(i,1) == "D" || cigar.substr(i,1) == "N")
		{
			non_match_counter++;
		}
	}
	return non_match_counter;
}

string empty_string = "";
string generate_string_of_Ms(int len)
{
	empty_string = "";
	for (int i = 0; i < len; i++)
	{
		empty_string = empty_string + "M";
	}
	return empty_string;
}

string temp_IN;
string generate_IsNs_seq(string long_align, string target_seq)
{
	temp_IN = "";
	for (int i = 0; i < long_align.size(); i++)
	{
		if(long_align.substr(i,1) == "I" || long_align.substr(i,1) == "N")
		{
			temp_IN = temp_IN + target_seq.substr(i,1);
		}
	}
	return temp_IN;
}


int main(int argc, char* argv[])
{
	int max_dist_from_cut_site = stoi(argv[5]);
	int dist_of_cut_site_from_start = 34;
	int dist_of_cut_site_from_end = 34;

	vector<cell_target_combo> cell_target_combo_vect;
	map<string,int> cell_target_pos_in_vect_map;

	cell_target_combo temp_c_t_combo_struct;

	ifstream alignment_file;
	alignment_file.open(argv[1]);
	string line;
	string curr_ref_seq = "";
	string curr_target_seq = "";
	string temp_seq_name;
	int pos_of_p;
	int umi_count;
	int total_umi_count;
	bool inputing_ref_seq = false;
	bool inputing_real_seq = false;

//open fasta alignment file; for each cell-target combo, save target seq, ref seq, seq_name (PCR_RT_targ), total number of UMIs, total UMI count (PCR amplicons).
	while(getline(alignment_file, line))
	{
		if (line.find("ref") != line.npos) //if reference sequence
		{
			temp_c_t_combo_struct.ref_seq = curr_ref_seq.substr(0, curr_ref_seq.size() - 20);
			temp_c_t_combo_struct.target_seq = curr_target_seq.substr(0, curr_target_seq.size() - 20); //I think the 20bp at the end are a primer binding site??
			temp_c_t_combo_struct.PCR_RT_Targ = temp_seq_name;
			temp_c_t_combo_struct.total_umi_count = total_umi_count;
			temp_c_t_combo_struct.umis = umi_count;
			if (curr_ref_seq != "")
			{
				cell_target_combo_vect.push_back(temp_c_t_combo_struct);
			}
			curr_ref_seq = "";
			curr_target_seq = "";
			temp_seq_name = "";
			umi_count = 0;
			temp_c_t_combo_struct = cell_target_combo();
			inputing_real_seq = false;
			inputing_ref_seq = true;

		}
		else if(line.substr(0,1) == ">") //not reference, but sequence name
		{
			//temp_seq_name = line.substr(1,line.size());
			//cout << line << endl;
			pos_of_p = line.find("p");
			temp_seq_name = line.substr(pos_of_p);
			umi_count = stoi(line.substr(line.find("_")+1,pos_of_p-(line.find("_")+2)));
			total_umi_count = stoi(line.substr(1,temp_seq_name.find("_") - 1));
			//cout << umi_count << endl;
			inputing_ref_seq = false;
			inputing_real_seq = true;
		}
		else if (inputing_ref_seq == true)
		{
			curr_ref_seq = curr_ref_seq + line;
		}
		else if (inputing_real_seq == true)
		{
			curr_target_seq = curr_target_seq + line;
		}
	}

// put last sequence into the vector
	temp_c_t_combo_struct.ref_seq = curr_ref_seq.substr(0, curr_ref_seq.size() - 20);
	temp_c_t_combo_struct.target_seq = curr_target_seq.substr(0, curr_target_seq.size() - 20);
	temp_c_t_combo_struct.PCR_RT_Targ = temp_seq_name;
	temp_c_t_combo_struct.total_umi_count = total_umi_count;
	temp_c_t_combo_struct.umis = umi_count;
	curr_ref_seq = "";
	curr_target_seq = "";
	temp_seq_name = "";
	umi_count = 0;
	cell_target_combo_vect.push_back(temp_c_t_combo_struct);


//for each cell-target combo: (1) convert alignment to MDIN code; (2) Generate "cigar" string from this code; (3) generate alignment which includes a single "I" for all length insertions (useful for when we'd like to plot all editing patterns; (4) Count insertion length
//Determine if edited: Some mismatches/dels should not be considered edits if they are not close enough to the cut site. If no edits exist within MAXDIST (often, 4) of cut site, then consider site unedited.
//If target is edited, remove any indels/mismatches which are discontinuous & not within 4 bases of "continuous" edit (function: correct_long_alignment_string); if indels/mismatches beyond 4 bases of last "continuous" edit, change them to matches.
	for(int i = 0; i < cell_target_combo_vect.size(); i++)
	{
		cell_target_combo_vect[i].long_alignment = make_alignment_string(cell_target_combo_vect[i].ref_seq, cell_target_combo_vect[i].target_seq);
		cell_target_combo_vect[i].cigar_string = generate_cigar_string(cell_target_combo_vect[i].long_alignment);
		cell_target_combo_vect[i].long_aligment_short_insertions = convert_alignment_to_short_insertions(cell_target_combo_vect[i].long_alignment);
		cell_target_combo_vect[i].num_Is = count(cell_target_combo_vect[i].long_alignment.begin(), cell_target_combo_vect[i].long_alignment.end(), 'I');
		if(cell_target_combo_vect[i].cigar_string == "70M")
		{
			cell_target_combo_vect[i].edited = false;
		}
		else
		{
			cell_target_combo_vect[i].edited = determine_if_edited_from_short_alignment(cell_target_combo_vect[i].long_aligment_short_insertions, dist_of_cut_site_from_start, max_dist_from_cut_site);
			if (cell_target_combo_vect[i].edited == false)
			{
				cell_target_combo_vect[i].cigar_string = "70M";
				cell_target_combo_vect[i].long_alignment = generate_string_of_Ms(70);
				cell_target_combo_vect[i].long_aligment_short_insertions = cell_target_combo_vect[i].long_alignment;
				cell_target_combo_vect[i].num_Is = 0;
			}
			if (cell_target_combo_vect[i].edited == true)
			{
				if(count_non_matches_in_cigar_string(cell_target_combo_vect[i].cigar_string) > 1)
				{
					cell_target_combo_vect[i].long_alignment = correct_long_alignment_string(cell_target_combo_vect[i].long_aligment_short_insertions, cell_target_combo_vect[i].long_alignment, dist_of_cut_site_from_start, max_dist_from_cut_site, 4);
					cell_target_combo_vect[i].cigar_string = generate_cigar_string(cell_target_combo_vect[i].long_alignment);
					cell_target_combo_vect[i].long_aligment_short_insertions = convert_alignment_to_short_insertions(cell_target_combo_vect[i].long_alignment);
					cell_target_combo_vect[i].num_Is = count(cell_target_combo_vect[i].long_alignment.begin(), cell_target_combo_vect[i].long_alignment.end(), 'I');
				}
			}
		}

		//if (cell_target_combo_vect[i].cigar_string != "70M" /*&& cell_target_combo_vect[i].edited == false*/)
		//if (cell_target_combo_vect[i].PCR_RT_Targ.find("H9_CTACGACGAG") != cell_target_combo_vect[i].PCR_RT_Targ.npos)
		//{
			/*
			cout << cell_target_combo_vect[i].PCR_RT_Targ << endl;
			cout << cell_target_combo_vect[i].ref_seq << endl;
			cout << cell_target_combo_vect[i].target_seq << endl;
			cout << cell_target_combo_vect[i].long_alignment << endl;
			cout << cell_target_combo_vect[i].cigar_string << endl;
			cout << "Matches_at_start=" << cell_target_combo_vect[i].num_matches_at_start << endl;
			cout << "Matches_at_end=" << cell_target_combo_vect[i].num_matches_at_end << endl;

			//cout << cell_target_combo_vect[i].cigar_string << endl;
			cout << cell_target_combo_vect[i].long_aligment_short_insertions << endl;
			*/
		//}



		//cout << cell_target_combo_vect[i].cigar_string << endl;
		//if (cell_target_combo_vect[i].long_alignment.find("N") != cell_target_combo_vect[i].long_alignment.npos)
		//{
		//	cout << cell_target_combo_vect[i].ref_seq << endl;
		//	cout << cell_target_combo_vect[i].target_seq << endl;
		//	cout << cell_target_combo_vect[i].long_alignment << endl;
		//	cout << cell_target_combo_vect[i].cigar_string << endl;
		//}
	}
	//cout << "Here we are!" << endl;

	//Make a cigar string which includes insertion information
	for(int i = 0; i < cell_target_combo_vect.size(); i++)
	{
		cell_target_combo_vect[i].IsNs = generate_IsNs_seq(cell_target_combo_vect[i].long_alignment, cell_target_combo_vect[i].target_seq);
		cell_target_combo_vect[i].cigar_plus_IsNs = cell_target_combo_vect[i].cigar_string + ":" + cell_target_combo_vect[i].IsNs;
	}

	//NEW STUFF
	//make a map of cell-target combos & their positions in the vector;
	//make a vector which contains the cell-target combos which appear more than once (these are the targets which contain more than one editing pattern in the dataset)
	map<string,vector<int> > cell_target_combo_positions_map;
	vector <string> combos_which_appear_more_than_once;

	for(int i = 0; i < cell_target_combo_vect.size(); i++)
	{
		if (cell_target_combo_positions_map.find(cell_target_combo_vect[i].PCR_RT_Targ) == cell_target_combo_positions_map.end())
		{
			cell_target_combo_positions_map[cell_target_combo_vect[i].PCR_RT_Targ].push_back(i);
		}
		else
		{
			cell_target_combo_positions_map[cell_target_combo_vect[i].PCR_RT_Targ].push_back(i);
			combos_which_appear_more_than_once.push_back(cell_target_combo_vect[i].PCR_RT_Targ);
		}
	}

	//Above, we corrected discontinuous editing patterns which are probably not due to CRISPR processes. This means that there are potentially cell_target combos which exist mutliple times in the vector but now have the same editing pattern. We'd like to combine these into one. The code below addresses this.
	
	string temp_pcr_rt_targ;
	string temp_cigar_string;
	int temp_num_umis;
	int temp_total_umis;
	int temp_pos_in_vect;
	int pos_in_vect_of_matching_seq;
	map<string,int> temp_cigar_to_position_map;
	for (int i = 0; i < combos_which_appear_more_than_once.size(); i++)
	{
		temp_cigar_to_position_map.clear();
		temp_pcr_rt_targ = combos_which_appear_more_than_once[i];
		//
		for (int j = 0; j < cell_target_combo_positions_map[temp_pcr_rt_targ].size(); j++)
		{
			temp_pos_in_vect = cell_target_combo_positions_map[temp_pcr_rt_targ][j];
			temp_cigar_string = cell_target_combo_vect[temp_pos_in_vect].cigar_plus_IsNs;
			temp_num_umis = cell_target_combo_vect[temp_pos_in_vect].umis;
			temp_total_umis = cell_target_combo_vect[temp_pos_in_vect].total_umi_count;

			if (temp_cigar_to_position_map.find(temp_cigar_string) == temp_cigar_to_position_map.end())
			{
				temp_cigar_to_position_map[temp_cigar_string] = temp_pos_in_vect;
				//temp_cigar_to_position_and_count_map[temp_cigar_string].push_back(temp_num_umis);
			}
			else
			{
				pos_in_vect_of_matching_seq = temp_cigar_to_position_map[temp_cigar_string];
				cell_target_combo_vect[pos_in_vect_of_matching_seq].umis+=temp_num_umis;
				cell_target_combo_vect[pos_in_vect_of_matching_seq].total_umi_count+=temp_total_umis;
				cell_target_combo_vect[temp_pos_in_vect].umis = 0;
				cell_target_combo_vect[temp_pos_in_vect].print_this = false;
			}
		}
	}


	map<string, int> short_align_count_map;

	string temp_short_alignment_plus_Is;
	for (int i = 0; i < cell_target_combo_vect.size(); i++)
	{
		temp_short_alignment_plus_Is = cell_target_combo_vect[i].long_aligment_short_insertions + to_string(cell_target_combo_vect[i].num_Is);
		if (short_align_count_map.find(temp_short_alignment_plus_Is) == short_align_count_map.end())
		{
			short_align_count_map[temp_short_alignment_plus_Is] = 1;
		}
		else
		{
			short_align_count_map[temp_short_alignment_plus_Is]++;
		}
	}

	map<string, int> edited_only_short_align_count_map;
	map<string, int> unedited_only_short_align_count_map;
	for (int i = 0; i < cell_target_combo_vect.size(); i++)
	{
		if(cell_target_combo_vect[i].edited == true)
		{
			temp_short_alignment_plus_Is = cell_target_combo_vect[i].long_aligment_short_insertions + to_string(cell_target_combo_vect[i].num_Is);
			if (edited_only_short_align_count_map.find(temp_short_alignment_plus_Is) == edited_only_short_align_count_map.end())
			{
				edited_only_short_align_count_map[temp_short_alignment_plus_Is] = 1;
			}
			else
			{
				edited_only_short_align_count_map[temp_short_alignment_plus_Is]++;
			}
		}
		else if (cell_target_combo_vect[i].edited == false)
		{
			temp_short_alignment_plus_Is = cell_target_combo_vect[i].long_aligment_short_insertions + to_string(cell_target_combo_vect[i].num_Is);
			if(unedited_only_short_align_count_map.find(temp_short_alignment_plus_Is) == unedited_only_short_align_count_map.end())
			{
				unedited_only_short_align_count_map[temp_short_alignment_plus_Is] = 1;
			}
			else
			{
				unedited_only_short_align_count_map[temp_short_alignment_plus_Is]++;
			}
		}
	}


	ofstream matrix_outfile;
	matrix_outfile.open(argv[2]);

	string sequence;
	int count;
	for (auto it = short_align_count_map.begin(); it != short_align_count_map.end(); it++)
	{
		sequence = it -> first;
		count = it -> second;
		matrix_outfile << count << "\t" << sequence.substr(70);
		for (int i = 0; i < 70; i++)
		{
			matrix_outfile << "\t" << sequence.substr(i,1);
		}
		matrix_outfile << endl;
	}


	ofstream edited_matrix_outfile;
	edited_matrix_outfile.open(argv[3]);
	for (auto it = edited_only_short_align_count_map.begin(); it != edited_only_short_align_count_map.end(); it++)
	{
		sequence = it -> first;
		count = it -> second;
		edited_matrix_outfile << count << "\t" << sequence.substr(70);
		for (int i = 0; i < 70; i++)
		{
			edited_matrix_outfile << "\t" << sequence.substr(i,1);
		}
		edited_matrix_outfile << endl;
	}

	ofstream unedited_matrix_outfile;
	unedited_matrix_outfile.open(argv[4]);

	for (auto it = unedited_only_short_align_count_map.begin(); it != unedited_only_short_align_count_map.end(); it++)
	{
		sequence = it -> first;
		count = it -> second;
		unedited_matrix_outfile << count << "\t" << sequence.substr(70);
		for (int i = 0; i < 70; i++)
		{
			unedited_matrix_outfile << "\t" << sequence.substr(i,1);
		}
		unedited_matrix_outfile << endl;
	}

	map<string,int> names_to_pos_in_vect_map_for_sorting;
	string temp_name_for_sorting;
	for (int i = 0; i < cell_target_combo_vect.size(); i++)
	{
		if(cell_target_combo_vect[i].umis != 0)
		{
			temp_name_for_sorting = cell_target_combo_vect[i].PCR_RT_Targ + cell_target_combo_vect[i].cigar_string + cell_target_combo_vect[i].IsNs;
			names_to_pos_in_vect_map_for_sorting[temp_name_for_sorting] = i;
		}
	}

	ofstream edited_unedited_outfile;
	edited_unedited_outfile.open(argv[6]);
	string temp_umi_and_name;
	string E_or_U;
	int temp_pos;
	for (auto it = names_to_pos_in_vect_map_for_sorting.begin(); it != names_to_pos_in_vect_map_for_sorting.end(); it++)
	{
		temp_pos = it -> second;
		temp_umi_and_name = "(" + to_string(cell_target_combo_vect[temp_pos].total_umi_count) + ")" + to_string(cell_target_combo_vect[temp_pos].umis) + "_" + cell_target_combo_vect[temp_pos].PCR_RT_Targ;
		if (cell_target_combo_vect[temp_pos].edited == true)
		{
			E_or_U = "E";
		}
		else
		{
			E_or_U = "U";
		}
		edited_unedited_outfile << temp_umi_and_name << "|" << E_or_U << "-" << cell_target_combo_vect[temp_pos].cigar_string << ":" << cell_target_combo_vect[temp_pos].IsNs << endl;
	}
/*
	for (int i = 0; i <  cell_target_combo_vect.size(); i++)
	{
		if (cell_target_combo_vect[i].umis != 0)
		{
			temp_umi_and_name = cell_target_combo_vect[i].PCR_RT_Targ.substr(cell_target_combo_vect[i].PCR_RT_Targ.find('_') + 1);
			if (cell_target_combo_vect[i].edited == true)
			{
				E_or_U = "E";
			}
			else
			{
				E_or_U = "U";
			}
			edited_unedited_outfile << temp_umi_and_name << "|" << E_or_U << "-" << cell_target_combo_vect[i].cigar_string << ":" << cell_target_combo_vect[i].IsNs << endl;
		}
	}
*/

	return 0;
}
