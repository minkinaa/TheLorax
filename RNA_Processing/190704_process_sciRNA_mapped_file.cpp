/*
 * 190704_process_sciRNA_mapped_file.cpp
 *
 *  Created on: Jul 4, 2019
 *      Author: anna
 *
 *      #read in list of all genes
 *
 *      argv1: all genes file
 *      argv2: main bed
 *      unmapped_intervals_outfile.open(argv[3]);
 *      multimapped_intervals_outfile.open(argv[4]);
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
	vector<int> gene_counts_vect;
};

string gene_tc;
string return_gene_to_count(string temp_exon_string, string temp_gene_string)
{
	if (temp_exon_string.find(",") == temp_exon_string.npos && temp_exon_string != ".")
	{
		gene_tc = temp_exon_string;
	}
	else if (temp_exon_string == ".")
	{
		if (temp_gene_string == ".")
		{
			gene_tc = "NOT_MAPPED";
			//cout << "got_not_mapped!!" << endl;
		}
		else if (temp_gene_string.find(",") == temp_gene_string.npos)
		{
			gene_tc = temp_gene_string;
		}
		else
		{
			gene_tc = "MAPS_TO_MULTIPLE_GENES";
		}
	}
	else if (temp_exon_string.find(",") != temp_exon_string.npos)
	{
		gene_tc = "MAPS_TO_MULTIPLE_EXONS";
	}

	return gene_tc;
}


int main(int argc, char* argv[])
{
	vector<cell> vector_of_cells;
	cell temp_cell;
	map<string,int> cell_id_to_vect_pos_map;

	map<string,int> gene_id_to_vect_pos_map;
	vector<string> all_genes_vect;

	map<string,int> gene_id_to_total_count;

	vector<string> cell_id_for_multimappers;
	vector<string> multimapper_string;

	map<string, int> chr_pos_to_count_map;

	ifstream(all_genes);
	all_genes.open(argv[1]);
	cout << "argv1=" << argv[1];
	string line;

	while(getline(all_genes, line))
	{
		cout << line << endl;
		all_genes_vect.push_back(line);
		gene_id_to_vect_pos_map[line] = all_genes_vect.size() - 1;
		gene_id_to_total_count[line] = 0;
	}
	all_genes.close();

	//fill temp cell
	for (int i = 0; i < all_genes_vect.size(); i++)
	{
		temp_cell.gene_counts_vect.push_back(0);
	}

	ofstream(unmapped_intervals_outfile);
	unmapped_intervals_outfile.open(argv[3]);

	ofstream(multimapped_intervals_outfile);
	multimapped_intervals_outfile.open(argv[4]);


	string temp_chr;
	string temp_start;
	string temp_end;
	string temp_cell_umi_id;
	string temp_plus_minus;
	string temp_exon_string;
	string temp_gene_string;
	string gene_to_count;

	string temp_cell_id;

	ifstream(main_bed);
	main_bed.open(argv[2]);
	stringstream ss;
	string elem;
	int counter = 0;

	int pos_in_vect_of_gene;
	int pos_of_cell_in_vector;
	bool not_mapped_or_multi_mapped = false;

	int line_counter = 0;
	string temp_chr_pos;

	while(getline(main_bed, line))
	{
		line_counter++;
		if (line_counter % 10000 == 0)
		{
			cout << line_counter << endl;
		}
		ss << line;
		//cout << line << endl;
		counter = 0;
		while(getline(ss, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_chr = elem;
				counter++;
			}
			else if (counter == 1)
			{
				temp_start = elem;
				counter++;
			}
			else if (counter == 2)
			{
				temp_end = elem;
				counter++;
			}
			else if (counter == 3)
			{
				temp_cell_umi_id = elem;
				counter++;
			}
			else if (counter == 4)
			{
				counter++;
			}
			else if (counter == 5)
			{
				temp_plus_minus = elem;
				counter++;
			}
			else if (counter == 6)
			{
				temp_exon_string = elem;
				counter++;
			}
			else if (counter == 7)
			{
				temp_gene_string = elem;
				counter++;
			}
			else
			{
				cout << "wait, this shouldn't happen" << endl;
			}
		}
		ss.clear();

		temp_cell_id = temp_cell_umi_id.substr(0,temp_cell_umi_id.find("_")) + "_" + temp_cell_umi_id.substr(temp_cell_umi_id.length() - 10, 10);
		if (cell_id_to_vect_pos_map.find(temp_cell_id) == cell_id_to_vect_pos_map.end())
		{
			temp_cell.cell_id = temp_cell_id;
			vector_of_cells.push_back(temp_cell);
			cout << vector_of_cells.size() << endl;
			temp_cell.cell_id = "";
			pos_of_cell_in_vector = vector_of_cells.size() - 1;
			cell_id_to_vect_pos_map[temp_cell_id] = pos_of_cell_in_vector;
			gene_to_count = return_gene_to_count(temp_exon_string, temp_gene_string);
			if (gene_to_count != "NOT_MAPPED" && gene_to_count != "MAPS_TO_MULTIPLE_GENES" &&  gene_to_count != "MAPS_TO_MULTIPLE_EXONS")
			{
				if (gene_id_to_vect_pos_map.find(gene_to_count) == gene_id_to_vect_pos_map.end())
				{
					cout << "Shoot, we are missing genes in our gene list!!!: " << gene_to_count << endl;
				}
				else
				{
					cout << "here" << endl;
					pos_in_vect_of_gene = gene_id_to_vect_pos_map[gene_to_count];
					vector_of_cells[pos_of_cell_in_vector].gene_counts_vect[pos_in_vect_of_gene]++;
					gene_id_to_total_count[gene_to_count]++;
				}
			}
			else
			{
				not_mapped_or_multi_mapped = true;
			}
		}
		else if (cell_id_to_vect_pos_map.find(temp_cell_id) != cell_id_to_vect_pos_map.end())
		{
			pos_of_cell_in_vector = cell_id_to_vect_pos_map[temp_cell_id];
			gene_to_count = return_gene_to_count(temp_exon_string, temp_gene_string);
			if (gene_to_count != "NOT_MAPPED" && gene_to_count != "MAPS_TO_MULTIPLE_GENES" &&  gene_to_count != "MAPS_TO_MULTIPLE_EXONS")
			{
				if (gene_id_to_vect_pos_map.find(gene_to_count) == gene_id_to_vect_pos_map.end())
				{
					cout << "Shoot, we are missing genes in our gene list!!!: " << gene_to_count << endl;
				}
				else
				{
					pos_in_vect_of_gene = gene_id_to_vect_pos_map[gene_to_count];
					vector_of_cells[pos_of_cell_in_vector].gene_counts_vect[pos_in_vect_of_gene]++;
					gene_id_to_total_count[gene_to_count]++;
				}
			}
			else
			{
				not_mapped_or_multi_mapped = true;
			}
		}
		if (not_mapped_or_multi_mapped == true)
		{
			if (gene_to_count == "NOT_MAPPED")
			{
				unmapped_intervals_outfile << temp_chr << "\t" << temp_start << "\t" << temp_end << "\t" << temp_cell_umi_id << "\t" << temp_plus_minus << "\t" << temp_exon_string << "\t" << temp_gene_string << endl;
			}
			else if (gene_to_count == "MAPS_TO_MULTIPLE_EXONS")
			{
				cell_id_for_multimappers.push_back(temp_cell_id);
				multimapper_string.push_back(temp_exon_string);
				multimapped_intervals_outfile << temp_chr << "\t" << temp_start << "\t" << temp_end << "\t" << temp_cell_umi_id << "\t" << temp_plus_minus << "\t" << temp_exon_string << "\t" << temp_gene_string << endl;

			}
			else if (gene_to_count == "MAPS_TO_MULTIPLE_GENES")
			{
				cell_id_for_multimappers.push_back(temp_cell_id);
				multimapper_string.push_back(temp_gene_string);
				multimapped_intervals_outfile << temp_chr << "\t" << temp_start << "\t" << temp_end << "\t" << temp_cell_umi_id << "\t" << temp_plus_minus << "\t" << temp_exon_string << "\t" << temp_gene_string << endl;
			}
			not_mapped_or_multi_mapped = false;
		}

		//fill map of insertion pos's;
		if (temp_plus_minus == "+")
		{
			temp_chr_pos = temp_chr + "_" + temp_start + "_" + temp_plus_minus;
		}
		else
		{
			temp_chr_pos = temp_chr + "_" + temp_end + "_" + temp_plus_minus;
		}

		if (chr_pos_to_count_map.find(temp_chr_pos) == chr_pos_to_count_map.end())
		{
			chr_pos_to_count_map[temp_chr_pos] = 1;
		}
		else
		{
			chr_pos_to_count_map[temp_chr_pos]++;
		}
	}

	cout << "done with initial pass, now dealing with multimappers" << endl;

	//deal with multimappers: this code gives the count to the gene which is most abundant across the dataset
	string string_of_genes;
	vector<string> multi_gene_temp_vector;
	string temp_gene;
	int temp_count;
	string temp_best_gene;
	int temp_highest_count = 0;

	int counter_of_unresolved_multimappers = 0;

	cout << "Number_of_multimappers" << cell_id_for_multimappers.size() << endl;

	for (int i = 0; i < cell_id_for_multimappers.size(); i++)
	{
		//cout << cell_id_for_multimappers[i] << endl;
		temp_highest_count = 0;
		temp_best_gene = "";
		string_of_genes = multimapper_string[i];
		//cout << string_of_genes << endl;
		ss << string_of_genes;
		while(getline(ss, elem, ','))
		{

			temp_gene = elem;
			//cout << elem << endl;
			temp_count = gene_id_to_total_count[temp_gene];
			if (temp_count > temp_highest_count)
			{
				temp_highest_count = temp_count;
				temp_best_gene = temp_gene;
			}
		}
		ss.clear();
		if (temp_highest_count > 0)
		{
			//cout << temp_highest_count << endl;
			pos_of_cell_in_vector = cell_id_to_vect_pos_map[cell_id_for_multimappers[i]];
			pos_in_vect_of_gene = gene_id_to_vect_pos_map[temp_best_gene];
			vector_of_cells[pos_of_cell_in_vector].gene_counts_vect[pos_in_vect_of_gene]++;
		}
		else
		{
			//cout << "this gene was zero" << endl;
			counter_of_unresolved_multimappers++;
		}
	}



	cout << "generating outfile!!" << endl;
	ofstream(cell_by_gene_outfile);
	cell_by_gene_outfile.open(argv[5]);

	int non_zero_counter = 0;

	for (int i = 0; i < vector_of_cells.size(); i++)
	{
		cell_by_gene_outfile <<  vector_of_cells[i].cell_id;
		for (int j = 0; j < vector_of_cells[i].gene_counts_vect.size(); j++)
		{
			cell_by_gene_outfile << "\t" << vector_of_cells[i].gene_counts_vect[j];
			if (vector_of_cells[i].gene_counts_vect[j] != 0)
			{

			}
		}
		cell_by_gene_outfile << endl;
	}

	//print map of counts per pos
	ofstream(counts_per_pos_outfile);
	counts_per_pos_outfile.open(argv[6]);

	for (auto it = chr_pos_to_count_map.begin(); it != chr_pos_to_count_map.end(); it++)
	{
		if (it -> second > 1)
		{
			counts_per_pos_outfile << it -> first << "\t" << it -> second << "\t" << endl;
		}
	}

/*
	//generate a matrix.mtx cell ranger file
	ofstream(matrix_cell_ranger);
	matrix_cell_ranger.open(argv[6]);

	matrix_cell_ranger << "%%MatrixMarket matrix coordinate real general" << endl << "%" << endl;
*/




	return 0;
}


