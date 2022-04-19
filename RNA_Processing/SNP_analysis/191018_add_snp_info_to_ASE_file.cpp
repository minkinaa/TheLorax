/*
 * 191018_add_snp_info_to_ASE_file.cpp
 *
 *  Created on: Oct 18, 2019
 *      Author: anna
 *
 *      argv1: snp file
 *      argv2: counts per base file
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

struct SNP{
	string chr_pos;
	string snp_id;
	string major_allele;
	string minor_allele;
};

int main(int argc, char* argv[])
{
	ifstream(snp_file);
	snp_file.open(argv[1]);

	SNP temp_SNP_struct;
	vector<SNP> SNP_vect;
	map<string, int> SNP_to_pos_in_vect_map;

	string line;
	stringstream temp;
	string elem;

	string temp_chr_pos_string;

	int line_counter = 0;
	int counter = 0;
	while(getline(snp_file, line))
	{
		line_counter++;
		if (line_counter %500000 == 0)
		{
			cout << line_counter << endl;
		}
		temp << line;
		counter = 0;
		temp_chr_pos_string = "";
		while(getline(temp, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_chr_pos_string = "chr" + elem + "_";
				counter++;
			}
			else if (counter == 1)
			{
				temp_chr_pos_string = temp_chr_pos_string + elem;
				temp_SNP_struct.chr_pos = temp_chr_pos_string;
				//cout << temp_chr_pos_string << endl;
				counter++;
			}
			else if (counter == 2)
			{
				temp_SNP_struct.snp_id = elem;
				counter++;
			}
			else if (counter == 3)
			{
				temp_SNP_struct.major_allele = elem;
				counter++;
			}
			else if (counter == 4)
			{
				temp_SNP_struct.minor_allele = elem;
				counter++;
			}
			else
			{
				cout << "Ooooops, this shouldn't happen!!" << endl;
			}
		}
		SNP_vect.push_back(temp_SNP_struct);
		SNP_to_pos_in_vect_map[SNP_vect[SNP_vect.size()-1].chr_pos] = SNP_vect.size()-1;
		temp.clear();
	}

	ifstream(base_count_file);
	base_count_file.open(argv[2]);

	cout << "ok, writing file" << endl;
	ofstream(outfile);
	outfile.open(argv[3]);

	string temp_snp;
	string temp_major;
	string temp_minor;

	getline(base_count_file, line);
	outfile << line << "\t" << "SNP" << "\t" << "major" << "\t" << "minor" << endl;


	line_counter = 0;
	while(getline(base_count_file, line))
	{
		if (line_counter %100000 == 0)
		{
			cout << line_counter << endl;
		}
		line_counter++;
		temp << line;
		counter = 0;
		temp_chr_pos_string = "";
		while(getline(temp, elem, '\t'))
		{
			if (counter == 0)
			{
				temp_chr_pos_string = elem;
				counter++;
			}
			else if (counter == 1)
			{
				temp_chr_pos_string = temp_chr_pos_string + "_" + elem;
				//cout << temp_chr_pos_string << endl;
				if (SNP_to_pos_in_vect_map.find(temp_chr_pos_string) != SNP_to_pos_in_vect_map.end())
				{
					temp_snp = SNP_vect[SNP_to_pos_in_vect_map[temp_chr_pos_string]].snp_id;
					temp_major = SNP_vect[SNP_to_pos_in_vect_map[temp_chr_pos_string]].major_allele;
					temp_minor = SNP_vect[SNP_to_pos_in_vect_map[temp_chr_pos_string]].minor_allele;
					outfile << line << "\t" << temp_snp << "\t" << temp_major << "\t" << temp_minor << endl;
				}
				else
				{
					//outfile << line << "\t" << "none" << "\t" << "none" << "\t" << "none" << endl;
				}
				counter++;
			}
		}
		temp.clear();
	}



	return 0;
}





