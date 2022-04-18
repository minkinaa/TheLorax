/*
 * 191203_CROPt_make_UMI_RT_BC_seq_output_file.cpp
 *
 *  Created on: Dec 3, 2019
 *      Author: anna
 *
 *      Purpose: Take in file of sorted sequence from each RT well, output file of structure:
 *      UMI_RTindex_GTbarcode_trimmed_sequence
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

int main(int argc, char* argv[])
{
	ifstream(sorted_seq_file);
	sorted_seq_file.open(argv[1]);
	string line;
	string temp_UMI;
	string temp_RT_index;
	int pos_of_GT_barcode;
	string temp_GTbarcode;
	string trimmed_seq;
	int check_seq_not_near_end;
	int whole_string_length;

	ofstream(parsed_output);
	parsed_output.open(argv[2]);

	int counter = 0;

	while(getline(sorted_seq_file, line))
	{
		counter++;
		if (counter % 1000 == 0){
			cout << counter << endl;
		}
		temp_UMI = line.substr(0,8);
		temp_RT_index = line.substr(8,10);
		whole_string_length = line.size();
		check_seq_not_near_end = line.find("ATTGGTCTTA");
		if (check_seq_not_near_end != line.npos && whole_string_length - check_seq_not_near_end > 22)
		{
			pos_of_GT_barcode = line.find("ATTGGTCTTA") + 10;
            temp_GTbarcode = line.substr(pos_of_GT_barcode, 10);
            trimmed_seq = line.substr(pos_of_GT_barcode + 10);
            parsed_output << temp_UMI << "_" << temp_RT_index << "_" << temp_GTbarcode << "_" <<  trimmed_seq << endl;
		}
		else if (line.find("CGTAAGGTCT") != line.npos)
		{
			check_seq_not_near_end = line.find("CGTAAGGTCT");
			if (whole_string_length - (check_seq_not_near_end - 20) > 22 && check_seq_not_near_end > 20)

			{
				pos_of_GT_barcode = line.find("CGTAAGGTCT") - 10;
				temp_GTbarcode = line.substr(pos_of_GT_barcode, 10);
				trimmed_seq = line.substr(pos_of_GT_barcode + 10);
				parsed_output << temp_UMI << "_" << temp_RT_index << "_" << temp_GTbarcode << "_" <<  trimmed_seq << endl;
			}
		}
		else
		{
			cout << "ATTGGTCTTA not in the right place!" << endl;
			cout << line << endl;
		}
	}

	return 0;
}
