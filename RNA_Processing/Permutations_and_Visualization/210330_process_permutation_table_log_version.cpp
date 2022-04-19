/*
 * 210330_process_permutation_table_log_version.cpp
 *
 *  Created on: Mar 30, 2021
 *      Author: anna
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
#include <math.h>

using namespace std;

double calculate_mean(vector <double> &vect_of_permutations){
	double temp_sum = 0;
	for (int i = 1; i < vect_of_permutations.size(); i++) //this ignores real mean
	{
		temp_sum = temp_sum + vect_of_permutations[i];
	}
	return(temp_sum/(double(vect_of_permutations.size()-1)));
}

double calculate_st_dev(vector <double> &vect_of_permutations, double temp_mean){
	double temp_squared_sum = 0;
	double temp_diff_sq;
	for(int i = 1; i < vect_of_permutations.size(); i++){
		temp_diff_sq = pow(abs(vect_of_permutations[i]-temp_mean),2);
		temp_squared_sum = temp_squared_sum + temp_diff_sq;
	}
	return(sqrt(temp_squared_sum/double(vect_of_permutations.size()-1)));
}

double calculate_median(vector <double> &vect_of_permutations){
	int size_of_vect = vect_of_permutations.size();
	if(size_of_vect%2 == 0){
		return((vect_of_permutations[size_of_vect/2] + vect_of_permutations[size_of_vect/2 + 1])/2);
	}
	else {
		size_of_vect =size_of_vect + 1;
		return(vect_of_permutations[size_of_vect/2]);
	}
}


int main(int argc, char* argv[]){

	ifstream(input_perm_file);
	input_perm_file.open(argv[1]);

	ofstream(output_file);
	output_file.open(argv[2]);

	output_file << "Grp1_col_grpNum\tGrp2_col_grpNum\tGene\tGrp1_Mean\tGrp2_Mean\tGrp1_Frac_NonZero\tGrp2_Frac_NonZero\tLog_Real_quotient\tpos_of_real\tPerm_Mean_LogRatio\tPerm_St_Dev_LogRatio\tPerm_Median_LogRatio" << endl;

	string line;
	string elem;
	stringstream ss;

	vector<string> first_seven_col_vect;
	vector<double> permutations_vect;
	vector<int> vect_of_pos_of_real_value;
	int new_pos_of_real;

	int counter;
	double temp_real_ratio;
	int vect_size;

	double temp_mean;
	double temp_st_dev;
	double temp_median;

	int temp_counter = 1;
	int line_counter = 0;

	getline(input_perm_file, line);

	bool has_pos_col = false;
	size_t pos_of_pos_column = line.find("pos_of_real");
	if (pos_of_pos_column != string::npos){
		has_pos_col = true;
		cout << "has extra column" << endl;
	}


	while(getline(input_perm_file, line)){
		//cout << line_counter << endl;
		line_counter++;
		if (line_counter %10000 == 0){
			cout << line_counter << endl;
		}
		ss << line;
		counter = 0;
		first_seven_col_vect.clear();
		permutations_vect.clear();
		//cout << "Cleared vector size: " << permutations_vect.size() << endl;
		vect_of_pos_of_real_value.clear();
		while(getline(ss, elem, '\t')){
			//cout << counter << endl;
			if (counter < 8)
			{
				first_seven_col_vect.push_back(elem);
				counter++;
			}
			else{
				permutations_vect.push_back(log(stod(elem)));
			}
		}
		ss.clear();
		if (has_pos_col == true){
			permutations_vect.pop_back();
		}
		//cout << "perm_vect_size: " << permutations_vect.size() << endl;
		temp_real_ratio = permutations_vect[0];
		//cout << temp_real_ratio << endl;
		temp_mean = calculate_mean(permutations_vect);
		temp_st_dev = calculate_st_dev(permutations_vect, temp_mean);

		sort(permutations_vect.begin(), permutations_vect.end());
		temp_median = calculate_median(permutations_vect);

		if (temp_counter == 1){
			for(int i = 0; i < permutations_vect.size(); i++){
				//cout << permutations_vect[i] << endl;
			}
			temp_counter++;
		}
		for (int i = 0; i < permutations_vect.size(); i++){
			if(permutations_vect[i] == temp_real_ratio){
				vect_of_pos_of_real_value.push_back(i+1);
			}
			//cout << "here3??" << endl;
			//else if (permutations_vect[i] > temp_real_ratio){
			//	break;
			//}
		}
		//cout << "Vector size: " << vect_of_pos_of_real_value.size() << endl;
		if (vect_of_pos_of_real_value.size() == 1){
			new_pos_of_real = vect_of_pos_of_real_value[0];
		}
		else
		{
			//cout << "here??" << endl;
			vect_size = vect_of_pos_of_real_value.size();
			new_pos_of_real = vect_of_pos_of_real_value[ceil(vect_size/2)];
			//cout << "here2??" << endl;
		}

		for (int i = 1; i < first_seven_col_vect.size(); i++){
			output_file << first_seven_col_vect[i] << "\t";
		}
		output_file << temp_real_ratio << "\t" << new_pos_of_real << "\t" << temp_mean << "\t" << temp_st_dev << "\t" << temp_median << endl;

	}

	return 0;
}



