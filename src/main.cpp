#include <iostream>
#include "Hamiltonian_functions_real.h"
#include "Optical_Matrix.h"
#include "KMesh.h"
#include "Sum.h"
#include "BandStructure.h"
#include <complex>
#include <memory>
#include <vector>
#include <array>
#include<utility>
#include <string>
using namespace std;
int main(int argc, char *argv[]){
	if (argc == 1){
		throw runtime_error("Error: Too few input arguments!");
	}	
	 cout << "\n";
	 cout << "             	  /\\     /\\ " << endl;
      	 cout << "		 /  \\~~~/  \\ " << endl;
     	 cout << "            	(   o   o   )" <<endl;
         cout << "	  	/     ^      \\ "  << endl;
   	 cout << "               /|             |\\ " << endl;
         cout << "              / |    _____     | \\ " << endl;
	 cout << "     	     /__|___|     |____|__\\ " << endl;
	 cout << "	    |______________________|" << endl;
	 cout << "| |    | |/ /  \\/  |/ __ \\|  __ \\|  ____| |" << endl;
 	 cout << "| |    | ' /| \\  / | |  | | |  | | |__  | |" << endl;
 	 cout << "| |    |  < | |\\/| | |  | | |  | |  __| | |" << endl;
         cout << "| |____| . \\| |  | | |__| | |__| | |____| |____ " << endl;
         cout << "|______|_|\\_\\_|  |_|\\____/|_____/|______|______|" << endl;
	 cout << "By: Ondřej Faiman " << endl;
	 for (int ii = 1; ii < argc; ii++){
	 	string arg = argv[ii];
		if (arg == "-h"){
			cout << endl;
			cout << "Help:" << endl;
			cout << "\tInclude all necessary information in the INPUT.txt file." << endl;
			cout << "\tRun 'export OMP_NUM_THREADS=N' to run the QI computation on N threads."<<endl;
			cout << "\tOptions:" << endl;
			cout << "\t -b: Write band structure to the BAND_LKMODEL.dat file. The kpath is given in the input file." << endl;
			cout << "\t -t: Write elements of the QI tensor to the QI_LKMODEL.dat file for the given range of frequencies." << endl;
			cout << "\t -d: Write also the decompostion of the QI tensor into band contributions to the QI_LKMODEL_DECOMP_.dat file for the given range of frequencies." << endl;
			return 1;
		}
	}
	int is_computed = 0;
	for (int ii = 1; ii < argc; ii++){
	 	string arg = argv[ii];
		if (arg == "-t"){
			if (is_computed == 1){
				cout << "Warning: -d option already used, skipping the -t option!\n";
		}
			else{
				is_computed++;
				cout << "Started calculating QI tensor!" << endl;
				make_qi_tensor(0);
				cout << "QI tensor written!" << endl;
			}
		}
		if (arg == "-d"){
			if (is_computed == 1){
				cout << "Warning: -t option already used, skipping the -d option!\n";
			}
			else{
				is_computed++;
				cout << "Started calculating QI tensor with decomposition!" << endl;
				make_qi_tensor(1);
				cout << "QI tensor written!" << endl;
			}
		}
		if (arg == "-b"){
		hamiltonian Ham("INPUT.txt");
		band_structure Band(Ham);
		Band.do_plot();
		cout << "Band structure written!" << endl;
		}
		if (arg != "-b" &&
			arg != "-h" &&
			arg != "-d" &&
			arg != "-t"){
			throw runtime_error("Error: Unknown argument "+arg+"!");
		}
	}
	cout << "Done!" << endl;
	return 0;
}

