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
			return 1;
		}
	}
	for (int ii = 1; ii < argc; ii++){
	 	string arg = argv[ii];
		if (arg == "-t"){
			cout << "Started calculating QI tensor!" << endl;
			make_qi_tensor();
			cout << "QI tensor written!" << endl;
		}
		if (arg == "-b"){
		hamiltonian Ham("INPUT.txt");
		band_structure Band(Ham);
		Band.do_plot();
		cout << "Band structure written!" << endl;
		}
		if (arg != "-b" &&
			arg != "-h" &&
			arg != "-t"){
			throw runtime_error("Error: Unknown argument "+arg+"!");
		}
	}
	cout << "Done!" << endl;
	return 0;
}

