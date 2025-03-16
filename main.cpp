#include <iostream>
#include "Hamiltonian_functions_real.h"
#include "BandStructure.h"
#include <complex>
#include <memory>
#include<utility>
using namespace std;
int main(){
 	hamiltonian Ham("INPUT.txt");
	band_structure Band(Ham);
	Band.do_plot();
}

