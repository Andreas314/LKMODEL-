#include <iostream>
#include "Hamiltonian_functions.h"
#include <complex>
#include <memory>
#include<utility>
using namespace std;
int main(){
 	hamiltonian Ham("INPUT.txt");
	Ham.diagonal_at_k_point(.25, .5, 1.0);
 	hamiltonian Ham_2(move(Ham));
	Ham_2.diagonal_at_k_point(.25, .5, 1.0);
	cout << Ham.get_value(0);
}

