#include <iostream>
#include <iomanip>
#include <string> 


using namespace std;

#include "libs/Simplex.h"


float M[8] = {1, 1, 1, 7,
			   -2, -5, 0, 0};
int vb[1] = {3};




int main() {
    SimplexTable table;
    float ans[4];
    
    table.InitParameters(3, 1, M, vb);
    
    table.PrintTable();
    cout << endl;
    table.SolveSimplex(ans);
    
    cout << "Sol. : ";
    for (int j = 0; j < 3; j++) {
	    cout << ans[j] << " ";
	}
	cout << endl << "com minimo: " << ans[3] << endl;
    
    return 0;
}
