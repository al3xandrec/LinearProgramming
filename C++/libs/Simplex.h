


#define XAN_SIMPLEX_H


const char separator    = ' ';
const int numWidth      = 15;

template<typename T> void printElement(T t) {
    cout << left << setw(numWidth) << setfill(separator) << t;
}



class SimplexTable {
    private:
        int restrictions, variables, *vb;
        float *table;
        
    public:
        void InitParameters(int var, int rest, float *M, int *VB) {
            variables = var;
            restrictions = rest;
            table = (float*)calloc((var+1)*(rest+1), sizeof(float));
            vb = (int*)calloc(rest, sizeof(int));
            
            for (int i = 0; i <= rest; i++) {
                for (int j = 0; j <= var; j++) {
                    table[j + (var+1)*i] = M[j + (var+1)*i];
                }
                vb[i] = VB[i];
            }
            
        }
        
        void PrintTable(void) {
        	// Headers
            for(int j = -2; j <= variables; j++) {
                if (j == -2) printElement("V.B.");
                else if (j == -1) printElement("Variaveis");
                else if (j < variables) {
                    int J = j+1;
                    string s = "x" + to_string(J);
                    printElement(s);
                }
                else if (j == variables) printElement("V.A.");
            }
            cout << endl;
            
            // Restrictions
            for(int i = 0; i < restrictions; i++) {
                for(int j = -2; j <= variables; j++) {
                    if (j == -2) printElement(vb[i]);
                    else if (j == -1) {
                        int I = i+1;
                        string s = "Cond." + to_string(I);
                        printElement(s);
                    }
                    else if (j <= variables) printElement(table[j + (variables+1)*i]);
                }
                cout << endl;
            }
            
            // Objective Funcion
            for(int j = -2; j <= variables; j++) {
                if (j == -2) printElement(" ");
                else if (j == -1) printElement("F.O.");
                else if (j < variables) printElement(table[j + (variables+1)*restrictions]);
                else if (j == variables) {
                    string s = "z + (" + to_string(table[j + (variables+1)*restrictions]) + ")";
                    printElement(s);
                }
            }
            cout << endl;
        }
        
        // Line operations and basic variables
        void SwapBV(int cond, int newBV) {
            vb[cond-1] = newBV;
        }
        
        void JacobiOP1(int mLine, float mValue) {
        	mLine--;
            for (int j = 0; j <= variables; j++) {
            	table[j + (variables+1)*mLine] = mValue*table[j + (variables+1)*mLine];
            }
        }
        
        void JacobiOP2(int tLine, int mLine, float mValue) {
        	tLine--;
        	mLine--;
            for (int j = 0; j <= variables; j++) {
            	table[j + (variables+1)*tLine] = table[j + (variables+1)*tLine] + mValue*table[j + (variables+1)*mLine];
            }
        }
        
        void SolveSimplex(float *ans) {
        	bool keepIterating = true;
        	
        	// Init Simplex
        	for (int j = 0; j < variables; j++) {
        		ans[j] = 0;
        	}
        	
        	for (int i = 0; i < restrictions; i++) {
        		ans[vb[i]-1] = table[variables + (variables+1)*i];
        	}
        	
        	while (keepIterating) {
        		
        		// Is there a negative element with the coeficients of the Objective Function?
        		int newBV = -1;
        		bool checkFO = false;
        		for(int j = 0; j < variables; j++) {
        			if (table[j + (variables+1)*restrictions] < 0 && !checkFO) {
        				checkFO = true;
        				newBV = j+1;
        			}
        		}
        		
        		keepIterating = checkFO;
        		if (checkFO) {
        			// Check which line the basic variable will be replaces
        			int lineForNewBV = -1;
        			float minRatio = -1;
        			for (int i = 0; i < restrictions; i++) {
        				float ratio;
						if (table[(newBV-1) + (variables+1)*i] != 0) {
							ratio = table[variables + (variables+1)*i]/table[(newBV-1) + (variables+1)*i];
						}
						else ratio = -1;
						
        				if ( ratio > 0 && ( minRatio < 0 || ratio < minRatio ) ) {
        					lineForNewBV = i;
        					minRatio = ratio;
        				}
        				
        				cout << "Cond." << i+1 << ": " << ratio << endl;
        			}
        			
        			cout << "Var: x" << newBV << " || Cond." << lineForNewBV+1 << endl << endl;
        			
        			// Line operations
        			SwapBV(lineForNewBV+1, newBV);
        			float inverseCoef = 1/table[(newBV-1) + (variables+1)*lineForNewBV];
        			JacobiOP1(lineForNewBV+1, inverseCoef);
        			
        			for (int i = 0; i <= restrictions; i++) {
        				if (i != lineForNewBV) {
        					JacobiOP2(i+1, lineForNewBV+1, (-1)*table[(newBV-1) + (variables+1)*i]);
        				}
        			}
        			
        			PrintTable();
        			cout << endl;
        		}
        		
        		// Store the new values in vector ans
        		for (int j = 0; j < variables; j++) {
	        		ans[j] = 0;
	        	}
	        	
	        	for (int i = 0; i < restrictions; i++) {
	        		ans[vb[i]-1] = table[variables + (variables+1)*i];
	        	}
	        	
	        	ans[variables] = (-1)*table[variables + (variables+1)*restrictions];
        	}
        }
};

