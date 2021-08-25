#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;

void dissimilarity_func(double* A1, double* A2, int n, int lambda, double* over, int l);

int main()
{
	using namespace std;
	
    std::ofstream out1;
    std::ofstream out2;
    std::ifstream inp;
    std::ifstream inp_fi;
    out1.open("D.dat");
    out2.open("D_k.dat");
    
    char buf[100], filename[100], symbol;
    int n_spins, n_measurements, n;
    double dissimilarity;
    int lambda, lambda0, diss_steps, diss_start;
    
    const char *keys[6] = { "n_spins", "n_measurements", "dissimilarity_steps",
        "starting_step", "lambda_0", "filename"};
    int key, str_comp;
    
    inp.open("inp.dat");
    for (int i = 0; i < 6; i++){
        inp >> buf;
        for (int j = 0; j < 6; j++){
            str_comp = strcmp(buf,keys[j]);
            if (str_comp == 0){
                key = j+1;
                break;
            }
        }
        
        switch (key)
        {
            case 1:
                inp >> n_spins;
                break;
            case 2:
                inp >> n_measurements;
                break;
            case 3:
                inp >> diss_steps;
                break;
            case 4:
                inp >> diss_start;
                break;
            case 5:
                inp >> lambda0;
                break;
            case 6:
                inp >> filename;
                break;
            default:
                cout << "Error: wrong key!" << endl;
                exit(0);
        }
    }
    inp.close();
    
    n = n_spins * n_measurements;
    
    int check = pow(lambda0, diss_steps);
    if ((n % check) != 0){
        cout << "Warning! The linear dimension of the input sequence must be divisible by k = lambda_0^{diss_steps}" << endl;
    }
    
    double *A1 = NULL;
    double *A2 = NULL;
    double *overA = NULL;

    A1 = (double*) malloc(n * sizeof(double));
    A2 = (double*) malloc(n * sizeof(double));
    overA = (double*) malloc(diss_steps * sizeof(double));
    
    srand(time(NULL));
    unsigned long long seed;
    seed=rand();
    
    inp_fi.open(filename);
        
    for (int i = 0; i < n; i++){
        inp_fi >> symbol;
        A1[i] = (double) symbol - '0';
        if (A1[i] < 0.5) {A1[i] = -1;}
    }
        
    inp_fi.close();
        
    for (int i = 0; i < diss_steps; i++){
        lambda = pow(lambda0, i + 1);
        if (i % 2 == 0){
            dissimilarity_func(A1,  A2,  n, lambda, overA, i);
        }
        else{
            dissimilarity_func(A2,  A1,  n, lambda, overA, i);
        }
    }
    
    dissimilarity = 0;
    for (int diss = diss_start; diss < diss_steps; diss++){
        dissimilarity += fabs(overA[diss]);
        out2 << diss << ' ' << fabs(overA[diss]) << endl;
    }
    out1 << dissimilarity << endl;
   
    out1.close();
    out2.close();
}

void dissimilarity_func(double* A1, double* A2, int n, int lambda, double* over, int l){
    
    double overlap2 = 0, overlap1 = 0, overlap3 = 0, ren_r, o1 = 0, o2 = 0, o3 = 0;
    double renorm;
    int n_bl = n/lambda;
    
    for (int n1_counti = 0; n1_counti < n_bl; n1_counti++){
        renorm = 0;
        for (int count_i=0; count_i < lambda; count_i++){
                renorm += A1[n1_counti * lambda + count_i];
        }
        ren_r = lambda;
        for (int count_i=0; count_i < lambda; count_i++){
                A2[n1_counti * lambda + count_i] = renorm/ren_r;
        }
    }
    
    for (int i = 0; i < n; i++){
        overlap1 += A1[i] * A1[i];
        overlap2 += A1[i] * A2[i];
        overlap3 += A2[i] * A2[i];
    }
    
    o1 = overlap1/n;
    o2 = overlap2/n;
    o3 = overlap3/n;
    over[l]= o2 - 0.5 * (o1 +  o3);
}
