#include <TROOT.h>
#include <iostream>

void FastSim(int numero, unsigned int seed, bool distr_z, bool distr_m, int m);
void Residui();
void Risoluzione();
void Efficienza();

void All() {
    gROOT->ProcessLine(".L FastSimulation.C");
    gROOT->ProcessLine(".L Residui.C");
    gROOT->ProcessLine(".L Risoluzione.C");
    gROOT->ProcessLine(".L Efficienza.C");

    int n,m;
    unsigned int seed;
    bool distr_z,distr_m;
    cout<<"Inserire numero scontri:     ";
    cin>>n;
    cout<<"Inserire seed TRandom:     ";
    cin>>seed;
    cout<<"Inserire distribuzione in z [1 = prof, 0 = uniforme]:     ";
    cin>>distr_z;
    cout<<"Inserire distribuzione molteplicità [1 = continua, 0 = costante]:     ";
    cin>>distr_m;
    if(distr_m == true){
        cout<<"Inserire distribuzione molteplicità continua [1 = prof, 0 = uniforme]:     ";
        cin>>m;      
    }else{
        cout<<"Inserire la molteplicità costante:     ";
        cin>>m;           
    }

    FastSim(n, seed, distr_z, distr_m, m);
    Residui();
    Risoluzione();
    Efficienza();

}