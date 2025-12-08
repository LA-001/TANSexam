#include "TROOT.h"
#include "TString.h"
#include <iostream>

void All() {

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

    TString cmd = Form("FastSim(%d, %u, %d, %d, %d);", n, seed, (bool)distr_z, (bool)distr_m, m);

    gROOT->ProcessLine(".L FastSimulation.C+");
    gROOT->ProcessLine(cmd);

    gROOT->ProcessLine(".L Residui.C+");
    gROOT->ProcessLine("Residui();");

    gROOT->ProcessLine(".L Risoluzione.C+");
    gROOT->ProcessLine("Risoluzione();");

    gROOT->ProcessLine(".L Efficienza.C+");
    gROOT->ProcessLine("Efficienza();");
}
