#include "Simulazione1.h"
#include <cmath>
#include <Riostream.h>
#include <TRandom3.h>
#include <TH1F.h>

ClassImp(Simulazione1);

Simulazione1::Simulazione1(): TObject(),
fRMSz(0.),
fRMSxy(0.),
fSeed(0),
fHm(NULL),
fEta(NULL)
{ 
    gRandom->SetSeed(fSeed);
    cout<<"DEFAULT CONSTRUCTOR - this:  "<<this<<endl;
}


Simulazione1::Simulazione1(TH1F *eta,TH1F *hm,unsigned int seed): TObject(),
fRMSz(53),
fRMSxy(0.1),
fSeed(seed),
fHm(hm),
fEta(eta)
{ 
    gRandom->SetSeed(fSeed);
    cout<<"STANDARD CONSTRUCTOR - this:  "<<this<<endl;
}


Simulazione1::~Simulazione1()
{
    cout<<"DESTRUCTOR - this:  "<<this<<endl;
}

int Simulazione1::Multiplicity(bool distr,int nScelto){
    int Multi = 0;
        //se è true e 0 si utilizza la distribuzione uniforme, se è true ma =! da 0 uso quella assegnata
    if(distr == true){
        if(nScelto == 0) Multi = floor(gRandom->Uniform(3,60));
        if(nScelto != 0) Multi = fHm->GetRandom();       
    }
    else{
        Multi = nScelto;
    }
    
return Multi;   
}

float Simulazione1::VertexSimXY(){
    float xy = gRandom->Gaus(0.,fRMSxy);

return xy;
}

float Simulazione1::VertexSimZ(bool distr){
    float z = 0;

    if(distr == true){
        z = gRandom->Gaus(0.,fRMSz);
    }else{
        z = gRandom->Uniform(-160.,160.);
    }
    
return z;
}

float Simulazione1::Azimut(){
    float azimut = gRandom->Rndm()*2.*M_PI;

return azimut;
}

float Simulazione1::Eta(){
    float eta = 0;

    eta = fEta->GetRandom();

return eta;
}

