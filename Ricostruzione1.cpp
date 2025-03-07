#include "Ricostruzione1.h"
#include <cmath>
#include <TMath.h>
#include <Riostream.h>
#include <TRandom3.h>
#include <vector>

ClassImp(Ricostruzione1);

Ricostruzione1::Ricostruzione1(): TObject(),
fRLayer1(0.),          
fRLayer2(0.),           
fSeed(0),
fRMSz(0.),              //millimetri
fRMSphi(0.)             //millimetri
{ 
    gRandom->SetSeed(fSeed);
    cout<<"DEFAULT CONSTRUCTOR - this:  "<<this<<endl;
}

Ricostruzione1::Ricostruzione1(unsigned int seed): TObject(),             
fRLayer1(40.),     
fRLayer2(70.),           
fSeed(seed),
fRMSz(0.12),              //millimetri
fRMSphi(0.03)             //millimetri
{ 
    gRandom->SetSeed(fSeed);
    cout<<"DEFAULT CONSTRUCTOR - this:  "<<this<<endl;
}

Ricostruzione1::~Ricostruzione1()
{
    cout<<"DESTRUCTOR - this:  "<<this<<endl;
}

float Ricostruzione1::SmearingPhi(float x, float y, int strato){
    float k = TMath::ATan2(y,x);
    float a = 0;
    if(k < 0) k += 2*M_PI;

    float PHIrec = 0;
    if(strato == 2){
        PHIrec = k + gRandom->Gaus(0,fRMSphi)/fRLayer1;
    }else{
        PHIrec = k + gRandom->Gaus(0,fRMSphi)/fRLayer2;
    }

    if(PHIrec > 2*M_PI) PHIrec += -2*M_PI;

    //cout<<k<<"      "<<PHIrec<<endl;

return PHIrec;
}

float Ricostruzione1::SmearingZ(float z){
    float zz = 0;
    do{
        zz = z + gRandom->Gaus(0,fRMSz);
    }while(zz > 135. || zz < -135.);

    //cout<<z<<"      "<<zz<<endl;

return zz;
}

vector<float> Ricostruzione1::Rumore(){
    vector<float> rumore(3);
    
    float r = gRandom->Rndm()*2;
    if(r < 1){
          rumore[0] = fRLayer1;
    }else rumore[0] = fRLayer2;

    rumore[1] = gRandom->Rndm()*2*M_PI;
    rumore[2] = gRandom->Uniform(-135.,135.);

return rumore;
}
