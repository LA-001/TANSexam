#include "Generazione.h"
#include <cmath>
#include <TMath.h>
#include <Riostream.h>
#include <TRandom3.h>
#include <TH1F.h>

ClassImp(Generazione);

Generazione::Generazione(): TObject(),
fRMSz(0.),
fRMSxy(0.),
fHm(NULL),
fEta(NULL)
{ 
    cout << "DEFAULT CONSTRUCTOR - this: " << this << endl;
}

Generazione::Generazione(TH1F *eta, TH1F *hm): TObject(),
fRMSz(53.),
fRMSxy(0.1),
fHm(hm),
fEta(eta)
{ 
    cout << "STANDARD CONSTRUCTOR - this: " << this << endl;
}

Generazione::~Generazione()
{
    cout << "DESTRUCTOR - this: " << this << endl;
}

int Generazione::Multiplicity(bool distr, int nScelto){
    int Multi = 0;

    if (distr == true) {
        if (nScelto == 0)
            Multi = gRandom->Uniform(2, 55);
        else
            Multi = fHm->GetRandom();       
    }
    else {
        Multi = nScelto;
    }

    return Multi;   
}

double Generazione::VertexSimXY(){
    double xy = gRandom->Gaus(0., fRMSxy);
    return xy;
}

double Generazione::VertexSimZ(bool distr){
    double z = 0;

    if (distr == true) {
        z = gRandom->Gaus(0., fRMSz);
    } else {
        z = gRandom->Uniform(-160., 160.);
    }

    return z;
}

double Generazione::Azimut(){
    double azimut = gRandom->Rndm() * 2. * M_PI;
    return azimut;
}

double Generazione::Tetha(){
    double eta = fEta->GetRandom();
    return 2. * TMath::ATan(TMath::Exp(-eta));
}
