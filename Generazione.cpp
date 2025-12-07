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
fGenMaxz(0.),
fGenMinz(0.),
fHm(nullptr),
fEta(nullptr)
{ 
    cout << "DEFAULT CONSTRUCTOR - this: " << this << endl;
}

Generazione::Generazione(TH1F *eta, TH1F *hm): TObject(),
fRMSz(53.),
fRMSxy(0.1),
fGenMaxz(160.),
fGenMinz(-160.),
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
            Multi = gRandom->Integer(53) + 3;
        else
            Multi = round(fHm->GetRandom());       
    }
    else {
        Multi = nScelto;
    }

    return Multi;   
}

double Generazione::VertexSimXY(){
    return gRandom->Gaus(0., fRMSxy);
}

double Generazione::VertexSimZ(bool distr){
    double z = 0;

    if (distr == true) {
        z = gRandom->Gaus(0., fRMSz);
    } else {
        z = gRandom->Uniform(fGenMinz, fGenMaxz);
    }

    return z;
}

double Generazione::Phi(){
    return gRandom->Rndm() * 2. * M_PI;
}

double Generazione::Theta(){
    double eta = fEta->GetRandom();
    return 2. * TMath::ATan(TMath::Exp(-eta));
}
