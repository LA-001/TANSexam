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
fHm(nullptr),
fEta(nullptr)
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

int Generazione::Multiplicity(bool distr, int nScelto) const{
    int Multi = 0;

    if (distr == true) {
        if (nScelto == 0)
            Multi = gRandom->Integer(53) + 3;        // genera uniformemente interi tra [0,52] + 3
        else
            Multi = round(fHm->GetRandom());       
    }
    else {
        Multi = nScelto;
    }

    return Multi;   
}

double Generazione::VertexSimXY() const{
    return gRandom->Gaus(0., fRMSxy);        // gaussiana generatrice della coordinata x e y
}

double Generazione::VertexSimZ(bool distr) const{
    double z = 0;

    if (distr == true) {
        z = gRandom->Gaus(0., fRMSz);        // gaussiana generatrice della coordinata z
    } else {
        z = gRandom->Uniform(-3*fRMSz, 3*fRMSz);        // genera uniformemente tra [-3sigma, 3sigma] con sigma la RMS della coordinata z
    }

    return z;
}

double Generazione::Phi() const{
    return gRandom->Rndm() * 2. * M_PI;       
}

double Generazione::Theta() const{
    double eta = fEta->GetRandom();
    return 2. * TMath::ATan(TMath::Exp(-eta));        // formula inversa della pseudorapidità eta
}
