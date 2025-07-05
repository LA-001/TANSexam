#include "Ricostruzione.h"
#include <cmath>
#include <TMath.h>
#include <Riostream.h>
#include <TRandom3.h>
#include <vector>
#include <TNtuple.h>
#include <Math/Polynomial.h> 
#include <complex> 

ClassImp(Ricostruzione);

Ricostruzione::Ricostruzione(): TObject(),             
fRLayer1(40.),            // millimetri
fRLayer2(70.),            // millimetri  
fRMSz(0.12),              // millimetri
fRMSphi(0.03),            // radianti (corretto)
fB(1.)                   // Tesla
{ 
    cout << "DEFAULT CONSTRUCTOR - this:  " << this << endl;
}

Ricostruzione::~Ricostruzione()
{
    cout << "DESTRUCTOR - this:  " << this << endl;
}

double Ricostruzione::SmearingPhi(double x, double y, int strato){
    double k = TMath::ATan2(y, x);
    if(k < 0) k += 2 * M_PI;

    double PHIrec = 0;
    if(strato == 2){
        PHIrec = k + gRandom->Gaus(0, fRMSphi) / fRLayer1;
    } else {
        PHIrec = k + gRandom->Gaus(0, fRMSphi) / fRLayer2;
    }

    if(PHIrec > 2 * M_PI) PHIrec -= 2 * M_PI;
    if(PHIrec < 0) PHIrec += 2 * M_PI;

    return PHIrec;
}

double Ricostruzione::SmearingZ(double z, double H){
    double zz = 0;
    do {
        zz = z + gRandom->Gaus(0, fRMSz);
    } while(zz > H / 2. || zz < -H / 2.);

    return zz;
}

double Ricostruzione::GenRandom(){
    return gRandom->Rndm();
}

void Ricostruzione::Rumore(TNtupleD *rec1, TNtupleD *rec2, int et, bool on){
    Double_t rumore[4];
   	
    if(on){
        Double_t r = gRandom->Rndm() * 2;
        rumore[0] = (r < 1) ? fRLayer1 : fRLayer2;
        rumore[1] = gRandom->Rndm() * 2 * M_PI;
        rumore[2] = gRandom->Uniform(-135., 135.);
        rumore[3] = static_cast<Double_t>(et);

        if(r < 1){
            rec1->Fill(rumore);
        } else {
            rec2->Fill(rumore);
        }
    }
}
