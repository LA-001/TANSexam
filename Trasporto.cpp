#include "Trasporto.h"
#include <TMath.h>
#include <Riostream.h>
#include <TRandom3.h>

ClassImp(Trasporto);

Trasporto::Trasporto(): TObject(),
fRPipe(30.),             
fRLayer1(40.),     
fRLayer2(70.),           
fHRiv(270.),
fRMSspace(0.001),         // 0.001 radianti
fRMSz(0.12),              // millimetri
fRMSphi(0.03)            // radianti
{ 
    cout<<"DEFAULT CONSTRUCTOR - this:  "<<this<<endl;
}

Trasporto::~Trasporto()
{
    cout<<"DESTRUCTOR - this:  "<<this<<endl;
}

void Trasporto::Condizione(double x,double y,double z,int strato, int evento){
    double k = TMath::Sqrt(x*x + y*y);
    if (strato == 1) {
        if (k >= (fRPipe-0.01) && k <= (fRPipe+0.01)) {
        } else {
            cout << "La condizione per la Beam Pipe è violata nell'evento: " << evento <<endl;
        }
    }

    if (strato == 2) {
        if (k >= (fRLayer1-0.1) && k <= (fRLayer1+0.1)) {
        } else {
            cout << "La condizione per il Layer 1 è violata nell'evento: " << evento <<", z = "<<z<<", raggio = "<<k<<endl;
        }
    }

    if (strato == 3) {
        if (k >= (fRLayer2-0.1) && k <= (fRLayer2+0.1)) {
        } else {
            cout << "La condizione per il Layer 2 è violata nell'evento: " << evento <<", z = "<<z<<", raggio = "<<k<<endl;
        }
    }
}

void Trasporto::EquazioneRetta(double punto[3], const double versori[3], double R){ 
    double pv00 = punto[0]*versori[0];
    double pv11 = punto[1]*versori[1];
    double vv00 = versori[0]*versori[0];
    double vv11 = versori[1]*versori[1];
    double pp00 = punto[0]*punto[0];
    double pp11 = punto[1]*punto[1];

    double delta = (pv00 + pv11)*(pv00 + pv11) - (vv00 + vv11)*(pp00 + pp11 - R*R);       

    if (delta < 0) {
        cout << "Errore: delta negativo, nessuna intersezione trovata." << endl;
        punto[0] = 1e4;
        punto[1] = 1e4;
        punto[2] = 1e4;
        return;
    }

    double t = ((-pv00 - pv11) + TMath::Sqrt(delta))/(vv00 + vv11);

    punto[0] = punto[0] + versori[0]*t;
    punto[1] = punto[1] + versori[1]*t;
    punto[2] = punto[2] + versori[2]*t;
    
}

void Trasporto::Scattering(double versore[3], bool on){
    if (!on) return;        \\ nel caso fosse false no si ha MCS, non cambia il vettore versore

    double T[3][3];

    double A = - versore[1]/TMath::Sqrt(versore[0]*versore[0] + versore[1]*versore[1]);
    double B =   versore[0]/TMath::Sqrt(versore[0]*versore[0] + versore[1]*versore[1]);
   
    T[0][0] = A;     T[0][1] = -versore[2]*B;                     T[0][2] = versore[0];
    T[1][0] = B;     T[1][1] = versore[2]*A;                      T[1][2] = versore[1];
    T[2][0] = 0.;    T[2][1] = (versore[0]*B - versore[1]*A);     T[2][2] = versore[2];

    double thp = fRMSspace;
    double php = gRandom->Rndm() * 2 * TMath::Pi();

    double u[3];
    u[0] = TMath::Sin(thp) * TMath::Cos(php);
    u[1] = TMath::Sin(thp) * TMath::Sin(php);
    u[2] = TMath::Cos(thp);

    for (Int_t i = 0; i < 3; i++) {
        versore[i] = 0.;
        for (Int_t j = 0; j < 3; j++) {
            versore[i] += T[i][j] * u[j];
        }
    }  

}

double Trasporto::SmearingPhi(const double x, const double y, const double R){
    double k = TMath::ATan2(y, x);
    if(k < 0) k += 2 * TMath::Pi();

    double PHIrec = 0;
    PHIrec = k + gRandom->Gaus(0, fRMSphi) / R;

    if(PHIrec > 2. * TMath::Pi()) PHIrec -= 2 * TMath::Pi();
    if(PHIrec < 0.) PHIrec += 2 * TMath::Pi();

    return PHIrec;
}

double Trasporto::SmearingZ(const double z){
    double zz = 0;
    zz = z + gRandom->Gaus(0, fRMSz);
    
    if(zz > fHRiv/2.) zz = fHRiv/2. - 0.0001;
    if(zz < -fHRiv/2.) zz = -fHRiv/2. + 0.0001;

    return zz;
}

double Trasporto::GenRandom(){
    return gRandom->Rndm();
}


void Trasporto::Rumore(Hit* xhitL1, Hit* xhitL2, TTree* T_hitL1, TTree* T_hitL2, int et, bool on) {
    if (!on) return;

    for(int i = 0; i < 2; i++){
     
        xhitL1->r = fRLayer1;
        xhitL1->phi = gRandom->Rndm() * 2 * TMath::Pi();
        xhitL1->z = gRandom->Uniform(-fHRiv/2., fHRiv/2.);
        xhitL1->etichetta = et;

        T_hitL1->Fill();

        xhitL2->r = fRLayer2;
        xhitL2->phi = gRandom->Rndm() * 2 * TMath::Pi();
        xhitL2->z = gRandom->Uniform(-fHRiv/2., fHRiv/2.);
        xhitL2->etichetta = et;

        T_hitL2->Fill();
    }
}

