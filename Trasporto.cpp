#include "Trasporto.h"
#include <TMath.h>
#include <cmath>
#include <Riostream.h>
#include <TRandom3.h>
#include <vector>

ClassImp(Trasporto);

Trasporto::Trasporto(): TObject(),
fRPipe(30.),             
fRLayer1(40.),     
fRLayer2(70.),           
fHRiv(350.),
fRMSspace(0.001)         // 0.001 radianti
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

void Trasporto::EquazioneRetta1(double punto[3], double phi, double tetha){

    double c1 = TMath::Sin(tetha)*TMath::Cos(phi);
    double c2 = TMath::Sin(tetha)*TMath::Sin(phi);
    double c3 = TMath::Cos(tetha);

    double delta = (punto[0]*c1 + punto[1]*c2)*(punto[0]*c1 + punto[1]*c2) 
                   - (c1*c1 + c2*c2)*(punto[0]*punto[0] + punto[1]*punto[1] - fRPipe*fRPipe);

    if (delta < 0) {
        cout << "Errore: delta negativo, nessuna intersezione trovata." << endl;
        punto[0] = 0;
        punto[1] = 0; 
        punto[2] = 0; 
        return;
    }

    double t = ((-punto[0]*c1 - punto[1]*c2) + TMath::Sqrt(delta))/(c1*c1 + c2*c2);

    punto[0] = punto[0] + c1*t;
    punto[1] = punto[1] + c2*t;
    punto[2] = punto[2] + c3*t; 
    
}

/******** Restituisce i punti di intersezione con i rispettivi layer dando in ingresso la nuova direzione (c1,c2,c3) ********/

void Trasporto::EquazioneRetta2(double punto[3], double versori[3], int strato){ 
    double delta = 0;

    if(strato == 2) 
        delta = (punto[0]*versori[0] + punto[1]*versori[1])*(punto[0]*versori[0] + punto[1]*versori[1]) 
                - (versori[0]*versori[0] + versori[1]*versori[1])*(punto[0]*punto[0] + punto[1]*punto[1] - fRLayer1*fRLayer1);

    else if(strato == 3) 
        delta = (punto[0]*versori[0] + punto[1]*versori[1])*(punto[0]*versori[0] + punto[1]*versori[1]) 
                - (versori[0]*versori[0] + versori[1]*versori[1])*(punto[0]*punto[0] + punto[1]*punto[1] - fRLayer2*fRLayer2);

    else{
        cout<<"Errore: codice identificativo del Layer non riconosciuto"<<endl;
    }

    if (delta < 0) {
        cout << "Errore: delta negativo, nessuna intersezione trovata." << endl;
        punto[0] = 0;
        punto[1] = 0; 
        punto[2] = 0; 
        return;
    }

    double t = ((-punto[0]*versori[0] - punto[1]*versori[1]) + TMath::Sqrt(delta))/(versori[0]*versori[0] + versori[1]*versori[1]);

    punto[0] = punto[0] + versori[0]*t;
    punto[1] = punto[1] + versori[1]*t;
    punto[2] = punto[2] + versori[2]*t;
    
}

/**************** Restituisce la direzione dopo lo scattering multiplo nel sistema di riferimento del laboratorio *************/

void Trasporto::Scattering1(double versori[3], double tetha, double azimut, bool on){

    double a = fRMSspace/TMath::Sqrt(2);      // tetha primo
    double b = gRandom->Rndm()*2*M_PI;        // azimut primo

    if(on == true){              
        versori[0] = -TMath::Sin(a)*TMath::Cos(b)*TMath::Sin(azimut) 
                     - TMath::Sin(a)*TMath::Sin(b)*TMath::Cos(azimut)*TMath::Cos(tetha) 
                     + TMath::Cos(a)*TMath::Sin(tetha)*TMath::Cos(azimut);

        versori[1] =  TMath::Sin(a)*TMath::Cos(b)*TMath::Cos(azimut) 
                     - TMath::Sin(a)*TMath::Sin(b)*TMath::Sin(azimut)*TMath::Cos(tetha) 
                     + TMath::Cos(a)*TMath::Sin(tetha)*TMath::Sin(azimut);

        versori[2] =  TMath::Sin(a)*TMath::Sin(b)*TMath::Sin(tetha)  
                     + TMath::Cos(a)*TMath::Cos(tetha);
    }else{
        versori[0] = sin(tetha)*cos(azimut);
        versori[1] = sin(tetha)*sin(azimut);
        versori[2] = cos(tetha);
    }
}

void Trasporto::Scattering2(double versori[3], bool on){

    double a = fRMSspace/TMath::Sqrt(2);      // tetha primo
    double b = gRandom->Rndm()*2*M_PI;        // azimut primo

    double radice = TMath::Sqrt(versori[0]*versori[0] + versori[1]*versori[1]);

    if(on == true){              
        versori[0] = -TMath::Sin(a)*TMath::Cos(b)*versori[1]/radice 
                     - TMath::Sin(a)*TMath::Sin(b)*versori[2]*versori[0]/radice 
                     + TMath::Cos(a)*versori[0];

        versori[1] =  TMath::Sin(a)*TMath::Cos(b)*versori[0]/radice 
                     - TMath::Sin(a)*TMath::Sin(b)*versori[1]*versori[2]/radice 
                     + TMath::Cos(a)*versori[1];

        versori[2] =  TMath::Sin(a)*TMath::Sin(b)*radice 
                     + TMath::Cos(a)*versori[2];
    }
}

double Trasporto::Phi(double x, double y){
    double k = TMath::ATan2(y,x);
    if(k < 0) k += 2*M_PI;
    return k;
}
