#include "Simulazione2.h"
#include <TMath.h>
#include <cmath>
#include <Riostream.h>
#include <TRandom3.h>
#include <vector>

ClassImp(Simulazione2);

Simulazione2::Simulazione2(): TObject(),
fRPipe(0.),            
fRLayer1(0.),          
fRLayer2(0.),           
fSeed(0),
fRMSspace(0.)
{ 
    gRandom->SetSeed(fSeed);
    cout<<"DEFAULT CONSTRUCTOR - this:  "<<this<<endl;
}

Simulazione2::Simulazione2(unsigned int seed): TObject(),
fRPipe(30),             
fRLayer1(40),     
fRLayer2(70),           
fSeed(seed),
fRMSspace(0.001)         //1 milliradianti
{ 
    gRandom->SetSeed(fSeed);
    cout<<"STANDARD CONSTRUCTOR - this:  "<<this<<endl;
}

Simulazione2::~Simulazione2()
{
    cout<<"DESTRUCTOR - this:  "<<this<<endl;
}

void Simulazione2::Condizioni(float x,float y,float z,int strato, int evento){
    float k = TMath::Sqrt(x*x + y*y);
     if (strato == 1) {
        if (k >= (fRPipe-0.1) && k <= (fRPipe+0.1)) {
        } else {
            cout << "La condizione per lo strato 1 è violata nell'evento: " << evento <<endl;
        }
    }

    if (strato == 2) {
        if (k >= (fRLayer1-0.1) && k <= (fRLayer1+0.1) && -135 <= z && z <= 135) {
        } else {
            cout << "La condizione per lo strato 2 è violata nell'evento: " << evento <<", z = "<<z<<", raggio = "<<k<<endl;
        }
    }

    if (strato == 3) {
        if (k >= (fRLayer2-0.1) && k <= (fRLayer2+0.1) && -135 <= z && z <= 135) {
        } else {
            cout << "La condizione per lo strato 3 è violata nell'evento: " << evento <<", z = "<<z<<", raggio = "<<k<<endl;
        }
    }
}

vector<float> Simulazione2::EquazioneRetta1(float x, float y, float z, float azimut, float tetha){
    vector<float> punto(3);

    float c1 = TMath::Sin(tetha)*TMath::Cos(azimut);
    float c2 = TMath::Sin(tetha)*TMath::Sin(azimut);
    float c3 = TMath::Cos(tetha);
    float delta = (x*c1 + y*c2)*(x*c1 + y*c2) - (c1*c1 + c2*c2)*(x*x + y*y - fRPipe*fRPipe);
    float t = ((-x*c1 - y*c2) + TMath::Sqrt(delta))/(c1*c1 + c2*c2);

    if (delta < 0) {
        cout << "Errore: delta negativo, nessuna intersezione trovata." << endl;
        punto[0] = punto[1] = punto[2] = 0; 
        return punto;
    }

    punto[0] = x + c1*t;
    punto[1] = y + c2*t;
    punto[2] = z + c3*t;
  
return punto; 
}

/********resistituisce i punti di intersezione con i rispettivi layer dando in ingresso la nuova direzione (c1,c2,c3)********/

vector<float> Simulazione2::EquazioneRetta2(float x, float y, float z, float c1, float c2, float c3, int strato){
    vector<float> punto(3);   
    float delta = 0;
    if(strato == 2) delta = (x*c1 + y*c2)*(x*c1 + y*c2) -(c1*c1 + c2*c2)*(x*x + y*y - fRLayer1*fRLayer1);
    if(strato == 3) delta = (x*c1 + y*c2)*(x*c1 + y*c2) -(c1*c1 + c2*c2)*(x*x + y*y - fRLayer2*fRLayer2);
    if(strato != 2 && strato != 3){
        cout<<"Errore: codice identificativo del Layer non riconosciuto"<<endl;
    }
    float t = ((-x*c1 - y*c2) + TMath::Sqrt(delta))/(c1*c1 + c2*c2);

    if (delta < 0) {
        cout << "Errore: delta negativo, nessuna intersezione trovata." << endl;
        punto[0] = punto[1] = punto[2] = 0; 
        return punto;
    }
    
    punto[0] = x + c1*t;
    punto[1] = y + c2*t;
    punto[2] = z + c3*t;
    
return punto; 
}

/****************mi restituisce la direzione dopo lo scattering multiplo nel SR del laboratorio*************/

vector<float> Simulazione2::Scattering1(float tetha, float azimut, bool on){
    float a = fRMSspace/TMath::Sqrt(2);      //tetha primo
    float b = gRandom->Rndm()*2*M_PI;                  //azimut primo
    vector<float> scatter(3);

    if(on == true){              
        scatter[0] = -TMath::Sin(a)*TMath::Cos(b)*TMath::Sin(azimut) - TMath::Sin(a)*TMath::Sin(b)*TMath::Cos(azimut)*TMath::Cos(tetha) + TMath::Cos(a)*TMath::Sin(tetha)*TMath::Cos(azimut);
        scatter[1] =  TMath::Sin(a)*TMath::Cos(b)*TMath::Cos(azimut) - TMath::Sin(a)*TMath::Sin(b)*TMath::Sin(azimut)*TMath::Cos(tetha) + TMath::Cos(a)*TMath::Sin(tetha)*TMath::Sin(azimut);
        scatter[2] =  TMath::Sin(a)*TMath::Sin(b)*TMath::Sin(tetha)  + TMath::Cos(a)*TMath::Cos(tetha);
    }else{
        scatter[0] = TMath::Sin(tetha)*TMath::Cos(azimut);
        scatter[1] = TMath::Sin(tetha)*TMath::Sin(azimut);
        scatter[2] = TMath::Cos(tetha);
    }

return scatter;
}

vector<float> Simulazione2::Scattering2(float c1, float c2, float c3, bool on){
    float a = fRMSspace/TMath::Sqrt(2);      //tetha primo
    float b = gRandom->Rndm()*2*M_PI;                  //azimut primo
    vector<float> scatter2(3);

    float radice = TMath::Sqrt(c1*c1 + c2*c2);

    if(on == true){              
        scatter2[0] = -TMath::Sin(a)*TMath::Cos(b)*c2/radice - TMath::Sin(a)*TMath::Sin(b)*c3*c1/radice + TMath::Cos(a)*c1;
        scatter2[1] =  TMath::Sin(a)*TMath::Cos(b)*c1/radice - TMath::Sin(a)*TMath::Sin(b)*c2*c3/radice + TMath::Cos(a)*c2;
        scatter2[2] =  TMath::Sin(a)*TMath::Sin(b)*radice + TMath::Cos(a)*c3;
    }else{
        scatter2[0] = c1;
        scatter2[1] = c2;
        scatter2[2] = c3;
    }

return scatter2;
}

float Simulazione2::Phi(float x, float y){
	float k = TMath::ATan2(y,x);
    if(k < 0) k += 2*M_PI;

return k;
}
