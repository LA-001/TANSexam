#include "MyEs.h"
#include <vector>
#include <TMath.h>
#include <Riostream.h>

ClassImp(MyEs);

MyEs::MyEs(): TObject()
{
     cout<<"DEFAULT CONSTRUCTOR - this:  "<<this<<endl;
}

MyEs::~MyEs()
{
    cout<<"DESTRUCTOR - this:  "<<this<<endl;
}

float MyEs::MediaVector(const vector<float>& vec){
    float somma = 0;
    int count = 0;

    for (float valore : vec) {
            somma += valore;
            count++;
    }

    if (count == 0) {
        //cout << "vettore vuoto" << endl;
        return 1650;
    }

return somma / count;
}

float MyEs::Intersezione(float r1, float z1, float r2, float z2){

    float m = (r2 - r1)/(z2 - z1);
    
return -r1/m + z1;
}

float MyEs::Intersezione2(float r1, float phi1, float z1, float r2, float phi2, float z2)
{
	float direzioneX = r2*TMath::Cos(phi2) - r1*TMath::Cos(phi1);
    float direzioneY = r2*TMath::Sin(phi2) - r1*TMath::Sin(phi1);
    float direzioneZ = z2 - z1;

    if (direzioneX == 0) direzioneX = 0.000001; 
    if (direzioneY == 0) direzioneY = 0.000001;

    float t1 = -r1*TMath::Sin(phi1) / direzioneY;
    float t2 = -r1*TMath::Cos(phi1) / direzioneX;

    float Zhit1 = z1 + t1 * direzioneZ;
    float Zhit2 = z1 + t2 * direzioneZ;

    //if(Zhit2 > 1000 && z2 != 10000) cout<<z1<<"      "<<z2<<"      "<<Zhit1<<"      "<<Zhit2<<endl;

    if(TMath::Abs(Zhit1) > 200 && TMath::Abs(Zhit2) < 200) return Zhit2;
    if(TMath::Abs(Zhit2) > 200 && TMath::Abs(Zhit1) < 200) return Zhit1;
	
    if(TMath::Abs(Zhit1 - Zhit2) < 80){
        return (Zhit1 + Zhit2)/2;
    }else{
        return 10000.;
    } 
}


vector<float> MyEs::IntervalloVtx(const vector<float>& vertice){

    vector<float> pippo;
    vector<float> pluto;

    unsigned int k = 1;
    float step = 2.5;

     for (float j = -180.; j <= (180. - step); j += 0.05) {	
        pippo.clear();

        for (unsigned int i = 0; i < vertice.size(); i++) {
          if (j <= vertice[i] && vertice[i] <= (j + step)) {
                   pippo.push_back(vertice[i]);
          }
        }

         if (pippo.size() > k) {
          pluto.clear();
          pluto = pippo;
          k = pluto.size();
          }
      }	


return pluto;
}
