#ifndef SIMULAZIONE_H
#define SIMULAZIONE_H

#include "TObject.h"
#include <TH1F.h>

class Simulazione1 : public TObject{

public:
Simulazione1();
Simulazione1(TH1F *eta,TH1F *hm,unsigned int seed=65539);
virtual ~Simulazione1();

float VertexSimXY();
float VertexSimZ(bool distr);
int Multiplicity(bool distr, int ndistr);
float Azimut();
float Eta();


private:
Simulazione1(const Simulazione1& source);
Simulazione1& operator=(const Simulazione1& source);

    //tutte le misure sono in millimetri

float fRMSz;              //rms della normale per z
float fRMSxy;             //rms della normale per x e y
int fSeed;
TH1F *fHm;
TH1F *fEta;

   ClassDef(Simulazione1,1);
};

#endif
