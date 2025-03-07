#ifndef RICOSTRUZIONE1_H
#define RICOSTRUZIONE1_H

#include "TObject.h"

class Ricostruzione1 : public TObject{

public:

Ricostruzione1();
float GetRLayer1() const {return fRLayer1;}
float GetRLayer2() const {return fRLayer2;}
Ricostruzione1(unsigned int seed=65539);
virtual ~Ricostruzione1();
float SmearingPhi(float x, float y, int strato);
float SmearingZ(float z);
vector<float> Rumore();
float GenRumore();

private:
Ricostruzione1(const Ricostruzione1& source);
Ricostruzione1& operator=(const Ricostruzione1& source);

float fRLayer1;         //raggio primo Layer
float fRLayer2;         //raggio secondo Layer
float fRMSz;            //RMS dello smearing di z
float fRMSphi;          //RMS dello smearing di phi
int fSeed;

  ClassDef(Ricostruzione1,1);
};

#endif
