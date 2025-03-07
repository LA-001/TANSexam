#ifndef SIMULAZIONE2_H
#define SIMULAZIONE2_H

#include "TObject.h"

class Simulazione2 : public TObject{

public:

Simulazione2();
Simulazione2(unsigned int seed=65539);
virtual ~Simulazione2();
void Condizioni(float x, float y, float z, int strato, int evento);
vector<float> EquazioneRetta1(float x, float y, float z, float azimut, float tetha);
vector<float> EquazioneRetta2(float x, float y, float z, float c1, float c2, float c3, int strato);
vector<float> Scattering1(float tetha, float azimut, bool on);
vector<float> Scattering2(float c1, float c2, float c3, bool on);
float Phi(float x, float y);

private:
Simulazione2(const Simulazione2& source);
Simulazione2& operator=(const Simulazione2& source);

float fRPipe;           //raggio Beam Pipe
float fRLayer1;         //raggio primo Layer
float fRLayer2;         //raggio secondo Layer
int fSeed;
float fRMSspace;

   ClassDef(Simulazione2,1);
};

#endif
