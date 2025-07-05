#ifndef RICOSTRUZIONE_H
#define RICOSTRUZIONE_H

#include <TNtupleD.h>

#include "TObject.h"

class Ricostruzione : public TObject{

public:

Ricostruzione();
double GetRLayer1() const {return fRLayer1;}
double GetRLayer2() const {return fRLayer2;}
virtual ~Ricostruzione();
double SmearingPhi(double x, double y, int strato);
double SmearingZ(double z, double H);
void Rumore(TNtupleD *rec1,TNtupleD *rec2, int et, bool on);
double GenRandom();

private:
Ricostruzione(const Ricostruzione& source);
Ricostruzione& operator=(const Ricostruzione& source);

double fRLayer1;         //raggio primo Layer
double fRLayer2;         //raggio secondo Layer
double fRMSz;            //RMS dello smearing di z
double fRMSphi;          //RMS dello smearing di phi
double fB;               //intensità campo magnetico parallelo al fascio

  ClassDef(Ricostruzione,1);
};

#endif
