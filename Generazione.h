#ifndef GENERAZIONE_H
#define GENERAZIONE_H

#include "TObject.h"
#include <TH1F.h>

class Generazione : public TObject {

public:
    Generazione();
    Generazione(TH1F *eta, TH1F *hm);
    virtual ~Generazione();

    double VertexSimXY();
    double VertexSimZ(bool distr);
    int Multiplicity(bool distr, int ndistr);
    double Azimut();
    double Tetha();

private:
    Generazione(const Generazione& source);
    Generazione& operator=(const Generazione& source);

    // tutte le misure sono in millimetri
    double fRMSz;     // rms della normale per z
    double fRMSxy;    // rms della normale per x e y
    TH1F *fHm;
    TH1F *fEta;

    ClassDef(Generazione,1);
};

#endif
