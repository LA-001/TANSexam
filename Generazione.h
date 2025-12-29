#ifndef GENERAZIONE_H
#define GENERAZIONE_H

#include <TObject.h>
#include <TH1F.h>

class Generazione : public TObject {

public:
    Generazione();
    Generazione(TH1F *eta, TH1F *hm);
    virtual ~Generazione();

    int Multiplicity(bool distr, int nScelto) const;
    double VertexSimXY() const;
    double VertexSimZ(bool distr) const;
    double Phi() const;
    double Theta() const;

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
