#ifndef TRASPORTO_H
#define TRASPORTO_H

#include "TObject.h"

class Trasporto : public TObject{

public:

    Trasporto();
    virtual ~Trasporto();

    void Condizione(double x, double y, double z, int strato, int evento);
    void EquazioneRetta1(double punto[3], double azimut, double tetha);
    void EquazioneRetta2(double punto[3], double versori[3], int strato);
    void Scattering1(double versori[3], double tetha, double azimut, bool on);
    void Scattering2(double versori[3], bool on);
    double Phi(double x, double y);
    double GetHRiv() const { return fHRiv; }

private:
    Trasporto(const Trasporto& source);
    Trasporto& operator=(const Trasporto& source);

    double fRPipe;           // raggio Beam Pipe
    double fRLayer1;         // raggio primo Layer
    double fRLayer2;         // raggio secondo Layer
    double fHRiv;
    double fRMSspace;

    ClassDef(Trasporto,1);
};

#endif
