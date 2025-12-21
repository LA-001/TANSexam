#ifndef TRASPORTO_H
#define TRASPORTO_H

#include "Defstruct.h"  
#include <TObject.h>
#include <TTree.h>
#include <TBranch.h>

class Trasporto : public TObject{

public:

    Trasporto();
    virtual ~Trasporto();

    void Condizione(double x, double y, double z, int strato, int evento);
    void EquazioneRetta(double punto[3], const double versori[3], double R);
    void Scattering(double versori[3], bool on);
    double GetHRiv() const { return fHRiv; }
    double GetRPipe() const { return fRPipe; }
    double GetRLayer1() const { return fRLayer1; }
    double GetRLayer2() const { return fRLayer2; }
    double SmearingPhi(const double x, const double y, const double R);
    double SmearingZ(const double z);
    void Rumore(Hit* xhit2, Hit* xhit3, TTree* hit2, TTree* hit3, int et, bool on);
    double GenRandom();


private:
    Trasporto(const Trasporto& source);
    Trasporto& operator=(const Trasporto& source);

    double fRPipe;           // raggio Beam Pipe
    double fRLayer1;         // raggio primo Layer
    double fRLayer2;         // raggio secondo Layer
    double fHRiv;            //lunghezza del rivelatore
    double fRMSspace;        //RMS spaziale del multiple scattering
    double fRMSz;            //RMS dello smearing di z
    double fRMSphi;          //RMS dello smearing di phi

    ClassDef(Trasporto,1);
};

#endif
