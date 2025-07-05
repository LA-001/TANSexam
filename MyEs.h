#ifndef MYES_H
#define MYES_H

#include "TObject.h"

class MyEs : public TObject{

public:
MyEs();
virtual ~MyEs();

double MediaVector(const vector<double>& vec1);
double Intersezione(double r1, double z1, double r2, double z2);
double Intersezione2(double r1, double phi1, double z1, double r2, double phi2, double z2);
vector<double> IntervalloVtx(const vector<double>& vertice, double min, double max);
vector<double> IntervalloVtx2(const std::vector<double>& vertice);

private:
MyEs(const MyEs& source);
MyEs& operator=(const MyEs& source);

    ClassDef(MyEs,1);
};
#endif

