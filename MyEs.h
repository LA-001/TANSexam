#ifndef MYES_H
#define MYES_H

#include "TObject.h"

class MyEs : public TObject{

public:
MyEs();
virtual ~MyEs();

float MediaVector(const vector<float>& vec);
float Intersezione(float r1, float z1, float r2, float z2);
float Intersezione2(float r1, float phi1, float z1, float r2, float phi2, float z2);
vector<float> IntervalloVtx(const vector<float>& vertice);

private:
MyEs(const MyEs& source);
MyEs& operator=(const MyEs& source);

    ClassDef(MyEs,1);
};
#endif
