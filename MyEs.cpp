#include "MyEs.h"
#include <vector>
#include <TMath.h>
#include <Riostream.h>
#include <TNtuple.h>
#include "TH1F.h"
#include "TCanvas.h"
#include "TSpectrum.h"

ClassImp(MyEs);

MyEs::MyEs(): TObject()
{
     cout<<"DEFAULT CONSTRUCTOR - this:  "<<this<<endl;
}

MyEs::~MyEs()
{
    cout<<"DESTRUCTOR - this:  "<<this<<endl;
}

double MyEs::MediaVector(const std::vector<double>& vec){
    double somma = 0;
    int count = 0;

    for (double valore : vec) {
        somma += valore;
        count++;
    }

    if (count == 0) {
        return 16500;
    }

    return somma / count;
}

double MyEs::Intersezione(double r1, double z1, double r2, double z2){

    double m = (r2 - r1)/(z2 - z1);
    
    return z1 - r1/m;
}

double MyEs::Intersezione2(double r1, double phi1, double z1, double r2, double phi2, double z2)
{
    double x1 = r1 * TMath::Cos(phi1);
    double y1 = r1 * TMath::Sin(phi1);
    double x2 = r2 * TMath::Cos(phi2);
    double y2 = r2 * TMath::Sin(phi2);

    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    double t_den = dx*dx + dy*dy;
    if (t_den < 1e-12) return 10000.;

    double t = -(x1*dx + y1*dy) / t_den;
    double z_proj = z1 + t * dz;

    if (std::abs(z_proj) > 200.0) return 10000.;

    return z_proj;
}

double MyEs::RunWind(const std::vector<double>& vertice){
    vector<double> pluto;
    if (vertice.empty()) {
        return 10000;
    }

    vector<double> v = vertice;
    sort(v.begin(), v.end());

    double W = 3.5;
    double mean = 0.;

    int n = (int)v.size();
    int j = 0;
    int best_i = 0, best_j = 0, best_count = 0;
    for (int i = 0; i < n; i++) {
        while (j < n && v[j] - v[i] <= W) j++;
        int count = j - i;
        if (count > best_count) {
            best_count = count;
            best_i = i;
            best_j = j;     //punta al primo valore esterno a destra della finestra
        }
    }

    if (best_count < 2) {
        return 10000;
    }

    pluto.assign(v.begin() + best_i, v.begin() + best_j);

    mean = MediaVector(pluto);

    return mean;
}
