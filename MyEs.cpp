#include "MyEs.h"
#include <vector>
#include <TMath.h>
#include <Riostream.h>
#include <TNtuple.h>
#include "TH1F.h"
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
    
    return -r1/m + z1;
}

double MyEs::Intersezione2(double r1, double phi1, double z1, double r2, double phi2, double z2)
{
    // Converti in coordinate cartesiane
    double x1 = r1 * TMath::Cos(phi1);
    double y1 = r1 * TMath::Sin(phi1);
    double x2 = r2 * TMath::Cos(phi2);
    double y2 = r2 * TMath::Sin(phi2);

    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    // Calcola t solo se direzione ≠ 0
    double t_den = dx*dx + dy*dy;
    if (t_den == 0) return 10000.;

    // Trova t che minimizza distanza da origine (cioè intersezione con asse r=0)
    double t = -(x1*dx + y1*dy) / t_den;

    // Ora ricava z(t)
    double z_inter = z1 + t * dz;

    if (TMath::Abs(z_inter) > 200) return 10000.;

    return z_inter;
}

vector<double> MyEs::IntervalloVtx(const std::vector<double>& vertice, double min, double max){
    std::vector<double> pippo;
    std::vector<double> pluto;

    unsigned int k = 2;
    double step = 2;

    for (double j = min ; j <= max ; j += 0.01) {	
        pippo.clear();

        for (unsigned int i = 0; i < vertice.size(); i++) {
            if (j <= vertice[i] && vertice[i] <= (j + step)) {
                pippo.push_back(vertice[i]);
            }
        }

        if (pippo.size() > k) {
            pluto.clear();
            pluto = pippo;
            k = pluto.size();
        }
    }	

    return pluto;
}

vector<double> MyEs::IntervalloVtx2(const std::vector<double>& vertice){

    vector<double> pluto;

    TH1D *hist = new TH1D("hist","hist",360,-180.,180.);
    for(unsigned int i = 0; i < vertice.size(); i++) hist->Fill(vertice[i]);

    int binMax = hist->GetMaximumBin();
    double xMax = hist->GetBinCenter(binMax);

    double min = xMax - 1.5;
    double max = xMax + 1.5;

    for(unsigned int i = 0; i < vertice.size(); i++){
        if(min <= vertice[i] && vertice[i] <= max) pluto.push_back(vertice[i]);
    }

    delete hist;

    return pluto;
}
