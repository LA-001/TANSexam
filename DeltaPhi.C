#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtupleD.h>
#include <TLeaf.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "Generazione.h"
#include "Trasporto.h"

using namespace std;

void DP(int numero, unsigned int seed, bool distr_z, bool distr_m, int m){

    TStopwatch timer;

    TFile *file = TFile::Open("fileRoot/kinem.root");
    if (!file || file->IsZombie()) {
        cout<<"Errore: impossibile aprire il file ROOT."<<endl;
        file->Close();
        delete file;
        return;
    }

    TH1F *hm = (TH1F*)file->Get("hm");
    TH1F *eta = (TH1F*)file->Get("heta2");
    if (!hm || !eta) {
        cout<<"Errore: istogrammi 'hm' o 'heta2' non trovati."<<endl;
        return;
    }

    gRandom->SetSeed(seed);

    Generazione *ptr = new Generazione(eta,hm);
    Trasporto *ptr2 = new Trasporto();

    TH1D *hist = new TH1D("hist","Delta Phi",200,0,0.006);
    hist->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hist->GetYaxis()->SetTitle("Conteggi");

    double H = ptr2->GetHRiv();

    timer.Start();

    for(int tot = 1; tot <= numero; tot++){	  

        double x0 = ptr->VertexSimXY();  
        double y0 = ptr->VertexSimXY();  
        double z0 = ptr->VertexSimZ(distr_z);   
        int molti = ptr->Multiplicity(distr_m, m); 

        double punto[3], versori[3];

        for(int i = 0; i < molti; i++){
            double phi = ptr->Azimut();
            double tetha = ptr->Tetha();

            versori[0] = TMath::Sin(tetha) * TMath::Cos(phi);
            versori[1] = TMath::Sin(tetha) * TMath::Sin(phi);
            versori[2] = TMath::Cos(tetha);

            punto[0] = x0;
            punto[1] = y0;
            punto[2] = z0;

            ptr2->EquazioneRetta(punto, versori, ptr2->GetRPipe());
            ptr2->Scattering(versori, true);
            ptr2->EquazioneRetta(punto, versori, ptr2->GetRLayer1());

            if(-H/2. <= punto[2] && punto[2] <= H/2.) {
                double Phi1 = ptr2->SmearingPhi(punto[0], punto[1], ptr2->GetRLayer1());

                ptr2->Scattering(versori, true);
                ptr2->EquazioneRetta(punto, versori, ptr2->GetRLayer2());

                if(-H/2. <= punto[2] && punto[2] <= H/2.) {
                    double Phi2 = ptr2->SmearingPhi(punto[0], punto[1], ptr2->GetRLayer2());

                    double deltaPhi = TMath::Abs(Phi1 - Phi2);

                    if (deltaPhi > M_PI) {
                        deltaPhi = 2 * M_PI - deltaPhi;
                    }  

                    hist->Fill(deltaPhi);
                }
            } 
        }
    }

    TCanvas *cPhi = new TCanvas("cPhi", "cPhi", 1200, 800);
    hist->DrawCopy();
    cPhi->SetGrid();
    cPhi->SaveAs("Plot/cPhi_3.png");

    timer.Stop();
    timer.Print();

    file->Close();
    delete file;

    delete ptr;
    delete ptr2;
}
