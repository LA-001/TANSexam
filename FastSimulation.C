#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TRandom3.h"
#include "TTree.h"

#include "Generazione.h"
#include "Trasporto.h"
#include "Defstruct.h"

using namespace std;

void FastSim(int numero, unsigned int seed, bool distr_z, bool distr_m, int m) {
    TStopwatch timer;

    TFile *file = TFile::Open("fileRoot/kinem.root");
    if (!file || file->IsZombie()) {
        cout << "Errore: impossibile aprire il file ROOT." << endl;
        return;
    }
    
    TH1F *hm = (TH1F*)file->Get("hm");
    TH1F *eta = (TH1F*)file->Get("heta2");
    if (!hm || !eta) {
        cout << "Errore: istogrammi 'hm' o 'heta2' non trovati." << endl;
        file->Close();
        return;
    }

    TFile fout("fileRoot/simulazione.root", "RECREATE");

    gRandom->SetSeed(seed);

    Generazione *ptr = new Generazione(eta, hm);
    Trasporto *ptr2 = new Trasporto();

    Hit xhitL1, xhitL2;
    Vrt xvrt;

    TTree *T_hitL1 = new TTree("T_hitL1","TTree hit Layer 1");
    T_hitL1->Branch("hitL1", &xhitL1, "r/D:phi/D:z/D:etichetta/I");

    TTree *T_hitL2 = new TTree("T_hitL2","TTree hit Layer 2");
    T_hitL2->Branch("hitL2", &xhitL2, "r/D:phi/D:z/D:etichetta/I");

    TTree *T_vrt = new TTree("T_vrt","TTree della VM");
    T_vrt->Branch("vrt", &xvrt, "x0/D:y0/D:z0/D:moltiplicita/I");

    double H = ptr2->GetHRiv();

    timer.Start();

    for (int tot = 1; tot <= numero; ++tot) {
        double x0 = ptr->VertexSimXY();
        double y0 = ptr->VertexSimXY();
        double z0 = ptr->VertexSimZ(distr_z);
        int molti = ptr->Multiplicity(distr_m, m);

        xvrt = {x0, y0, z0, molti};
        T_vrt->Fill();

        double punto[3], versori[3];

        for (int i = 0; i < molti; ++i) {
            double phi = ptr->Phi();
            double theta = ptr->Theta();

            versori[0] = TMath::Sin(theta) * TMath::Cos(phi);
            versori[1] = TMath::Sin(theta) * TMath::Sin(phi);
            versori[2] = TMath::Cos(theta);

            punto[0] = x0;
            punto[1] = y0;
            punto[2] = z0;

            ptr2->EquazioneRetta(punto, versori, ptr2->GetRPipe());
            ptr2->Scattering(versori, true);
            ptr2->EquazioneRetta(punto, versori, ptr2->GetRLayer1());

            if (-H/2. <= punto[2] && punto[2] <= H/2.) {
                xhitL1.r = ptr2->GetRLayer1();
                xhitL1.phi = ptr2->SmearingPhi(punto[0], punto[1], ptr2->GetRLayer1());
                xhitL1.z = ptr2->SmearingZ(punto[2]);
                xhitL1.etichetta = tot;

                T_hitL1->Fill();

                ptr2->Scattering(versori, true);
                ptr2->EquazioneRetta(punto, versori, ptr2->GetRLayer2());

                if (-H/2. <= punto[2] && punto[2] <= H/2.) {
                    xhitL2.r = ptr2->GetRLayer2();
                    xhitL2.phi = ptr2->SmearingPhi(punto[0], punto[1], ptr2->GetRLayer2());
                    xhitL2.z = ptr2->SmearingZ(punto[2]);
                    xhitL2.etichetta = tot;
                    
                    T_hitL2->Fill();
                }
            }
        }

        ptr2->Rumore(&xhitL1, &xhitL2, T_hitL1, T_hitL2, tot, true);
    }

    timer.Stop();
    timer.Print();

    T_hitL1->Write();
    T_hitL2->Write();
    T_vrt->Write();

    fout.Purge();
    fout.Close();

    file->Close();
    delete file;
    delete ptr;
    delete ptr2;
}
