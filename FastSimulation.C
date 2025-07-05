#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtupleD.h>
#include <TLeaf.h>
#include <TStopwatch.h>
#include "TRandom.h"

#include "Generazione.h"
#include "Trasporto.h"
#include "Ricostruzione.h"

using namespace std;

void FastSim(int numero, unsigned int seed, bool distr_z, bool distr_m, int m){

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

    TFile fout("fileRoot/simulazione.root","RECREATE");

    gRandom->SetSeed(seed);

    Generazione *ptr = new Generazione(eta,hm);
    Trasporto *ptr2 = new Trasporto();
    Ricostruzione *ptr3 = new Ricostruzione();

    TNtupleD *nt = new TNtupleD("nt","generazione vertice","x0:y0:z0:phi:tetha");
    double xnt[5]; 
    TNtupleD *hit2 = new TNtupleD("hit2","hit layer2","r2:phi2:z2:etichetta");
    double xhit2[4];
    TNtupleD *hit3 = new TNtupleD("hit3","hit layer3","r3:phi3:z3:etichetta");
    double xhit3[4];
    TNtupleD *vrt = new TNtupleD("vrt","verità montecarlo","x0:y0:z0:moltiplicita");
    double xvrt[4]; 

    double x0, y0, z0, molti, phi, tetha;
    double punto[3], versori[3];
    double H = ptr2->GetHRiv();

    timer.Start();

    for(int tot = 1; tot<= numero; tot++){	  

        x0 = ptr->VertexSimXY();  
        y0 = ptr->VertexSimXY();  
        z0 = ptr->VertexSimZ(distr_z);   
        molti = ptr->Multiplicity(distr_m, m); 

        for (int i = 0; i < molti; i++) {
            xnt[0] = x0; 
            xnt[1] = y0; 
            xnt[2] = z0; 
            xnt[3] = ptr->Azimut();
            xnt[4] = ptr->Tetha();
            
            nt->Fill(xnt);
        }

        xvrt[0] = x0;
        xvrt[1] = y0; 
        xvrt[2] = z0;
        xvrt[3] = molti;
        vrt->Fill(xvrt);

        for(int ev = 0; ev < nt->GetEntries(); ev++){
            nt->GetEntry(ev);

            punto[0] = nt->GetLeaf("x0")->GetValue();
            punto[1] = nt->GetLeaf("y0")->GetValue();
            punto[2] = nt->GetLeaf("z0")->GetValue();
            phi = nt->GetLeaf("phi")->GetValue();
            tetha = nt->GetLeaf("tetha")->GetValue();

            ptr2->EquazioneRetta1(punto, phi, tetha);
            ptr2->Scattering1(versori, tetha, phi, false);
            ptr2->EquazioneRetta2(punto, versori, 2);

            if(-H/2. <= punto[2] && punto[2] <= H/2.){
                xhit2[0] = TMath::Sqrt(punto[0]*punto[0] + punto[1]*punto[1]);
                xhit2[1] = ptr3->SmearingPhi(punto[0], punto[1], 2);
                xhit2[2] = ptr3->SmearingZ(punto[2], H);
                xhit2[3] = tot;

                ptr2->Scattering2(versori, false);
                ptr2->EquazioneRetta2(punto, versori, 3);

                if(-H/2. <= punto[2] && punto[2] <= H/2.){
                    xhit3[0] = TMath::Sqrt(punto[0]*punto[0] + punto[1]*punto[1]);
                    xhit3[1] = ptr3->SmearingPhi(punto[0], punto[1], 2);
                    xhit3[2] = ptr3->SmearingZ(punto[2], H);
                    xhit3[3] = tot;
                }

                hit2->Fill(xhit2);
                if(-H/2. <= punto[2] && punto[2] <= H/2.) hit3->Fill(xhit3);
            } 
        }

        if(ptr3->GenRandom() < 0.0001) ptr3->Rumore(hit2, hit3, tot, true);

        nt->Reset();
    }

    timer.Stop();
    timer.Print();

    hit2->Write();
    hit3->Write();
    vrt->Write();
    fout.Purge();
    fout.Close();

    file->Close();
    delete file;

    delete ptr;
    delete ptr2;
    delete ptr3;
}
