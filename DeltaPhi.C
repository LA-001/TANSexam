#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TLeaf.h>
#include <TStopwatch.h>
#include "TRandom.h"

#include "Generazione.h"
#include "Trasporto.h"
#include "Ricostruzione.h"

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
    Ricostruzione *ptr3 = new Ricostruzione();

    TH1F *hist = new TH1F("hist","Delta Phi",200,0,0.01);
    TNtuple *nt = new TNtuple("nt","generazione vertice","x0:y0:z0:phi:tetha");
    float xnt[5]; 
    TNtuple *hit2 = new TNtuple("hit2","hit layer2","r2:phi2:z2:etichetta");
    float xhit2[4];
    TNtuple *hit3 = new TNtuple("hit3","hit layer3","r3:phi3:z3:etichetta");
    float xhit3[4]; 


    float x0, y0, z0, molti, phi, tetha;
    float punto[3],versori[3];
    float H = ptr2->GetHRiv();
    float deltaPhi;

    timer.Start();

//**************Generazione dei vertici e delle molteplicità**********************************

    for(int tot = 1; tot<= numero; tot++){	  

        x0 = ptr->VertexSimXY();  
        y0 = ptr->VertexSimXY();  
        z0 = ptr->VertexSimZ(distr_z);   
        molti = ptr->Multiplicity(distr_m,m); 

        for (int i = 0; i < molti; i++) {
                xnt[0] = x0; 
                xnt[1] = y0; 
                xnt[2] = z0; 
                xnt[3] = ptr->Azimut();
                xnt[4] = ptr->Tetha();
                
                nt->Fill(xnt);
                
        }

//**************Trasporto delle particelle generate*******************************************

        for(int ev=0; ev<nt->GetEntries(); ev++){
            nt->GetEntry(ev);

            punto[0] = nt->GetLeaf("x0")->GetValue();
            punto[1] = nt->GetLeaf("y0")->GetValue();
            punto[2] = nt->GetLeaf("z0")->GetValue();
            phi = nt->GetLeaf("phi")->GetValue();
            tetha = nt->GetLeaf("tetha")->GetValue();

            ptr2->EquazioneRetta1(punto, phi, tetha);  //trasporto da generazione a Beam Pipe

            ptr2->Scattering1(versori, tetha, phi, true);  //multiple scattering su Beam Pipe

            ptr2->EquazioneRetta2(punto, versori, 2);  //trasporto da Beam Pipe a Layer 1

            if(-H/2. <= punto[2] && punto[2] <= H/2.){
            
                xhit2[0] = TMath::Sqrt(punto[0]*punto[0] + punto[1]*punto[1]);
                xhit2[1] = ptr3->SmearingPhi(punto[0],punto[1],2);
                xhit2[2] = ptr3->SmearingZ(punto[2],H);
                xhit2[3] = tot;

                ptr2->Scattering2(versori, true);  //multiple scattering sul Layer 1

                ptr2->EquazioneRetta2(punto, versori, 3);  //trasporto dal Layer 1 a Layer 2

                if(-H/2. <= punto[2] && punto[2] <= H/2.){

                    xhit3[0] = TMath::Sqrt(punto[0]*punto[0] + punto[1]*punto[1]);
                    xhit3[1] = ptr3->SmearingPhi(punto[0],punto[1],2);
                    xhit3[2] = ptr3->SmearingZ(punto[2],H);
                    xhit3[3] = tot;

                    deltaPhi = TMath::Abs(xhit2[1] - xhit3[1]);

                    if (deltaPhi > M_PI) {
                        deltaPhi = 2* M_PI - deltaPhi;
                    }  

                    //cout<<deltaPhi<<endl;

                    hist->Fill(deltaPhi);

                }
            } 
        }


        nt->Reset();
    }

hist->DrawCopy();

timer.Stop();
timer.Print();

file->Close();
delete file;

delete ptr;
delete ptr2;
delete ptr3;
}
