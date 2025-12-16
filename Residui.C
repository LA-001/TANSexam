#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TStyle.h"

#include <vector>
#include <map>

#include "MyEs.h"
#include "Defstruct.h"

using namespace std;

void Residui() {

    TStopwatch timer;

    MyEs *ptr = new MyEs();

    TFile *fin = TFile::Open("fileRoot/simulazione.root");
    if (!fin || fin->IsZombie()) {
        cout << "Errore: impossibile aprire il file ROOT." << endl;
        fin->Close();
        delete fin;
        return;
    }

    TFile fout("fileRoot/residui.root", "RECREATE");

    Vrt mc;
    Hit hitL1, hitL2;
    Ric xvtx;

    TTree *T_vtx = new TTree("T_vtx","TTree vertici ricostruiti e veri");
    T_vtx->Branch("vtx", &xvtx, "z0/D:zRic/D:molti/I");

    TTree *T_MC = (TTree*)fin->Get("T_vrt");
    TTree *T_rec1 = (TTree*)fin->Get("T_hitL1");
    TTree *T_rec2 = (TTree*)fin->Get("T_hitL2");
    T_MC->SetBranchAddress("vrt", &mc);
    T_rec1->SetBranchAddress("hitL1", &hitL1);
    T_rec2->SetBranchAddress("hitL2", &hitL2);

    int lab = -1;
    double Phi_VM = 0.004, z_max = 170., z_min = -170.;

    vector<double> vertice;
    vector<double> intervalloVer;
    map<int, vector<Hit>> hitsByLabel;      //map<chiave,valore> nomeMappa, uso come chiave l'etichetta int e gli aggiungo il vettore associato

    timer.Start();

    for (long i = 0; i < T_rec2->GetEntries(); ++i) {
        T_rec2->GetEntry(i);
        hitsByLabel[hitL2.etichetta].push_back(hitL2);        //nomeMappa[chiave].push_back(vettore), è come se mi creasse una tabella con una colonna le chiavi e un'altra i vettori
    }

    TCanvas *cv[3];
    TH1D *hist[3];
    char nome[10];

    hist[0] = new TH1D("hist0", "Residui",  300, -600, 600);
    hist[1] = new TH1D("hist1", "Residui con molteplicita' 6", 101, -600, 600);
    hist[2] = new TH1D("hist2", "Residui con molteplicita' compresa tra [45,55]", 101, -500, 500);

    for(int i = 0; i <= 2; i++){
        hist[i]->GetXaxis()->SetTitle("Residui [#mum]");
        hist[i]->GetYaxis()->SetTitle("Conteggi");
        hist[i]->SetLineColor(kBlack);
    }

//-----------------------------------------------------------------------------------------------------------------

    for(long ev = 0; ev < T_rec1->GetEntries(); ++ev){
        T_rec1->GetEntry(ev);

        if (lab != hitL1.etichetta){
            double media = ptr->RunWind(vertice);

            if (z_min <= media && media <= z_max) {
                T_MC->GetEntry(hitL1.etichetta - 2);

                hist[0]->Fill((media - mc.z0) * 1000);
                if (mc.moltiplicita == 6) hist[1]->Fill((media - mc.z0) * 1000);
                if (mc.moltiplicita > 44 && mc.moltiplicita < 56) hist[2]->Fill((media - mc.z0) * 1000);

                xvtx = {mc.z0, media, mc.moltiplicita};
                T_vtx->Fill();
            }

            vertice.clear();
            lab = hitL1.etichetta;
        }

        double phi_L1 = hitL1.phi;

        map<int, vector<Hit>>::iterator it = hitsByLabel.find(hitL1.etichetta);
        if (it != hitsByLabel.end()) {          //se .find() non ha riscontro di etichetta rilascia un iteratore speciale .end(), se così fosse non avrei nulla da cercare e dunque salta l'if
            for (const Hit& h : it->second) {   //Per ogni elemento (Hit) contenuto nel vettore it->second, crea un riferimento costante chiamato h e fai qualcosa con esso
                double phi_L2 = h.phi;

                double deltaPhi = TMath::Abs(phi_L1 - phi_L2);
                if (deltaPhi > TMath::Pi()) deltaPhi = 2 * TMath::Pi() - deltaPhi;

                if (deltaPhi <= Phi_VM) {
                    double inter = ptr->Intersezione(hitL1.r, hitL1.z, h.r, h.z);
                    if (inter >= z_min && inter <= z_max) {
                        vertice.push_back(inter);
                    }
                }
            }
        }
    }

    timer.Stop();

//-----------------------------------------------------------------------------------------------------------------

    gStyle->SetOptStat(11);
    gStyle->SetOptFit(1);

    for(int i = 0; i < 3; i++){
        snprintf(nome, 5, "c%i", i + 1);
        cv[i] = new TCanvas(nome, nome, 1200, 800);
        cv[i]->SetGrid();
        cv[i]->SetLeftMargin(0.13);
        hist[i]->Fit("gaus");
        hist[i]->DrawCopy();
        hist[i]->Write();
        cv[i]->SaveAs(Form("Plot/%s.png", nome));
    }

    T_vtx->Write();
    fout.Purge();
    fout.Close();

    fin->Close();
    delete fin;
    delete ptr;

    timer.Print();
}
