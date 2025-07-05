#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtupleD.h>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStopwatch.h>
#include <TStyle.h>

#include "MyEs.h"

using namespace std;

void Residui(){
	
	TStopwatch timer;

    MyEs *ptr = new MyEs();

	TFile *fin = TFile::Open("fileRoot/simulazione.root");
    if (!fin || fin->IsZombie()) {
        cout<<"Errore: impossibile aprire il file ROOT."<<endl;
        fin->Close();
        delete fin;
        return;
    }

	TFile fout("fileRoot/output.root", "RECREATE");

    double z0, molti, labelMC;
	int et = 1;
	double r2, phi2, z2, r3, phi3, z3, label1, label2;
    double inter, deltaPhi = 0, A, B, C, D;
    int k = 1;
    int lablab = 1;
    double media; 
    int inizio = 0;
    vector<double> vertice;
    vector<double> intervalloVer;

    TNtupleD *MC = (TNtupleD*)fin->Get("vrt");
	TNtupleD *rec1 = (TNtupleD*)fin->Get("hit2");
    TNtupleD *rec2 = (TNtupleD*)fin->Get("hit3");
	
    MC->SetBranchAddress("z0", &z0);
    MC->SetBranchAddress("moltiplicita", &molti);
    
	rec1->SetBranchAddress("r2", &r2);
	rec1->SetBranchAddress("phi2", &phi2);
	rec1->SetBranchAddress("z2", &z2);
    rec1->SetBranchAddress("etichetta", &label1);

	rec2->SetBranchAddress("r3", &r3);
	rec2->SetBranchAddress("phi3", &phi3);
	rec2->SetBranchAddress("z3", &z3);
	rec2->SetBranchAddress("etichetta", &label2);

    TCanvas *cv[3];
    TH1F *hist[3];
    char nome[5];

    hist[0] = new TH1F("hist0", "Residui",  300, -600, 600);
    hist[1] = new TH1F("hist1", "Residui con molteplicita' 6", 100, -600, 600);
    hist[2] = new TH1F("hist2", "Residui con molteplicita' compresa tra [45,55]", 100, -500, 500);

    for(int i = 0; i <= 2; i++){
        hist[i]->GetXaxis()->SetTitle("Residui [#mum]");
        hist[i]->GetYaxis()->SetTitle("Conteggi");
        hist[i]->SetLineColor(kBlack);
    }

	TNtupleD *vtx = new TNtupleD("vtx", "vertici ricostruiti e generati", "z0:molti:zRic");
    double xvtx[3];
    TNtupleD *ver = new TNtupleD("ver", "vertici ricostruiti da mediare", "zric:etichetta");
    double xver[2];

	timer.Start();

    for(long ev = 0; ev < rec1->GetEntries(); ev++){
        rec1->GetEntry(ev);
        
		A = r2;
        B = phi2;            
        C = z2;

        if(label1 != lablab){
            lablab = label1;
            inizio = k;

            if(vertice.size() > 2){
                //sort(vertice.begin(), vertice.end());
                //intervalloVer = ptr->IntervalloVtx(vertice,vertice[0],vertice[vertice.size()-1]);
                intervalloVer = ptr->IntervalloVtx2(vertice);
        	    media = ptr->MediaVector(intervalloVer);

                if(-170. <= media && media <= 170.){
                    MC->GetEntry(label1 - 2);
                            
        	        hist[0]->Fill((media - z0) * 1000);
				    if (molti == 6) hist[1]->Fill((media - z0) * 1000);
        		    if (molti > 44 && molti < 56) hist[2]->Fill((media - z0) * 1000);
                    
                    xvtx[0] = z0;
				    xvtx[1] = molti;
				    xvtx[2] = media;
                            
			        vtx->Fill(xvtx);
                }
                vertice.clear(); 
            }
        }

        for (long ev2 = inizio; ev2 < rec2->GetEntries(); ev2++) {
            rec2->GetEntry(ev2);

            if(label2 != label1){
                k = ev2;
                break;
            }
                            
            D = phi3;

            deltaPhi = TMath::Abs(B - D);

            if (deltaPhi > TMath::Pi()) {
                deltaPhi = 2 * TMath::Pi() - deltaPhi;
            } 

            if (deltaPhi <= 0.0018) {
                inter = ptr->Intersezione(A, C, r3, z3);
                if (inter >= -165. && inter <= 165.) {
                    vertice.push_back(inter);
                }
            }
        }
    }
    
    timer.Stop();
	timer.Print();

	gStyle->SetOptStat(11);
	gStyle->SetOptFit(1);
 
    for(int i = 0; i < 3; i++){
        snprintf(nome, 5, "c%i", i + 1);
        cv[i] = new TCanvas(nome, nome, 1200, 800);
        cv[i]->SetGrid();
        cv[i]->SetLeftMargin(0.1);
        hist[i]->Fit("gaus");
        hist[i]->DrawCopy();
        hist[i]->Write();
        cv[i]->SaveAs(Form("Plot/%s.png", nome));
    }

    vtx->Write();
    fout.Purge();
    fout.Close();

    fin->Close();
    delete fin;

    delete ptr; 
}
