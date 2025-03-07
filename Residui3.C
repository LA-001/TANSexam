#include "MyEs.h"
#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtuple.h>
#include <vector>
#include <TLeaf.h>
#include <algorithm>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraphErrors.h>

using namespace std;

void Residui3(){
	
	TStopwatch timer;

    MyEs *ptr = new MyEs();

	TFile *fin = TFile::Open("fileRoot/ricostruzione_tree.root");
	TFile *fout = new TFile("fileRoot/output.root", "RECREATE");

    typedef struct{
        float r,phi,z;
        int lab;
    } HIT;

    static HIT point1;
    static HIT point2;

    float z0, molti, labelMC;
	int et = 1;

    TNtuple *MC = (TNtuple*)fin->Get("MC");
	TTree *tree1 = (TTree*)fin->Get("T1");
    TTree *tree2 = (TTree*)fin->Get("T2");
	
    MC->SetBranchAddress("z0", &z0);
    MC->SetBranchAddress("molti", &molti);
    MC->SetBranchAddress("etichetta",&labelMC);

	TBranch *b1 = tree1->GetBranch("HitL1");
    TBranch *b2 = tree2->GetBranch("HitL2");
    
    b1->SetAddress(&point1.r);
    b2->SetAddress(&point2.r);

    float inter, deltaPhi = 0,A,B,C,D;
    vector<float> vertice;
    vector<float> intervalloVer;

    TH1F *hist = new TH1F("hist", "residui [mm]",  300, -0.7, 0.7);
    TH1F *hist2 = new TH1F("hist2", "molteplicita' 4", 100, -1., 1.);
    TH1F *hist3 = new TH1F("hist3", "molteplicita' compresa [45,55]", 100, -0.5, 0.5);

	TNtuple *vtx = new TNtuple("vtx", "vertici ricostruiti e generati", "z0:molti:zRic:moltiRic");
    float xvtx[4];

	int labb = 2;
    float diff = 0;

    float media,a = -56666; //numero random impossibile

	timer.Start();

    for(int ev = 0; ev < tree1->GetEntries(); ev++){
        tree1->GetEvent(ev);
        
		A = point1.r;
        B = point1.phi;            
        C = point1.z; 

        for (int ev2 = 0; ev2 < tree2->GetEntries(); ev2++) {
            tree2->GetEvent(ev2);
                            
            if(point1.lab == point2.lab){
                D = point2.phi;

                deltaPhi = TMath::Abs(B - D);

                if (deltaPhi > TMath::Pi()) {
                    deltaPhi = 2 * TMath::Pi() - deltaPhi;
                }

                if (deltaPhi <= 0.002) {
                    inter = ptr->Intersezione(A, C, point2.r, point2.z);
                    if (inter >= -180. && inter <= 180.) {
                        vertice.push_back(inter);
                    }
                }
            }
        }
        
        if(point1.lab != labb - 1){
        sort(vertice.begin(), vertice.end());
        
                    if(vertice.size() != 0){
					
                    intervalloVer = ptr->IntervalloVtx(vertice);
        			media = ptr->MediaVector(intervalloVer);

                        if(-180. <= media && media <= 180.){
                            MC->GetEntry(point1.lab - 2);
                            
        			        hist->Fill(media - z0);
					        if (molti == 4) hist2->Fill(media - z0);
        			        if (molti > 44 && molti < 56) hist3->Fill(media - z0);
                    
                            xvtx[0] = z0;
					        xvtx[1] = molti;
					        xvtx[2] = media;
					        xvtx[3] = intervalloVer.size();
                            
					        vtx->Fill(xvtx);

                            cout<<z0<<"     "<<media<<"     "<<molti<<"     "<<xvtx[3]<<endl;
                        }
                    }   //chiude if(vertice.size() != 0) 
        }
        if(point1.lab == labb){
            labb++;
            vertice.clear();
        }         
    }
        
    timer.Stop();
	timer.Print();
 
	TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
    hist->Draw();
	c1->Update();
    TCanvas *c2 = new TCanvas("c2", "c2", 700, 500);
    hist2->Draw();
	c2->Update();
    TCanvas *c3 = new TCanvas("c3", "c3", 700, 500);
    hist3->Draw();
	c3->Update();

    vtx->Write();

delete ptr; 
}
