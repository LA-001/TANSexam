#include "Ricostruzione1.h"
#include <Riostream.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtuple.h>
#include <vector>
#include <TLeaf.h>
#include <algorithm>
#include <limits.h>
#include <TCanvas.h>

using namespace std;

void Ricostruzione(unsigned int seed){

TFile fin("fileRoot/simulazione.root");
TFile fout("fileRoot/ricostruzione.root","recreate");

Ricostruzione1 *ptr3 = new Ricostruzione1(seed);

float x2,y2,z2,x3,y3,z3,label;
float z0,molti,random;
vector<float> rumore(3);

TNtuple *hit = (TNtuple*)fin.Get("hit");
hit->SetBranchAddress("x2",&x2);
hit->SetBranchAddress("y2",&y2);
hit->SetBranchAddress("z2",&z2);
hit->SetBranchAddress("x3",&x3);
hit->SetBranchAddress("y3",&y3);
hit->SetBranchAddress("z3",&z3);
hit->SetBranchAddress("etichetta",&label);

TNtuple *vrt = (TNtuple*)fin.Get("vrt");
vrt->SetBranchAddress("z0", &z0);
vrt->SetBranchAddress("moltiplicita", &molti);

TNtuple *rec1 = new TNtuple("rec1","hit registrati L1","r2:phi2:z2:etichetta");
float xrec1[4];
TNtuple *rec2 = new TNtuple("rec2","hit registrati L2","r3:phi3:z3:etichetta");
float xrec2[4];
TNtuple *MC = new TNtuple("MC","verità montecarlo","z0:molti:etichetta");
float xMC[3];

int evento = -3, a = -3, et = 0;

	for(int ev=0; ev < hit->GetEntries(); ev++){
		hit->GetEntry(ev);

        //cambiamento in coordinate cilindriche applicando lo smearing della ricostruzione
		    xrec1[0] = TMath::Sqrt(x2*x2 + y2*y2);
            xrec1[1] = ptr3->SmearingPhi(x2,y2,2);
            xrec1[2] = ptr3->SmearingZ(z2);

            xrec2[0] = TMath::Sqrt(x3*x3 + y3*y3);
            xrec2[1] = ptr3->SmearingPhi(x3,y3,3);

        if(z3 != 10000.){
            xrec2[2] = ptr3->SmearingZ(z3);
        }

        if(label != a){
        et++;
        a = label;
        }

		xrec1[3] = et;
        xrec2[3] = et;

		rec1->Fill(xrec1);
        rec2->Fill(xrec2);
        
        if(label != evento){

            //rumore
            random = ptr3->GenRumore();
            if(random < 0.0001){
                rumore = ptr3->Rumore();
                if(rumore[0] == ptr3->GetRLayer1()){
                    xrec1[0] = rumore[0];
                    xrec1[1] = rumore[2];
                    xrec1[2] = rumore[2];
                    xrec1[3] = et;                    
                    
                    rec1->Fill(xrec1);
                }else{
                    xrec2[0] = rumore[0];
                    xrec2[1] = rumore[2];
                    xrec2[2] = rumore[2];
                    xrec2[3] = et;                    
                    
                    rec2->Fill(xrec2);
                }
            }           

            //verità montecarlo            
            vrt->GetEntry(label-1);
            xMC[0] = z0;
            xMC[1] = molti;
            
            MC->Fill(xMC);
            evento = label;
        }       
	}

fout.Write();
fout.Close();

delete ptr3;
}
