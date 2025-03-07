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
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>

using namespace std;

void Ricostruzione(unsigned int seed){

TFile fin("fileRoot/simulazione.root");
TFile fout("fileRoot/ricostruzione_tree.root","recreate");

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

TTree *tree1 = new TTree("T1", "TTree del L1");
TTree *tree2 = new TTree("T2", "TTree del L2");

TClonesArray *ptrhits1 = new TClonesArray("hits1",100);
TClonesArray &hitt1 = *ptrhits1;
TClonesArray *ptrhits2 = new TClonesArray("hits2",100);
TClonesArray &hitt2 = *ptrhits2;

typedef struct{
    float r,phi,z;
    int lab;
} HIT;

static HIT point1;
static HIT point2;

//dichiaro i 2 branches

//tree1->Branch("HitL1",&point1.r, "r/F:phi:z:lab/I");
tree1->Branch("HitL1",&ptrhits1);
//tree2->Branch("HitL2",&point2.r, "r/F:phi:z:lab/I");
tree2->Branch("HitL2",&ptrhits2);

TNtuple *MC = new TNtuple("MC","verità montecarlo","z0:molti:etichetta");
float xMC[3];

int evento = -3, a = -3, et = 0;

	for(int ev=0; ev < hit->GetEntries(); ev++){
		hit->GetEntry(ev);

        //cambiamento in coordinate cilindriche applicando lo smearing della ricostruzione
        new(hitt1[ev]) HIT(TMath::Sqrt(x2*x2 + y2*y2),ptr3->SmearingPhi(x2,y2,2),ptr3->SmearingZ(z2),et);		

        if(z3 != 10000.){
            new(hitt2[ev]) HIT(TMath::Sqrt(x3*x3 + y3*y3),ptr3->SmearingPhi(x3,y3,3),ptr3->SmearingZ(z3),et);
        }

        if(label != a){
        et++;
        a = label;
        }

		tree1->Fill();
        tree2->Fill();
        
        if(label != evento){

            //rumore
            /*random = ptr3->GenRumore();
            if(random < 0.0001){
                rumore = ptr3->Rumore();
                if(rumore[0] == ptr3->GetRLayer1()){
                    point1.r = rumore[0];
                    point1.phi = rumore[2];
                    point1.z = rumore[2];
                    point1.lab = et;                    
                    
                    tree1->Fill();
                }else{
                    point2.r = rumore[0];
                    point2.phi = rumore[2];
                    point2.z = rumore[2];
                    point2.lab = et;                    
                    
                    tree2->Fill();
                }
            } */          

            //verità montecarlo            
            vrt->GetEntry(label-1);
            xMC[0] = z0;
            xMC[1] = molti;
            
            MC->Fill(xMC);
            evento = label;
        }       
	}

cout<<tree1->GetEntries()<<endl;
cout<<tree2->GetEntries()<<endl;

fout.Write();
fout.Close();

delete ptr3;
}
