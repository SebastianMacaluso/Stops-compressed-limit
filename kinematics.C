 //
// check parton-level kinematics for MadGraph runs
//
//  specialized to study parton-level direct stop pair production in l+jets (though some results are adapted to dilep by 
//  pretending quarks are leptons/neutrinos), or a single leptonic stop decay in its rest frame.  in the latter case, 
//  many of the "hadronic side" variables simply never get set.
//

#include<fstream>
#include<iostream>
#include<vector>
#include<complex>
#include<algorithm>
#include <cmath>

#include "LHEF.h"  // Les Houches event format headers

#include "TROOT.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMinuit.h"
#include <TStyle.h>


TH1F *hjet1pt = new TH1F("hjet1pt","MT",20,0,300);
TH1F *hjet2pt = new TH1F("hjet2pt","pT jets",20,0,300);
TH1F *hjet3pt = new TH1F("hjet3pt","MT",20,0,300);
TH1F *hjet4pt = new TH1F("hjet4pt","pT jets",20,0,300);
TH1F *hjet5pt = new TH1F("hjet5pt","pT jets",30,0,3.15);
TH1F *hjet6pt = new TH1F("hjet6pt","dPhiMin",30,0,3.15);
TH1F *hjet7pt = new TH1F("hjet7pt","pT jets",150,0,500);









using namespace std;

typedef complex<double> dcomplex;

bool smear = false;

int lumi = 300;

double ptjcut = 20;

const double pi = 3.14159265359;      

// event file info

//const int nSims = 1;  string dataDir = ".";  string fileNames[nSims] = {"multijets400_7j_Brock"};

//const int nSims = 5;  string dataDir = ".";  string fileNames[nSims] = {"multijets_60_80_btag_cuts","multijets_60_80_btag_nocuts","multijets_60_80_nobtag_nocuts","multijets400_7j_Brock","multijets400_7j_Seb"};

//const int nSims = 2;  string dataDir = ".";  string fileNames[nSims] = {"t1t1~j_250_75","tt~j"};

//const int nSims = 1;  string dataDir = ".";  string fileNames[nSims] = {"t1t1~j_250_75"};

const int nSims = 10;  string dataDir = ".";  string fileNames[nSims] = {"tt~_incl","bbjj_ccjj","bbjz_ccjz","bbjw_ccjw","tt~w_incl","tt~z_incl","t1t1j_250_75","t1t1jj_250_75","t1t1jj_250_75_MINJETS_ON","t1t1jj_250_75_MINJETS_OFF"};


//const int nSims = 2;  string dataDir = ".";  string fileNames[nSims] = {"t1t1j","t1t1j"};

//const int nSims = 14;  string dataDir = ".";  string fileNames[nSims] = {"tt~_incl","tt~z_incl","tt~w_incl","bbjj_ccjj","bbjz_ccjz","bbjw_ccjw","t1t1~j_no_matching_273_100","t1t1~_225_50_incl","t1t1~_250_75_incl","t1t1~_275_100_incl","t1t1~_300_125_incl","t1t1~_350_175_incl","t1t1~_400_225_incl","t1t1~_448_270_incl"};

//const int nSims = 10;  string dataDir = ".";  string fileNames[nSims] = {"tt~","tt~j","tt~jj","tt~l","tt~lj","tt~ljj","t1t1~_250_75","t1t1~j_250_75","t1t1~_448_270","t1t1~j_448_270"};


//const int nSims = 10;  string dataDir = ".";  string fileNames[nSims] = {"tt~","tt~j","tt~jj","tt~l","tt~lj","tt~ljj","t1t1~_250_75","t1t1~j_250_75","t1t1~_448_270","multijets400_Brock"};

//const int nSims = 13;  string dataDir = ".";  string fileNames[nSims] = {"tt~","tt~j","tt~jj","tt~l","tt~lj","tt~ljj","t1t1~_183_10","t1t1~j_183_10","t1t1~_250_75","t1t1~j_250_75","t1t1~_448_270","t1t1~j_448_270","multijets"};


//const int nSims = 13;  string dataDir = ".";  string fileNames[nSims] = {"tt~","tt~j","tt~jj","tt~l","tt~lj","tt~ljj","t1t1~_183_10","t1t1~j_183_10","t1t1~_250_75","t1t1~j_250_75","t1t1~_448_270","t1t1~j_448_270","multijets400_7j_Brock"};






int treeColors[14] = {5,1,4,7,2,6,3,8,9,10,11,12,13,14};  
//int treeStyles[10] = {1,1,2,2,1,2,3,9,1,1};

//int treeColors[10] = {8,2,1,4,4,51,6,3,2,1};  
int treeStyles[14] = {1,1,1,1,1,1, 1,1,1,1,1,1,1,1};



// Event trees
TTree * trees[nSims];

// Tree variables
int    sim;               // sim #
int    event;             // event #
double weight;
// global event variables


double mt1;
double mt2;
double allpartpt;
double tquark1pt;
double tquark2pt;
double ISRpt;
double MHTpT;
double METpT;
double MR;
double MT_MHT;
double MT_MET;
double HT;
double MET_sqrtHT;
double MET_sumHTMET;
double MHT_sqrtHT;
double MHT_sumHTMHT;
double dPhiISRMET;
double dPhiISRMHT;
double dPhiMinMHT;
double dPhiMinMET;
double dPhiMax;
double dPhiJetMHT;
double dPhiMin2;
double dPhiMax2;
double dPhiJetMHT2;
double dPhiMin3;
double dPhiMax3;
double dPhiJetMHT3;

double dPhiMinISR;

double jet1pt;
double jet2pt;
double jet3pt;
double jet4pt;
double jet5pt;
double jet6pt;
double jet7pt;




struct particle {
  TLorentzVector v;
  int id;
  int idp;  // parent id
  int idgp; // grandparent id
  int ii;   // particle #
  int ip;   // parent particle #
  int igp;  // grandparent particle #
};

// random # generator
TRandom * randGen = new TRandom(12345);

// function for comparing Pt's of two 4-vectors (for sorting lists)
bool comparePt(particle p1, particle p2) { return (p1.v.Pt() > p2.v.Pt()); }

// function for getting charge sign of particles
double getCharge(int id)
{
  if (id == 0)  return 0;
  int absid = abs(id);
  int signid = id/absid;
  if (absid > 1000000) absid = absid%1000000;
  if (absid == 1 || absid == 3 || absid == 5)     return -signid*1./3.;
  if (absid == 2 || absid == 4 || absid == 6)     return  signid*2./3.;
  if (absid == 11 || absid == 13 || absid == 15)  return -signid;
  if (absid == 24 || absid == 37)                 return  signid;
  return 0;
}


// compare plots of variables from different trees
TH1F * comparisonHistos[nSims+1];
char lumiTitle[300];
int nBins =20;
const char * drawOption = "e";
bool logable = true;  // flag in how to set plot minima.  false => min=0, not suitable for log plotting!
bool normalize = false;   // normalize the histos to unity?
bool noTitle = false;
bool stack = false;  // draw as histo stack
bool fillin = false;   // fill in histos with simColors
bool bold = true;  // make the histos bold
bool legend = true;  // draw legend
bool sumBG = false;   // add up all "backgrounds", assumed to be everything after the first sample
bool quiet = false;  // suppress printout
TLegend * comparisonLegend = new TLegend(0.2,0.8,1.0,1.0);
//TLegend * comparisonLegend = new TLegend(0.65,0.9,1.0,1.0);
//TLegend * comparisonLegend = new TLegend(0.65,0.7,1.0,1.0);
THStack * comparisonStack = new THStack("comparisonStack","");
//void compare(string varKey, string cut = "1", double xmin = -1e12, double xmax = -1e12);
void compare(string varKey, string cut = "1", double xmin = -1e12, double xmax = -1e12);

// compare plots of different variables from the same tree
// uses some of the same global switches as the compare function
const int nVarsMax = 30;
TH1F * varComparisonHistos[nVarsMax];
string varKeys[nVarsMax];
TLegend * varComparisonLegend = new TLegend(0.52,0.57,0.87,0.87);
bool ignoreZeroBin = false;  // for normalization, ignores the zero bin in varCompare...for easily comparing variables multiplied by conditionals
void varCompare(int nVars, int tree = 0, string cut = "1", double xmin = -1e12, double xmax = -1e12);


// additional functions and supporting global variables
TLorentzVector lPlusT,lMinusT;  // transverse-plane Lorentz vectors of leptons
double PTmiss[2] = {0,0};  // missing transverse momentum (x,y)
void compute_mT2(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);  // compute mT2 of W subsystems
const double etaPhiSmear = 0.025;  // eta/phi smearing estimate from CMS particle-flow, for ~30 GeV jets (gets better at higher pt)
double jetResolution(double);
void printMeansAndMedians()
{
  double q[1];
  double prob[1] = {0.5};
  for (int i = 0; i < nSims; i++) {
    comparisonHistos[i]->GetQuantiles(1,q,prob);
    cout << "mean " << comparisonHistos[i]->GetMean() << ", median " << q[0] << endl;
  }
}








/*b tagging function.. takes in as arguments a collection of bjets as a vector of particle structs.. using a flat distribution function, a random number is chosen between 0 and 1.. the bjet is kept only if this random number falls below 0.75 MQP*/

//TF1 *flatfunc = new TF1("flatfunc","1.",0,1); //Create flat distribution function
/*
vector<particle> btageff(vector<particle> bOG) {
  vector<particle> effbtag; //Create a vector of particle structs that will be used to contain the surviving bjets MQP
  // vector<particle> failbtag;
 
  for (int i=0; i<(int)bOG.size(); i++) {

    double beff = flatfunc->GetRandom(); //Get random number between 0 and 1
    //cout << beff << "value of the flat function" << endl;
    if (beff<0.75) effbtag.push_back(bOG[i]); //If beff is below 0.75, keep the bjet
    else parjetsansb.push_back(bOG[i]);
  }

  //delete flatfunc;
  return effbtag;
  //return failbtag;

}
*/
/*jet smearing function.. takes in as arguments a collection of jets as a vector of particle structs as well as the average number of pile up events per crossing mu.. only mu values of 0, 50, or 140 will work.. otherwise an error message will come up and the original unsmeared jets will be returned.. MQP*/
/*
vector<particle> jetsmear(vector<particle> jetOG, int mu) {
  vector<particle> smearjets; //Create a vector of particle structs that will be used to contain the smeared jets MQP
  vector<double> param; //Create a vector of doubles that will be used to select the right smearing function depending on mu MQP

  if (mu==0) { param.push_back(0.37*0.37); param.push_back(0.0); param.push_back(0.02*0.02); } //Pushing back parameters for jet pT dependent smearing MQP
  if (mu==50) { param.push_back(0.55*0.55); param.push_back(8.2*8.2); param.push_back(0.02*0.02); }
  if (mu==140) { param.push_back(0.65*0.65); param.push_back(15.0*15.0); param.push_back(0.02*0.02); }

  if ((int)param.size()!=3) { cout << "invalid choice of mu, returning unsmeared jets" << endl; return jetOG; }

  if ((int)param.size()==3) {
    for (int i=0; i<(int)jetOG.size(); i++) { //Looping over jets MQP
      double ptj=jetOG[i].v.Pt(); 
      double fracres=sqrt(((param[0])/(ptj)) + ((param[1])/(ptj*ptj)) + param[2]); //Calculate jet pT dependent fractional smearing MQP
      double smearfac=randGen->Gaus(1,fracres); //Generate random number from Gaussian with fractional jet smearing MQP
      //double smearfac2=randGen->Gaus(1,fracres); //Generate random number from Gaussian with fractional jet smearing MQP


 TLorentzVector smearj(smearfac*jetOG[i].v); //Smearing the jet MQP
      particle smearp;  smearp.v = smearj;  smearp.id = jetOG[i].id;  smearp.idp = jetOG[i].idp;  smearp.idgp = jetOG[i].idgp,  smearp.ii = jetOG[i].ii;  smearp.ip = jetOG[i].ip;  smearp.igp = jetOG[i].igp; //Creating the particle struct for the smeared jet MQP
      smearjets.push_back(smearp); //Pushing back the smeared jets MQP


    }

    return smearjets;
  }
}
*/











//////////////////////////////////////////////////////////////////////////////
// MAIN function

void kinematics_june() {

  sprintf(lumiTitle,"# events (%i fb^{-1})",lumi);


 //btag test stuff
  
 // int bnumOG=0;
 // int bnumafter=0;

  for (int simIt = 0; simIt < nSims; simIt++) {
  
    // Define the tree
    char treeName[300];
    sprintf(treeName,"tree_%i",simIt);
    trees[simIt] = new TTree(treeName,"a tree");

    // Initialize the tree
    trees[simIt]->Branch("sim", &sim, "sim/I");
    trees[simIt]->Branch("event", &event, "event/I");
    trees[simIt]->Branch("weight", &weight, "weight/D");

    // global event variables

   
    trees[simIt]->Branch("allpartpt", &allpartpt, "allpartpt/D");
    trees[simIt]->Branch("mt1", &mt1, "mt1/D");
    trees[simIt]->Branch("mt2", &mt2, "mt2/D");
    trees[simIt]->Branch("MR", &MR, "MR/D");
    trees[simIt]->Branch("MT_MHT", &MT_MHT, "MT_MHT/D");
    trees[simIt]->Branch("MT_MET", &MT_MET, "MT_MET/D");
    trees[simIt]->Branch("ISRpt", &ISRpt, "ISRpt/D");
    trees[simIt]->Branch("METpT", &METpT, "METpT/D");
    trees[simIt]->Branch("MHTpT", &MHTpT, "MHTpT/D");
    trees[simIt]->Branch("HT", &HT, "HT/D");
    trees[simIt]->Branch("MET_sqrtHT", &MET_sqrtHT, "MET_sqrtHT/D");
    trees[simIt]->Branch("MET_sumHTMET", &MET_sumHTMET, "MET_sumHTMET/D");
    trees[simIt]->Branch("MHT_sqrtHT", &MHT_sqrtHT, "MHT_sqrtHT/D");
    trees[simIt]->Branch("MHT_sumHTMHT", &MHT_sumHTMHT, "MHT_sumHTMHT/D");

    trees[simIt]->Branch("dPhiISRMHT", &dPhiISRMHT, "dPhiISRMHT/D");
    trees[simIt]->Branch("dPhiISRMET", &dPhiISRMET, "dPhiISRMET/D");
    trees[simIt]->Branch("dPhiMinMHT", &dPhiMinMHT, "dPhiMinMHT/D");  
    trees[simIt]->Branch("dPhiMinMET", &dPhiMinMET, "dPhiMinMET/D");  
    trees[simIt]->Branch("dPhiMax", &dPhiMax, "dPhiMax/D");  
    trees[simIt]->Branch("dPhiJetMHT", &dPhiJetMHT, "dPhiJetMHT/D");    
    trees[simIt]->Branch("dPhiMin2", &dPhiMin2, "dPhiMin2/D");  
    trees[simIt]->Branch("dPhiMax2", &dPhiMax2, "dPhiMax2/D");  
    trees[simIt]->Branch("dPhiJetMHT2", &dPhiJetMHT2, "dPhiJetMHT2/D"); 
    trees[simIt]->Branch("dPhiMin3", &dPhiMin3, "dPhiMin3/D");  
    trees[simIt]->Branch("dPhiMax3", &dPhiMax3, "dPhiMax3/D");  
    trees[simIt]->Branch("dPhiJetMHT3", &dPhiJetMHT3, "dPhiJetMHT3/D");

    trees[simIt]->Branch("dPhiMinISR", &dPhiMinISR, "dPhiMinISR/D");  

    trees[simIt]->Branch("jet1pt", &jet1pt, "jet1pt/D");
    trees[simIt]->Branch("jet2pt", &jet2pt, "jet2pt/D");
    trees[simIt]->Branch("jet3pt", &jet3pt, "jet3pt/D");
    trees[simIt]->Branch("jet4pt", &jet4pt, "jet4pt/D");
    trees[simIt]->Branch("jet5pt", &jet5pt, "jet5pt/D");
    trees[simIt]->Branch("jet6pt", &jet6pt, "jet6pt/D");
    trees[simIt]->Branch("jet7pt", &jet7pt, "jet7pt/D");





    //trees[simIt]->Branch("", &, "/D");





    sim = simIt;

    // Open a stream connected to an event file:
    ifstream ifs;
    ifs.open((dataDir+"/"+fileNames[simIt]+".lhe").c_str());
    cout << endl << "///////////////////////////////////////////////////////////" << endl;
    cout << "//   " << dataDir+"/"+fileNames[simIt] << endl;
    if (!ifs.good()) { cout << "BAD LHE FILE" << endl;  continue; }  

    // Create the Reader object:
    LHEF::Reader reader(ifs);
    
    // Now loop over all events:
    long ieve = 0;
    bool done = false;
    done = !reader.readEvent();
    while ( !done ) {

      ++ieve;
      event = ieve;
      
      if (ieve%1000==0)  {
 	cout << "Event #" << ieve << endl;
      }

   
     
      int npart = reader.hepeup.NUP;

      weight = reader.hepeup.XWGTUP*1000*lumi;  // 1000 => pb->fb

      TLorentzVector  initial4vector;
      TLorentzVector  gen4vector;

      int lepCharge = 0;

      std::vector<TLorentzVector>ISRsum;
      std::vector<TLorentzVector>ISRsum1;
      std::vector<TLorentzVector>ISRsum2;
      
      vector<particle> ISR;

      vector<particle> metpar;
      vector<particle> sumetpar;
      vector<particle> parjetall;
      vector<particle> parjetsansb;
      vector<particle> parjetb;
      vector<particle> allpart;
      vector<particle> parlepton;
      vector<particle> pardetected;
      vector<particle> vp_MET;

      particle metpar_smear;
      particle METall_smear;

      vector<double> DRminFinal;
      vector<int> JindexFinal;
      vector<int> JindexFinala;
      vector<int> JindexFinalb;

      double DRmin = 999;
      double DRmin2 = 999;
     
      double indexJa = -9;
      double indexJa2 = -9;
      double indexJb = -9;
      double indexJb2 = -9;

      TLorentzVector tquark1;
      TLorentzVector tquark2;

      vector<particle> stuffMET;




      // Loop over particles in the event
      for (int i = 0; i < npart; i++) {
	vector<double> fiveVector = reader.hepeup.PUP[i];
	int id = reader.hepeup.IDUP[i];
	int stat = reader.hepeup.ISTUP[i];
	int ip = reader.hepeup.MOTHUP[i].first - 1;
	int idp = reader.hepeup.IDUP[ip];
	int igp = reader.hepeup.MOTHUP[ip].first - 1;
	int idgp = 0;

	if (ip >= 0 && igp >= 0)  idgp = reader.hepeup.IDUP[igp];
	
	TLorentzVector vect(fiveVector[0],fiveVector[1],fiveVector[2],fiveVector[3]);
	particle part;  part.v = vect;  part.id = id;  part.idp = idp;  part.idgp = idgp,  part.ip = ip;  part.igp = igp;
	if (stat == -1) {
	  initial4vector += vect;  
	  continue; 
	}   

	
	if ( stat==1 && (abs(id)==4|| abs(id)==5 || abs(id)==7 || abs(id)==11 || abs(id)==13 || abs(id)==15 )) {
	  pardetected.push_back(part); //Grabbing all detected particles SM
	}
	
	
	if ( stat==1  && (abs(id)==4|| abs(id)==5 ||abs(id)==7)) {
	  parjetall.push_back(part); //Grabbing all jets SM
	}

	if ( stat==1 && (abs(id)==11 || abs(id)==13 ||abs(id)==15 )) {
	  parlepton.push_back(part); //Grabbing all leptons SM
	}

	if ( (abs(id)==4|| abs(id)==7 ) && stat==1 ) {
	  parjetsansb.push_back(part);//Grabbing all the jets except the b ones SM
	}


	if ( stat==1 && abs(id)==5) {
	  parjetb.push_back(part); //Grabbing all bjets MQP
	}
      
	  if ( stat==1 && (abs(id)==12 || abs(id)==14 || abs(id)==16 || abs(id)==1000022 || abs(id)==1000023 || abs(id)==1000025 || abs(id)==1000035)) {
	    vp_MET.push_back(part); //Grabbing all neutrinos and neutralinos SM
	  }

	/*
	if (stat==1 &&( abs(id)==1000022 ||  abs(id)==12 || abs(id)==14 || abs(id)==16 ))
	  {
	    //METstuff=part;
	   stuffMET.push_back(part);
	  }
	*/ 

	
	if ( stat==1 && ((abs(id)<20 && abs(id)>0)|| abs(id)==1000022)) {
	  allpart.push_back(part); //Grabbing all particles  SM
	 
	}

      }  



      ///////////////////////////////////////////////////////////////////////////////////////////////////////////

      //Order particles from greater to lower Pt
      sort(parjetall.begin(), parjetall.end(), comparePt);
      sort(parjetb.begin(), parjetb.end(), comparePt);
     
      double mw=80.4;
      double mtop=173.1;
      
      //if((int)parjetall.size()<7)  cout << "fewer than 7 jets-> trouble" << endl;

      if((int)parjetsansb.size()<1)  cout << " jet trouble" << endl;

      if((int)parjetsansb.size()<1 && (int)parjetb.size()<1)  cout << " jet and b trouble" << endl;

      if((int)allpart.size()<1)  cout << " all part trouble" << endl;

      if((int)parjetsansb.size()<1 && (int)parjetb.size()>0)  cout << " no jets,  but b jets, trouble" << endl;
      
      if((int)parjetsansb.size()<1 && (int)allpart.size()>0)  cout << " no jets, but allpart  trouble" << endl;
      
      if((int)parjetsansb.size()<1 && (int)stuffMET.size()>0)  cout << " no jets but MET trouble" << endl;

      /*
      if (parlepton.size()>1){
      cout << parlepton.size() << "# of leptons " << endl;
      
      for (int i=0; i<(int)parlepton.size(); i++)
	{
	  
	  cout << parlepton[i].v.Pt() << "lepton pt" << endl;
	  
	}
      
      }

      */


      //TLorentzVectors
      TLorentzVector allparticles;
      
      for (int i=0; i<(int)allpart.size(); i++)
	{
	  
	  allparticles=allparticles+allpart[i].v;
	  
	}
      
      
      //Buckets o Tops
      
      vector<particle> parjetall_smear_no_ptjcut;
      vector<particle> pardetected_smear;
      vector<particle> parjetb_effed;
      vector<particle> parjetb_smear_no_ptjcut;
      vector<particle> parjetsansb_smear_no_ptjcut;
      
 
      /*
      if (simIt<2){parjetb_effed=parjetb; }
      

      if (simIt>1){
	// cout << (int)parjetb.size() << "bjets before" << endl; 
	//cout << (int)parjetsansb.size() << "non bjets before" << endl;

      
      for (int i=0; i<(int)parjetb.size(); i++) {
	
	double beff = flatfunc->GetRandom(); //Get random number between 0 and 1
	//cout << beff << "value of the flat function" << endl;
	if (beff<0.75) parjetb_effed.push_back(parjetb[i]); //If beff is below 0.75, keep the bjet
	else parjetsansb.push_back(parjetb[i]);//If beff>0.75 trow it to parjetsansb 
      }
      
      //cout << (int)parjetsansb.size() << "non bjets AFTER" << endl;
      // cout << (int)parjetb_effed.size() << "bjets  AFTER" << endl;

      //cout << "---------------------------------------------------------------------" << endl;
     
       }
*/
      ///////////////////////////// NO SMEARING AND NO B-TAGGING EFFICIENCY (ALREADY DONE WITH THE TRIGGER CODE) //////////////////////////////////////
      
      
      //Smeared jets and b-jets that I will use to do the ISR identification
      parjetb_effed=parjetb;
      parjetall_smear_no_ptjcut=parjetall; //Smearing the jets with mu=50 MQP
      pardetected_smear=pardetected; //Smearing the jets and leptons (electrons and muons) with mu=50 MQP. Is the lepton smearring the same as the jet smearring?
      // vector<particle> parjetb_effed=btageff(parjetb); //Applyi
     
      parjetb_smear_no_ptjcut=parjetb_effed; //Smearing the b jets with mu=50 SM
      parjetsansb_smear_no_ptjcut=parjetsansb; //Smearing all the jets except the b ones with mu=50 SM
      
      /*
      if (simIt>2){
      
      //Smeared jets and b-jets that I will use to do the ISR identification
      parjetall_smear_no_ptjcut=jetsmear(parjetall,50); //Smearing the jets with mu=50 MQP
      pardetected_smear=jetsmear(pardetected,50); //Smearing the jets and leptons (electrons and muons) with mu=50 MQP. Is the lepton smearring the same as the jet smearring?
      // vector<particle> parjetb_effed=btageff(parjetb); //Applyi
     
      parjetb_smear_no_ptjcut=jetsmear(parjetb_effed,50); //Smearing the b jets with mu=50 SM
      parjetsansb_smear_no_ptjcut=jetsmear(parjetsansb,50); //Smearing all the jets except the b ones with mu=50 SM
      }
      */


     
      

      // bnumOG=bnumOG+(int)parjetb.size();
      //bnumafter=bnumafter+(int)parjetb_effed.size();
      
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

      //I apply a ptjcut on jets after smearing

      vector<particle> parjetall_smear;
      
      for (int m=0; m<(int)parjetall_smear_no_ptjcut.size(); m++)
	{
	  if (parjetall_smear_no_ptjcut[m].v.Pt()> ptjcut) {parjetall_smear.push_back(parjetall_smear_no_ptjcut[m]);}
	}  
      
      

      vector<particle> parjetb_smear;
      
      for (int m=0; m<(int)parjetb_smear_no_ptjcut.size(); m++)
	{
	  if (parjetb_smear_no_ptjcut[m].v.Pt()> ptjcut) {parjetb_smear.push_back(parjetb_smear_no_ptjcut[m]);}
	}    

      
      vector<particle> parjetsansb_smear;
      
      for (int m=0; m<(int)parjetsansb_smear_no_ptjcut.size(); m++)
	{
	  if (parjetsansb_smear_no_ptjcut[m].v.Pt()> ptjcut) {parjetsansb_smear.push_back(parjetsansb_smear_no_ptjcut[m]);}
	}   
      
      /////////////////////////////////////////      SM          //////////////////////////////////////////////////////////////////////////


      //We calculate the missing HT
      TLorentzVector MHT;
      
      for (int i=0; i<(int)parjetall_smear.size(); i++)
	{
	  
	  MHT=MHT-parjetall_smear[i].v;
	  
	}
      particle MHTpar;  MHTpar.v = MHT;   //Creating the new particle struct SM

      ////////////////////////////////////////        SM       ///////////////////////////////////////////
      
      //We find the true MET as the sum of neutrinos and neutralinos

      TLorentzVector MET;
      
      for (int i=0; i<(int)vp_MET.size(); i++)
	{
	  
	  MET=MET+vp_MET[i].v;
	  
	}
      particle METpar;  METpar.v = MET;   //Creating the new particle struct SM
      //////////////////////////////////////////////////////////////////////////////////////////////
      // We order the jets from greater to lower pT

	
      sort(parjetall.begin(), parjetall.end(), comparePt);	    
      sort(parjetall_smear.begin(), parjetall_smear.end(), comparePt);
      sort(parjetb_smear.begin(), parjetb_smear.end(), comparePt);
      sort(parjetb.begin(), parjetb.end(), comparePt);
      sort(parjetsansb_smear.begin(), parjetsansb_smear.end(), comparePt);
      
   
      ///////////////////////////////////////////////////////////////////////////////////////////////

    

      
      vector<particle> ISRdef_smear;
      vector<particle> ISRdef;

      /*
      vector<particle> ISRb1hope_smear;
      vector<particle> ISRb2hope_smear;

      int numISRdefnot=0;
      */

      //if ((int)parjetb_smear.size()!=2) {cout<<"less than 2 b tags"<< endl; continue;}




 //ISR jet identification (for smeared jets)
      if ((int)parjetb_smear.size()>=2 && (int)parjetsansb_smear.size()>0) {
	for (int i=0; i<(int)parjetsansb_smear.size(); i++) {
	  if ((parjetsansb_smear[i].v+parjetb_smear[0].v).Mag()>(mtop*1.2) && (parjetsansb_smear[i].v+parjetb_smear[1].v).Mag()>(mtop*1.2)) ISRdef_smear.push_back(parjetsansb_smear[i]);

	}
      }

	
    

      sort(ISRdef_smear.begin(), ISRdef_smear.end(), comparePt);
      // sort(ISRdef.begin(), ISRdef.end(), comparePt);

     


      /////////////////////////////////////        SM      /////////////////////////////////////////////////////////////


      // Top quark 1 and 2  reconstruction  
      if ((int)ISRdef_smear.size()>0 && (int)parjetb_smear.size()>=2 && (int)parjetsansb_smear.size()>0 && (int)parjetall_smear.size()>6 && (int)parlepton.size()==0 ) {



	
	jet1pt=parjetall_smear[0].v.Pt();
	jet2pt=parjetall_smear[1].v.Pt();
	jet3pt=parjetall_smear[2].v.Pt();
	jet4pt=parjetall_smear[3].v.Pt();
	jet5pt=parjetall_smear[4].v.Pt();
	jet6pt=parjetall_smear[5].v.Pt();
	jet7pt=parjetall_smear[6].v.Pt();


	//	cout << jet7pt << "pt of the 7th hardest jet" << endl;

	//cout << parjetall[6].v.Pt() << "pt of the 7th hardest jet" << endl;

	//////////////////////////////////            SM        ////////////////////////////////////////////////////////
	//---------------------------------------------------------------------------------------------------------

	//I define the MHT variables

	MHTpT=MHT.Pt(); 

	//------------------------------------------------------------------------------------------------------------------
	//I define HT, MHT_sqrtHT and HT_sumHTMHT
	double ht=0;
	for (int i=0; i<(int)parjetall_smear.size(); i++){
	       ht=ht+fabs(parjetall_smear[i].v.Pt());
	       //cout << parjetall_smear[i].v.Pt() << "abs(PT)" << endl;
	     }

	       HT=ht;
	       //cout << HT << "Total abs(PT)" << endl;

	    MHT_sqrtHT=MHTpT/sqrt(HT);
	    MHT_sumHTMHT=MHTpT/(HT+MHTpT);

	    //  cout << METpt << "MHT" << endl;
	    //cout << MHT_sqrtHT << "MHT/sqrt(HT)" << endl;
	    //cout << MHT_sumHTMHT << " MHT/(HT+MHT)" << endl;
	//-------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------
	//I define the MET variables

	METpT=MET.Pt(); 

	//------------------------------------------------------------------------------------------------------------------
	//I define  MET_sqrtHT and HT_sumHTMET


	    MET_sqrtHT=METpT/sqrt(HT);
	    MET_sumHTMET=METpT/(HT+METpT);


	    ////////////////////////             SM               ///////////////////////////////////////////////////
	//-----------------------------------------------------------------------------------------------------------

	//I loop over the 3 hardest jets and find the closest one in Delta Phi to the MET and then plot Delta Phi between them

	double DPhi_minMET= 999;
	double indexJ_minMET= -9; 
	for (int j=0; j<3; j++)
	  {
	    //cout << j << "J" << endl;   
	    if (fabs(parjetall_smear[j].v.DeltaPhi(MET)) < DPhi_minMET) 
	      {
		DPhi_minMET =fabs(parjetall_smear[j].v.DeltaPhi(MET));
		indexJ_minMET = j;
	      }
	  }

	//if (indexJ_min==0){dPhiMinISR=fabs(parjetall_smear[indexJ_min].v.DeltaPhi(MHT));}


	dPhiMinMET=fabs(parjetall_smear[indexJ_minMET].v.DeltaPhi(MET));










	//------------------------------------------------------------------------------------------------------------------------
	//I loop over the 3 hardest jets and find the closest one in Delta Phi to the MHT and then plot Delta Phi between them

	double DPhi_min= 999;
	double indexJ_min= -9; 
	for (int j=0; j<3; j++)
	  {
	    //cout << j << "J" << endl;   
	    if (fabs(parjetall_smear[j].v.DeltaPhi(MHT)) < DPhi_min) 
	      {
		DPhi_min =fabs(parjetall_smear[j].v.DeltaPhi(MHT));
		indexJ_min = j;
	      }
	  }

	if (indexJ_min==0){dPhiMinISR=fabs(parjetall_smear[indexJ_min].v.DeltaPhi(MHT));}


	dPhiMinMHT=fabs(parjetall_smear[indexJ_min].v.DeltaPhi(MHT));

	//cout << indexJ << "indexJ" << endl;
	//cout << dPhiMin << "Delta Phi Min" << endl;

	//I loop over the 3 hardest jets and find the 2nd  closest one in Delta Phi to the MHT and then plot Delta Phi between them

	double DPhi_min2= 999;
	double indexJ_min2= -9; 
	for (int m=0; m<3; m++)
	  {
	    //cout << m << "m" << endl;   
	    if (parjetall_smear[m].v!=parjetall_smear[indexJ_min].v && fabs(parjetall_smear[m].v.DeltaPhi(MHT)) < DPhi_min2) 
	      {
		DPhi_min2 =fabs(parjetall_smear[m].v.DeltaPhi(MHT));
		indexJ_min2 = m;
	      }
	  }

	dPhiMin2=fabs(parjetall_smear[indexJ_min2].v.DeltaPhi(MHT));


	//I find the 3rd closest jet in Delta Phi to the MET (among the 3 hardest) and then plot Delta Phi between them

	double indexJ_min3=-9;
	for(int m=0; m<3; m++)
	  {
	    if (parjetall_smear[m].v!=parjetall_smear[indexJ_min].v && parjetall_smear[m].v!=parjetall_smear[indexJ_min2].v)
	      {indexJ_min3=m;}
	  }

	double DPhi_min3=fabs(parjetall_smear[indexJ_min3].v.DeltaPhi(MHT));
	dPhiMin3=fabs(parjetall_smear[indexJ_min3].v.DeltaPhi(MHT));

	////------------------------------------------------------------------------------------------------
	//I loop over the 3 hardest jets and find the furthest one in Delta Phi to the MHT and then plot Delta Phi between them

	double DPhi_max= -10;
	double indexJ_max= -9; 
	for (int j=0; j<3; j++)
	  {
	    //cout << j << "J" << endl;   
	    if (fabs(parjetall_smear[j].v.DeltaPhi(MHT)) > DPhi_max) 
	      {
		DPhi_max =fabs(parjetall_smear[j].v.DeltaPhi(MHT));
		indexJ_max = j;
	      }
	  }

	dPhiMax=fabs(parjetall_smear[indexJ_max].v.DeltaPhi(MHT));

	//I loop over the 3 hardest jets and find 2nd the furthest one in Delta Phi to the MET and then plot Delta Phi between them

	double DPhi_max2=-10;
	double indexJ_max2=-9;

	for (int m=0; m<3; m++)
	  {
	    //cout << m << "M" << endl;   
	    if (parjetall_smear[m].v!=parjetall_smear[indexJ_max].v && fabs(parjetall_smear[m].v.DeltaPhi(MHT)) > DPhi_max2) 
	      {
		DPhi_max2 =fabs(parjetall_smear[m].v.DeltaPhi(MHT));
		indexJ_max2 = m;
	      }
	  }

	dPhiMax2=fabs(parjetall_smear[indexJ_max2].v.DeltaPhi(MHT));

	//I find the furthest jet in Delta Phi to the MHT (among the 3 hardest) and then plot Delta Phi between them

       	double indexJ_max3=-9;
	for(int m=0; m<3; m++)
	  {
	    if (parjetall_smear[m].v!=parjetall_smear[indexJ_max].v && parjetall_smear[m].v!=parjetall_smear[indexJ_max2].v)
	      {indexJ_max3=m;}
	  }

	double DPhi_max3=fabs(parjetall_smear[indexJ_max3].v.DeltaPhi(MHT));
	dPhiMax3=fabs(parjetall_smear[indexJ_max3].v.DeltaPhi(MHT));

	//------------------------------------------------------------------------------------------------------
	//I define a variable where I plot either the closest or furthest jet to the MHT, picking between the one that is  closer to 0 (for the closest) or closer to Pi (for the furthest)

	if (fabs(DPhi_max - pi)< fabs(DPhi_min)) {dPhiJetMHT=fabs(parjetall_smear[indexJ_max].v.DeltaPhi(MHT));}
	else {dPhiJetMHT=fabs(parjetall_smear[indexJ_min].v.DeltaPhi(MHT));}
	
	
	if (fabs(DPhi_max2 - pi)< fabs(DPhi_min2)) {dPhiJetMHT2=fabs(parjetall_smear[indexJ_max2].v.DeltaPhi(MHT));}
	else {dPhiJetMHT2=fabs(parjetall_smear[indexJ_min2].v.DeltaPhi(MHT));}
	
  
	if (fabs(DPhi_max3 - pi)< fabs(DPhi_min3)) {dPhiJetMHT3=fabs(parjetall_smear[indexJ_max3].v.DeltaPhi(MHT));}
	else {dPhiJetMHT3=fabs(parjetall_smear[indexJ_min3].v.DeltaPhi(MHT));}
	
 

	//////////////////////////////////               SM                    //////////////////////////////////////////////////////////

	// I find the top quark 1 and 2 candidates
	
	//I pick the 2 b-tagged jets with highest pt and add the remaining b-tagged jets to the collection of jets
      	for (int m=2; m<(int)parjetb_smear.size(); m++){parjetsansb_smear.push_back(parjetb_smear[m]);
	  // cout << parjetb_smear[m].v.Pt() << "extra b jets pt" << endl;
}
      
      
      
	for (int m=0; m<(int)parjetsansb_smear.size(); m++)
	  { for (int k=0; k<(int)parjetsansb_smear.size(); k++)
	      if (parjetsansb_smear[m].v!=ISRdef_smear[0].v && parjetsansb_smear[k].v!=parjetsansb_smear[m].v && parjetsansb_smear[k].v!=ISRdef_smear[0].v && parjetsansb_smear[k].v.DeltaR(parjetsansb_smear[m].v) < DRmin) 
		{
		  DRmin =  parjetsansb_smear[k].v.DeltaR(parjetsansb_smear[m].v);
		  indexJa = m;
		  indexJb = k;
		  
		}
	  }
      
	
	DRminFinal.push_back(DRmin);  
	JindexFinala.push_back(indexJa);
	JindexFinalb.push_back(indexJb);



	for (int m=0; m<(int)parjetsansb_smear.size(); m++)
	  { for (int k=0; k<(int)parjetsansb_smear.size(); k++)
	      if ( parjetsansb_smear[JindexFinala[0]].v!=parjetsansb_smear[k].v && parjetsansb_smear[JindexFinalb[0]].v!=parjetsansb_smear[k].v && 
parjetsansb_smear[JindexFinala[0]].v!=parjetsansb_smear[m].v && parjetsansb_smear[JindexFinalb[0]].v!=parjetsansb_smear[m].v &&
		   parjetsansb_smear[m].v!=ISRdef_smear[0].v &&/* k!=m && */ parjetsansb_smear[k].v!=parjetsansb_smear[m].v &&  parjetsansb_smear[k].v!=ISRdef_smear[0].v && parjetsansb_smear[k].v.DeltaR(parjetsansb_smear[m].v) < DRmin2) 
		{
		  DRmin2 =  parjetsansb_smear[k].v.DeltaR(parjetsansb_smear[m].v);
		  indexJa2 = m;
		  indexJb2 = k;
		  
		}
	  }

	DRminFinal.push_back(DRmin2);  
	JindexFinala.push_back(indexJa2);
	JindexFinalb.push_back(indexJb2);

	
	//	for (int m=0; m<(int)JindexFinala.size(); m++){ cout << JindexFinala[m] << "jet 1 index" << endl;}
	//	cout <<  "------------------------------------" << endl;
	

	double DRwJetb1=(parjetsansb_smear[JindexFinala[0]].v+parjetsansb_smear[JindexFinalb[0]].v).DeltaR(parjetb_smear[0].v);
	double DRwJetb2=(parjetsansb_smear[JindexFinala[0]].v+parjetsansb_smear[JindexFinalb[0]].v).DeltaR(parjetb_smear[1].v);
	  
	if (DRwJetb1< DRwJetb2)
	  {

	    tquark1=parjetsansb_smear[JindexFinala[0]].v+parjetsansb_smear[JindexFinalb[0]].v+parjetb_smear[0].v;
	    mt1=tquark1.Mag();

	    tquark2=parjetsansb_smear[JindexFinala[1]].v+parjetsansb_smear[JindexFinalb[1]].v+parjetb_smear[1].v;
	    mt2=tquark2.Mag();
	  }
	else 
	  {
	    tquark1=parjetsansb_smear[JindexFinala[0]].v+parjetsansb_smear[JindexFinalb[0]].v+parjetb_smear[1].v;
	    mt1=tquark1.Mag();
	    
	    tquark2=parjetsansb_smear[JindexFinala[1]].v+parjetsansb_smear[JindexFinalb[1]].v+parjetb_smear[0].v;
	    mt2=tquark2.Mag();

	  }

	/*
	if (mt2<0){
	  cout <<  tquark2.E()<< "energy" << endl;
	  cout <<  tquark2.Px()<< "px" << endl;
	  cout <<  tquark2.Py()<< "py" << endl;
	  cout <<  tquark2.Pz()<< "pz" << endl;
	}
	*/

	/*
	//I find the angle between the MET and each of the jets that I use to reconstruct the top quarks
	dPhiJa1MET=fabs(parjetsansb_smear[JindexFinala[0]].v.DeltaPhi(METall));

	dPhiJa2MET=fabs(parjetsansb_smear[JindexFinala[1]].v.DeltaPhi(METall));
	
	dPhiJb1MET=fabs(parjetsansb_smear[JindexFinalb[0]].v.DeltaPhi(METall));
	
	dPhiJb2MET=fabs(parjetsansb_smear[JindexFinalb[1]].v.DeltaPhi(METall));

	*/



	//	cout << tquark1.E() << "tquark1 energy" << endl;
	
        //Absolute values of the momentum of tquark1 and tquark2



       ////////////////////////////              SM                 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      
       //I define Razor variables using the top reconstructed candidates

	MR=2*sqrt(pow(((tquark1.E())*tquark2.Pz()-(tquark2.E())*tquark1.Pz()),2)/(pow((tquark1.Pz()-tquark2.Pz()),2)-pow(((tquark1.E())-(tquark2.E())),2))); // Razor mass of eq. (5) of 1006.2727


      
	MT_MET=sqrt(0.5*((sqrt(pow(MET.Px(),2)+pow(MET.Py(),2)))*((sqrt(pow(tquark1.Px(),2)+pow(tquark1.Py(),2)))+(sqrt(pow(tquark2.Px(),2)+pow(tquark2.Py(),2))))-MET.Px()*(tquark1.Px()+tquark2.Px())-MET.Py()*(tquark1.Py()+tquark2.Py())));//Razor variable of eq. (6) of 1310.4827 for MET

      	
	MT_MHT=sqrt(0.5*((sqrt(pow(MHT.Px(),2)+pow(MHT.Py(),2)))*((sqrt(pow(tquark1.Px(),2)+pow(tquark1.Py(),2)))+(sqrt(pow(tquark2.Px(),2)+pow(tquark2.Py(),2))))-MHT.Px()*(tquark1.Px()+tquark2.Px())-MHT.Py()*(tquark1.Py()+tquark2.Py())));//Razor variable of eq. (6) of 1310.4827 for MHT

	//////////////////////////////////////////////////////////////////////////////////////////


	/*
      dPhijet1MET=fabs(parjetall_smear[0].v.DeltaPhi(METall));    
      dPhijet2MET=fabs(parjetall_smear[1].v.DeltaPhi(METall));    
      dPhijet3MET=fabs(parjetall_smear[2].v.DeltaPhi(METall));    
	*/
      
      //if(dPhijet1MET>0.628319 && dPhijet2MET>0.628319 && dPhijet3MET>0.628319 ){
      
      ISRpt=ISRdef_smear[0].v.Pt();
      
      
      
      //ISRpx=ISRdef_smear[0].v.Px();
      
      //ISRpx0=ISRdef[0].v.Px();
      
      // ISRpy=ISRdef_smear[0].v.Py();
      /*
      jet1pt=parjetall_smear[0].v.Pt();
      jet2pt=parjetall_smear[1].v.Pt();
      jet3pt=parjetall_smear[2].v.Pt();
      */
      /*
	cout <<(jet1pt) << " jet1pt ";//
	cout <<(jet2pt) << " jet2pt ";//
	cout <<(jet3pt) << " jet3pt ";//
      */
      

     
      
      //Et=sumetpar[0].v.E();
      
      
      allpartpt=allparticles.Pt();
      


      //dPhi btw ISR smear and the MHT or MET
      dPhiISRMHT=fabs(ISRdef_smear[0].v.DeltaPhi(MHT));
      dPhiISRMET=fabs(ISRdef_smear[0].v.DeltaPhi(MET));

      // dPhiMETall_smearMETall=fabs(METall_smear.v.DeltaPhi(METall));

	
      // Fill tree
      trees[simIt]->Fill();
      /*
      hjet1pt->Fill(jet1pt);
      hjet2pt->Fill(jet2pt);
      hjet3pt->Fill(jet3pt);
      hjet4pt->Fill(jet4pt);
      hjet5pt->Fill(jet5pt);
      hjet6pt->Fill(jet6pt);
      hjet7pt->Fill(jet7pt);
      

      if ( ISRpt>350 ){
      if (simIt==0){
      hjet1pt->Fill(MT);
      //hjet2pt->Fill(dPhiMin2);
      //hjet3pt->Fill(dPhiMin3);
      }
      }

   if ( ISRpt<80){
      if (simIt==0){
      hjet3pt->Fill(MT);
      //hjet2pt->Fill(dPhiMin2);
      //hjet3pt->Fill(dPhiMin3);
      }
   }
      */
      /*
      if (simIt==1){
 
      hjet3pt->Fill(MT);
      //hjet5pt->Fill(dPhiMin2);
      //hjet6pt->Fill(dPhiMin3);
      }
      }
      */
      /*
      if ( ISRpt<150 && ISRpt>50 && dPhiISRMET>3){
      if (simIt==0){
      hjet2pt->Fill(MT);
      //hjet2pt->Fill(dPhiMin2);
      //hjet3pt->Fill(dPhiMin3);
      }
      if (simIt==1){
 
      hjet4pt->Fill(MT);
      //hjet5pt->Fill(dPhiMin2);
      // hjet6pt->Fill(dPhiMin3);
      }
      }
      */
 /*
	if (ISRpt>0){
	if (simIt==1){ dPhiISRMet1->Fill(dPhiMETISRsmear);}

	if (simIt==1){ dPhiISRMet4->Fill(dPhiISRdef_smearMETall_smear);}

	
	}

	  */



	 }

      // }
      
  

  

    




	
     

      done = !reader.readEvent();
 
      //if (ieve > 10000)  break;///////




    }

  }


  //trees[0]->StartViewer();
  
  /*
  dPhiISRMet1->Scale(7.095*1000.*300./46000.);

  
  dPhiISRMet4->Scale(7.095*1000.*300./46000.);
  */
  /*
  jet1pt->Scale(671.513*1000.*300/97355);
  jet2pt->Scale(671.513*1000.*300/97355);
  jet3pt->Scale(671.513*1000.*300/97355);
  jet4pt->Scale(671.513*1000.*300/97355);
  jet5pt->Scale(671.513*1000.*300/97355);
  jet6pt->Scale(671.513*1000.*300/97355);
  jet7pt->Scale(671.513*1000.*300/97355);
  */

  // TCanvas* c1 = new TCanvas ("c1", "Test Plots", 600, 600);
  //c1->SetFillColor(kWhite);
  // hjet7pt->Draw();
  //hjet7pt->SetLineColor(1);

  //  hjet6pt->Draw();
  //hjet6pt->SetLineColor(2);
  
  // hjet5pt->Draw("same");
  //hjet5pt->SetLineColor(2);
  
  // hjet3pt->Draw();
  //hjet3pt->SetLineColor(1);
  
  //hjet1pt->Draw("same");
  //hjet1pt->SetLineColor(2);

  // hjet3pt->Draw("same");
  //hjet3pt->SetLineColor(3);



  // hjet4pt->Draw("same");
  //hjet4pt->SetLineColor(4);  


  

  

  /*
  TCanvas* c2 = new TCanvas ("c2", "Test Plots", 600,600);
  //c1->SetFillColor(kWhite);
  dPhiISRMet->Draw();
  //dPhiISRMet4j5->SetLineColor(kRed);
  //dPhiISRMet4j5->Draw("same");
  
*/

}








///////////////////////////////////////////////////////////////////////////////
// Additional functions

// calculate fractional energy smearing for a (sub)jet
double jetResolution(double Ejet)
{

  // Resolution constants fudged from CMS PF 09-001, fig 9.  This overshoots a bit below 50 GeV, undershoots a bit at 100, and overshoots 
  // a bit again asymptotically in the 100's.  ("A bit" here means 10~20%, so nothing to be too concerned about.)  The overshoot becomes
  // most obvious at the lowest pt's on CMS's plot, ~15 GeV, where they get 15% and I get 21%.
  const double a = 0;    // electronics/physics noise (we technically automatically incorporate physics noise, so this is hopefully conservative)
  const double b = 0.8;   // sampling fluctuations
  const double c = 0.04;  // nonlinearities/nonuniformities
  const double Emin = 20;

  if (Ejet == 0)  return 0;
  if (Ejet < Emin)  return sqrt(a*a/Emin/Emin+b*b/Emin+c*c); // regulate low-pt
  return sqrt(a*a/Ejet/Ejet + b*b/Ejet + c*c);
}


// Minuit function for calculating mt2 configuration of WW system
void compute_mT2(Int_t &npar, Double_t * gin, Double_t & f, Double_t * par, Int_t flag)
{

  // 2+1 vector of neutrino associated with lepton
  double vPlus_px = par[0];
  double vPlus_py = par[1];
  TLorentzVector vPlusT(vPlus_px,vPlus_py,0,0);  vPlusT.SetE(vPlusT.Pt());

  // 2+1 vector of neutrino associated with antilepton
  double vMinus_px = PTmiss[0]-vPlus_px;
  double vMinus_py = PTmiss[1]-vPlus_py;
  TLorentzVector vMinusT(vMinus_px,vMinus_py,0,0);  vMinusT.SetE(vMinusT.Pt());

  // transverse masses of systems a and b, combined with their candidate neutrinos
  double mTb = (lPlusT+vPlusT).M();
  double mTa = (lMinusT+vMinusT).M();

  f = max(mTb,mTa);
  //cout << f << endl;

  return;
}







// compare plots of variables from different trees
// 
// * 29 May 2010
//   . tried to export to a separate ".h" file for portability, but failed.  main can find the function, but CINT can't.
//     this seems to be a general problem with anything defined via header (simple global variables, structs, classes, functions).  
//     there must be a workaround to this that I don't know of...
//   . made tempHisto fully local
//   . made definitions of comparisonHistos local
//   . fixed normalized histos...bins are now normalized wrt the full variable range instead of display range
//   . can now display negative ranges without problems
//   . changed output to reflect # weighted and unweighted events, inside the requested display range and total passing cuts
//   . add background summation option
//   . assorted other minor changes to comments, etc
// * 10 Aug 2010
//   . first version of varCompare, largely mimicking functionality of 29 May 2010 version of compare
// * 10 Sept 2010
//   . added option for line styles
// * 26 Sept 2010
//   . added option for bold lines
// * 10 Dec 2010
//   . fixed minor bug:  check for 0 integral before attempting to normalize
// * 29 Dec 2010
//   . varCompare: added "ignoreZeroBin" option, which ignores the 0-bin for normalizing...for comparing shapes of variables multiplying conditionals
// * 18 Jan 2011
//   . added "quiet" option
// * 10 Feb 2011
//   . switched from char* to c++ string (way fewer compiler warnings now!  also much more graceful handling of varKeys)
// * 24 Aug 2011
//   . sumw2 activated to get correct bin errors, "HIST" added to draw option to auto-suppress drawing error bars (unless we add "e" to drawOption)
// * 29 Dec 2011
//   . propagated "HIST" modification to stack drawing.  undid a hack where the first histo in a stack is filled with grey.
//
void compare(string varKey, string cut, double xmin, double xmax)
{
  
  int nTrees = nSims;  // in case a different global label for # different physics processes
  TTree ** allTrees = trees;
  string * physics = fileNames;
  
  // first, fill temporary histo to find variable range and overall normalization
  // start off with some appropriately small range and a large # of bins (tempHisto is set rebinnable)
  string weightedCut = "weight*("+cut+")";
  int nPass[nTrees];          // # of tree events passing cut, to monitor statistics
  double nWeightedPass[nTrees];   // weighted # events passing
  double nWeightedPassBG = 0;     // "background" total
  double varMin = 1e12;
  double varMax = -1e12;
  for (int tr = 0; tr < nTrees; tr++) {
    TH1F * tempHisto = new TH1F("tempHisto","",10000,0,0.1);
    tempHisto->SetBit(TH1F::kCanRebin);
    allTrees[tr]->Draw((varKey+" >> "+tempHisto->GetName()).c_str(),weightedCut.c_str());
    nPass[tr] = tempHisto->GetEntries();
    nWeightedPass[tr] = tempHisto->Integral();
    if (tr > 0) nWeightedPassBG += nWeightedPass[tr];
    for (int binIt = 1; binIt <= tempHisto->GetNbinsX(); binIt++) {
      double binContent  = tempHisto->GetBinContent(binIt);
      double binLowEdge  = tempHisto->GetBinLowEdge(binIt);
      double binHighEdge = tempHisto->GetBinWidth(binIt) + binLowEdge;
      if (binContent == 0)  continue;
      if (binLowEdge  < varMin)  varMin = binLowEdge;
      if (binHighEdge > varMax)  varMax = binHighEdge;
    }
    delete tempHisto;
  }
  
  // add some padding at the edges
  double varRange = varMax-varMin;
  varMin -= varRange*0.00;
  varMax += varRange*0.001;  // ROOT sometimes chops off events near the end

  // or reset to user binning
  if (xmin > -1e10) varMin = xmin;
  if (xmax > -1e10) varMax = xmax;

  // prepare comparison histos (including one extra for summed "backgrounds")
  for (int tr = 0; tr < nTrees+1; tr++) {
    if (comparisonHistos[tr] == NULL) {  // define the comparison histos if this is the first time
      char histoName[100];
      sprintf(histoName,"comparisonHisto%i",tr);
      comparisonHistos[tr] = new TH1F(histoName,"",nBins,0,1);  // range will be recomputed later
      //if (tr==1)  comparisonHistos[tr]->SetLineColor(17);
      comparisonHistos[tr]->SetMarkerColor(treeColors[tr]);
      comparisonHistos[tr]->SetStats(false);
      comparisonHistos[tr]->Sumw2();
    }
    comparisonHistos[tr]->Reset();
    comparisonHistos[tr]->SetBins(nBins,varMin,varMax);
    comparisonHistos[tr]->SetLineWidth(1);
    if (bold) comparisonHistos[tr]->SetLineWidth(2);
    if (!sumBG)                       comparisonHistos[tr]->SetLineStyle(treeStyles[tr]);
    else if (tr > 0 && tr < nTrees)   comparisonHistos[tr]->SetLineStyle(2);  // individual BG components are dashed
  }
     











  comparisonHistos[0]->SetLineColor(1);
  
  comparisonHistos[1]->SetLineColor(2); 
  
  comparisonHistos[2]->SetLineColor(4);
  comparisonHistos[3]->SetLineColor(6);
  comparisonHistos[4]->SetLineColor(28);
  comparisonHistos[5]->SetLineColor(12);
  


  comparisonHistos[6]->SetLineColor(3);
  comparisonHistos[7]->SetLineColor(5);
  
  comparisonHistos[8]->SetLineColor(46);
  comparisonHistos[9]->SetLineColor(7);
  //comparisonHistos[10]->SetLineColor(7);
  
  //comparisonHistos[11]->SetLineColor(9);
  //comparisonHistos[12]->SetLineColor(8);
  //comparisonHistos[13]->SetLineColor(2);


    /*
  
  comparisonHistos[0]->SetLineColor(1);
  
  comparisonHistos[1]->SetLineColor(2); 
  
  comparisonHistos[2]->SetLineColor(4);
  comparisonHistos[3]->SetLineColor(46);
  comparisonHistos[4]->SetLineColor(28);
  comparisonHistos[5]->SetLineColor(12);
  


  comparisonHistos[6]->SetLineColor(3);
  comparisonHistos[7]->SetLineColor(4);
  
  comparisonHistos[8]->SetLineColor(5);
  comparisonHistos[9]->SetLineColor(6);
  comparisonHistos[10]->SetLineColor(7);
  
  comparisonHistos[11]->SetLineColor(9);
  comparisonHistos[12]->SetLineColor(8);
  comparisonHistos[13]->SetLineColor(2);

    
  comparisonHistos[6]->SetLineWidth(3);
  comparisonHistos[7]->SetLineWidth(3);
  comparisonHistos[8]->SetLineWidth(3);
  comparisonHistos[9]->SetLineWidth(3);
  comparisonHistos[10]->SetLineWidth(3);
  comparisonHistos[11]->SetLineWidth(3);
  comparisonHistos[12]->SetLineWidth(3);
  comparisonHistos[13]->SetLineWidth(3);

 comparisonHistos[6]->SetLineStyle(2);
  comparisonHistos[7]->SetLineStyle(2);
  comparisonHistos[8]->SetLineStyle(2);
  comparisonHistos[9]->SetLineStyle(2);
  comparisonHistos[10]->SetLineStyle(2);
  comparisonHistos[11]->SetLineStyle(2);
  comparisonHistos[12]->SetLineStyle(2);
  comparisonHistos[13]->SetLineStyle(2);

    */

  //comparisonHistos[6]->SetLineWidth(1);
  //comparisonHistos[7]->SetLineWidth(1);
  //comparisonHistos[8]->SetLineWidth(1);
  //comparisonHistos[9]->SetLineWidth(1);
  //comparisonHistos[10]->SetLineWidth(1);
  //comparisonHistos[11]->SetLineWidth(1);
  // comparisonHistos[12]->SetLineWidth(1);
  //comparisonHistos[13]->SetLineWidth(1);

  //comparisonHistos[0]->SetLineWidth(4);
  //comparisonHistos[1]->SetLineWidth(4);
  //comparisonHistos[2]->SetLineWidth(4);
  //comparisonHistos[3]->SetLineWidth(4);
  //comparisonHistos[4]->SetLineWidth(4);
  //comparisonHistos[5]->SetLineWidth(4);
 
  /*
  comparisonHistos[0]->SetLineStyle(7);
  comparisonHistos[1]->SetLineStyle(3);
  comparisonHistos[2]->SetLineStyle(3);
  comparisonHistos[3]->SetLineStyle(9);
  comparisonHistos[4]->SetLineStyle(3);
  comparisonHistos[5]->SetLineStyle(3);
  */



  
  // fill the comparison histos
  double maxBinContent = 0;
  double minBinContent = 1e12;
  int nDisplayed[nTrees];
  double nWeightedDisplayed[nTrees];
  for (int tr = 0; tr < nTrees; tr++) {
    allTrees[tr]->Draw((varKey+" >> "+comparisonHistos[tr]->GetName()).c_str(),weightedCut.c_str());
    nDisplayed[tr] = comparisonHistos[tr]->GetEntries();
    nWeightedDisplayed[tr] = comparisonHistos[tr]->Integral();
    //if (tr > 0)  comparisonHistos[nTrees]->Add(comparisonHistos[tr]);  // summed "backgrounds"
    if (!quiet) {
      cout << nWeightedDisplayed[tr] << " / " << nWeightedPass[tr] << " weighted    " 
	   << nDisplayed[tr]         << " / " << nPass[tr]         << " actual      ";
      if (nWeightedPass[tr] > 0)  cout << "(" <<  nWeightedDisplayed[tr]/nWeightedPass[tr] << ")" << endl;
      else                        cout << endl;
    }
    if (normalize && nWeightedPass[tr] > 0)  comparisonHistos[tr]->Scale(1./(nWeightedPass[tr]));  // if requested, rescale all histos to unity
    for (int binIt = 1; binIt <= comparisonHistos[tr]->GetNbinsX(); binIt++) {
      if (comparisonHistos[tr]->GetBinContent(binIt) < minBinContent &&
	  comparisonHistos[tr]->GetBinContent(binIt) > 0)
	minBinContent = comparisonHistos[tr]->GetBinContent(binIt);
      if (comparisonHistos[tr]->GetBinContent(binIt) > maxBinContent)
	maxBinContent = comparisonHistos[tr]->GetBinContent(binIt);
    }
  }


 
  //--------------------------------------------------------------------------------------------------------


  
  
//Find bin integral(To check if the cross sections are being read correctly)
/*  
   cout << "------------------------------------------------------------------------" << endl;

  cout <<nWeightedDisplayed[0] << " tt~ inclusive bin integral ";//  
  cout <<nWeightedDisplayed[1] << " tt~l inclusive bin integral ";//
  cout <<nWeightedDisplayed[2] << " bbjj-ccjj bin integral ";//
  cout <<nWeightedDisplayed[3] << " bbjz-ccjz bin integral ";//  
  cout <<nWeightedDisplayed[4] << " bbjw-ccjw bin integral ";//
  cout <<nWeightedDisplayed[5] << " t1t1~_183_10 inclusive bin integral ";//
  cout <<nWeightedDisplayed[6] << " t1t1~_250_75 inclusive bin integral ";//
  cout <<nWeightedDisplayed[7] << " t1t1~_448_270 inclusive bin integral ";// 
*/ 
 /*
 cout <<nWeightedDisplayed[8] << " t1t1~_250_75 bin integral ";//
  cout <<nWeightedDisplayed[9] << " t1t1~j_250_75 bin integral ";//
  cout <<nWeightedDisplayed[10] << " t1t1~_448_270 bin integral ";//
  cout <<nWeightedDisplayed[11] << " t1t1~j_448_270 bin integral ";//
  cout <<nWeightedDisplayed[12] << " multijets bin integral ";//
  */
  cout << "------------------------------------------------------------------------" << endl;


 //Ratios signal/Background

 /*
  cout <<nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0] << " tt~jj + tt~j + tt~  bin integral ";//
  cout <<nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3] << " tt~ljj + tt~lj + tt~l  bin integral ";//
  cout <<nWeightedDisplayed[6]+nWeightedDisplayed[7] << " t1t1~ + t1t1~j 183_10  bin integral ";//
  cout <<nWeightedDisplayed[8]+nWeightedDisplayed[9] << " t1t1~ + t1t1~j 250_75  bin integral ";//
  cout <<nWeightedDisplayed[10]+nWeightedDisplayed[11] << " t1t1~ + t1t1~j 448_270  bin integral ";//

 */
  cout << "------------------------------------------------------------------------" << endl;


  cout <<(nWeightedDisplayed[6])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_t1t1j/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;
  

 cout <<(nWeightedDisplayed[7])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_t1t1jj/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;

  /*
 cout <<(nWeightedDisplayed[8])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_250_75/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;


 cout <<(nWeightedDisplayed[9])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_275_100/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;



 cout <<(nWeightedDisplayed[10])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_300_125/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;



 cout <<(nWeightedDisplayed[11])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_350_175/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;


 cout <<(nWeightedDisplayed[12])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_400_225/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;
 

 cout <<(nWeightedDisplayed[13])/(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_448_270/(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;


*/
  cout << "------------------------------------------------------------------------" << endl;
  cout << "------------------------------------------------------------------------" << endl;
  
 





  cout <<(nWeightedDisplayed[6])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_t1t1j/sqrt(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;

  cout <<(nWeightedDisplayed[7])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_t1t1jj/sqrt(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;
  /*
  cout <<(nWeightedDisplayed[8])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_250_75/sqrt(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;

  cout <<(nWeightedDisplayed[9])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_275_100/sqrt(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;

  cout <<(nWeightedDisplayed[10])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_300_125/sqrt(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;

  cout <<(nWeightedDisplayed[11])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_350_175/sqrt(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;

  cout <<(nWeightedDisplayed[12])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_400_225/sqrt(bg) ";//
  cout << "------------------------------------------------------------------------" << endl;

  

  cout <<(nWeightedDisplayed[13])/sqrt(nWeightedDisplayed[5]+nWeightedDisplayed[4]+nWeightedDisplayed[3]+nWeightedDisplayed[2]+nWeightedDisplayed[1]+nWeightedDisplayed[0]) << " Ratio sig_448_270/sqrt(bg) ";//
*/
  cout << "------------------------------------------------------------------------" << endl;
  
  cout << "------------------------------------------------------------------------" << endl;

   


 



  // normalize summed "backgrounds" if requested, and get max/min bin contents
  if (normalize && nTrees > 1 && nWeightedPassBG > 0) {
    comparisonHistos[nTrees]->Scale(1./nWeightedPassBG);
  }
  if (sumBG) {
    for (int binIt = 1; binIt <= comparisonHistos[nTrees]->GetNbinsX(); binIt++) {
      if (comparisonHistos[nTrees]->GetBinContent(binIt) < minBinContent &&
	  comparisonHistos[nTrees]->GetBinContent(binIt) > 0)
	minBinContent = comparisonHistos[nTrees]->GetBinContent(binIt);
      if (comparisonHistos[nTrees]->GetBinContent(binIt) > maxBinContent)
	maxBinContent = comparisonHistos[nTrees]->GetBinContent(binIt);
    }
  }
  


  // adjust vertical range, etc, and draw
  comparisonLegend->Clear();
  for (int tr = 0; tr < nTrees+1; tr++) {  
    if (!sumBG && tr == nTrees)  continue;
    if (logable)   comparisonHistos[tr]->SetMinimum(minBinContent*0.5);
    else           comparisonHistos[tr]->SetMinimum(0);
    if (logable)   comparisonHistos[tr]->SetMaximum(maxBinContent*2);
    else           comparisonHistos[tr]->SetMaximum(maxBinContent*1.05);
    comparisonHistos[tr]->GetXaxis()->SetTitle(varKey.c_str());
    if (normalize)   comparisonHistos[tr]->GetYaxis()->SetTitle("normalized rate");
    //comparisonHistos[tr]->GetYaxis()->SetRange(0,0.7);
    // comparisonHistos[tr]->GetYaxis()->SetRangeUser(0, 1);
    else             comparisonHistos[tr]->GetYaxis()->SetTitle(lumiTitle);
    comparisonHistos[tr]->GetYaxis()->SetTitleOffset(1.2);
    comparisonHistos[tr]->GetYaxis()->SetLabelSize(0.035);
    // if (normalize)  comparisonHistos[tr]->Scale(1./comparisonHistos[tr]->Integral()); 
    // if (normalize)  comparisonHistos[tr]->GetYaxis()->SetRangeUser(0,0.1);
    string histoTitle = varKey;
    if (!(atoi(cut.c_str()) > 0))  histoTitle += " ("+cut+")";
    comparisonHistos[tr]->SetTitle(histoTitle.c_str());
    if (noTitle)  comparisonHistos[tr]->SetTitle("");
    string drawOptionFull = drawOption;
    drawOptionFull = "HIST " + drawOptionFull;
    if (tr != 0)  drawOptionFull += " same";
    if (fillin) {
      comparisonHistos[tr]->SetFillStyle(1001);
      comparisonHistos[tr]->SetFillColor(treeColors[tr]);
      if (treeColors[tr]==1)  comparisonHistos[tr]->SetFillColor(17);
      //if (tr == 0)  comparisonHistos[tr]->SetFillColor(17);
    }
    else  {
      comparisonHistos[tr]->SetFillStyle(0);
      comparisonHistos[tr]->SetFillColor(10);
    }
    if (!stack)  comparisonHistos[tr]->Draw(drawOptionFull.c_str());
    if (tr < nTrees){
      comparisonLegend->AddEntry(comparisonHistos[tr],physics[tr].c_str());
    }
    else if (tr == nTrees)
      comparisonLegend->AddEntry(comparisonHistos[nTrees],"summed backgrounds");
  }
  
  /*
  if (stack) {
    comparisonStack = new THStack("comparisonStack","");  // memory leak?  what memory leak?  (gotta keep redefining stack cuz it doesn't register changes to histos)
    for (int tr = nTrees-1; tr >= 0; tr--) {
       comparisonStack->Add(comparisonHistos[tr]);
    }
    comparisonStack->Draw("HIST");
    comparisonStack->GetXaxis()->SetTitle( comparisonHistos[0]->GetXaxis()->GetTitle() );
    comparisonStack->GetYaxis()->SetTitle( comparisonHistos[0]->GetYaxis()->GetTitle() );
    comparisonStack->GetXaxis()->SetTitleOffset(1.2);
    comparisonStack->GetYaxis()->SetLabelSize(0.035);
  }
  */
  if (legend) {
    comparisonLegend->Draw();
    comparisonLegend->SetBorderSize(0);
    comparisonLegend->SetFillColor(10);
  }

  
  return;
}


// compare plots of different variables from the same tree
// 
//   see "version information" above compare
//
void varCompare(int nVars, int tree, string cut, double xmin, double xmax)
{

  TTree ** allTrees = trees;

  if (!quiet) {
    for (int i = 0; i < nVars; i++) {
      cout << varKeys[i] << endl;
    }
  }

  // first, fill temporary histo to find variable range and overall normalization
  // start off with some appropriately small range and a large # of bins (tempHisto is set rebinnable)
  string weightedCut = "weight*("+cut+")";
  int nPass[nVarsMax];          // # of tree events passing cut, to monitor statistics
  double nWeightedPass[nVarsMax];   // weighted # events passing
  double varMin = 1e12;
  double varMax = -1e12;
  for (int varIt = 0; varIt < nVars; varIt++) {
    TH1F * tempHisto = new TH1F("tempHisto","",10000,0,0.1);
    tempHisto->SetBit(TH1F::kCanRebin);
    allTrees[tree]->Draw((varKeys[varIt]+" >> "+tempHisto->GetName()).c_str(),weightedCut.c_str());
    if (ignoreZeroBin)  tempHisto->SetBinContent(tempHisto->FindBin(0.),0);
    nPass[varIt] = tempHisto->GetEntries();
    nWeightedPass[varIt] = tempHisto->Integral();
    for (int binIt = 1; binIt <= tempHisto->GetNbinsX(); binIt++) {
      double binContent  = tempHisto->GetBinContent(binIt);
      double binLowEdge  = tempHisto->GetBinLowEdge(binIt);
      double binHighEdge = tempHisto->GetBinWidth(binIt) + binLowEdge;
      if (binContent == 0)  continue;
      if (binLowEdge  < varMin)  varMin = binLowEdge;
      if (binHighEdge > varMax)  varMax = binHighEdge;
    }
    delete tempHisto;
  }

  // add some padding at the edges
  double varRange = varMax-varMin;
  varMin -= varRange*0.00;
  varMax += varRange*0.001;

  // or reset to user binning
  if (xmin > -1e10) varMin = xmin;
  if (xmax > -1e10) varMax = xmax;
  
  // prepare comparison histos
  for (int varIt = 0; varIt < nVars; varIt++) {
    if (varComparisonHistos[varIt] == NULL) {  // define the comparison histos if this is the first time
      char histoName[100];
      sprintf(histoName,"varComparisonHisto%i",varIt);
      varComparisonHistos[varIt] = new TH1F(histoName,"",nBins,0,1);  // range will be recomputed later
      varComparisonHistos[varIt]->SetLineColor(treeColors[varIt]);
      varComparisonHistos[varIt]->SetMarkerColor(treeColors[varIt]);
      varComparisonHistos[varIt]->SetLineStyle(treeStyles[varIt]);
      varComparisonHistos[varIt]->SetStats(false);
      varComparisonHistos[varIt]->Sumw2();
    }
    varComparisonHistos[varIt]->Reset();
    varComparisonHistos[varIt]->SetBins(nBins,varMin,varMax);
    varComparisonHistos[varIt]->SetLineWidth(1);
    if (bold) varComparisonHistos[varIt]->SetLineWidth(2);
  }

  // fill the comparison histos
  double maxBinContent = 0;
  double minBinContent = 1e12;
  int nDisplayed[nVarsMax];
  double nWeightedDisplayed[nVarsMax];
  for (int varIt = 0; varIt < nVars; varIt++) {
    allTrees[tree]->Draw((varKeys[varIt]+" >> "+varComparisonHistos[varIt]->GetName()).c_str(),weightedCut.c_str());
    nDisplayed[varIt] = varComparisonHistos[varIt]->GetEntries();
    nWeightedDisplayed[varIt] = varComparisonHistos[varIt]->Integral();
    if (!quiet) {
      cout << nWeightedDisplayed[varIt] << " / " << nWeightedPass[varIt] << " weighted    " 
	   << nDisplayed[varIt]         << " / " << nPass[varIt]         << " actual      ";
      if (nWeightedPass[varIt] > 0)  cout << "(" <<  nWeightedDisplayed[varIt]/nWeightedPass[varIt] << ")" << endl;
      else                           cout << endl;
    }
    if (normalize && nWeightedPass[varIt] > 0)  varComparisonHistos[varIt]->Scale(1./(nWeightedPass[varIt]));  // if requested, rescale all histos to unity
    for (int binIt = 1; binIt <= varComparisonHistos[varIt]->GetNbinsX(); binIt++) {
      if (varComparisonHistos[varIt]->GetBinContent(binIt) < minBinContent &&
	  varComparisonHistos[varIt]->GetBinContent(binIt) > 0)
	minBinContent = varComparisonHistos[varIt]->GetBinContent(binIt);
      if (varComparisonHistos[varIt]->GetBinContent(binIt) > maxBinContent)
	maxBinContent = varComparisonHistos[varIt]->GetBinContent(binIt);
    }
  }

  // adjust vertical range, etc, and draw
  varComparisonLegend->Clear();
  for (int varIt = 0; varIt < nVars; varIt++) {  
    if (logable)   varComparisonHistos[varIt]->SetMinimum(minBinContent*0.5);
    else           varComparisonHistos[varIt]->SetMinimum(0);
    if (logable)   varComparisonHistos[varIt]->SetMaximum(maxBinContent*2);
    else           varComparisonHistos[varIt]->SetMaximum(maxBinContent*1.05);
    varComparisonHistos[varIt]->GetXaxis()->SetTitle("");
    if (normalize)   varComparisonHistos[varIt]->GetYaxis()->SetTitle("normalized rate");
    else             varComparisonHistos[varIt]->GetYaxis()->SetTitle(lumiTitle);
    varComparisonHistos[varIt]->GetYaxis()->SetTitleOffset(1.2);
    varComparisonHistos[varIt]->GetYaxis()->SetLabelSize(0.035);
    string histoTitle = varKeys[varIt];
    if (!(atoi(cut.c_str()) > 0))  histoTitle += " ("+cut+")";
    varComparisonHistos[varIt]->SetTitle(histoTitle.c_str());
    if (noTitle)  varComparisonHistos[varIt]->SetTitle("");
    string drawOptionFull = drawOption;
    drawOptionFull = "HIST " + drawOptionFull;
    if (varIt != 0)  drawOptionFull += " same";
    varComparisonHistos[varIt]->SetFillStyle(0);
    varComparisonHistos[varIt]->SetFillColor(10);
    varComparisonHistos[varIt]->Draw(drawOptionFull.c_str());
    varComparisonLegend->AddEntry(varComparisonHistos[varIt],varKeys[varIt].c_str());
  }

  if (legend) {
    varComparisonLegend->Draw();
    varComparisonLegend->SetBorderSize(0);
    varComparisonLegend->SetFillColor(10);
  }

  return;

}








