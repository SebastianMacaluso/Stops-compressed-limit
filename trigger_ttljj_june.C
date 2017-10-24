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

using namespace std;
typedef complex<double> dcomplex;

bool smear = false;
bool normalize = false;


const int jetID = 7;              // ID of a light jet

const double ptjcut = 20;

// event file info

const int SimFirst=1;
const int SimLast=48;
const int nSims = SimLast-SimFirst+1; 

char const * samplename="jets.tt~ljj";
char const * outputDir = ".";
char const * inputDir = "/het/p1/macaluso/multi_jobs/tt~ljj__pt1_450_noID_eff_results";  

int eventrun[nSims]; int eventpass[nSims]; double xsecfile[nSims]; double sigmafile[nSims];

// global event variables

struct lheobject {
  TLorentzVector v;
  int id;
  int stat;
};

// random # generator

TRandom * randGen = new TRandom(12345);

// function for comparing Pt's of two 4-vectors (for sorting lists)
bool comparePt(lheobject p1, lheobject p2) { return (p1.v.Pt() > p2.v.Pt()); }

TF1 *flatfunc = new TF1("flatfunc","1.",0,1); //Create flat distribution function


/*
//b tagging function.. takes in as arguments a collection of bjets as a vector of lheobject structs.. using a flat distribution function, a random number is chosen between 0 and 1.. the bjet is kept only if this random number falls below 0.75 MQP

vector<lheobject> btageff(vector<lheobject> bOG) {
  vector<lheobject> effbtag; //Create a vector of lheobject structs that will be used to contain the surviving bjets MQP

  for (int i=0; i<(int)bOG.size(); i++) {
    double beff = flatfunc->GetRandom(); //Get random number between 0 and 1
    
    if (beff<0.75) effbtag.push_back(bOG[i]); //If beff is below 0.75, keep the bjet
  }

  return effbtag;
}
*/


/*jet smearing function.. takes in as arguments a collection of jets as a vector of lheobject structs as well as the average number of pile up events per crossing mu.. only mu values of 0, 50, or 140 will work.. otherwise an error message will come up and the original unsmeared jets will be returned.. MQP*/

vector<lheobject> jetsmear(vector<lheobject> jetOG, int mu) {
  vector<lheobject> smearjets; //Create a vector of lheobject structs that will be used to contain the smeared jets MQP
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
      lheobject smearp;  smearp.v = smearj; smearp.id = jetOG[i].id; smearp.stat = jetOG[i].stat; //Creating the lheobject struct for the smeared jet MQP
      smearjets.push_back(smearp); //Pushing back the smeared jets MQP
    }

    return smearjets;
  }
}

//////////////////////////////////////////////////////////////////////////////
// MAIN function

void trigger_ttljj_june() {

  //Checking true # files for reweight
  
  int nSimsTrue=0;
  
  for (int i = 0; i < nSims; i++) {
    char inputcheck[200];
    int filenum=i+SimFirst;
    sprintf(inputcheck,"%s/%s-%d.lhe",inputDir,samplename,filenum);
    
    ifstream ifs;
    ifs.open(inputcheck);
    if (ifs.good()) nSimsTrue++;
    if (!ifs.good()) { cout << "Missing file here" << endl; }




    if (ifs.good()) {
      
      LHEF::Reader reader(ifs);
      
      sigmafile[i]=reader.heprup.XSECUP[0];
      
      
      if (sigmafile[i]==0)  {
	cout << "nSimstrue before discarding empty file " << nSimsTrue << endl;
	nSimsTrue--;
	cout << "nSimstrue after discarding empty file " << nSimsTrue << endl;
      }
      
    }


  
  }

  //LHE writer

  char outputfile[200];
  sprintf(outputfile,"%s/%s_%d_%d_btag_cuts_noID_eff.lhe",outputDir,samplename,SimFirst,SimLast);
  LHEF::Writer lhe_writer(outputfile);

  for (int simIt = 0; simIt < nSims; simIt++) {
    
    //LHE reader

    char inputfile[200];
    int filenum=simIt+SimFirst;
    sprintf(inputfile,"%s/%s-%d.lhe",inputDir,samplename,filenum);
    
    ifstream ifs;
    ifs.open(inputfile);
    cout << endl << "///////////////////////////////////////////////////////////" << endl;
    cout << "//   " << inputfile << endl;
    if (!ifs.good()) { cout << "BAD LHE FILE" << endl;  continue; }  


    LHEF::HEPRUP lhe_header;
    //lhe_header.resize(1);
    //lhe_header.IDBMUP.first = totalGen;
    //lhe_header.XSECUP[0] = sigma;
    lhe_writer.heprup = lhe_header;
    lhe_writer.init();  // write empty header



    LHEF::Reader reader(ifs);
    
    xsecfile[simIt]=reader.heprup.XSECUP[0];

    //Now loop over all events:

    long ieverun = 0;
    long ievepass = 0;
    bool done = false;
    done = !reader.readEvent();
    while ( !done ) {
      ++ieverun;
      
      if (ieverun%1000==0)  {
 	cout << "Event #" << ieverun << endl;
      }

      int npart = reader.hepeup.NUP;
      double weightperfile = reader.hepeup.XWGTUP;
      //cout << "weight per file " << weightperfile << endl;
      double weight = weightperfile/(double)nSimsTrue;

      vector<lheobject> vp_lightjets;
      vector<lheobject> vp_bjets;
      vector<lheobject> vp_leps;
      vector<lheobject> vp_MET;
      TLorentzVector  initial4vector;
      
      // Loop over particles in the event

      for (int i = 0; i < npart; i++) {
	vector<double> fiveVector = reader.hepeup.PUP[i];
	
	int id = reader.hepeup.IDUP[i];
	int stat = reader.hepeup.ISTUP[i];
	
	TLorentzVector vect(fiveVector[0],fiveVector[1],fiveVector[2],fiveVector[3]);
	
	lheobject part;  part.v = vect;  part.id = id;  part.stat = stat;
	
	if (stat == -1) { 
	  initial4vector += vect;  
	  continue; 
	}
	
	if ( stat==1 && (abs(id)==7 || abs(id)==4)) {
	  vp_lightjets.push_back(part); //Grabbing all jets MQP
	}

	if ( stat==1 && abs(id)==5) {
	  vp_bjets.push_back(part); //Grabbing all jets MQP
	}

	if ( stat==1 && (abs(id)==11 || abs(id)==13 || abs(id)==15)) {
	  vp_leps.push_back(part); //Grabbing all leptons MQP
	}
	if ( stat==1 && (abs(id)==12 || abs(id)==14 || abs(id)==16 || abs(id)==1000022 || abs(id)==1000023 || abs(id)==1000025 || abs(id)==1000035)) {
	  vp_MET.push_back(part); //Grabbing all neutrinos and neutralinos SM
	}
      }


      ////////////////////////////    SM ///////////////////////////////////////////
      // b tagging efficiency

      // cout << (int)vp_bjets.size() << "bjets before" << endl; 
      // cout << (int)vp_lightjets.size() << "non bjets before" << endl;
      
      vector<lheobject> vp_btagjets;
      
      for (int i=0; i<(int)vp_bjets.size(); i++) {
	
	double beff = flatfunc->GetRandom(); //Get random number between 0 and 1
	//cout << beff << "value of the flat function" << endl;
	if (beff<0.7) vp_btagjets.push_back(vp_bjets[i]); //If beff is below 0.7, keep the bjet
	else {
	  // cout << (int)vp_bjets[i].id << "b jet ID before" << endl;
   
	  vp_bjets[i].id=jetID; 
	  // cout << (int)vp_bjets[i].id << "b jet ID AFTER" << endl;
	  vp_lightjets.push_back(vp_bjets[i]);//If beff>0.7 we move it to vp_lightjets
	}
	//	cout << (int)vp_bjets[i].id << "b jet ID FINAL" << endl;
      }
      
      //cout << (int)vp_lightjets.size() << "non bjets AFTER" << endl;
      //cout << (int)vp_btagjets.size() << "bjets  AFTER" << endl;

      //cout << "---------------------------------------------------------------------" << endl;
     
      //----------------------------------------------------------------------------------------------------

      
      //vector<lheobject> vp_btagjets=btageff(vp_bjets); //Applying btag efficiency
      //vector<lheobject> vp_btagjets=vp_bjets; // NOT Applying btag efficiency
      vector<lheobject> vp_smearbtagjets_no_ptjcut=jetsmear(vp_btagjets,50); //Smearing tagged bjets
      vector<lheobject> vp_smearlightjets_no_ptjcut=jetsmear(vp_lightjets,50); //Smearing the jets with mu=50 MQP

      vector<lheobject> vp_smearjets_no_ptjcut=vp_smearlightjets_no_ptjcut;
      vp_smearjets_no_ptjcut.insert(vp_smearjets_no_ptjcut.end(),vp_smearbtagjets_no_ptjcut.begin(),vp_smearbtagjets_no_ptjcut.end()); //Combining all jets


      ///////////////////////////////////////           SM                /////////////////////////////////////////////////
      //I apply a ptjcut on jets after smearing


      vector<lheobject> vp_smearbtagjets;
      
      for (int m=0; m<(int)vp_smearbtagjets_no_ptjcut.size(); m++)
	{
	  if (vp_smearbtagjets_no_ptjcut[m].v.Pt()> ptjcut) {vp_smearbtagjets.push_back(vp_smearbtagjets_no_ptjcut[m]);}
	}    

      
      vector<lheobject> vp_smearlightjets;
      
      for (int m=0; m<(int)vp_smearlightjets_no_ptjcut.size(); m++)
	{
	  if (vp_smearlightjets_no_ptjcut[m].v.Pt()> ptjcut) {vp_smearlightjets.push_back(vp_smearlightjets_no_ptjcut[m]);}
	}   

      
      vector<lheobject> vp_smearjets;
      
      for (int m=0; m<(int)vp_smearjets_no_ptjcut.size(); m++)
	{
	  if (vp_smearjets_no_ptjcut[m].v.Pt()> ptjcut) {vp_smearjets.push_back(vp_smearjets_no_ptjcut[m]);}
	}  
      
      
      
      /////////////////////////////////////////      SM          //////////////////////////////////////////////////////////////////////////




      sort(vp_smearbtagjets.begin(), vp_smearbtagjets.end(), comparePt);
      sort(vp_smearlightjets.begin(), vp_smearlightjets.end(), comparePt);
      sort(vp_smearjets.begin(), vp_smearjets.end(), comparePt);

      TLorentzVector l_MH;
      double HT=0;

      for (int i=0; i<(int)vp_smearjets.size(); i++) {
	l_MH=l_MH-vp_smearjets[i].v;
	HT=HT+vp_smearjets[i].v.Pt();
      }

      TLorentzVector l_MHT; l_MHT.SetPx(l_MH.Px()); l_MHT.SetPy(l_MH.Py()); l_MHT.SetE(sqrt(l_MH.Px()*l_MH.Px() + l_MH.Py()*l_MH.Py())); 
      lheobject p_MHT; p_MHT.v = l_MHT;  p_MHT.id = 99;  p_MHT.stat = 1;
 
      //ISR identification

      vector<lheobject> vp_ISRcand;
      double mtop=173.;

      if ((int)vp_smearbtagjets.size()>1) {
	for (int i=0; i<(int)vp_smearlightjets.size(); i++) {
	  double mbj1=(vp_smearlightjets[i].v+vp_smearbtagjets[0].v).Mag();
	  double mbj2=(vp_smearlightjets[i].v+vp_smearbtagjets[1].v).Mag();

	  if (mbj1>(1.2*mtop) && mbj2>(1.2*mtop)) {
	    vp_ISRcand.push_back(vp_smearlightjets[i]);
	  }
	}
      }

      sort(vp_ISRcand.begin(), vp_ISRcand.end(), comparePt);

      //dPhiMin Calculation
      
      double dPhiMin_calc= 999;

      for (int j=0; j<3; j++) {
	if ((int)vp_smearjets.size()>j) {
	  if (fabs(vp_smearjets[j].v.DeltaPhi(p_MHT.v)) < dPhiMin_calc) {
	    dPhiMin_calc=fabs(vp_smearjets[j].v.DeltaPhi(p_MHT.v));
	  }
	}
      }

      double dPhiMin=0;

      if (dPhiMin_calc!=999) dPhiMin=dPhiMin_calc;
      
      // placing cuts

      if ((int)vp_ISRcand.size()>0 && /* (int)vp_smearjets.size()>6 &&*/ (int)vp_leps.size()==0 && dPhiMin>0.5 && vp_ISRcand[0].v.Pt()>400 && p_MHT.v.Pt()/sqrt(HT)>3) {
      
	++ievepass;
	// write LHE output
	vector<lheobject> allParticles;
	allParticles.insert(allParticles.end(),vp_MET.begin(),vp_MET.end());
	allParticles.insert(allParticles.end(),vp_smearjets.begin(),vp_smearjets.end());
	allParticles.insert(allParticles.end(),vp_leps.begin(),vp_leps.end());
	//allParticles.push_back(p_MHT); //Push back MHT vector into particle list for lhe insertion MQP
	
	LHEF::HEPEUP lhe_event;
	lhe_event.resize(allParticles.size());
	lhe_event.XWGTUP = weight;
	
	for (int it = 0; it < (int)allParticles.size(); it++) {
	  lheobject localPart = allParticles[it];
	  
	  lhe_event.ISTUP[it] = localPart.stat;
	  lhe_event.IDUP[it] = localPart.id;
	  lhe_event.PUP[it][0] = localPart.v.Px();
	  lhe_event.PUP[it][1] = localPart.v.Py();
	  lhe_event.PUP[it][2] = localPart.v.Pz();
	  lhe_event.PUP[it][3] = localPart.v.E();
	  lhe_event.PUP[it][4] = localPart.v.Mag();
	}
	
	lhe_writer.hepeup = lhe_event;
	lhe_writer.writeEvent();
		 }

      done = !reader.readEvent();
       }

    eventrun[simIt]=ieverun;
    eventpass[simIt]=ievepass;
  }

  int totalrun=0;
  int totalpass=0;
  double xsecsum=0;

  for (int i=0; i<nSims; i++) {
    cout << "Total Events in File# " << i+SimFirst << " : " << eventrun[i] << endl;
    cout << "Events Passing in File# " << i+SimFirst << " : " << eventpass[i] << endl;
    cout << "///////////////////////////////" << endl;

    totalrun=totalrun+eventrun[i];
    totalpass=totalpass+eventpass[i];
    xsecsum=xsecsum+xsecfile[i];
  }

  cout << "///////////////////////////////" << endl;
  cout << "///////////////////////////////" << endl;
  cout << "Total Events in This Batch: " << totalrun << endl;
  cout << "Events Passing in This Batch: " << totalpass << endl;
  //cout << "Average cross section before reweighting in This Batch: " << xsecsum << endl;
  cout << "Average cross section after reweighting in This Batch: " << xsecsum/(nSimsTrue) << endl;
  cout << "Total number of files in This Batch: " << nSims << endl;
  cout << "Total files that are not empty " << nSimsTrue << endl;
  

  /*
  double xsecavg=xsecsum/(double)nSims;
  double xsecfinal=xsecavg*((double)totalpass/(double)totalrun);

  LHEF::HEPRUP lhe_header;
  lhe_header.resize(1);
  lhe_header.IDBMUP.first = totalpass;
  lhe_header.XSECUP[0] = xsecfinal;
  lhe_writer.heprup = lhe_header;
  lhe_writer.init();  // write empty header
  */
}
