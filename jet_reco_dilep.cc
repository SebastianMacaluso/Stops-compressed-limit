// jet_reco.cc    (BAT   August 2009)
//
//
// Jet reconstruction (using fastjet) for studying dileptonic ttbar events
//
// NOTE:  This version has been modified to read in the MC@NLO data from Can, originally generated for our stop project.
//        In particular, it understands when to look for negative-weight events.
//        It also either snips off parton-level info, or pretends that parton-level objects are final-state and ignores the real final-state entirely.
//



//----------------------------------------------------------------------
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "LHEF.h" // Les Houches event format headers
#include<iostream> // needed for io
#include<fstream>  // file io
#include<sstream>  // needed for internal io
#include<vector> 
#include<string>


///////     SM //////////////////////////////////////////////////////////

#include<complex>
#include<algorithm>
#include <cmath>
#include <math.h>



///////////////////////////////////////////////////////////////////////






using namespace std;

bool onlyPartons = false;  // only output parton-level and pretend that they're final-state...still does jet reco (kind of a hack)


// variables controlling cuts (should not be weaker than any generator-level cuts)
const double ptlcut = 5;         // pt cut on individual hard leptons
const double ptjcut = 15;         // pt cut on jets
const double etalcut = 2.5;       // max eta for leptons
const double etajcut = 5.0;       // max eta for jets
const double etabcut = 2.5;       // max eta for b-tag
const double isoFrac = 0.90;      // fractional cone pt threshold for lepton iso:  ptl/(ptl+pt"hadrons"-in-cone))
const double isoR = 0.4;          // iso cone size
const double Rjet = 0.4;          // jet clustering radius
const double ReJet = 0.4;        // electro-jet delta R

//-------------------------------            SM        ------------------------------

const double eR = 0.2;           //  Delta R to see if the object is considered to be an electron  SM    
const double etacutmu = 2.4;       // max eta for muons SM
const double etacute = 2.47;        // max eta for electrons SM
const double ptcutlep = 10;         // pt cut on leptons SM

//------------------------------------------------------------------------------------

const int jetID = 7;              // ID of a light jet

const double pi = 3.14159265359;


// determine the charge of a particle based on id
int charge(int);

// calculate absolute phi distance
double Dphi(const fastjet::PseudoJet &, const fastjet::PseudoJet &);

// calculate eta-phi distance
double DR(const fastjet::PseudoJet &, const fastjet::PseudoJet &);

// b-tagging and c-tagging
void bcTag(vector<fastjet::PseudoJet> &, const vector<fastjet::PseudoJet> &, double);

// master b-tag function, used with a global switch to select which method to use
//bool isAbjet(const fastjet::PseudoJet &, const fastjet::ClusterSequence &, double, double);

// b-tag #1: tagged if a b-hadron was "nearby"
//bool nearbyB(const fastjet::PseudoJet &, double);

// b-tag #2: uniquely assign b-hadrons to nearest jets (with distance cutoff), avoid multiple jets tagged for same b
//void assignBs(const vector<fastjet::PseudoJet> &, vector<double> &, vector<double> &, vector<int> &, double);

// lists of the b-hadrons in the event
vector<fastjet::PseudoJet> b_hadrons;

// b tagging stuff
const double searchDist = 1.0;     // distance scale used in b-tagging methods that search for nearby ancestral b-hadrons, measured wrt R
const double matchDist = 1.0;      // distance scale used in b-tagging methods that match ancestral b-hadrons to closest jet, measured wrt R
const double btagMethod = 2;       // 1: look for nearby ancestral b-hadrons 
                                   // 2: match ancestral b-hadrons to closest jet
//const bool bQuarkTag = false;      // use b-quarks instead of b-hadrons for method #2 (mainly for parton-level shower analysis)


fastjet::PseudoJet blankJet(0,0,0,0);





//---------------------------        SM        --------------------------------------------
double rand01()//random number between 0 and 1
      { return (double)rand()/RAND_MAX;}




//---------------------------        SM        --------------------------------------------
double rand02()//random number between 0 and 1
      { return (double)rand()/RAND_MAX;}


//----------------------------       SM        ----------------------------------------------


// function for comparing Delta R's  (for sorting lists) from lower to geater delta R
bool compareDeltaR(double R1, double R2) { return (R1 < R2); }


//----------------------------       SM        ----------------------------------------------



// function for comparing Pt's of two 4-vectors (for sorting lists)
bool comparePt(fastjet::PseudoJet p1, fastjet::PseudoJet p2) { return (p1.pt() > p2.pt()); }

///-----------------------------------------------





/// MAIN
int main (int argc, char ** argv) {
  


  // input directory stuff
  char * inputDir = "/het/p1/macaluso/topjet_fastjet/";
  char * fileTag = "ht800";
  double reweight = 1.0;  // SET EQUAL TO BR FOR SELECTIVELY DECAYED SAMPLES

  // override with command-line input
  if (argc == 4) {
    inputDir = argv[1];
    fileTag = argv[2];
    reweight = atof(argv[3]);
  }
  if (argc == 3) {
    inputDir = argv[1];
    fileTag = argv[2];
  }

  
  cout << "input from : " << inputDir << endl;
  cout << "file tag   : " << fileTag << endl;
  cout << "reweighting: " << reweight << endl;
  char filename[200];

  srand(11232698);


  // output files of reconstructed events
  //ofstream recoData;
  //sprintf(filename,"%s/jets.%s.dat",inputDir,fileTag);
  //recoData.open(filename);

  // LHE file of reconstructed events
  sprintf(filename,"%s/jets.%s.lhe",inputDir,fileTag);
  LHEF::Writer lhe_writer(filename);

  // variables characterizing each data file
  double sigma;
  int nGen, nPass;

  // event variables
  int eventNum,nGod,nLep,nB,nPart;

  // particle variables
  int id;
  float px, py, pz, e;
  int bAncestor;

  // variables characterizing the full data set
  int totalGen = 0;
  int nLoaded = 0;
  int nFinal = 0;
  double sigmaSum = 0.;

  // counters for event failure
  int nTooFewLeptons = 0;
  int nTooFewIsoLeptons = 0;
  int nTooManyIsoLeptons = 0;
  int nFailOppositeCharge = 0;
  int nTooFewJets = 0;
  int nTooManyClosedJets_e = 0; //SM
  int nTooManyClosedJets_mu = 0; //SM




  // first check if there are multiple files and if we need to apply reweightings
  int nFiles = 0;
  bool multipleFiles = true;  // default: look for multiple files labelled 01, 02, etc
  string datfilenames[100];
  while (1) {

    nFiles++;
    
    ifstream datfileTemp;
    if (multipleFiles) {
      // if (nFiles < 10)
      //sprintf(filename,"%s/output.%s.0%i.dat",inputDir,fileTag,nFiles);
      //else
	sprintf(filename,"%s/output.%s.%i.dat",inputDir,fileTag,nFiles);
    }
    datfileTemp.open(filename);
    if(!datfileTemp.good() && nFiles == 1)  {
      multipleFiles = false;
      sprintf(filename,"%s/output.%s.dat",inputDir,fileTag);
    }
    datfilenames[nFiles-1] = string(filename);
    
    ifstream sigfile;
    if (!multipleFiles)
      sprintf(filename,"%s/sigma.%s.dat",inputDir,fileTag);
    // else if (nFiles < 10)
    //sprintf(filename,"%s/sigma.%s.0%i.dat",inputDir,fileTag,nFiles);
    else
      sprintf(filename,"%s/sigma.%s.%i.dat",inputDir,fileTag,nFiles);
    sigfile.open(filename);
    if (!sigfile.good())  break;
    sigfile >> sigma >> nGen >> nPass;
    cout << "Checking sigma file #" << nFiles << endl;
    cout << nGen << "   generated, " << nPass << " passed" << endl << endl;
    totalGen += nGen;
    sigmaSum += sigma;
    sigfile.close();

    if (!multipleFiles) break;
  }
  if (multipleFiles) nFiles--; // overcounted on last file check
  sigma = sigmaSum/1000./nFiles*reweight;   // assuming files contain unweighted events (1000 for fb->pb from my old output conventions)
  double weight = sigma/totalGen;
  bool usePerEventWeights = (sigmaSum < 0);  // if sigma is negative, use as a flag to look for per-event weights instead


  LHEF::HEPRUP lhe_header;
  lhe_header.resize(1);
  lhe_header.IDBMUP.first = totalGen;
  lhe_header.XSECUP[0] = sigma;
  lhe_writer.heprup = lhe_header;
  lhe_writer.init();  // write empty header

  // read input files
  for (int fileIt = 0; fileIt < nFiles; fileIt++) {

    ifstream datfile;
    datfile.open(&datfilenames[fileIt][0]);
    cout << datfilenames[fileIt] << endl;
    if(!datfile.good())  {
      cout << "FAILING TO FIND:  " << datfilenames[fileIt] << endl;
      break;
    }
    cout << endl << "------------------------------------------------------" << endl;
    cout << "Open file #" << fileIt+1 << endl;

    // iterate over events in the data file
    while (1) {
      
      b_hadrons.clear();

      vector<fastjet::PseudoJet> godParticles;
      vector<int>                godStatuses;
      vector<fastjet::PseudoJet> hadrons;
      vector<fastjet::PseudoJet> allLeptons;
      vector<fastjet::PseudoJet> failedLeptons;

      vector<fastjet::PseudoJet> electrons; //SM
      vector<fastjet::PseudoJet> recelectrons; //SM
      vector<fastjet::PseudoJet> recelectrons2; //SM
      vector<fastjet::PseudoJet> muons; //SM
      vector<fastjet::PseudoJet> taus; //SM

      vector<fastjet::PseudoJet> initmuons; //SM

      vector<fastjet::PseudoJet> ehadrons; //SM
      fastjet::PseudoJet E_ehadrons; //SM
      vector<fastjet::PseudoJet> failed_electrons; //SM
     

      vector<fastjet::PseudoJet> sparticles;
      vector<fastjet::PseudoJet> neutrinos;
      vector<fastjet::PseudoJet> isoLeptons;
      vector<fastjet::PseudoJet> finalLeptons;
      vector<fastjet::PseudoJet> jets;

     
      vector<fastjet::PseudoJet> recelectrons3;  //SM
      vector<fastjet::PseudoJet> recmuons; //SM
      vector<fastjet::PseudoJet> jets2; //SM
      vector<fastjet::PseudoJet> jets3; //SM
      






      
      // load God-level info
      if (usePerEventWeights)  { datfile >> eventNum >> nGod >> weight;  weight *= reweight/nFiles; }  //  *** ASSUMES SAME # EVENTS PER FILE
      else  datfile >> eventNum >> nGod;
      if (!datfile.good()) break;
      string line;
      getline(datfile,line); // increment line
      for (int i = 0; i < nGod; i++) {
	int status = 2;
	getline(datfile,line);
	//datfile >> id >> px >> py >> pz >> e;
	int hasStatus = sscanf(line.c_str(),"%*d %*f %*f %*f %*f %*d")+1;
	if (hasStatus)  sscanf(line.c_str(),"%d %f %f %f %f %d",&id,&px,&py,&pz,&e,&status);
	else            sscanf(line.c_str(),"%d %f %f %f %f",&id,&px,&py,&pz,&e);
	godParticles.push_back(fastjet::PseudoJet(px,py,pz,e));
	godParticles.back().set_user_index(id);
	if (status == 2 && abs(id) != 6 && abs(id) >= 1 && abs(id) <= 16)  status = 1;  // correct mistaken sign of "secondary hard particle" codes from MC@NLO translation <----- and now turn them into "final-state" particles!!
	godStatuses.push_back(status);
      }
      if (!onlyPartons) {
	godParticles.clear();
	godStatuses.clear();
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // load b hadrons
      datfile >> id >> nB;
      //cout << endl << "b-hadrons" << endl;//BAT
      for (int i = 0; i < nB; i++) {
	datfile >> id >> px >> py >> pz >> e;
	b_hadrons.push_back(fastjet::PseudoJet(px,py,pz,e));
	b_hadrons.back().set_user_index(id);
	//cout << "  (" << b_hadrons.back().rap() << "," << b_hadrons.back().phi() << ")  pt=" << b_hadrons.back().perp() << endl; //BAT
      }
      // load everything else
      datfile >> id >> nPart;
      for (int i = 0; i < nPart; i++) {
        datfile >> id >> px >> py >> pz >> e;
	fastjet::PseudoJet part(px,py,pz,e);
	part.set_user_index(id);
	if (abs(id) > 1000000)
	  sparticles.push_back(part);
	else if (abs(id) == 12 || abs(id) == 14 || abs(id) == 16)
	  neutrinos.push_back(part);
	//	else if ((abs(id) == 11 || abs(id) == 13) && part.perp() > ptlcut && abs(part.eta()) < etalcut)
	  // note:  leptons that can't be identified because of too little pt or too-high eta get tossed in as "hadrons" below
	//	  allLeptons.push_back(part);
	else if (abs(id)==11 && part.perp() > ptcutlep && abs(part.eta()) < etacute ) electrons.push_back(part);  //SM
	else if (abs(id)==13 && part.perp() > ptcutlep && abs(part.eta()) < etacutmu ) initmuons.push_back(part);  //SM


	else
	  hadrons.push_back(part);
      }
      
      //-------------------------------------------------------------------------------------------------     

      //We follow most of the physics objet reconstruction from 1406.1122
      

      ////////////          SM        /////////////////////////////////////////////
      // We cut on the electrons depending on the efficiency.
      
      if(electrons.size()>0)
	{
	  for (int eIt=0; eIt<(int)electrons.size(); eIt++)
	    {
	      double eeff= rand01();
	      if (eeff<0.1) hadrons.insert(hadrons.end(),electrons.begin(),electrons.end()); //If eeff is below 0.1, dump it to hadrons
	      else recelectrons.push_back(electrons[eIt]); //If eeff > 0.1 keep the electron
	    }
	}
      
      
      ////////////          SM        ///////////////////////////////////////////
      //We sum the 4 momentum vectors of all the hadrons for a cone within DeltaR<0.2 around the electron. If the energy of the final vector (E_ehadrons) is more than 10% of the energy of the electron, then we throw out the electron.
      
      if(recelectrons.size()>0 && hadrons.size()>0)
	{
	  for (int eit2=0; eit2<(int)recelectrons.size(); eit2++)
	    {
	      
	      for (int hit=0; hit<(int)hadrons.size(); hit++)
		{
		  if (DR(recelectrons[eit2],hadrons[hit]) < eR) ehadrons.push_back(hadrons[hit]);
		}
	     	     
	      for (int Eit=0; Eit<(int)ehadrons.size(); Eit++)   { E_ehadrons=E_ehadrons+ehadrons[Eit];}
	      if ((E_ehadrons.E()/recelectrons[eit2].E())<0.1) recelectrons2.push_back(recelectrons[eit2]);
	      else failed_electrons.push_back(recelectrons[eit2]);
	      
	    }
	}

   
      /////////////       SM         ///////////////////////////////////////////////////////////
      //We dump the failed electrons to hadrons.
     
      if (failed_electrons.size()>0)
	{      
	  for (int eit3=0; eit3<(int)failed_electrons.size(); eit3++)
	    {
	      hadrons.insert(hadrons.end(),failed_electrons.begin(),failed_electrons.end());
	      
	    } 
	}


///////////////-------------------------------//////////////----------------------------------------////////////////////////////
        
        //-------------------------------------------------------------------------------------------------
        
        //We follow most of the physics objet reconstruction from 1406.1122
        
        
        ////////////          SM        /////////////////////////////////////////////
        // We cut on the muons depending on the efficiency.
        
        if(initmuons.size()>0)
        {
            for (int eIt=0; eIt<(int)initmuons.size(); eIt++)
            {
                double mmff= rand02();
                if (mmff<0.05) hadrons.insert(hadrons.end(),initmuons.begin(),initmuons.end()); //If mmff is below 0.05, dump it to hadrons
                else muons.push_back(initmuons[eIt]); //If eeff > 0.05 keep the muon
            }
        }
        



      ///////////////////////////////////////////////////////////////////////////////////////////////


      nLoaded++;
      
      if (nLoaded % 1000 == 0) {
	cout << eventNum << endl;
	//break;
      }
      

      if (!onlyPartons) {
	

	/*
	// search for isolated leptons
	for (int lepIt = 0; lepIt < allLeptons.size(); lepIt++) {
	  double lpt = allLeptons[lepIt].perp();
	  double hpt = 0;
	  for (int partIt = 0; partIt < hadrons.size(); partIt++) {
	    if (DR(allLeptons[lepIt],hadrons[partIt]) < isoR)  hpt += hadrons[partIt].perp();
	  }
	  if      (lpt/(lpt+hpt) > isoFrac)   isoLeptons.push_back(allLeptons[lepIt]);
	  else                                failedLeptons.push_back(allLeptons[lepIt]);
	}
	
      */
	
	// take the leading two leptons, check opposite-charged
	/*if (allLeptons.size() < 2) {
	  nTooFewLeptons++;
	  continue;
	}
	if (isoLeptons.size() < 2) {
	  nTooFewIsoLeptons++;
	  continue;
	}
	if (isoLeptons.size() > 2) {
	nTooManyIsoLeptons++;
	  continue;
	  }*/



      //	isoLeptons = sorted_by_pt(isoLeptons);
      //	finalLeptons=isoLeptons;





	/*finalLeptons.push_back(isoLeptons[0]);
	finalLeptons.push_back(isoLeptons[1]);
	if ( charge(finalLeptons[0].user_index()) != -charge(finalLeptons[1].user_index()) ) {
	  nFailOppositeCharge++;
	  continue;
	}
	for (int lepIt = 2; lepIt < isoLeptons.size();  lepIt++)  failedLeptons.push_back(isoLeptons[lepIt]);

	*/
	
	// append failed leptons to the "hadron" list
	//if (failedLeptons.size() > 0)  hadrons.insert(hadrons.end(),failedLeptons.begin(),failedLeptons.end());

	
	// construct jet algorithm
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rjet, recomb_scheme, strategy);
	
	// run the jet clustering
	fastjet::ClusterSequence clust_seq(hadrons,jet_def);
	
	// extract the inclusive jets
	vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptjcut);
	
	// apply cuts
	//cout << "jets" << endl; //BAT
	for (int jetIt = 0; jetIt < inclusive_jets.size(); jetIt++) 
	  {
	  fastjet::PseudoJet part = inclusive_jets[jetIt];

	  if (fabs(part.eta()) > etajcut || part.perp() < ptjcut)
	    { 
	    //cout << "  fail (" << part.rap() << "," << part.phi() << ")  pt=" << part.perp() << endl;   
	    continue; 
	    }
	  part.set_user_index(jetID);  
	  jets.push_back(part);
	  }
	jets = sorted_by_pt(jets);
	if (jets.size() < 7)
	  {
	  nTooFewJets++;
	  continue;
	  }
	///----------------------------------------------------------------------------------------------------------


	


	//////////////////////        SM         //////////////////////////////////////////////////////////////////////
	// We compare delta R between the electrons and the jets. First we find the DRmin=min Delta R(jet,lepton). If DRmin<0.4 then we add the lepton to that jet. If there are more than 1 lepton close to the same jet, then we add all of them to that jet. After this, we obtain the final jets and leptons.
	
	double DRmine = 999;
	double indexJe = -9;
	vector<double> DRminFinale;
	vector<int> JindexFinale;
	vector<int> eIndex;
		
	if(recelectrons2.size()>0)
	  {
	    for (int muit4=0; muit4<(int)recelectrons2.size(); muit4++)
	      {
		for (int jit3=0; jit3<(int)jets.size(); jit3++)
		  {
		   
		    if (DR(recelectrons2[muit4],jets[jit3]) < DRmine) 
		      {
			DRmine = DR(recelectrons2[muit4],jets[jit3]);
			indexJe = jit3;
		      }
		  }
	
		if (DRmine< ReJet) 
		  {
		    //cout << DRmin << " Delta R min" << endl;
		    
		    DRminFinale.push_back(DRmine);
		    
		    eIndex.push_back(muit4);
		    JindexFinale.push_back(indexJe);

		  }
		else recelectrons3.push_back(recelectrons2[muit4]);
	      }
	  }

	// cout << JindexFinale.size() << " # of electrons close to a jet" << endl;
	
	
	for (int m=0; m<(int)jets.size(); m++) 
	  {	
	    jets2.push_back(jets[m]);
	    // cout << jets2[m].px() << " jet2 px before adding electrons" << endl;
	  }


	//	cout << jets2.size() << " # jets2" << endl;
	//cout << JindexFinale.size() << " # electrons close to a jet" << endl;
	


	for (int m=0; m<(int)JindexFinale.size(); m++) 
	  {
	    jets2[JindexFinale[m]]=jets2[JindexFinale[m]]+recelectrons2[eIndex[m]];
	    // cout << recelectrons2[eIndex[m]].px() << "   px electrons close to a jet" << endl;
	  }


	/*
	for (int m=0; m<(int)jets2.size(); m++) 
	  {
	    cout << jets2[m].px() << " jet2 after adding electrons  px" << endl;
	  }
	
	
	cout <<  " ---------------------------------------------------" << endl;

	*/

	//////////////////////        SM         //////////////////////////////////////////////////////////////////////
	// We compare delta R between the electrons and the jets. First we find the DRmin=min Delta R(jet,lepton). If DRmin<0.4 then we add the lepton to that jet. If there are more than 1 lepton close to the same jet, then we add all of them to that jet. After this, we obtain the final jets and leptons.

	
	double DRmin = 999;
	double indexJ = -9;
	vector<double> DRminFinal;
	vector<int> JindexFinal;
	vector<int> muonIndex;
		
	if(muons.size()>0)
	  {
	    for (int muit4=0; muit4<(int)muons.size(); muit4++)
	      {
		for (int jit3=0; jit3<(int)jets2.size(); jit3++)
		  {
		   
		    if (DR(muons[muit4],jets2[jit3]) < DRmin) 
		      {
			DRmin = DR(muons[muit4],jets2[jit3]);
			indexJ = jit3;
		      }
		  }
	
		if (DRmin< ReJet) 
		  {
		    //cout << DRmin << " Delta R min" << endl;
	   
		    DRminFinal.push_back(DRmin);
		    
		    muonIndex.push_back(muit4);
		    JindexFinal.push_back(indexJ);

		  }
		else recmuons.push_back(muons[muit4]);
	      }

	  }
	    

	//cout << JindexFinal.size() << " # of muons close to a jet" << endl;
	// if (JindexFinal.size()>2)  cout << JindexFinal.size() << " # of muons close to a jet" << endl;

	
	for (int m=0; m<(int)jets2.size(); m++) 
	  {	
	    jets3.push_back(jets2[m]);
	    // cout << jets3[m].px() << " jet3 px before adding muons" << endl;
	  }


	//cout << jets3.size() << " # jets 3" << endl;
	//cout << JindexFinal.size() << " # muons close to a jet" << endl;



	for (int m=0; m<(int)JindexFinal.size(); m++) 
	  {
	    jets3[JindexFinal[m]]=jets3[JindexFinal[m]]+muons[muonIndex[m]];
	    //cout << muons[muonIndex[m]].px() << "   px muons close to a jet" << endl;
	  }


	/*
	for (int m=0; m<(int)jets3.size(); m++) 
	  {
	    cout << jets3[m].px() << " jet 3 after adding muons  px" << endl;
	  }

	
	cout <<  " ---------------------------------------------------" << endl;

	*/

	
	for (int mu=0; mu<(int)jets2.size(); mu++) 
	      {
		sort(jets2.begin(), jets2.end(), comparePt);
		sort(jets3.begin(), jets3.end(), comparePt);
		
		//cout << jets3[mu].pt() << " jets3  pt" << endl;
		
	      }
	
	//	cout <<  " ---------------------------------------------------" << endl;
	  
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////



	// apply b-tagging

	bcTag(jets3,b_hadrons,Rjet);

	/*
	vector<double> ptBs;
	vector<double> DRBs;
	vector<int> idBs;
	assignBs(jets,ptBs,DRBs,idBs,Rjet);
	for (int jetIt = 0; jetIt < jets.size(); jetIt++) {
	  if ( isAbjet(jets[jetIt],clust_seq,ptBs[jetIt],Rjet) ) {
	    jets[jetIt].set_user_index(idBs[jetIt]);
	  }
	  //cout << "  " <<  jets[jetIt].user_index() <<   "    (" << jets[jetIt].rap() << "," << jets[jetIt].phi() << ")  pt=" << jets[jetIt].perp() << endl; //BAT
	}
	*/

      }

      nFinal++;
      

      // write LHE output
      vector<fastjet::PseudoJet> allParticles;
      allParticles.insert(allParticles.end(),godParticles.begin(),godParticles.end());
      if (!onlyPartons) {
	allParticles.insert(allParticles.end(),sparticles.begin(),sparticles.end());
	allParticles.insert(allParticles.end(),neutrinos.begin(),neutrinos.end());
	//	allParticles.insert(allParticles.end(),finalLeptons.begin(),finalLeptons.end());
		allParticles.insert(allParticles.end(),recelectrons3.begin(),recelectrons3.end());    ///SM
		allParticles.insert(allParticles.end(),recmuons.begin(),recmuons.end());    ///SM
	//allParticles.insert(allParticles.end(),jets3.begin(),jets3.end());
	//allParticles.insert(allParticles.end(),recelectrons2.begin(),recelectrons2.end());
	//allParticles.insert(allParticles.end(),muons.begin(),muons.end());
	allParticles.insert(allParticles.end(),jets3.begin(),jets3.end());
      }
      LHEF::HEPEUP lhe_event;
      lhe_event.resize(allParticles.size());
      lhe_event.XWGTUP = weight;
      for (int it = 0; it < allParticles.size(); it++) {
	fastjet::PseudoJet localPart = allParticles[it];
	lhe_event.ISTUP[it] = 1;
	if (it < godParticles.size())  lhe_event.ISTUP[it] = godStatuses[it];  // particle from hard event
	lhe_event.IDUP[it] = localPart.user_index();
	for (int mu = 0; mu < 4; mu++)  lhe_event.PUP[it][mu] = localPart.four_mom()[mu];
	lhe_event.PUP[it][4] = localPart.m();
      }
      lhe_writer.hepeup = lhe_event;
      lhe_writer.writeEvent();
      

    } // end this event, go back and read next event

    datfile.close();

  } // end iteration over data files



  cout << endl << "-------------------------------------------------------" << endl;
  cout << totalGen << " generated" << endl;
  cout << nLoaded << " pass PYTHIA filters" << endl;
  cout << "-------------------------------------------" << endl;
  cout << nTooFewLeptons << " less than two lepton candidates" << endl;
  cout << nTooFewIsoLeptons << " less than two iso leptons" << endl;
  cout << nTooManyIsoLeptons << " more than two iso leptons" << endl;
  cout << nFailOppositeCharge << " fail opposite-charge dileptons" << endl;
  cout << nTooFewJets << " fail >=7 jets" << endl;
  cout << "-------------------------------------------" << endl;
  cout << nFinal << " pass" << endl;

}







////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions




//------------------------------------------------
// determine the charge of a particle based on id
//
//    *** the entire set of visible final state particles is
//         11 (e), 13 (mu), 22 (photon), 211 (pi^+/-), 321 (K^+/-), 130 (K_L)
//         310 (K_S)  (actually, I've never see this, but in principle, at large boost...)
//         2212 (p), 2112 (n) 
//
int charge(int id)
{
  if (id == -11 || id == -13 || id ==  211 || id ==  321 || id ==  2212)  return  1;
  if (id ==  11 || id ==  13 || id == -211 || id == -321 || id == -2212)  return -1;
  return 0;
}


//-----------------------------------------------
// calculate |\Delta phi| between two 4-vectors
//
double Dphi(const fastjet::PseudoJet & jet1, const fastjet::PseudoJet & jet2)
{
  double dphi = jet2.phi_std() - jet1.phi_std();
  if (dphi > pi)  dphi -= 2*pi;
  if (dphi < -pi) dphi += 2*pi;
  return fabs(dphi);
}


//----------------------------------------------------------------------
/// calculate \Delta R between two objects
double DR(const fastjet::PseudoJet & jet1, const fastjet::PseudoJet & jet2) {

  double dy = jet2.rap() - jet1.rap();
  double dphi = jet2.phi_std() - jet1.phi_std();

  if (dphi > pi)  dphi -= 2*pi;
  if (dphi < -pi) dphi += 2*pi;
  
  return sqrt(dy*dy + dphi*dphi);
}


//-----------------------------------------------------------------------------------
// determine the signed heavy quark ID (bottom or charm) inside of a heavy hadron
int HFID(int id)
{
  // note that I only store HF hadrons immediately before the heavy quark decays,
  // so don't worry about, e.g., radially-excited mesons which quickly release
  // pions/photons.  also, I don't independently store charms that came from bottom decay.
  if (id == 0)  return 0;
  int absid = abs(id);
  if (absid <= 10)  return id;  // bare quark
  int digit2 = ( absid-1000 *(absid/1000 ) - (absid%100 ) ) / 100 ;
  int digit3 = ( absid-10000*(absid/10000) - (absid%1000) ) / 1000;
  int flavor = max(digit2,digit3) * id/abs(id);
  if (absid < 1000 && abs(flavor) == 5)  flavor *= -1;  // meson codes are backward for down-type (e.g., 521 contains anti-b)
  return flavor;
}

//-----------------------------------
// b-tagging and c-tagging
//
void bcTag(vector<fastjet::PseudoJet> & jets, const vector<fastjet::PseudoJet> & bcHadrons, double Rtag)
{
  // first, find the unique closest jet to each b hadron
  vector<int> jetIndices;
  for (unsigned int itb = 0; itb < bcHadrons.size(); itb++) {
    double DRmin = 100000.;
    int indexMin = -1;
    for (unsigned int itj = 0; itj < jets.size(); itj++) {
      if (jets[itj].perp() < 1.0)  continue;  // don't try to match to "blank" jets
      double DRlocal = DR(bcHadrons[itb],jets[itj]);
      if (DRlocal < DRmin)  {  DRmin = DRlocal;  indexMin = itj; }
    }
    if (DRmin < Rtag)   jetIndices.push_back(indexMin);
    else                jetIndices.push_back(-1);
  }

  // next, check for hard b/c hadrons matched to each jet (pick off the hardest flavor)
  for (unsigned int itj = 0; itj < jets.size(); itj++) {
    double ptLocal;
    double ptMax = 0;
    int idLead = 7;  // 1 = unflavored (default)
    for (unsigned int itb = 0; itb < jetIndices.size(); itb++) {
      if (jetIndices[itb] == (int)itj) {
	ptLocal = bcHadrons[itb].perp();
	if (ptLocal > ptMax)  { ptMax = ptLocal; idLead = HFID(bcHadrons[itb].user_index()); }//Labels the jet with the b jet label (5)
      }
    }
    jets[itj].set_user_index(idLead);//labels the jet with 7
  }
}




//----------------------------------------------------------------------
/// master b-tag function, used with a global switch to select which method to use
/*
bool isAbjet(const fastjet::PseudoJet & theJet, const fastjet::ClusterSequence & theClusterSequence, double ptb, double R) {

  if (fabs(theJet.eta()) > etabcut)  return false;   // cutoff of LHC tracking

  if (btagMethod == 1)
    return nearbyB(theJet,R);

  if (btagMethod == 2)
    return (ptb > 0);

  return false;
}



//----------------------------------------------------------------------
/// b-tag #1:  determine whether this jet is a b-jet by searching for nearby b hadrons
bool nearbyB(const fastjet::PseudoJet & theJet, double R) {

  double Rmax = R*matchDist;

  for (unsigned int it = 0; it < b_hadrons.size(); it++)
    if ( DR(theJet,b_hadrons[it]) < Rmax)
      return true;

  return false;
}


//----------------------------------------------------------------------
/// btag #2:  match jets uniquely to b-hadrons.  determine the (scalar-summed) b hadron pt's within each jet (0 => no b's) and DR of leading b hadron
void assignBs(const vector<fastjet::PseudoJet> & theJets, vector<double> & ptbs, vector<double> & DRbs, vector<int> & idbs, double R) {

  double Rmax = R*searchDist;

  vector<fastjet::PseudoJet> b_particles = b_hadrons;
  //if (bQuarkTag) b_particles = b_quarks;

  // first, find the closest jet to each b hadron
  vector<int> jetIndices;
  vector<double> jetDRs;
  for (unsigned int itb = 0; itb < b_particles.size(); itb++) {
    double DRmin = 100000.;
    int indexMin = -1;
    for (unsigned int itj = 0; itj < theJets.size(); itj++) {
      double DRlocal = DR(b_particles[itb],theJets[itj]);
      if (DRlocal < DRmin)  {  DRmin = DRlocal;  indexMin = itj; }
    }
    if (DRmin < Rmax)  {  jetIndices.push_back(indexMin);  jetDRs.push_back(DRmin); }
    else               {  jetIndices.push_back(-1);        jetDRs.push_back(-1);    }
  }

  // next, (scalar) sum up the b hadron pt within each jet and pay attention to the DR of the leader.
  for (unsigned int itj = 0; itj < theJets.size(); itj++) {
    double ptblocal = 0.;
    double ptblocalmax = 0;
    double DRlead = -1;
    int idLead = 0;
    for (unsigned int itb = 0; itb < jetIndices.size(); itb++) {
      if (jetIndices[itb] == itj) {
	ptblocal += b_particles[itb].perp();
	if (ptblocal > ptblocalmax)  { ptblocalmax = ptblocal; DRlead = jetDRs[itb]; idLead = b_particles[itb].user_index(); }
      }
    }
    ptbs.push_back(ptblocal);
    DRbs.push_back(DRlead);
    if (idLead != 0) {
      if (abs(idLead) > 1000)  idLead =  5*abs(idLead)/idLead;
      else                     idLead = -5*abs(idLead)/idLead;  // mesons codes are backward...e.g., 521 contains an anti-b    
    }
    idbs.push_back(idLead);
  }
}
*/


