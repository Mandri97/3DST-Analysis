using namespace std;

// Histograms {{{

// Energy deposit of neutrons (signal and background)
auto hist_sig_edep	     = new TH1F("hist_sig_edep",        "Neutron energy deposit;Energy [MeV];Events", 35, 0, 35);
auto hist_bkg_edep	     = new TH1F("hist_bkg_edep",        "Neutron energy deposit;Energy [MeV];Events", 35, 0, 35);

// Time-of-Flight of neutrons (signal and background)
auto hist_sig_time	     = new TH1F("hist_sig_time",        "Time distribution; Time [ns];Events", 25, 0, 25);
auto hist_bkg_time	     = new TH1F("hist_bkg_time",        "Time distribution; Time [ns]; Events", 25, 0, 25);

auto hist_sig_eRec	     = new TH1F("hist_sig_eRec", "      Neutron reconstructed energy; Energy [GeV]; Events", 20, 0, 1);
auto hist_bkg_eRec	     = new TH1F("hist_bkg_eRec",        "Neutron reconstructed energy; Energy [GeV]; Events", 20, 0, 1);
auto hist_sig_eRecPrimary    = new TH1F("hist_sig_eRecPrimary", "Neutron reconstructed energy; Energy [GeV]; Events", 20, 0, 1);

auto hist_sig_SmearTrue	     = new TH1F("hist_sig_SmearTrue",   "(E_smear - E_true)/E_true; Energy [Gev]; Events", 40, -1, 1);
auto hist_bkg_SmearTrue	     = new TH1F("hist_bkg_SmearTrue",   "(E_smear - E_true)/E_true; Energy [Gev]; Events", 40, -1, 1);

auto hist_sig_MeanSmear	     = new TH1F("hist_sig_MeanSmear",   "Mean (E_smear - E_true)/E_true; Lever Arm [cm]; Energy [Gev]", 20, 0, 200);
auto hist_bkg_MeanSmear	     = new TH1F("hist_bkg_MeanSmear",   "Mean (E_smear - E_true)/E_true; Lever Arm [cm]; Energy [Gev]", 20, 0, 200);

auto hist_sig_StdSmear	     = new TH1F("hist_sig_StdSmear",    "Std (E_smear - E_true)/E_true; Lever Arm [cm]; Energy [Gev]", 20, 0, 200);
auto hist_bkg_StdSmear	     = new TH1F("hist_bkg_StdSmear",    "Std (E_smear - E_true)/E_true; Lever Arm [cm]; Energy [Gev]", 20, 0, 200);

auto hist_sig_SmearLvr	     = new TH2F("hist_sig_SmearLvr",    "(E_smear - E_true)/E_true vs. Lever Arm; Lever Arm [cm]; Energy [Gev]", 20, 0, 200, 40, -1, 1);
auto hist_bkg_SmearLvr	     = new TH2F("hist_bkg_SmearLvr",    "(E_smear - E_true)/E_true vs. Lever Arm; Lever Arm [cm]; Energy [Gev]", 20, 0, 200, 40, -1, 1); 
auto hist_sig_Lvr	     = new TH1F("hist_sig_Lvr",         "Lever Arm; Length [cm]; Events", 20, 0, 200);
auto hist_bkg_Lvr	     = new TH1F("hist_bkg_Lvr",         "Lever Arm; Length [cm]; Events", 20, 0, 200);

auto hist_sig_LvrTime	     = new TH2F("hist_sig_LvrTime",     "Lever Arm vs. Time; Lever Arm [cm]; Time [ns]", 300, 0, 200, 40, 0, 20);
auto hist_bkg_LvrTime	     = new TH2F("hist_bkg_LvrTime",     "Lever Arm vs. Time; Lever Arm [cm]; Time [ns]", 300, 0, 200, 40, 0, 20);

auto hist_bkg_CYL	     = new TH1F("hist_bkg_CYL",         "Background from Cylender; Energy [MeV]; Events", 35, 0, 35);
auto hist_bkg_3DS	     = new TH1F("hist_bkg_3DS",         "Background from 3DST; Energy [MeV]; Events", 35, 0, 35);
auto hist_bkg_Roc	     = new TH1F("hist_bkg_Roc",         "Background from Rock; Energy [MeV]; Events", 35, 0, 35);

// Purity, Resolution
auto hist_purity_ArmTime     = new TH2F("hist_purity_ArmTime",   "Purity; Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
auto hist_resolution_ArmTime = new TH2F("hist_resolution_ArmTime","Standard Deviation Energy; Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);

// Number of cube fired per event and hit in a cube per event
auto hist_cubeFire	             = new TH1F("hist_cubeFire", "Cube fired per event", 100, 0, 100);
auto hist_hitMultiplicity	     = new TH1F("hist_hitMultiplicity", "Hit Multiplicity", 25, 0, 25);

// Neutrino vertex of neutron background
auto hist_bkg_vtx_neutrino = new TH3F("hist_bkg_vtx_neutrino", "Background neutrino interaction vertex; X[cm];Y[cm];Z[cm]",
					100, -500, 500, 1000, -500, 500, 2100, -1500, 600);
auto hist_bkg_vtx_yx       = new TH2F("hist_bkg_vtx_yx", "Background neutrino interaction vertex; Y[cm]; X[cm]",
					1000, -500, 500, 1000, -500, 500);
auto hist_bkg_vtx_zy       = new TH2F("hist_bkg_vtx_zy", "Background neutrino interaction vertex; Z[cm]; Y[cm]",
					2100, -1500, 600, 1000, -500, 500);
auto hist_bkg_vtx_zx       = new TH2F("hist_bkg_vtx_zx", "Background neutrino interaction vertex; Z[cm]; X[cm]",
					2100, -1500, 600, 100, -500, 500);

// Separete time, edep and lever arm into primary and secondary neutron
// if it has a ParentID == 0 || ParentID == -1 --> Primary
//						   Secondary otherwise
auto hist_sigTime_primary = new TH1F("hist_sig_time_noPi", "Time distribution; Time [ns];Events", 25, 0, 25);
auto hist_sigLvr_primary  = new TH1F("hist_sig_Lvr_noPi",  "Lever Arm; Track Length [cm]; Events", 20, 0, 200);
auto hist_sigEdep_primary = new TH1F("hist_sig_edep_noPi", "Neutron energy deposit;Energy [MeV];Events", 35, 0, 35);

auto hist_sigTime_secondary = new TH1F("hist_sigTime_secondary", "Time distribution; Time [ns];Events", 25, 0, 25);
auto hist_sigLvr_secondary  = new TH1F("hist_sigLvr_secondary",  "Lever Arm; Track Length [cm]; Events", 20, 0, 200);
auto hist_sigEdep_secondary = new TH1F("hist_sigEdep_secondary", "Neutron energy deposit;Energy [MeV];Events", 35, 0, 35);

auto hist_bkgTime_primary = new TH1F("hist_bkg_time_noPi", "Time distribution; Time [ns];Events", 25, 0, 25);
auto hist_bkgLvr_primary  = new TH1F("hist_bkg_Lvr_noPi",  "Lever Arm; Track Length [cm]; Events", 20, 0, 200);
auto hist_bkgEdep_primary = new TH1F("hist_bkg_edep_noPi", "Neutron energy deposit;Energy [MeV];Events", 35, 0, 35);

auto hist_bkgTime_secondary = new TH1F("hist_bkgTime_secondary", "Time distribution; Time [ns];Events", 25, 0, 25);
auto hist_bkgLvr_secondary  =  new TH1F("hist_bkgLvr_secondary", "Lever Arm; Track Length [cm]; Events", 20, 0, 200);
auto hist_bkgEdep_secondary = new TH1F("hist_bkgEdep_secondary", "Neutron energy deposit;Energy [MeV];Events", 35, 0, 35);

// }}}

// count
long int bkg = 0,
	 sig = 0,
	 nbCC  = 0;

int en001 = 0;

// TODO: Change this according to the energy deposit threshold
float energyHitCut = 0.5;

// Structure that will contains neutron hit information {{{
struct Hit_t {
    float timeWindow,           // time windows of the hit
	  timeSmear,		// smear time
	  energyDeposit,        // energy deposited by the neutron
	  trackLength,	        // lever arm
	  trueRec,		// true reconstructed energy
	  smearRec,		
	  vtxSignal[3],		// neutrino vertex position of the neutron
	  vtxTime;		// neutrino  vertex time

    int bkgLoc,			// neutrino vertex position
	neutronParentID,	// Where the neutron come from 
	neutronParentPdg;	// PDG of neutron parent

    bool isTherePion50,		// Is there a pion with KE > 50 MeV in FS particles
	 isThereProton300;      // Is there a proton with KE > 300 MeV in FS particles
};
// }}}

// Statistics class {{{
// This class calculate the Mean and RMS on the fly.
// Will add the link later
class Statistics_t{
    private:
	int n;
	float oldMean, newMean, oldStd, newStd;
	vector<float> values;

    public:
	Statistics_t(): n(0) {}

	void Push(float x){
	    n++;

	    if (n == 1){
		oldMean = newMean = x;
		oldStd = 0.0;
	    }else{
		newMean = oldMean + (x - oldMean)/n;
		newStd  = oldStd  + (x - oldMean) * (x - newMean);

		oldMean = newMean;
		oldStd  = newStd;
	    }
	}

	float Mean() const{
	    return (n > 0) ? newMean : 0.0;
	}

	float Variance() const{
	    return (n > 1) ? newStd / n  : 0.0;
	}

	float Std() const{
	    return pow(Variance(), 0.5);
	}

	int GetN() const{
	    return n;
	}
};
//}}}

// Statistical variable per Bin {{{
// Contains Hit sorted by Lever arm and Time respecticely per bin
vector<Hit_t> sig_Neutron_Lvr[20], sig_Neutron_Time[25];
vector<Hit_t> bkg_Neutron_Lvr[20], bkg_Neutron_Time[25];

// Testing variables
Statistics_t  sig_Stats_Lvr[20], bkg_Stats_Lvr[20];
vector<float> sig_Value_Lvr[20], bkg_Value_Lvr[20], sig_Value_Time[25], bkg_Value_Time[25];
map<string, int> countCubeHit;
Statistics_t sig_Stats_TrueE[30];
//}}}


// Function {{{
float recEnergy(Float_t l, Float_t t){
    return pow((723 * l)/t, 2) * 0.000000001;
}

float kineticEnergy(Float_t px, Float_t py, Float_t pz, Float_t e){
    float mass = pow(e*e - px*px - py*py - pz*pz, 0.5);
    return e - mass;
}

// Neutron origin
// TODO: Change this according to the geometry
int neutronOrigin(float x, float y, float z){
    if (abs(x) < 400 && abs(y) < 400 && z < -400 && z > -1400)     return  1; // Cylender
    else if (z < -1400 || z > 400 || abs(x) > 400 || abs(y) > 400) return -1; // rock
    else if (abs(x) < 170 && abs(y) < 170 && z > -100 && z < 150)  return -10; // 3DST + TPC
    else if (abs(x) < 250 && abs(y) < 250 && z < 200)              return   2; // ECAL + Magnet

    return -10;
}

// }}}

void analyze(string file){
    auto _file = new TFile(TString(file));
    auto tree = (TTree*)_file->Get("tree");

    // no tree
    if (tree == NULL){
	_file->Close();
	return ;
    }

    int nevents = tree->GetEntries();

    Float_t t_fsPx[1000], t_fsPy[1000], t_fsPz[1000], t_fsE[1000];

    Float_t t_neutronHitX[1000],      t_neutronHitY[1000],     t_neutronHitZ[1000], 
	    t_neutronHitE[1000],      t_neutronHitPDG[1000],   t_neutronHitT[1000],
	    t_neutronHitSmearT[1000], t_neutronParentID[1000], t_neutronParentPdg[1000];

    Float_t t_vtx[3], t_vtxTime, t_PionKE;

    Int_t   PDG, 
	    t_nFS,
	    t_fsPdg[1000];

    tree->SetBranchAddress("vtx",	       &t_vtx);
    tree->SetBranchAddress("nFS",	       &t_nFS);
    tree->SetBranchAddress("fsE",	       &t_fsE);
    tree->SetBranchAddress("fsPx",	       &t_fsPx);
    tree->SetBranchAddress("fsPy",	       &t_fsPy);
    tree->SetBranchAddress("fsPz",	       &t_fsPz);
    tree->SetBranchAddress("piKE",	       &t_PionKE);
    tree->SetBranchAddress("fsPdg",	       &t_fsPdg);
    tree->SetBranchAddress("vtxTime",	       &t_vtxTime);
    tree->SetBranchAddress("neutronHitX",      &t_neutronHitX);
    tree->SetBranchAddress("neutronHitY",      &t_neutronHitY);
    tree->SetBranchAddress("neutronHitZ",      &t_neutronHitZ);
    tree->SetBranchAddress("neutronHitE",      &t_neutronHitE);
    tree->SetBranchAddress("neutronHitT",      &t_neutronHitT);
    tree->SetBranchAddress("neutronHitPDG",    &t_neutronHitPDG);
    tree->SetBranchAddress("neutronParentId",  &t_neutronParentID);
    tree->SetBranchAddress("neutronHitSmearT", &t_neutronHitSmearT);
    tree->SetBranchAddress("neutronParentPDG", &t_neutronParentPdg);

    // Flag that indicate if a signal or background was detected in a spill 
    bool  is_Sig = false, 
	  is_Bkg = false;

    // Variable for the selected (true) signal and background that is detected
    Hit_t sig_earliestHit, bkg_earliestHit;

    // Initialization
    sig_earliestHit.timeWindow = 1000;
    bkg_earliestHit.timeWindow = 1000;
   
    // Loop over all events
    for (int ievent = 0; ievent < nevents; ievent++){
	tree->GetEntry(ievent);
    
	bool isTherePionAbove50    = false;
	bool isThereProtonAbove300 = false;

	// only study event in Fiducial Volume
	// TODO: Change this according to the geometry
	if (abs(t_vtx[0]) < 50 && abs(t_vtx[1]) < 50 && abs(t_vtx[2]) < 50){

	    // Flag to a CC event
	    bool is_CC = false;

	    // Searching for a muon or electron/positron
	    for (int inFS = 0; inFS < t_nFS; inFS++){
		if (abs(t_fsPdg[inFS]) == 11 || abs(t_fsPdg[inFS]) == 13){
		    is_CC = true;
		    break;
		}
	    }

	    // skip to the next event if there is no CC event
	    if (!is_CC) continue;
	    nbCC++;

	    // Check if there is a pion with KE > 50 MeV in FS particles
	    isTherePionAbove50 = t_PionKE > 50;

	    // Check if there is a proton with KE > 300 MeV in FS particles
	    for (int inFS = 0; inFS < t_nFS; inFS++){
		if (t_fsPdg[inFS] == 2212){
		    float ke = kineticEnergy(t_fsPx[inFS], t_fsPy[inFS], t_fsPz[inFS], t_fsE[inFS]);
		    isThereProtonAbove300 =  ke > 200;
		    break;
		}
	    }

	    int   nbCubeFired  = 0;
	    map<string, Hit_t> hitPerCube;

	    /* SIGNAL: Neutron Information
	     ---------------------------
	     Look for neutron hit induced by a CC event in the FV across all the cube in the 3DST 
	     Then select the earliest activated cube as a signal
	     */

	    for (int ineutronHit = 0; ineutronHit < 100; ineutronHit++){
		// Look for neutron hit across all the 3DST 
		// TODO: Change this according to the geometry
		if (abs(t_neutronHitX[ineutronHit]) < 120 && abs(t_neutronHitY[ineutronHit]) < 120 &&
		    abs(t_neutronHitZ[ineutronHit]) < 100){
		    
		    // calculate the lever arm for signal
		    float trackLength = pow(
			pow(t_vtx[0] - t_neutronHitX[ineutronHit], 2) + 
			pow(t_vtx[1] - t_neutronHitY[ineutronHit], 2) +
			pow(t_vtx[2] - t_neutronHitZ[ineutronHit], 2), 0.5);

		    // calculate the signal window
		    float signalWindow	    = t_neutronHitT[ineutronHit]      - t_vtxTime;
		    float signalWindowSmear = t_neutronHitSmearT[ineutronHit] - t_vtxTime;
		    
		    // Fix a bug from edep-sim
		    if (signalWindow      == 1) signalWindow      = 0.5;
		    if (signalWindowSmear == 1) signalWindowSmear = 0.5;

		    // Affect the neutron hit to Hit_t structure
		    if (signalWindow > 0){
			Hit_t temp;

			// basic kinetic variabel
			temp.trueRec	   = recEnergy(trackLength, signalWindow);
			temp.smearRec	   = recEnergy(trackLength, signalWindowSmear);
			temp.timeSmear     = signalWindowSmear;
			temp.timeWindow    = signalWindow;
			temp.trackLength   = trackLength;
			temp.energyDeposit = t_neutronHitE[ineutronHit];

			temp.vtxSignal[0]  = t_vtx[0];
			temp.vtxSignal[1]  = t_vtx[1];
			temp.vtxSignal[2]  = t_vtx[2];

			temp.isTherePion50    = isTherePionAbove50;
			temp.isThereProton300 = isThereProtonAbove300;
			temp.neutronParentID  = t_neutronParentID[ineutronHit];
			temp.neutronParentPdg = t_neutronParentPdg[ineutronHit];

			temp.vtxTime = t_vtxTime;

			/*
			 The key here is the position of the cube.
			 Be aware that t_neutronHitX, t_neutronHitY and t_neutronHitZ
			 must all be in cm.
			 */
			string key = string(Form("%d_%d_%d", (int) t_neutronHitX[ineutronHit],
							     (int) t_neutronHitY[ineutronHit],
							     (int) t_neutronHitZ[ineutronHit]));

			// hit multiplicity per cube
			auto findKey_mul = countCubeHit.find(key);
			if (findKey_mul != countCubeHit.end()) countCubeHit.at(key) += 1;
			else countCubeHit[key] = 1;

			/*
			+ If the cube has already been activated by a neutron
			    - See which neutron hit the cube first
			    - Affect the earliest neutron hit to the cube
			    - Sum up the energy deposit in the cube
			+ Affect the neutron hit to the cube otherwise
			*/
			auto findKey_hitCubeEvent = hitPerCube.find(key);
			if (findKey_hitCubeEvent != hitPerCube.end()){
			    if (hitPerCube.at(key).timeWindow < temp.timeWindow)
				hitPerCube.at(key).energyDeposit += temp.energyDeposit;
			    else{
				auto tempEnergy = hitPerCube.at(key).energyDeposit;
				hitPerCube.at(key) = temp;
				hitPerCube.at(key).energyDeposit += tempEnergy;
			    }
			}else{
			    hitPerCube[key] = temp;
			    nbCubeFired++;
			} 
		
			// Raise signal flag
			is_Sig = true;
		    }
		}
	    }

	    // Cube fire per event
	    if (nbCubeFired) hist_cubeFire->Fill(nbCubeFired);

	    // Searching for the first activated cube in the 3DST with an energy deposit > 0.5 MeV
	    for (auto hit : hitPerCube){
		if (hit.second.energyDeposit < 0.01) en001++;

		if (sig_earliestHit.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut) 
		    sig_earliestHit = hit.second;
	    }
	}
    }

    // Skip to the next spill if no signal were found
    if (!is_Sig){
	_file->Close();
	return;
    }

    /*
    BACKGROUND: Neutron information
    -------------------------------
    + Search for neutron in the whole 3DST induced by neutrino in:
      1. NC event in the out-of-FV region
      2. All event outside the 3DST (CC + NC)
    + Then check if the background happens inside the signal window
    + If so, select the background, signal otherwise
    */

    for (int ievent = 0; ievent < nevents; ievent++){
	tree->GetEntry(ievent);

	/*
	`outFV_inseide3DST`: The region between 3DST and out-of FV
	`out3DST`: Region outside the whole 3DST
	*/
	bool outFV_inside3DST = false,
	     out3DST          = false;

	// TODO: Change this according to the size of 3DST and FV	
	if (abs(t_vtx[0]) < 120 && abs(t_vtx[1]) < 120 && abs(t_vtx[2]) < 100){
	    if (abs(t_vtx[0]) > 50 || abs(t_vtx[1]) > 50 || abs(t_vtx[2]) > 50) outFV_inside3DST = true;
	}

	
	// TODO: Change this according to the geometry
	if (abs(t_vtx[0]) > 120 || abs(t_vtx[1]) > 120 || abs(t_vtx[2]) > 100) out3DST = true;

	/*
	`NC`: Neutral Current event
	`CC`: Charge Current event
	*/
	bool NC    = false,
	     CC    = false;

	// look for NC events
	if (outFV_inside3DST){
	    for (int inFS = 0; inFS < t_nFS; inFS++){
		// look for muon or electron
		if (abs(t_fsPdg[inFS]) == 11 || abs(t_fsPdg[inFS]) == 13){
		    CC = true;
		    break;
		}
	    }

	    // no muon or electron found
	    if (!CC){
		//look for pion
		for (int inFS = 0; inFS < t_nFS; inFS++){
		    if (abs(t_fsPdg[inFS]) == 211 || abs(t_fsPdg[inFS]) == 111){
			break;
			NC = true;
		    }
		}
	    }
	}

	// Search
	if ((outFV_inside3DST && NC) || out3DST){

	    map<string, Hit_t> hitPerCube;

	    for (int ineutronHit = 0; ineutronHit < 100; ineutronHit++){
		// TODO: Change this according to the geometry
		if (abs(t_neutronHitX[ineutronHit]) < 120 && abs(t_neutronHitY[ineutronHit]) < 120 &&
		    abs(t_neutronHitZ[ineutronHit]) < 100){

		    // Refering the background from FV signal vertex
		    float trackLength = pow(
			    pow(sig_earliestHit.vtxSignal[0] - t_neutronHitX[ineutronHit], 2) + 
			    pow(sig_earliestHit.vtxSignal[1] - t_neutronHitY[ineutronHit], 2) +
			    pow(sig_earliestHit.vtxSignal[2] - t_neutronHitZ[ineutronHit], 2), 0.5);			    

		    float backgroundWindow	= t_neutronHitT[ineutronHit]      - t_vtxTime;
		    float backgroundWindowSmear = t_neutronHitSmearT[ineutronHit] - t_vtxTime;

		    // fix bug from edep-sim
		    if (backgroundWindow      == 1) backgroundWindow      = 0.5;
		    if (backgroundWindowSmear == 1) backgroundWindowSmear = 0.5;

		    // Affect the neutron hit to Hit_t structure
		    if (backgroundWindow > 0){
			Hit_t temp;

			temp.timeWindow    = backgroundWindow;
			temp.timeSmear     = backgroundWindowSmear;
			temp.energyDeposit = t_neutronHitE[ineutronHit];
			temp.trackLength   = trackLength;
			temp.trueRec	   = recEnergy(trackLength, backgroundWindow);
			temp.smearRec	   = recEnergy(trackLength, backgroundWindowSmear);

			// Where the neutron comes from ?
			temp.bkgLoc = neutronOrigin(t_vtx[0], t_vtx[1], t_vtx[2]);

			temp.vtxSignal[0] = t_vtx[0];
			temp.vtxSignal[1] = t_vtx[1];
			temp.vtxSignal[2] = t_vtx[2];

			string key = string(Form("%d_%d_%d", (int) t_neutronHitX[ineutronHit],
							     (int) t_neutronHitY[ineutronHit],
							     (int) t_neutronHitZ[ineutronHit]));
			auto findKey = hitPerCube.find(key);

			/*
			+ If the cube has already been activated by a neutron
			    - See which neutron hit the cube first
			    - Affect the earliest neutron hit to the cube
			    - Sum up the energy deposit in the cube
			+ Affect the neutron hit to the cube otherwise
			*/
			if (findKey != hitPerCube.end()){
			    if (hitPerCube.at(key).timeWindow < temp.timeWindow)
				hitPerCube.at(key).energyDeposit += temp.energyDeposit;
			    else{
				auto tempEnergy = hitPerCube.at(key).energyDeposit;
				hitPerCube.at(key) = temp;
				hitPerCube.at(key).energyDeposit += tempEnergy;
			    }
			}else{
			    hitPerCube[key] = temp;
			} 

		    }
		}
	    }

	    // looking for the earliest hit
	    for (auto hit : hitPerCube){
		if (bkg_earliestHit.timeWindow > hit.second.timeWindow && hit.second.energyDeposit > energyHitCut)
		    bkg_earliestHit = hit.second;
	    }
	}
    }

    // Activate background flag is the background is happening inside the signal window
    if (sig_earliestHit.timeWindow > bkg_earliestHit.timeWindow) is_Bkg = true;

// Fill histogram {{{
    if (is_Bkg){
	// Background vertex
	hist_bkg_vtx_neutrino->Fill(bkg_earliestHit.vtxSignal[0], bkg_earliestHit.vtxSignal[1], bkg_earliestHit.vtxSignal[2]);
	hist_bkg_vtx_yx->Fill(bkg_earliestHit.vtxSignal[1], bkg_earliestHit.vtxSignal[0]);
	hist_bkg_vtx_zy->Fill(bkg_earliestHit.vtxSignal[2], bkg_earliestHit.vtxSignal[1]);
	hist_bkg_vtx_zx->Fill(bkg_earliestHit.vtxSignal[2], bkg_earliestHit.vtxSignal[0]);

	hist_bkg_Lvr     ->Fill(bkg_earliestHit.trackLength);
	hist_bkg_time    ->Fill(bkg_earliestHit.timeWindow);
	hist_bkg_eRec    ->Fill(bkg_earliestHit.trueRec);
	hist_bkg_edep    ->Fill(bkg_earliestHit.energyDeposit);
	hist_bkg_LvrTime ->Fill(bkg_earliestHit.trackLength, bkg_earliestHit.timeWindow);
	
	float bkg_SmearTrue = (bkg_earliestHit.smearRec - bkg_earliestHit.trueRec) / bkg_earliestHit.trueRec;

	hist_bkg_SmearTrue ->Fill(bkg_SmearTrue);
	hist_bkg_SmearLvr  ->Fill(bkg_earliestHit.trackLength, bkg_SmearTrue);

	// Sort background hit by lever arm (10 cm bin width) and time window
	if (bkg_earliestHit.trackLength < 200 && bkg_earliestHit.timeWindow < 25){
	    int iTime = (int) bkg_earliestHit.timeWindow;
	    int iArm  = (int) bkg_earliestHit.trackLength / 10;

	    bkg_Neutron_Lvr[iArm].push_back(bkg_earliestHit);
	    bkg_Neutron_Time[iTime].push_back(bkg_earliestHit);
	}

	// Where the background hit comes from
	if (bkg_earliestHit.bkgLoc ==  2) hist_bkg_3DS->Fill(bkg_earliestHit.energyDeposit);
	if (bkg_earliestHit.bkgLoc ==  1) hist_bkg_CYL->Fill(bkg_earliestHit.energyDeposit);
	if (bkg_earliestHit.bkgLoc == -1) hist_bkg_Roc->Fill(bkg_earliestHit.energyDeposit);

	if (bkg_earliestHit.trackLength < 170){
	    int index = (int) bkg_earliestHit.trackLength / 10; // ---> index in [0, 20[ which correspond to the 10 bin width

	    if (abs(bkg_SmearTrue) < 1) bkg_Stats_Lvr[index].Push(bkg_SmearTrue); 
	}

	// separation bkgnal 1D
	if (bkg_earliestHit.neutronParentID == -1 || bkg_earliestHit.neutronParentID == 0){
	    hist_bkgTime_primary->Fill(bkg_earliestHit.timeWindow);
	    hist_bkgLvr_primary->Fill(bkg_earliestHit.trackLength);
	    hist_bkgEdep_primary->Fill(bkg_earliestHit.energyDeposit);
	}else{
	    hist_bkgTime_secondary->Fill(bkg_earliestHit.timeWindow);
	    hist_bkgLvr_secondary->Fill(bkg_earliestHit.trackLength);
	    hist_bkgEdep_secondary->Fill(bkg_earliestHit.energyDeposit);
	}

	bkg++;
    }else if (is_Sig){
	hist_sig_Lvr     ->Fill(sig_earliestHit.trackLength);
	hist_sig_time    ->Fill(sig_earliestHit.timeWindow);
	hist_sig_eRec    ->Fill(sig_earliestHit.trueRec);
	hist_sig_edep    ->Fill(sig_earliestHit.energyDeposit);
	hist_sig_LvrTime ->Fill(sig_earliestHit.trackLength, sig_earliestHit.timeWindow);

	float sig_SmearTrue = (sig_earliestHit.smearRec - sig_earliestHit.trueRec) / sig_earliestHit.trueRec;
	hist_sig_SmearTrue ->Fill(sig_SmearTrue);
	hist_sig_SmearLvr  ->Fill(sig_earliestHit.trackLength, sig_SmearTrue);
	if (sig_earliestHit.trackLength < 200 && abs(sig_SmearTrue) < 1){
	    int index = (int) sig_earliestHit.trackLength / 10; 

	    sig_Stats_Lvr[index].Push(sig_SmearTrue); 
	    sig_Value_Lvr[index].push_back(sig_SmearTrue);
	}

	// Sort signal hit by lever arm (10 cm bin width) and time window
	if (sig_earliestHit.trackLength < 200 && sig_earliestHit.timeWindow < 25){
	    int iTime = (int) sig_earliestHit.timeWindow;
	    int iArm  = (int) sig_earliestHit.trackLength / 10;

	    sig_Neutron_Lvr[iArm].push_back(sig_earliestHit);
	    sig_Neutron_Time[iTime].push_back(sig_earliestHit);
	}

	// resolution of primary neutron along true energy
	if (sig_earliestHit.neutronParentID == 0 || sig_earliestHit.neutronParentID == -1){
	    hist_sig_eRecPrimary->Fill(sig_earliestHit.trueRec);

	    // calculating index with 0.33 MeV bin width (0 - 1 GeV)
	    int index = (int) ((sig_earliestHit.trueRec * 10) / 0.33);
	    float resolution = (sig_earliestHit.smearRec - sig_earliestHit.trueRec) / sig_earliestHit.trueRec;

	    if (abs(resolution) < 1) sig_Stats_TrueE[index].Push(resolution);
	}

	if (sig_earliestHit.timeWindow < 25){
	    int index = (int) sig_earliestHit.timeWindow;

	   if (abs(sig_SmearTrue) < 1) sig_Value_Time[index].push_back(sig_SmearTrue);
	}

	// separation of primary and secondary neutron signal
	if (sig_earliestHit.neutronParentID == -1 || sig_earliestHit.neutronParentID == 0){
	    hist_sigTime_primary->Fill(sig_earliestHit.timeWindow);
	    hist_sigLvr_primary->Fill(sig_earliestHit.trackLength);
	    hist_sigEdep_primary->Fill(sig_earliestHit.energyDeposit);
	}else{
	    hist_sigTime_secondary->Fill(sig_earliestHit.timeWindow);
	    hist_sigLvr_secondary->Fill(sig_earliestHit.trackLength);
	    hist_sigEdep_secondary->Fill(sig_earliestHit.energyDeposit);
	}

	sig++;
    }
// }}}

    _file->Close();
}

void neutron(string filename){
    // Read all files
    ifstream input(filename);
    
    auto outName = filename.substr(0, 4);

    string file;
    // analyze each spill separately
    while (input >> file) analyze(file);

    gStyle->SetOptStat(kFALSE);
    gROOT->SetBatch(kTRUE);


// Customize histogram {{{
    /**************************************************************************************************************/
    auto can10 = new TCanvas("can10", "can", 900, 600);
    hist_sig_Lvr->SetFillStyle(3004);
    hist_sig_Lvr->SetLineWidth(2);
    hist_sig_Lvr->SetLineColor(kAzure + 3);
    hist_sig_Lvr->SetFillColor(kAzure + 4);
    hist_sig_Lvr->Draw();

    hist_bkg_Lvr->SetFillStyle(3005);
    hist_bkg_Lvr->SetLineWidth(2);
    hist_bkg_Lvr->SetFillColor(kTeal - 6);
    hist_bkg_Lvr->SetLineColor(kTeal - 5);
    hist_bkg_Lvr->Draw("same");

    auto legend10 = new TLegend(0.75, 0.8, 0.88, 0.87);
    legend10->AddEntry(hist_sig_Lvr, "Signal");
    legend10->AddEntry(hist_bkg_Lvr, "Background");
    legend10->Draw();
    can10->SaveAs(Form("neutron%s_%.1f_NC_CC_leverArm.root", outName.c_str(), energyHitCut));

    /**************************************************************************************************************/

    auto can1 = new TCanvas("can1", "can", 900, 600);
    hist_sig_time->SetFillStyle(3004);
    hist_sig_time->SetLineWidth(2);
    hist_sig_time->SetLineColor(kAzure + 3);
    hist_sig_time->SetFillColor(kAzure + 4);
    hist_sig_time->Draw();

    hist_bkg_time->SetFillStyle(3005);
    hist_bkg_time->SetLineWidth(2);
    hist_bkg_time->SetFillColor(kTeal - 6);
    hist_bkg_time->SetLineColor(kTeal - 5);
    hist_bkg_time->Draw("same");

    auto legend1 = new TLegend(0.75, 0.8, 0.88, 0.87);
    legend1->AddEntry(hist_sig_time, "Signal");
    legend1->AddEntry(hist_bkg_time, "Background");
    legend1->Draw();
    can1->SaveAs(Form("neutron%s_%.1f_NC_CC_time.root", outName.c_str(), energyHitCut));    

    /**************************************************************************************************************/

    auto can2 = new TCanvas("can2", "can", 900, 600);
    hist_sig_edep->SetFillStyle(3004);
    hist_sig_edep->SetLineWidth(2);
    hist_sig_edep->SetLineColor(kAzure + 3);
    hist_sig_edep->SetFillColor(kAzure + 4);
    hist_sig_edep->Draw();

    hist_bkg_edep->SetFillStyle(3005);
    hist_bkg_edep->SetLineWidth(2);
    hist_bkg_edep->SetFillColor(kTeal - 6);
    hist_bkg_edep->SetLineColor(kTeal - 5);
    hist_bkg_edep->Draw("same");

    auto legend2 = new TLegend(0.75, 0.8, 0.88, 0.87);
    legend2->AddEntry(hist_sig_edep, "Signal");
    legend2->AddEntry(hist_bkg_edep, "Background");
    legend2->Draw();
    can2->SaveAs(Form("neutron%s_%.1f_NC_CC_earliest.root", outName.c_str(), energyHitCut));
    
    /**************************************************************************************************************/

    auto can3 = new TCanvas("can3", "can", 900, 600);
    hist_sig_eRec->SetFillStyle(3004);
    hist_sig_eRec->SetLineWidth(2);
    hist_sig_eRec->SetLineColor(kAzure + 3);
    hist_sig_eRec->SetFillColor(kAzure + 4);
    hist_sig_eRec->Draw();

    hist_bkg_eRec->SetFillStyle(3005);
    hist_bkg_eRec->SetLineWidth(2);
    hist_bkg_eRec->SetFillColor(kTeal - 6);
    hist_bkg_eRec->SetLineColor(kTeal - 5);
    hist_bkg_eRec->Draw("same");

    auto legend3 = new TLegend(0.75, 0.8, 0.88, 0.87);
    legend3->AddEntry(hist_sig_eRec, "Signal");
    legend3->AddEntry(hist_bkg_eRec, "Background");
    legend3->Draw();
    can3->SaveAs(Form("neutron%s_%.1f_NC_CC_recEnergy.root", outName.c_str(), energyHitCut));

    /**************************************************************************************************************/

    auto can7 = new TCanvas("can7", "can", 900, 600);
    hist_sig_LvrTime->Draw("colz");
    can7->SaveAs(Form("neutron%s_%.1f_NC_CC_sig_2LeverTime.root", outName.c_str(), energyHitCut));

    auto can8 = new TCanvas("can8", "can", 900, 600);
    hist_bkg_LvrTime->Draw("colz");
    can8->SaveAs(Form("neutron%s_%.1f_NC_CC_bkg_2LeverTime.root", outName.c_str(), energyHitCut));

    /**************************************************************************************************************/

    auto can4 = new TCanvas("can4", "can", 900, 600);
    // Background
    hist_bkg_3DS->SetFillStyle(3005);
    hist_bkg_3DS->SetLineWidth(2);
    hist_bkg_3DS->SetFillColor(kTeal - 6);
    hist_bkg_3DS->SetLineColor(kTeal - 5);
    hist_bkg_3DS->Draw();
    can4->SaveAs(Form("neutron%s_%.1f_NC_CC_hit3DSTSystem.root", outName.c_str(), energyHitCut));
    
    cout << "Bkg from 3DST System: " << hist_bkg_3DS->GetEntries() << endl;

    auto can6 = new TCanvas("can6", "can", 900, 600);
    hist_bkg_Roc->SetFillStyle(3005);
    hist_bkg_Roc->SetLineWidth(2);
    hist_bkg_Roc->SetFillColor(kTeal - 6);
    hist_bkg_Roc->SetLineColor(kTeal - 5);
    hist_bkg_Roc->Draw();
    can6->SaveAs(Form("neutron%s_%.1f_NC_CC_hitRoc.root", outName.c_str(), energyHitCut));

    cout << "Bkg from rock: " << hist_bkg_Roc->GetEntries() << endl;

    auto can9 = new TCanvas("can9", "can", 900, 600);
    // Background
    hist_bkg_CYL->SetFillStyle(3005);
    hist_bkg_CYL->SetLineWidth(2);
    hist_bkg_CYL->SetFillColor(kTeal - 6);
    hist_bkg_CYL->SetLineColor(kTeal - 5);
    hist_bkg_CYL->Draw();
    can9->SaveAs(Form("neutron%s_%.1f_NC_CC_hitCylinder.root", outName.c_str(), energyHitCut));

    cout << "Bkg from cylinder: " << hist_bkg_CYL->GetEntries() << endl;
    /**************************************************************************************************************/

    auto can12 = new TCanvas("can12", "can", 900, 600);
    hist_sig_SmearTrue->SetFillStyle(3004);
    hist_sig_SmearTrue->SetLineWidth(2);
    hist_sig_SmearTrue->SetLineColor(kAzure + 3);
    hist_sig_SmearTrue->SetFillColor(kAzure + 4);
    hist_sig_SmearTrue->Draw();
    can12->SaveAs(Form("neutron%s_%.1f_NC_CC_diff_SmearTrue_rec.root", outName.c_str(), energyHitCut));

    /**************************************************************************************************************/

    auto can18 = new TCanvas("can18", "can", 900, 600);
    hist_sig_SmearLvr->SetFillStyle(3004);
    hist_sig_SmearLvr->SetLineWidth(2);
    hist_sig_SmearLvr->SetLineColor(kAzure + 3);
    hist_sig_SmearLvr->SetFillColor(kAzure + 4);
    hist_sig_SmearLvr->Draw("colz");
    can18->SaveAs(Form("neutron%s_%.1f_NC_CC_sig_2Lvr_vs_EnergyRatio.root", outName.c_str(), energyHitCut));

    auto can19 = new TCanvas("can19", "can", 900, 600);
    hist_bkg_SmearLvr->SetFillStyle(3005);
    hist_bkg_SmearLvr->SetLineWidth(2);
    hist_bkg_SmearLvr->SetLineColor(kTeal - 6);
    hist_bkg_SmearLvr->SetFillColor(kTeal - 5);
    hist_bkg_SmearLvr->Draw("colz");
    can19->SaveAs(Form("neutron%s_%.1f_NC_CC_bkg_2Lvr_vs_EnergyRatio.root", outName.c_str(), energyHitCut));

    /**************************************************************************************************************/

    for (int i = 0; i < 20; i++){
	hist_sig_MeanSmear->SetBinContent(i + 1, sig_Stats_Lvr[i].Mean());
	hist_bkg_MeanSmear->SetBinContent(i + 1, bkg_Stats_Lvr[i].Mean());

	hist_sig_StdSmear->SetBinContent(i + 1, sig_Stats_Lvr[i].Std());
	hist_bkg_StdSmear->SetBinContent(i + 1, bkg_Stats_Lvr[i].Std());
    }

    auto can20 = new TCanvas("can20", "can", 900, 600);
    hist_sig_MeanSmear->SetFillStyle(3004);
    hist_sig_MeanSmear->SetLineWidth(2);
    hist_sig_MeanSmear->SetLineColor(kAzure + 3);
    hist_sig_MeanSmear->SetFillColor(kAzure + 4);
    hist_sig_MeanSmear->Draw();
    can20->SaveAs(Form("neutron%s_%.1f_NC_CC_Mean_EnergyRatio.root", outName.c_str(), energyHitCut));

    auto h1 = new TH1F("h1" ,"h1", 30, 0, 1);
    auto h2 = new TH1F("h2" ,"h2", 30, 0, 1);

    for (int k = 0; k < 30; k++){
	h1->SetBinContent(k + 1, sig_Stats_TrueE[k].Mean());
	h2->SetBinContent(k + 1, sig_Stats_TrueE[k].Std());
    }
    
    /**************************************************************************************************************/

    auto can21 = new TCanvas("can21", "can", 900, 600);
    hist_sig_StdSmear->SetFillStyle(3004);
    hist_sig_StdSmear->SetLineWidth(2);
    hist_sig_StdSmear->SetLineColor(kAzure + 3);
    hist_sig_StdSmear->SetFillColor(kAzure + 4);
    hist_sig_StdSmear->Draw();
    can21->SaveAs(Form("neutron%s_%.1f_NC_CC_Std_EnergyRatio.root", outName.c_str(), energyHitCut));
// }}}

    /**************************************************************************************************************/

    auto hist_resolutionPrimary_ArmTime = new TH2F("hist_resolutionPrimary_ArmTime", "Purity; Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);

    // Calculate the resolution as function of lever arm and time
    for (int iArm = 0; iArm < 20; iArm++){
	for (int iTime = 0; iTime < 25; iTime++){
	    // Calculate the mean and std for each grid
	    Statistics_t temp, tempPrimary;

	    // loop around (Hit_t) Sig_Neutron_Hit to find the value for each grid
	    for (int ivec = 0; ivec < sig_Neutron_Lvr[iArm].size(); ivec++){
		Hit_t tempHit = sig_Neutron_Lvr[iArm].at(ivec);

		if (iTime == (int) tempHit.timeWindow){
		    float sig_SmearTrue = (tempHit.smearRec - tempHit.trueRec) / tempHit.trueRec;

		    if (abs(sig_SmearTrue) < 1){
			temp.Push(sig_SmearTrue);

			if (tempHit.neutronParentID == 0 || tempHit.neutronParentID == -1)
			    tempPrimary.Push(sig_SmearTrue);
		    }
		}
	    }

	    // Fill histogram for a specific grid
	    hist_resolution_ArmTime  ->SetBinContent(iArm + 1, iTime + 1, temp.Std());
	    hist_resolutionPrimary_ArmTime->SetBinContent(iArm + 1, iTime + 1, tempPrimary.Std());
	}
    }
    
    /**************************************************************************************************************/

// Fill histogram {{{
    for (auto hit : countCubeHit) hist_hitMultiplicity->Fill(hit.second);
   
    auto can26 = new TCanvas("can26", "can", 900, 600);
    hist_hitMultiplicity->Draw();
    can26->SaveAs(Form("neutron%s_%.1f_NC_CC_hitMultiplicity.root", outName.c_str(), energyHitCut));

    auto can27 = new TCanvas("can27", "can", 900, 600);
    hist_cubeFire->Draw();
    can27->SaveAs(Form("neutron%s_%.1f_NC_CC_CubeFiredPerevent.root", outName.c_str(), energyHitCut));
// }}}

    // Calculate the purity
    // Loop around (Hit_t) sig_Neutron_Lvr and evaluate the purity for each grid
    for (int iLvr = 0; iLvr < 20; iLvr++){
	for (int iTime = 0; iTime < 25; iTime++){

	    int nbBkg = 0, nbSig = 0;

	    // Count signal
	    for (auto vec : sig_Neutron_Lvr[iLvr]){
		if (iTime == (int) vec.timeWindow) nbSig++;
	    }

	    // Count background
	    for (auto vec : bkg_Neutron_Lvr[iLvr]){
		if (iTime == (int) vec.timeWindow) nbBkg++;
	    }

	    if (nbSig == 0 && nbBkg == 0) continue;

	    float purity = (float)nbSig/(nbSig + nbBkg);
	    
	    // Fill histogram
	    hist_purity_ArmTime->SetBinContent(iLvr + 1, iTime + 1, purity);
	}
    }

    auto hist_purity_Primary_PionCut    = new TH2F("hist_purity_Primary_PionCut",    "Purity; Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
    auto hist_purity_Primary_ProtonCut  = new TH2F("hist_purity_Primary_ProtonCut",  "Purity; Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
    auto hist_purity_Primary_noCut      = new TH2F("hist_purity_Primary_noCut",      "Purity; Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
    auto hist_purity_Primary_withAllCut = new TH2F("hist_purity_Primary_withAllCut", "Purity; Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);

    auto hist_Parent_Pdg_PionCut    = new TH1F("hist_Parent_Pdg_PionCut",    "Secondary particle Parent PDG", 3000, 0, 3000);
    auto hist_Parent_Pdg_ProtonCut  = new TH1F("hist_Parent_Pdg_ProtonCut",  "Secondary particle Parent PDG", 3000, 0, 3000);
    auto hist_Parent_Pdg_noCut      = new TH1F("hist_Parent_Pdg_noCut",      "Secondary particle Parent PDG", 3000, 0, 3000);
    auto hist_Parent_Pdg_withAllCut = new TH1F("hist_Parent_Pdg_withAllCut", "Secondary particle Parent PDG", 3000, 0, 3000);    

    // Primary purity with pion and proton cut
    for (int iLvr = 0; iLvr < 20; iLvr++){
	for (int iTime = 0; iTime < 25; iTime++){

	    int nbSigPrimary_PionCut    = 0, nbSigSecondary_PionCut    = 0;
	    int nbSigPrimary_ProtonCut  = 0, nbSigSecondary_ProtonCut  = 0;
	    int nbSigPrimary_noCut      = 0, nbSigSecondary_noCut      = 0;
	    int nbSigPrimary_withAllCut = 0, nbSigSecondary_withAllCut = 0;

	    for (auto vec : sig_Neutron_Lvr[iLvr]){
		if (iTime == (int) vec.timeWindow){
		    // primary neutron
		    if (vec.neutronParentID <= 0){
			if (vec.isTherePion50 == false) nbSigPrimary_PionCut++;
			if (vec.isThereProton300 == false) nbSigPrimary_ProtonCut++;
			if (vec.isTherePion50 == false && vec.isThereProton300 == false) nbSigPrimary_withAllCut++;
			nbSigPrimary_noCut++;
		    }

		    // secondary neutron
		    if (vec.neutronParentID > 0){
			if (vec.isTherePion50 == false){
			    nbSigSecondary_PionCut++;
			    hist_Parent_Pdg_PionCut->Fill(vec.neutronParentPdg);
			}
			if (vec.isThereProton300 == false){
			    nbSigSecondary_ProtonCut++;
			    hist_Parent_Pdg_ProtonCut->Fill(vec.neutronParentPdg);
			}
			if (vec.isTherePion50 == false && vec.isThereProton300 == false){
			    nbSigSecondary_withAllCut++;
			    hist_Parent_Pdg_withAllCut->Fill(vec.neutronParentPdg);
			}
			nbSigSecondary_noCut++;
			hist_Parent_Pdg_noCut->Fill(vec.neutronParentPdg);
		    }
		}
	    }

	    // no cut
	    if (nbSigPrimary_noCut != 0 && nbSigPrimary_noCut != nbSigSecondary_noCut){
		float purity = (float) (nbSigPrimary_noCut - nbSigSecondary_noCut) / nbSigPrimary_noCut;	
		if (abs(purity) <= 1) 
		    hist_purity_Primary_noCut->SetBinContent(iLvr + 1, iTime + 1, purity);
	    }

	    // Pion cut
	    if (nbSigPrimary_PionCut != 0 && nbSigPrimary_PionCut != nbSigSecondary_PionCut){
		float purity = (float) (nbSigPrimary_PionCut - nbSigSecondary_PionCut) / nbSigPrimary_PionCut;	
		if (abs(purity) <= 1) {
		    hist_purity_Primary_PionCut->SetBinContent(iLvr + 1, iTime + 1, purity);

		}
	    }

	    // Proton cut
	    if (nbSigPrimary_ProtonCut != 0 && nbSigPrimary_ProtonCut != nbSigSecondary_ProtonCut){
		float purity = (float) (nbSigPrimary_ProtonCut - nbSigSecondary_ProtonCut) / nbSigPrimary_ProtonCut;	
		if (abs(purity) <= 1){
		    hist_purity_Primary_ProtonCut->SetBinContent(iLvr + 1, iTime + 1, purity);
		}
	    }

	    // with all cut
	    if (nbSigPrimary_withAllCut != 0 && nbSigPrimary_withAllCut != nbSigSecondary_withAllCut){
		float purity = (float) (nbSigPrimary_withAllCut - nbSigSecondary_withAllCut) / nbSigPrimary_withAllCut;	
		if (abs(purity) <= 1){
		    hist_purity_Primary_withAllCut->SetBinContent(iLvr + 1, iTime + 1, purity);
		}
	    }
	}
    }

    // Fill histogram {{{
    auto can22 = new TCanvas("can22", "can", 900, 600);

    hist_sigTime_primary->SetFillColor(kAzure + 3);
    hist_sigLvr_primary->SetFillColor(kAzure + 3);
    hist_sigEdep_primary->SetFillColor(kAzure + 3);

    hist_sigTime_secondary->SetFillColor(kAzure + 6);
    hist_sigLvr_secondary->SetFillColor(kAzure + 6);
    hist_sigEdep_secondary->SetFillColor(kAzure + 6);

    hist_sigTime_secondary->SetLineColor(kAzure - 4);
    hist_sigLvr_secondary->SetLineColor(kAzure - 4);
    hist_sigEdep_secondary->SetLineColor(kAzure - 4);

    hist_bkgTime_primary->SetFillColor(kTeal + 3);
    hist_bkgLvr_primary->SetFillColor(kTeal + 3);
    hist_bkgEdep_primary->SetFillColor(kTeal + 3);

    hist_bkgTime_secondary->SetFillColor(kTeal + 6);
    hist_bkgLvr_secondary->SetFillColor(kTeal + 6);
    hist_bkgEdep_secondary->SetFillColor(kTeal + 6);

    hist_bkgTime_secondary->SetLineColor(kTeal - 4);
    hist_bkgLvr_secondary->SetLineColor(kTeal - 4);
    hist_bkgEdep_secondary->SetLineColor(kTeal - 4);
    //}}}

    auto fffii = new TFile("primaryPurity.root", "RECREATE");
    hist_purity_Primary_PionCut->Write();
    hist_purity_Primary_ProtonCut->Write();
    hist_purity_Primary_noCut->Write();
    hist_purity_Primary_withAllCut->Write();

    hist_Parent_Pdg_PionCut->Write();
    hist_Parent_Pdg_ProtonCut->Write();
    hist_Parent_Pdg_noCut->Write();
    hist_Parent_Pdg_withAllCut->Write();

    hist_resolutionPrimary_ArmTime->Write();

    h1->Write();
    h2->Write();

    hist_sig_eRec->Write();
    hist_sig_eRecPrimary->Write();

    fffii->Close();

    auto fi = new TFile("hist_2D_purity.root", "RECREATE");
    hist_purity_ArmTime->Write();
    hist_resolution_ArmTime->Write();

    hist_sigTime_secondary->Write();
    hist_sigLvr_secondary->Write();
    hist_sigEdep_secondary->Write();

    hist_bkgTime_primary->Write();
    hist_bkgLvr_primary->Write();
    hist_bkgEdep_primary->Write();

    hist_bkgTime_secondary->Write();
    hist_bkgLvr_secondary->Write();
    hist_bkgEdep_secondary->Write();

    hist_bkg_time->Write();
    hist_bkg_edep->Write();
    hist_bkg_Lvr->Write();
    
    hist_bkg_vtx_yx->Write();
    hist_bkg_vtx_zy->Write();
    hist_bkg_vtx_zx->Write();

    fi->Close();

    auto _vtx = new TFile("hist_bkg_vtx.root", "RECREATE");
    hist_bkg_vtx_neutrino->Write();
    _vtx->Close();
}
