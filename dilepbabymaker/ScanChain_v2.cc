#include "ScanChain_v2.h"
#include "LeptonScaleFactors.h"
#include "fakerate_v3.h"
#include "puWeight.h"
#include "rooutil/rooutil.cc"

using namespace std;

//##############################################################################################################
babyMaker_v2::babyMaker_v2() : coreBtagSFFastSim(true), babymode(kWWWBaby), eventlist("eventlist.txt")
{
}

//##############################################################################################################
babyMaker_v2::~babyMaker_v2()
{
}

//##############################################################################################################
int babyMaker_v2::ProcessCMS4(TString filepaths, int max_events, int idx, bool verbose)
{
    // Create TChain to process
    TChain* chain = RooUtil::FileUtil::createTChain("Events", filepaths);

    // Initialize Looper
    looper.init(chain, &cms3, max_events);

    // Initializer job index
    job_index = idx;

    // Run the loop
    ScanChain_v2(verbose);

    // Exit
    return 0;
}

//##############################################################################################################
std::once_flag flag_init; // a flag for std::call_once
void babyMaker_v2::ScanChain_v2(bool verbose)
{
    while (looper.nextEvent())
    {
        try
        {
            if (verbose)
                cout << "[verbose] Processed " << looper.getNEventsProcessed() << " out of " << looper.getTChain()->GetEntries() << endl;

            // Some of the sample specific things need to be set prior to processing the very first event
            // Those are called inside Init, they have explicit "call_once" feature
            // The reason this is inside the loop is because it first need to access the file name before
            std::call_once(flag_init, &babyMaker_v2::Init, this);

            // Now process the baby ntuples
            Process();

            // For Gen-level weights the total weights need to be saved so must be called without skipping events (i.e. before isPass() or making sure i fill everything)
            FillGenWeights();

        }
        catch (const std::ios_base::failure& e)
        {
            if (isData())
            {
                std::cout << "Found bad event in data" << std::endl;
                FATALERROR(__FUNCTION__);
                exit(2);
            }
            tx->clear(); // clear the TTree of any residual stuff
            if (!looper.handleBadEvent())
                break;
        }
    }

    looper.printStatus();

    SaveOutput();
}

//##############################################################################################################
void babyMaker_v2::Init()
{
    // Set year via "GlobalConfig gconf"
    SetYear();

    // Provide which file it is and whether it is fast sim or not to JEC to determine which file to load
    coreJec.setJECFor(looper.getCurrentFileName(), isFastSim());

    // Provide which year it is to determine which JER config files to load
    coreJer.setJERFor(gconf.year);

    // Output root file
    CreateOutput();

    // Signal or BSM samples need some specifics things to be set prior to starting the processing
    AddOutput();

    // Set lepton ID configuration via "GlobalConfig gconf"
    SetLeptonID();
}


//##############################################################################################################
void babyMaker_v2::CreateOutput()
{
    ofile = new TFile(Form("output_%d.root", job_index), "recreate");
    t = new TTree("t", "All events");
    tx = new RooUtil::TTreeX(t);

    tx->createBranch<Int_t>("run");
    tx->createBranch<Int_t>("lumi");
    tx->createBranch<unsigned long long>("evt");
    tx->createBranch<int>("isData");
    tx->createBranch<float>("evt_scale1fb");
    tx->createBranch<float>("xsec_br");
    tx->createBranch<int>("evt_passgoodrunlist");
    tx->createBranch<TString>("CMS4path");
    tx->createBranch<int>("CMS4index");

    // Below are things related to event weights which must be done before skipping events
    tx->createBranch<float>("weight_fr_r1_f1");
    tx->createBranch<float>("weight_fr_r1_f2");
    tx->createBranch<float>("weight_fr_r1_f0p5");
    tx->createBranch<float>("weight_fr_r2_f1");
    tx->createBranch<float>("weight_fr_r2_f2");
    tx->createBranch<float>("weight_fr_r2_f0p5");
    tx->createBranch<float>("weight_fr_r0p5_f1");
    tx->createBranch<float>("weight_fr_r0p5_f2");
    tx->createBranch<float>("weight_fr_r0p5_f0p5");
    tx->createBranch<float>("weight_pdf_up");
    tx->createBranch<float>("weight_pdf_down");
    tx->createBranch<float>("weight_alphas_down");
    tx->createBranch<float>("weight_alphas_up");
    tx->createBranch<float>("weight_isr");
    tx->createBranch<float>("weight_isr_up");
    tx->createBranch<float>("weight_isr_down");

    h_neventsinfile = new TH1F("h_neventsinfile", "", 15, 0, 15);
    h_neventsinfile->SetBinContent(1, looper.getTChain()->GetEntries()); // this is the bin with value = 0
}

//##############################################################################################################
void babyMaker_v2::AddOutput()
{

    // The main branches are created
    AddBabyOutput();


    // Some BSM or SM signals related output (e.g. truth level study)
    AddTruthStudyOutput();
}

//##############################################################################################################
void babyMaker_v2::AddBabyOutput()
{
    // Based on the baby mode now preselect and if it doesn't pass return
    switch (babymode)
    {
        case kWWWBaby:   AddWWWBabyOutput(); return; break;
        case kFRBaby:    AddWWWBabyOutput(); return; break;
        case kOSBaby:    AddWWWBabyOutput(); return; break;
        case kTnPBaby:   AddTnPBabyOutput(); return; break;
        case kAllBaby:   AddWWWBabyOutput(); return; break;
        case kPOGBaby:   AddWWWBabyOutput(); return; break;
        case kLooseBaby: AddWWWBabyOutput(); return; break;
        default: return;
    }
}


//##############################################################################################################
void babyMaker_v2::AddTruthStudyOutput()
{
    if (isDoublyChargedHiggs())
    {
        AddDoublyChargedHiggsOutput();
    }

    if (isWprime() || isVH())
    {
        AddWprimeOutput();
    }

    if (isWHSUSY())
    {
        AddWHsusyOutput();
    }

    if (isWWW() || isVH() || isWHSUSY())
    {
        AddWWWSignalOutput();
    }
}

//##############################################################################################################
void babyMaker_v2::AddWWWBabyOutput()
{
    tx->createBranch<int>("HLT_DoubleMu");
    tx->createBranch<int>("HLT_DoubleEl");
    tx->createBranch<int>("HLT_DoubleEl_DZ");
    tx->createBranch<int>("HLT_DoubleEl_DZ_2");
    tx->createBranch<int>("HLT_MuEG");
    tx->createBranch<int>("HLT_SingleEl8");
    tx->createBranch<int>("HLT_SingleEl17");
    tx->createBranch<int>("HLT_SingleIsoEl8");
    tx->createBranch<int>("HLT_SingleIsoEl17");
    tx->createBranch<int>("HLT_SingleIsoEl23");
    tx->createBranch<int>("HLT_SingleIsoMu8");
    tx->createBranch<int>("HLT_SingleIsoMu17");
    tx->createBranch<int>("HLT_PFMET140_PFMHT140_IDTight");

    tx->createBranch<int>("mc_HLT_DoubleMu");
    tx->createBranch<int>("mc_HLT_DoubleEl");
    tx->createBranch<int>("mc_HLT_DoubleEl_DZ");
    tx->createBranch<int>("mc_HLT_DoubleEl_DZ_2");
    tx->createBranch<int>("mc_HLT_MuEG");
    tx->createBranch<int>("mc_HLT_SingleEl8");
    tx->createBranch<int>("mc_HLT_SingleEl17");
    tx->createBranch<int>("mc_HLT_SingleIsoEl8");
    tx->createBranch<int>("mc_HLT_SingleIsoEl17");
    tx->createBranch<int>("mc_HLT_SingleIsoEl23");
    tx->createBranch<int>("mc_HLT_SingleIsoMu8");
    tx->createBranch<int>("mc_HLT_SingleIsoMu17");
    tx->createBranch<int>("mc_HLT_PFMET140_PFMHT140_IDTight");

    tx->createBranch<int>("pass_duplicate_ee_em_mm");
    tx->createBranch<int>("pass_duplicate_mm_em_ee");

    tx->createBranch<int>("is2016");
    tx->createBranch<int>("is2017");

    // Until something weird about 2016 is resolved
    tx->createBranch<int>("HLT_MuEG_2016");
    tx->createBranch<int>("mc_HLT_MuEG_2016");
    tx->createBranch<int>("pass_duplicate_ee_em2016_mm");
    tx->createBranch<int>("pass_duplicate_mm_em2016_ee");

    tx->createBranch<int>("passTrigger");

    tx->createBranch<vector<LorentzVector>>("lep_p4");
    tx->createBranch<vector<float>>("lep_pt");
    tx->createBranch<vector<float>>("lep_eta");
    tx->createBranch<vector<float>>("lep_phi");
    tx->createBranch<vector<float>>("lep_coneCorrPt");
    tx->createBranch<vector<float>>("lep_ip3d");
    tx->createBranch<vector<float>>("lep_ip3derr");
    tx->createBranch<vector<int>>("lep_isTriggerSafe_v1");
    tx->createBranch<vector<int>>("lep_lostHits");
    tx->createBranch<vector<int>>("lep_convVeto");
    tx->createBranch<vector<int>>("lep_motherIdSS");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_3l_fo");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_3l_tight");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_fo");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_tight");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_veto");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_fo_noiso");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_tight_noiso");
    tx->createBranch<vector<int>>("lep_pass_VVV_cutbased_veto_noiso");
    tx->createBranch<vector<int>>("lep_pass_POG_veto");
    tx->createBranch<vector<int>>("lep_pass_POG_loose");
    tx->createBranch<vector<int>>("lep_pass_POG_medium");
    tx->createBranch<vector<int>>("lep_pass_POG_tight");
    tx->createBranch<vector<int>>("lep_pdgId");
    tx->createBranch<vector<float>>("lep_dxy");
    tx->createBranch<vector<float>>("lep_dz");
    tx->createBranch<vector<float>>("lep_pterr");
    tx->createBranch<vector<float>>("lep_relIso04DB");
    tx->createBranch<vector<float>>("lep_relIso03EA");
    tx->createBranch<vector<float>>("lep_relIso03EALep");
    tx->createBranch<vector<float>>("lep_relIso03EAv2");
    tx->createBranch<vector<float>>("lep_relIso04EAv2");
    tx->createBranch<vector<float>>("lep_relIso03EAv2Lep");
    tx->createBranch<vector<int>>("lep_tightCharge");
    tx->createBranch<vector<float>>("lep_trk_pt");
    tx->createBranch<vector<int>>("lep_charge");
    tx->createBranch<vector<float>>("lep_etaSC");
    tx->createBranch<vector<float>>("lep_MVA");
    tx->createBranch<vector<int>>("lep_isMediumPOG");
    tx->createBranch<vector<int>>("lep_isTightPOG");
    tx->createBranch<vector<int>>("lep_isFromW");
    tx->createBranch<vector<int>>("lep_isFromZ");
    tx->createBranch<vector<int>>("lep_isFromB");
    tx->createBranch<vector<int>>("lep_isFromC");
    tx->createBranch<vector<int>>("lep_isFromL");
    tx->createBranch<vector<int>>("lep_isFromLF");
    tx->createBranch<vector<int>>("lep_genPart_index");
    tx->createBranch<vector<float>>("lep_r9");
    tx->createBranch<vector<int>>("lep_nlayers");

    tx->createBranch<float>("el_pt");
    tx->createBranch<float>("el_eta");
    tx->createBranch<float>("el_phi");
    tx->createBranch<float>("el_relIso03EA");
    tx->createBranch<float>("el_relIso03EALep");
    tx->createBranch<float>("el_ip3d");
    tx->createBranch<float>("mu_pt");
    tx->createBranch<float>("mu_eta");
    tx->createBranch<float>("mu_phi");
    tx->createBranch<float>("mu_relIso04DB");
    tx->createBranch<float>("mu_relIso03EA");
    tx->createBranch<float>("mu_relIso03EALep");
    tx->createBranch<float>("mu_ip3d");

    tx->createBranch<float>("lbnt_pt");
    tx->createBranch<float>("lbnt_coneCorrPt");
    tx->createBranch<float>("lbnt_abseta");
    tx->createBranch<float>("lbnt_pdgId");
    tx->createBranch<float>("lbnt_el_pt");
    tx->createBranch<float>("lbnt_el_coneCorrPt");
    tx->createBranch<float>("lbnt_el_abseta");
    tx->createBranch<float>("lbnt_mu_pt");
    tx->createBranch<float>("lbnt_mu_coneCorrPt");
    tx->createBranch<float>("lbnt_mu_abseta");

    tx->createBranch<vector<LorentzVector>>("jets_p4");
    tx->createBranch<vector<LorentzVector>>("jets_up_p4");
    tx->createBranch<vector<LorentzVector>>("jets_dn_p4");
    tx->createBranch<vector<float>>("jets_csv");
    tx->createBranch<vector<float>>("jets_up_csv");
    tx->createBranch<vector<float>>("jets_dn_csv");
    tx->createBranch<vector<float>>("jets_jer_csv");
    tx->createBranch<vector<float>>("jets_jerup_csv");
    tx->createBranch<vector<float>>("jets_jerdn_csv");
    tx->createBranch<vector<LorentzVector>>("jets_jer_p4");
    tx->createBranch<vector<LorentzVector>>("jets_jerup_p4");
    tx->createBranch<vector<LorentzVector>>("jets_jerdn_p4");

    tx->createBranch<vector<LorentzVector>>("jets30_p4");
    tx->createBranch<vector<LorentzVector>>("jets30_up_p4");
    tx->createBranch<vector<LorentzVector>>("jets30_dn_p4");
    tx->createBranch<vector<LorentzVector>>("jets30_jer_p4");
    tx->createBranch<vector<LorentzVector>>("jets30_jerup_p4");
    tx->createBranch<vector<LorentzVector>>("jets30_jerdn_p4");

    tx->createBranch<vector<LorentzVector>>("ak8jets_p4");
    tx->createBranch<vector<float>>("ak8jets_softdropMass");
    tx->createBranch<vector<float>>("ak8jets_prunedMass");
    tx->createBranch<vector<float>>("ak8jets_trimmedMass");
    tx->createBranch<vector<float>>("ak8jets_mass");
    tx->createBranch<vector<float>>("ak8jets_nJettinessTau1");
    tx->createBranch<vector<float>>("ak8jets_nJettinessTau2");
    tx->createBranch<vector<float>>("ak8jets_softdropPuppiSubjet1");
    tx->createBranch<vector<float>>("ak8jets_softdropPuppiSubjet2");
    tx->createBranch<vector<float>>("ak8jets_puppi_softdropMass");
    tx->createBranch<vector<float>>("ak8jets_puppi_nJettinessTau1");
    tx->createBranch<vector<float>>("ak8jets_puppi_nJettinessTau2");
    tx->createBranch<vector<float>>("ak8jets_puppi_eta");
    tx->createBranch<vector<float>>("ak8jets_puppi_phi");
    tx->createBranch<vector<float>>("ak8jets_puppi_pt");
    tx->createBranch<vector<float>>("ak8jets_puppi_mass");
    tx->createBranch<vector<float>>("ak8jets_deep_rawdisc_h4q");
    tx->createBranch<vector<float>>("ak8jets_deep_rawdisc_hbb");
    tx->createBranch<vector<float>>("ak8jets_deep_rawdisc_top");
    tx->createBranch<vector<float>>("ak8jets_deep_rawdisc_qcd");
    tx->createBranch<vector<float>>("ak8jets_deep_rawdisc_w");
    tx->createBranch<vector<float>>("ak8jets_deep_rawdisc_z");
    tx->createBranch<vector<float>>("ak8jets_deep_rawdisc_zbb");

    tx->createBranch<float>("met_pt");
    tx->createBranch<float>("met_phi");
    tx->createBranch<float>("met_up_pt");
    tx->createBranch<float>("met_up_phi");
    tx->createBranch<float>("met_dn_pt");
    tx->createBranch<float>("met_dn_phi");
    tx->createBranch<float>("met_gen_pt");
    tx->createBranch<float>("met_gen_phi");
    tx->createBranch<float>("met_jer_pt");
    tx->createBranch<float>("met_jerup_pt");
    tx->createBranch<float>("met_jerdn_pt");
    tx->createBranch<float>("met_jer_phi");
    tx->createBranch<float>("met_jerup_phi");
    tx->createBranch<float>("met_jerdn_phi");

    tx->createBranch<int>("firstgoodvertex");
    tx->createBranch<int>("nTrueInt");
    tx->createBranch<int>("nVert");

    tx->createBranch<int>("nisoTrack_mt2_cleaned_VVV_cutbased_veto");

    tx->createBranch<float>("weight_btagsf");
    tx->createBranch<float>("weight_btagsf_heavy_DN");
    tx->createBranch<float>("weight_btagsf_heavy_UP");
    tx->createBranch<float>("weight_btagsf_light_DN");
    tx->createBranch<float>("weight_btagsf_light_UP");

    tx->createBranch<float>("gen_ht");

    tx->createBranch<vector<LorentzVector>>("genPart_p4");
    tx->createBranch<vector<int>>("genPart_motherId");
    tx->createBranch<vector<int>>("genPart_pdgId");
    tx->createBranch<vector<int>>("genPart_charge");
    tx->createBranch<vector<int>>("genPart_status");
    tx->createBranch<int>("ngenLep");
    tx->createBranch<int>("ngenLepFromTau");
    tx->createBranch<int>("ngenLepFromBoson");

    tx->createBranch<int>("Flag_AllEventFilters");
    tx->createBranch<int>("Flag_EcalDeadCellTriggerPrimitiveFilter");
    tx->createBranch<int>("Flag_HBHEIsoNoiseFilter");
    tx->createBranch<int>("Flag_HBHENoiseFilter");
    tx->createBranch<int>("Flag_badChargedCandidateFilter");
    tx->createBranch<int>("Flag_badMuonFilter");
    tx->createBranch<int>("Flag_badMuonFilterv2");
    tx->createBranch<int>("Flag_badChargedCandidateFilterv2");
    tx->createBranch<int>("Flag_eeBadScFilter");
    tx->createBranch<int>("Flag_ecalBadCalibFilter");
    tx->createBranch<int>("Flag_globalTightHalo2016");
    tx->createBranch<int>("Flag_goodVertices");
    tx->createBranch<int>("Flag_ecalLaserCorrFilter");
    tx->createBranch<int>("Flag_hcalLaserEventFilter");
    tx->createBranch<int>("Flag_trackingFailureFilter");
    tx->createBranch<int>("Flag_CSCTightHaloFilter");
    tx->createBranch<int>("Flag_CSCTightHalo2015Filter");
    tx->createBranch<int>("Flag_badMuons");
    tx->createBranch<int>("Flag_duplicateMuons");
    tx->createBranch<int>("Flag_noBadMuons");

    tx->createBranch<int>("fastsimfilt");

    // Summary variables
    tx->createBranch<int>("nVlep");
    tx->createBranch<int>("nTlep");
    tx->createBranch<int>("nTlepSS");
    tx->createBranch<int>("nLlep");
    tx->createBranch<int>("nLlep3L");
    tx->createBranch<int>("nTlep3L");

    tx->createBranch<int>("nSFOS");
    tx->createBranch<int>("nSFOSinZ");

    tx->createBranch<int>("nj");
    tx->createBranch<int>("nj_up");
    tx->createBranch<int>("nj_dn");
    tx->createBranch<int>("nj_jer");
    tx->createBranch<int>("nj_jerup");
    tx->createBranch<int>("nj_jerdn");

    tx->createBranch<int>("nj30");
    tx->createBranch<int>("nj30_up");
    tx->createBranch<int>("nj30_dn");
    tx->createBranch<int>("nj30_jer");
    tx->createBranch<int>("nj30_jerup");
    tx->createBranch<int>("nj30_jerdn");

    tx->createBranch<int>("nb");
    tx->createBranch<int>("nb_up");
    tx->createBranch<int>("nb_dn");
    tx->createBranch<int>("nb_jer");
    tx->createBranch<int>("nb_jerup");
    tx->createBranch<int>("nb_jerdn");

    tx->createBranch<float>("Ml0j0");
    tx->createBranch<float>("Ml0j0_up");
    tx->createBranch<float>("Ml0j0_dn");
    tx->createBranch<float>("Ml0j0_jer");
    tx->createBranch<float>("Ml0j0_jerup");
    tx->createBranch<float>("Ml0j0_jerdn");

    tx->createBranch<float>("Ml0j1");
    tx->createBranch<float>("Ml0j1_up");
    tx->createBranch<float>("Ml0j1_dn");
    tx->createBranch<float>("Ml0j1_jer");
    tx->createBranch<float>("Ml0j1_jerup");
    tx->createBranch<float>("Ml0j1_jerdn");

    tx->createBranch<float>("Ml1j0");
    tx->createBranch<float>("Ml1j0_up");
    tx->createBranch<float>("Ml1j0_dn");
    tx->createBranch<float>("Ml1j0_jer");
    tx->createBranch<float>("Ml1j0_jerup");
    tx->createBranch<float>("Ml1j0_jerdn");

    tx->createBranch<float>("Ml1j1");
    tx->createBranch<float>("Ml1j1_up");
    tx->createBranch<float>("Ml1j1_dn");
    tx->createBranch<float>("Ml1j1_jer");
    tx->createBranch<float>("Ml1j1_jerup");
    tx->createBranch<float>("Ml1j1_jerdn");

    tx->createBranch<float>("MinMlj");
    tx->createBranch<float>("MinMlj_up");
    tx->createBranch<float>("MinMlj_dn");
    tx->createBranch<float>("MinMlj_jer");
    tx->createBranch<float>("MinMlj_jerup");
    tx->createBranch<float>("MinMlj_jerdn");

    tx->createBranch<float>("SumMinMlj01");
    tx->createBranch<float>("SumMinMlj01_up");
    tx->createBranch<float>("SumMinMlj01_dn");
    tx->createBranch<float>("SumMinMlj01_jer");
    tx->createBranch<float>("SumMinMlj01_jerup");
    tx->createBranch<float>("SumMinMlj01_jerdn");

    tx->createBranch<float>("MaxMlj");
    tx->createBranch<float>("MaxMlj_up");
    tx->createBranch<float>("MaxMlj_dn");
    tx->createBranch<float>("MaxMlj_jer");
    tx->createBranch<float>("MaxMlj_jerup");
    tx->createBranch<float>("MaxMlj_jerdn");

    tx->createBranch<float>("SumMlj");
    tx->createBranch<float>("SumMlj_up");
    tx->createBranch<float>("SumMlj_dn");
    tx->createBranch<float>("SumMlj_jer");
    tx->createBranch<float>("SumMlj_jerup");
    tx->createBranch<float>("SumMlj_jerdn");

    tx->createBranch<float>("Ml0jj");
    tx->createBranch<float>("Ml0jj_up");
    tx->createBranch<float>("Ml0jj_dn");
    tx->createBranch<float>("Ml0jj_jer");
    tx->createBranch<float>("Ml0jj_jerup");
    tx->createBranch<float>("Ml0jj_jerdn");

    tx->createBranch<float>("Ml1jj");
    tx->createBranch<float>("Ml1jj_up");
    tx->createBranch<float>("Ml1jj_dn");
    tx->createBranch<float>("Ml1jj_jer");
    tx->createBranch<float>("Ml1jj_jerup");
    tx->createBranch<float>("Ml1jj_jerdn");

    tx->createBranch<float>("MinMljj");
    tx->createBranch<float>("MinMljj_up");
    tx->createBranch<float>("MinMljj_dn");
    tx->createBranch<float>("MinMljj_jer");
    tx->createBranch<float>("MinMljj_jerup");
    tx->createBranch<float>("MinMljj_jerdn");

    tx->createBranch<float>("MaxMljj");
    tx->createBranch<float>("MaxMljj_up");
    tx->createBranch<float>("MaxMljj_dn");
    tx->createBranch<float>("MaxMljj_jer");
    tx->createBranch<float>("MaxMljj_jerup");
    tx->createBranch<float>("MaxMljj_jerdn");

    tx->createBranch<float>("SumMljj");
    tx->createBranch<float>("SumMljj_up");
    tx->createBranch<float>("SumMljj_dn");
    tx->createBranch<float>("SumMljj_jer");
    tx->createBranch<float>("SumMljj_jerup");
    tx->createBranch<float>("SumMljj_jerdn");

    tx->createBranch<float>("Mjj");
    tx->createBranch<float>("Mjj_up");
    tx->createBranch<float>("Mjj_dn");
    tx->createBranch<float>("Mjj_jer");
    tx->createBranch<float>("Mjj_jerup");
    tx->createBranch<float>("Mjj_jerdn");

    tx->createBranch<float>("DRjj");
    tx->createBranch<float>("DRjj_up");
    tx->createBranch<float>("DRjj_dn");
    tx->createBranch<float>("DRjj_jer");
    tx->createBranch<float>("DRjj_jerup");
    tx->createBranch<float>("DRjj_jerdn");

    tx->createBranch<LV>("jet0_wtag_p4");
    tx->createBranch<LV>("jet0_wtag_p4_up");
    tx->createBranch<LV>("jet0_wtag_p4_dn");
    tx->createBranch<LV>("jet0_wtag_p4_jer");
    tx->createBranch<LV>("jet0_wtag_p4_jerup");
    tx->createBranch<LV>("jet0_wtag_p4_jerdn");

    tx->createBranch<LV>("jet1_wtag_p4");
    tx->createBranch<LV>("jet1_wtag_p4_up");
    tx->createBranch<LV>("jet1_wtag_p4_dn");
    tx->createBranch<LV>("jet1_wtag_p4_jer");
    tx->createBranch<LV>("jet1_wtag_p4_jerup");
    tx->createBranch<LV>("jet1_wtag_p4_jerdn");

    tx->createBranch<float>("MjjDR1");
    tx->createBranch<float>("MjjDR1_up");
    tx->createBranch<float>("MjjDR1_dn");
    tx->createBranch<float>("MjjDR1_jer");
    tx->createBranch<float>("MjjDR1_jerup");
    tx->createBranch<float>("MjjDR1_jerdn");

    tx->createBranch<float>("DRjjDR1");
    tx->createBranch<float>("DRjjDR1_up");
    tx->createBranch<float>("DRjjDR1_dn");
    tx->createBranch<float>("DRjjDR1_jer");
    tx->createBranch<float>("DRjjDR1_jerup");
    tx->createBranch<float>("DRjjDR1_jerdn");

    tx->createBranch<LV>("jet0_wtag_p4_DR1");
    tx->createBranch<LV>("jet0_wtag_p4_DR1_up");
    tx->createBranch<LV>("jet0_wtag_p4_DR1_dn");
    tx->createBranch<LV>("jet0_wtag_p4_DR1_jer");
    tx->createBranch<LV>("jet0_wtag_p4_DR1_jerup");
    tx->createBranch<LV>("jet0_wtag_p4_DR1_jerdn");

    tx->createBranch<LV>("jet1_wtag_p4_DR1");
    tx->createBranch<LV>("jet1_wtag_p4_DR1_up");
    tx->createBranch<LV>("jet1_wtag_p4_DR1_dn");
    tx->createBranch<LV>("jet1_wtag_p4_DR1_jer");
    tx->createBranch<LV>("jet1_wtag_p4_DR1_jerup");
    tx->createBranch<LV>("jet1_wtag_p4_DR1_jerdn");

    tx->createBranch<float>("MjjVBF");
    tx->createBranch<float>("MjjVBF_up");
    tx->createBranch<float>("MjjVBF_dn");
    tx->createBranch<float>("MjjVBF_jer");
    tx->createBranch<float>("MjjVBF_jerup");
    tx->createBranch<float>("MjjVBF_jerdn");

    tx->createBranch<float>("DetajjVBF");
    tx->createBranch<float>("DetajjVBF_up");
    tx->createBranch<float>("DetajjVBF_dn");
    tx->createBranch<float>("DetajjVBF_jer");
    tx->createBranch<float>("DetajjVBF_jerup");
    tx->createBranch<float>("DetajjVBF_jerdn");

    tx->createBranch<float>("MjjL");
    tx->createBranch<float>("MjjL_up");
    tx->createBranch<float>("MjjL_dn");
    tx->createBranch<float>("MjjL_jer");
    tx->createBranch<float>("MjjL_jerup");
    tx->createBranch<float>("MjjL_jerdn");

    tx->createBranch<float>("DetajjL");
    tx->createBranch<float>("DetajjL_up");
    tx->createBranch<float>("DetajjL_dn");
    tx->createBranch<float>("DetajjL_jer");
    tx->createBranch<float>("DetajjL_jerup");
    tx->createBranch<float>("DetajjL_jerdn");

    tx->createBranch<float>("MllSS");
    tx->createBranch<float>("MeeSS");
    tx->createBranch<float>("Mll3L");
    tx->createBranch<float>("Mee3L");
    tx->createBranch<float>("Mll3L1");
    tx->createBranch<float>("M3l");
    tx->createBranch<float>("Pt3l");
    tx->createBranch<float>("M01");
    tx->createBranch<float>("M02");
    tx->createBranch<float>("M12");
    tx->createBranch<int>("isSFOS01");
    tx->createBranch<int>("isSFOS02");
    tx->createBranch<int>("isSFOS12");

    tx->createBranch<float>("DPhi3lMET");
    tx->createBranch<float>("DPhi3lMET_up");
    tx->createBranch<float>("DPhi3lMET_dn");
    tx->createBranch<float>("DPhi3lMET_jer");
    tx->createBranch<float>("DPhi3lMET_jerup");
    tx->createBranch<float>("DPhi3lMET_jerdn");
    tx->createBranch<float>("DPhi3lMET_gen");
    tx->createBranch<float>("MTmax");
    tx->createBranch<float>("MTmax_up");
    tx->createBranch<float>("MTmax_dn");
    tx->createBranch<float>("MTmax_jer");
    tx->createBranch<float>("MTmax_jerup");
    tx->createBranch<float>("MTmax_jerdn");
    tx->createBranch<float>("MTmax_gen");
    tx->createBranch<float>("MTmin");
    tx->createBranch<float>("MTmin_up");
    tx->createBranch<float>("MTmin_dn");
    tx->createBranch<float>("MTmin_jer");
    tx->createBranch<float>("MTmin_jerup");
    tx->createBranch<float>("MTmin_jerdn");
    tx->createBranch<float>("MTmin_gen");

    tx->createBranch<float>("MT3rd");
    tx->createBranch<float>("MT3rd_up");
    tx->createBranch<float>("MT3rd_dn");
    tx->createBranch<float>("MT3rd_jer");
    tx->createBranch<float>("MT3rd_jerup");
    tx->createBranch<float>("MT3rd_jerdn");
    tx->createBranch<float>("MT3rd_gen");
    tx->createBranch<float>("MTmax3L");
    tx->createBranch<float>("MTmax3L_up");
    tx->createBranch<float>("MTmax3L_dn");
    tx->createBranch<float>("MTmax3L_jer");
    tx->createBranch<float>("MTmax3L_jerup");
    tx->createBranch<float>("MTmax3L_jerdn");
    tx->createBranch<float>("MTmax3L_gen");

    tx->createBranch<int>("passSSee");
    tx->createBranch<int>("passSSem");
    tx->createBranch<int>("passSSmm");

    tx->createBranch<int>("lep_idx0_SS");
    tx->createBranch<int>("lep_idx1_SS");

    tx->createBranch<TString>("bkgtype");
    tx->createBranch<int>("vetophoton");

    tx->createBranch<float>("purewgt");
    tx->createBranch<float>("purewgt_up");
    tx->createBranch<float>("purewgt_dn");

    tx->createBranch<float>("ffwgt");
    tx->createBranch<float>("ffwgt_up");
    tx->createBranch<float>("ffwgt_dn");
    tx->createBranch<float>("ffwgt_el_up");
    tx->createBranch<float>("ffwgt_el_dn");
    tx->createBranch<float>("ffwgt_mu_up");
    tx->createBranch<float>("ffwgt_mu_dn");
    tx->createBranch<float>("ffwgt_closure_up");
    tx->createBranch<float>("ffwgt_closure_dn");
    tx->createBranch<float>("ffwgt_closure_el_up");
    tx->createBranch<float>("ffwgt_closure_el_dn");
    tx->createBranch<float>("ffwgt_closure_mu_up");
    tx->createBranch<float>("ffwgt_closure_mu_dn");
    tx->createBranch<float>("ffwgt_full_up");
    tx->createBranch<float>("ffwgt_full_dn");

    tx->createBranch<float>("ffwgtqcd");
    tx->createBranch<float>("ffwgtqcd_up");
    tx->createBranch<float>("ffwgtqcd_dn");

    tx->createBranch<float>("lepsf");
    tx->createBranch<float>("lepsf_up");
    tx->createBranch<float>("lepsf_dn");

    tx->createBranch<float>("trigeff");
    tx->createBranch<float>("trigeff_up");
    tx->createBranch<float>("trigeff_dn");

    tx->createBranch<float>("trigsf");
    tx->createBranch<float>("trigsf_up");
    tx->createBranch<float>("trigsf_dn");

    tx->createBranch<float>("musmear_sf");

    tx->clear();
}

//##############################################################################################################
void babyMaker_v2::AddPOGBabyOutput()
{
}

//##############################################################################################################
void babyMaker_v2::AddTnPBabyOutput()
{

    tx->createBranch<LV>("dilep_p4");
    tx->createBranch<float>("dilep_mass");
    tx->createBranch<LV>("p4");
    tx->createBranch<bool>("conv_vtx_flag");
    tx->createBranch<bool>("evt_isRealData");
    tx->createBranch<bool>("isPF");

    tx->createBranch<int>("id");

    tx->createBranch<float>("AbsTrkIso");
    tx->createBranch<float>("TrkAn04");
    tx->createBranch<float>("RelIso03EA");
    tx->createBranch<float>("RelIso03EAv2");
    tx->createBranch<float>("RelIso03EAv2wlep");
    tx->createBranch<float>("miniiso");
    tx->createBranch<float>("miniisoDB");

    tx->createBranch<float>("dZ");
    tx->createBranch<float>("dxyPV");
    tx->createBranch<float>("dxyPV_err");
    tx->createBranch<float>("ip3d");
    tx->createBranch<float>("ip3derr");

    tx->createBranch<bool>("passes_HAD_veto_noiso_v3");
    tx->createBranch<bool>("passes_HAD_veto_v3");
    tx->createBranch<bool>("passes_SS_tight_noiso_v3");
    tx->createBranch<bool>("passes_SS_tight_v3");
    tx->createBranch<bool>("passes_VVV_tight");
    tx->createBranch<bool>("passes_VVV_tight_noiso");
    tx->createBranch<bool>("passes_VVV_3l_tight");
    tx->createBranch<bool>("passes_VVV_3l_tight_noiso");
    tx->createBranch<bool>("passes_VVV_fo");
    tx->createBranch<bool>("passes_VVV_fo_noiso");
    tx->createBranch<bool>("passes_VVV_3l_fo");
    tx->createBranch<bool>("passes_VVV_3l_fo_noiso");
    tx->createBranch<bool>("passes_VVV_veto");
    tx->createBranch<bool>("passes_VVV_veto_noiso");
    tx->createBranch<bool>("passes_POG_mediumID");
    tx->createBranch<bool>("passes_POG_looseID");
    tx->createBranch<bool>("passes_POG_MVA_HZZ");
    tx->createBranch<bool>("passes_POG_MVA_wp80");
    tx->createBranch<bool>("passes_POG_MVA_wp90");

    tx->createBranch<LV>("mc_motherp4");
    tx->createBranch<LV>("mc_p4");
    tx->createBranch<int>("mc_motherid");
    tx->createBranch<int>("motherID");

    tx->createBranch<bool>("threeChargeAgree");
    tx->createBranch<float>("dEtaIn");
    tx->createBranch<float>("dPhiIn");
    tx->createBranch<float>("eOverPIn");
    tx->createBranch<float>("ecalEnergy");
    tx->createBranch<float>("etaSC");
    tx->createBranch<float>("etaSCwidth");
    tx->createBranch<float>("phiSCwidth");
    tx->createBranch<float>("hOverE");
    tx->createBranch<float>("sigmaIEtaIEta_full5x5");
    tx->createBranch<float>("mva");
    tx->createBranch<float>("mva_25ns");
    tx->createBranch<float>("mva_2017");

    tx->createBranch<float>("ptErr");
    tx->createBranch<float>("trk_pt");

    tx->createBranch<int>("evt_event");
    tx->createBranch<int>("evt_lumiBlock");
    tx->createBranch<int>("evt_run");
    tx->createBranch<int>("nvtx");

    tx->createBranch<LV>("tag_p4");
    tx->createBranch<float>("tag_RelIso03EA");
    tx->createBranch<float>("tag_ecalEnergy");
    tx->createBranch<float>("tag_mva_25ns");
    tx->createBranch<int>("tag_mc_motherid");
    tx->createBranch<int>("tag_charge");
    tx->createBranch<int>("tag_HLT_Ele27_eta2p1_WPTight_Gsf");
    tx->createBranch<int>("tag_HLT_Ele32_eta2p1_WPTight_Gsf");
    tx->createBranch<int>("tag_HLT_Ele32_WPTight_Gsf");
    tx->createBranch<int>("tag_HLT_Ele35_WPTight_Gsf");
    tx->createBranch<int>("tag_HLT_IsoMu24");
    tx->createBranch<int>("tag_HLT_IsoMu27");
    tx->createBranch<int>("tag_HLT_IsoTkMu24");

    tx->createBranch<int>("charge");
    tx->createBranch<int>("ckf_charge");

    tx->createBranch<int>("exp_innerlayers");
    tx->createBranch<int>("gfit_validSTAHits");
    tx->createBranch<int>("gsf_validHits");

    tx->createBranch<int>("idx");
    tx->createBranch<int>("jet_close_lep_idx");
    tx->createBranch<int>("mc_id");
    tx->createBranch<int>("pid_PFMuon");
    tx->createBranch<int>("sccharge");
    tx->createBranch<int>("validHits");
    tx->createBranch<int>("validPixelHits");

    tx->clear();
}

//##############################################################################################################
void babyMaker_v2::AddDoublyChargedHiggsOutput()
{
    tx->createBranch<int>("iswhchannel");
    tx->createBranch<LV>("w_p4");
    tx->createBranch<LV>("q0_p4");
    tx->createBranch<LV>("q1_p4");
    tx->createBranch<LV>("h_p4");
    tx->createBranch<LV>("w0_p4");
    tx->createBranch<LV>("w1_p4");
    tx->createBranch<LV>("l0_p4");
    tx->createBranch<LV>("l1_p4");
    tx->createBranch<LV>("v0_p4");
    tx->createBranch<LV>("v1_p4");

    tx->clear();
}

//##############################################################################################################
void babyMaker_v2::AddWprimeOutput()
{
    tx->createBranch<int>("higgsdecay");
    tx->createBranch<int>("nlep");
    tx->createBranch<int>("nquark");
    tx->createBranch<int>("wa_id");
    tx->createBranch<LV>("higgs_p4");
    tx->createBranch<vector<LV>>("decay_p4");
    tx->createBranch<vector<int>>("decay_id");
    tx->createBranch<vector<int>>("decay_isstar");
    tx->createBranch<vector<LV>>("lepton_p4");
    tx->createBranch<vector<int>>("lepton_id");
    tx->createBranch<vector<int>>("lepton_isstar");
    tx->createBranch<vector<LV>>("quark_p4");
    tx->createBranch<vector<int>>("quark_id");
    tx->createBranch<vector<int>>("quark_isstar");
    tx->createBranch<float>("lqq_max_dr");
    tx->createBranch<float>("lq0_dr");
    tx->createBranch<float>("lq1_dr");
    tx->createBranch<float>("qq_dr");

    tx->createBranch<vector<LV>>("boosted0_decay_p4");
    tx->createBranch<vector<int>>("boosted0_decay_id");
    tx->createBranch<vector<int>>("boosted0_decay_isstar");
    tx->createBranch<vector<float>>("boosted0_decay_h_dr");
    tx->createBranch<vector<float>>("boosted0_decay_h_deta");
    tx->createBranch<vector<float>>("boosted0_decay_h_dphi");
    tx->createBranch<vector<float>>("boosted0_decay_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted0_decay_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted0_lepton_p4");
    tx->createBranch<vector<int>>("boosted0_lepton_id");
    tx->createBranch<vector<int>>("boosted0_lepton_isstar");
    tx->createBranch<vector<float>>("boosted0_lepton_h_dr");
    tx->createBranch<vector<float>>("boosted0_lepton_h_deta");
    tx->createBranch<vector<float>>("boosted0_lepton_h_dphi");
    tx->createBranch<vector<float>>("boosted0_lepton_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted0_lepton_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted0_quark_p4");
    tx->createBranch<vector<int>>("boosted0_quark_id");
    tx->createBranch<vector<int>>("boosted0_quark_isstar");
    tx->createBranch<vector<float>>("boosted0_quark_h_dr");
    tx->createBranch<vector<float>>("boosted0_quark_h_deta");
    tx->createBranch<vector<float>>("boosted0_quark_h_dphi");
    tx->createBranch<vector<float>>("boosted0_quark_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted0_quark_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted250_decay_p4");
    tx->createBranch<vector<int>>("boosted250_decay_id");
    tx->createBranch<vector<int>>("boosted250_decay_isstar");
    tx->createBranch<vector<float>>("boosted250_decay_h_dr");
    tx->createBranch<vector<float>>("boosted250_decay_h_deta");
    tx->createBranch<vector<float>>("boosted250_decay_h_dphi");
    tx->createBranch<vector<float>>("boosted250_decay_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted250_decay_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted250_lepton_p4");
    tx->createBranch<vector<int>>("boosted250_lepton_id");
    tx->createBranch<vector<int>>("boosted250_lepton_isstar");
    tx->createBranch<vector<float>>("boosted250_lepton_h_dr");
    tx->createBranch<vector<float>>("boosted250_lepton_h_deta");
    tx->createBranch<vector<float>>("boosted250_lepton_h_dphi");
    tx->createBranch<vector<float>>("boosted250_lepton_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted250_lepton_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted250_quark_p4");
    tx->createBranch<vector<int>>("boosted250_quark_id");
    tx->createBranch<vector<int>>("boosted250_quark_isstar");
    tx->createBranch<vector<float>>("boosted250_quark_h_dr");
    tx->createBranch<vector<float>>("boosted250_quark_h_deta");
    tx->createBranch<vector<float>>("boosted250_quark_h_dphi");
    tx->createBranch<vector<float>>("boosted250_quark_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted250_quark_h_dphi_rotated");

    tx->createBranch<float>("boosted250_lqq_max_dr");
    tx->createBranch<float>("boosted250_lq0_dr");
    tx->createBranch<float>("boosted250_lq1_dr");
    tx->createBranch<float>("boosted250_qq_dr");

    tx->createBranch<vector<LV>>("boosted500_decay_p4");
    tx->createBranch<vector<int>>("boosted500_decay_id");
    tx->createBranch<vector<int>>("boosted500_decay_isstar");
    tx->createBranch<vector<float>>("boosted500_decay_h_dr");
    tx->createBranch<vector<float>>("boosted500_decay_h_deta");
    tx->createBranch<vector<float>>("boosted500_decay_h_dphi");
    tx->createBranch<vector<float>>("boosted500_decay_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted500_decay_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted500_lepton_p4");
    tx->createBranch<vector<int>>("boosted500_lepton_id");
    tx->createBranch<vector<int>>("boosted500_lepton_isstar");
    tx->createBranch<vector<float>>("boosted500_lepton_h_dr");
    tx->createBranch<vector<float>>("boosted500_lepton_h_deta");
    tx->createBranch<vector<float>>("boosted500_lepton_h_dphi");
    tx->createBranch<vector<float>>("boosted500_lepton_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted500_lepton_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted500_quark_p4");
    tx->createBranch<vector<int>>("boosted500_quark_id");
    tx->createBranch<vector<int>>("boosted500_quark_isstar");
    tx->createBranch<vector<float>>("boosted500_quark_h_dr");
    tx->createBranch<vector<float>>("boosted500_quark_h_deta");
    tx->createBranch<vector<float>>("boosted500_quark_h_dphi");
    tx->createBranch<vector<float>>("boosted500_quark_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted500_quark_h_dphi_rotated");

    tx->createBranch<float>("boosted500_lqq_max_dr");
    tx->createBranch<float>("boosted500_lq0_dr");
    tx->createBranch<float>("boosted500_lq1_dr");
    tx->createBranch<float>("boosted500_qq_dr");

    tx->createBranch<vector<LV>>("boosted1000_decay_p4");
    tx->createBranch<vector<int>>("boosted1000_decay_id");
    tx->createBranch<vector<int>>("boosted1000_decay_isstar");
    tx->createBranch<vector<float>>("boosted1000_decay_h_dr");
    tx->createBranch<vector<float>>("boosted1000_decay_h_deta");
    tx->createBranch<vector<float>>("boosted1000_decay_h_dphi");
    tx->createBranch<vector<float>>("boosted1000_decay_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted1000_decay_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted1000_lepton_p4");
    tx->createBranch<vector<int>>("boosted1000_lepton_id");
    tx->createBranch<vector<int>>("boosted1000_lepton_isstar");
    tx->createBranch<vector<float>>("boosted1000_lepton_h_dr");
    tx->createBranch<vector<float>>("boosted1000_lepton_h_deta");
    tx->createBranch<vector<float>>("boosted1000_lepton_h_dphi");
    tx->createBranch<vector<float>>("boosted1000_lepton_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted1000_lepton_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted1000_quark_p4");
    tx->createBranch<vector<int>>("boosted1000_quark_id");
    tx->createBranch<vector<int>>("boosted1000_quark_isstar");
    tx->createBranch<vector<float>>("boosted1000_quark_h_dr");
    tx->createBranch<vector<float>>("boosted1000_quark_h_deta");
    tx->createBranch<vector<float>>("boosted1000_quark_h_dphi");
    tx->createBranch<vector<float>>("boosted1000_quark_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted1000_quark_h_dphi_rotated");

    tx->createBranch<float>("boosted1000_lqq_max_dr");
    tx->createBranch<float>("boosted1000_lq0_dr");
    tx->createBranch<float>("boosted1000_lq1_dr");
    tx->createBranch<float>("boosted1000_qq_dr");

    tx->createBranch<vector<LV>>("boosted1500_decay_p4");
    tx->createBranch<vector<int>>("boosted1500_decay_id");
    tx->createBranch<vector<int>>("boosted1500_decay_isstar");
    tx->createBranch<vector<float>>("boosted1500_decay_h_dr");
    tx->createBranch<vector<float>>("boosted1500_decay_h_deta");
    tx->createBranch<vector<float>>("boosted1500_decay_h_dphi");
    tx->createBranch<vector<float>>("boosted1500_decay_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted1500_decay_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted1500_lepton_p4");
    tx->createBranch<vector<int>>("boosted1500_lepton_id");
    tx->createBranch<vector<int>>("boosted1500_lepton_isstar");
    tx->createBranch<vector<float>>("boosted1500_lepton_h_dr");
    tx->createBranch<vector<float>>("boosted1500_lepton_h_deta");
    tx->createBranch<vector<float>>("boosted1500_lepton_h_dphi");
    tx->createBranch<vector<float>>("boosted1500_lepton_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted1500_lepton_h_dphi_rotated");

    tx->createBranch<vector<LV>>("boosted1500_quark_p4");
    tx->createBranch<vector<int>>("boosted1500_quark_id");
    tx->createBranch<vector<int>>("boosted1500_quark_isstar");
    tx->createBranch<vector<float>>("boosted1500_quark_h_dr");
    tx->createBranch<vector<float>>("boosted1500_quark_h_deta");
    tx->createBranch<vector<float>>("boosted1500_quark_h_dphi");
    tx->createBranch<vector<float>>("boosted1500_quark_h_deta_rotated");
    tx->createBranch<vector<float>>("boosted1500_quark_h_dphi_rotated");

    tx->createBranch<float>("boosted1500_lqq_max_dr");
    tx->createBranch<float>("boosted1500_lq0_dr");
    tx->createBranch<float>("boosted1500_lq1_dr");
    tx->createBranch<float>("boosted1500_qq_dr");

    tx->clear();

}

//##############################################################################################################
void babyMaker_v2::AddWWWSignalOutput()
{
    tx->createBranch<int>("iswhwww");
    tx->createBranch<int>("www_channel"); // 0 == all hadronic, 1 == 1lep , 2 == 2SS, 3 == 3lep, 4 == 2OS
    tx->createBranch<int>("has_tau"); // 0 == all hadronic, 1 == 1lep , 2 == 2SS, 3 == 3lep, 4 == 2OS

    tx->createBranch<vector<LV>>("w_p4");
    tx->createBranch<vector<int>>("w_islep");
    tx->createBranch<vector<int>>("w_isstar");
    tx->createBranch<vector<int>>("w_isH");

    tx->createBranch<vector<LV>>("l_p4");
    tx->createBranch<vector<float>>("l_w_pt");
    tx->createBranch<vector<float>>("l_w_eta");
    tx->createBranch<vector<float>>("l_w_phi");
    tx->createBranch<vector<float>>("l_w_mass");
    tx->createBranch<vector<int>>("l_w_id");
    tx->createBranch<vector<int>>("l_isstar");
    tx->createBranch<vector<int>>("l_isH");
    tx->createBranch<vector<int>>("l_istau");

    tx->createBranch<vector<LV>>("q_p4");
    tx->createBranch<vector<float>>("q_w_pt");
    tx->createBranch<vector<float>>("q_w_eta");
    tx->createBranch<vector<float>>("q_w_phi");
    tx->createBranch<vector<float>>("q_w_mass");
    tx->createBranch<vector<int>>("q_w_id");
    tx->createBranch<vector<int>>("q_isstar");
    tx->createBranch<vector<int>>("q_isH");

    tx->createBranch<float>("dRllSS");
    tx->createBranch<float>("dRqqSS");

    tx->createBranch<float>("DPhill_higgs");
    tx->createBranch<float>("Mll_higgs");
    tx->createBranch<float>("MT_higgs");

    tx->clear();

}

//##############################################################################################################
void babyMaker_v2::AddWHsusyOutput()
{
    h_nevents_SMS = new TH3F("h_counterSMS", "h_counterSMS", 37, 99, 1024, 19, -1, 474, 35, 0.5, 35.5); //15000 bins!
    h_nevents_SMS->Sumw2();
    h_nrawevents_SMS = new TH2F("histNEvts", "h_histNEvts", 37, 99, 1024, 19, -1, 474); //x=mStop, y=mLSP
    h_nrawevents_SMS->Sumw2();
    fxsec = new TFile("xsec_susy_13tev.root", "READ");
    if (fxsec -> IsZombie()) {
        std::cout << "Somehow xsec_stop_13TeV.root is corrupted. Exit..." << std::endl;
        exit(3);
    }
    hxsec = (TH1D *) fxsec -> Get("h_xsec_c1n2");
    tx->createBranch<float>("chimass");
    tx->createBranch<float>("lspmass");

    tx->clear();

}

//##############################################################################################################
void babyMaker_v2::SetYear()
{
    TString filename = looper.getCurrentFileName();
    if (
            filename.Contains("Run2016")
            || filename.Contains("Moriond17")
            || filename.Contains("RunIISummer16")
            || filename.Contains("run2_data2016")
            || filename.Contains("run2_moriond17")
            || filename.Contains("V00-00-02_2017Sep27")
       )
    {
        gconf.year = 2016;
        return;
    }
    if (
            filename.Contains("Run2017")
            || filename.Contains("RunIIFall17")
            || filename.Contains("_mc2017_")
            || filename.Contains("run2_mc2017")
       )
    {
        gconf.year = 2017;
        return;
    }
    if (
            filename.Contains("Run2018")
            || filename.Contains("RunIISpring18")
            || filename.Contains("RunIISummer18")
            || filename.Contains("run2_mc2018")
       )
    {
        gconf.year = 2018;
        return;
    }
    std::cout << "ERROR: did not recognize which year" << std::endl;
    FATALERROR(__FUNCTION__);
}

//##############################################################################################################
void babyMaker_v2::SetLeptonID()
{
    // Set lepton ID based on baby mode
    switch (babymode)
    {
        case kWWWBaby:   SetWWWAnalysisLeptonID(); break;
        case kFRBaby:    SetWWWAnalysisLeptonID(); break;
        case kOSBaby:    SetWWWAnalysisLeptonID(); break;
        case kTnPBaby:   SetWWWAnalysisLeptonID(); break;
        case kAllBaby:   SetWWWAnalysisLeptonID(); break;
        case kPOGBaby:   SetPOGAnalysisLeptonID(); break;
        case kLooseBaby: SetWWWAnalysisLeptonID(); break;
        default: return;
    }
}

//##############################################################################################################
void babyMaker_v2::SetWWWAnalysisLeptonID()
{
    // Setting general configuration
    gconf.wwwcfg["vetoid"] = "lep_pass_VVV_cutbased_veto";
    gconf.wwwcfg["looseid"] = "lep_pass_VVV_cutbased_fo";
    gconf.wwwcfg["tightid"] = "lep_pass_VVV_cutbased_tight";
    gconf.wwwcfg["3llooseid"] = "lep_pass_VVV_cutbased_3l_fo";
    gconf.wwwcfg["3ltightid"] = "lep_pass_VVV_cutbased_3l_tight";

    if (gconf.year == 2016)
    {
        gconf.ea_version = 2;
        //_________________________________
        // Isolation configuration
        // Same-sign muons
        gconf.mu_reliso_veto      = 0.4;
        gconf.mu_reliso_fo        = 0.4;
        gconf.mu_reliso_tight     = 0.03;
        gconf.mu_addlep_veto      = true; // false is what i think it should be at least when we get to do reanalysis
        gconf.mu_addlep_fo        = true;
        gconf.mu_addlep_tight     = true;
        // Same-sign electrons
        gconf.el_reliso_veto      = 0.4;
        gconf.el_reliso_fo        = 0.4;
        gconf.el_reliso_tight     = 0.03;
        gconf.el_addlep_veto      = true; // false is what i think it should be at least when we get to do reanalysis
        gconf.el_addlep_fo        = true;
        gconf.el_addlep_tight     = true;
        // Three-lepton muons (Shares same veto as same-sign)
        gconf.mu_reliso_3l_fo     = 0.4;
        gconf.mu_reliso_3l_tight  = 0.07;
        gconf.mu_addlep_3l_fo     = true;
        gconf.mu_addlep_3l_tight  = true;
        // Three-lepton electrons (Shares same veto as same-sign)
        gconf.el_reliso_3l_fo     = 0.4;
        gconf.el_reliso_3l_tight  = 0.05;
        gconf.el_addlep_3l_fo     = true;
        gconf.el_addlep_3l_tight  = true;
    }
    else if (gconf.year == 2017)
    {
        gconf.ea_version = 4;
        //_________________________________
        // Isolation configuration
        // Same-sign muons
        gconf.mu_reliso_veto      = 0.4;
        gconf.mu_reliso_fo        = 0.4;
        gconf.mu_reliso_tight     = 0.03;
        gconf.mu_addlep_veto      = false;
        gconf.mu_addlep_fo        = true;
        gconf.mu_addlep_tight     = true;
        // Same-sign electrons
        gconf.el_reliso_veto      = 0.4;
        gconf.el_reliso_fo        = 0.4;
        gconf.el_reliso_tight     = 0.03;
        gconf.el_addlep_veto      = false;
        gconf.el_addlep_fo        = true;
        gconf.el_addlep_tight     = true;
        // Three-lepton muons (Shares same veto as same-sign)
        gconf.mu_reliso_3l_fo     = 0.4;
        gconf.mu_reliso_3l_tight  = 0.07;
        gconf.mu_addlep_3l_fo     = true;
        gconf.mu_addlep_3l_tight  = true;
        // Three-lepton electrons (Shares same veto as same-sign)
        gconf.el_reliso_3l_fo     = 0.4;
        gconf.el_reliso_3l_tight  = 0.05;
        gconf.el_addlep_3l_fo     = true;
        gconf.el_addlep_3l_tight  = true;
    }
    else
    {
        std::cout << "year not recognized! gconf.year = " << gconf.year << std::endl;
        FATALERROR(__FUNCTION__);
    }
}

//##############################################################################################################
void babyMaker_v2::SetPOGAnalysisLeptonID()
{
    gconf.wwwcfg["vetoid"] = "lep_pass_POG_veto";
    gconf.wwwcfg["looseid"] = "lep_pass_POG_loose";
    gconf.wwwcfg["tightid"] = "lep_pass_POG_tight";
    gconf.wwwcfg["3llooseid"] = "lep_pass_POG_loose";
    gconf.wwwcfg["3ltightid"] = "lep_pass_POG_tight";

    if (gconf.year == 2016 || gconf.year == 2017)
    {
        gconf.ea_version = 2;
        //_________________________________
        // Isolation configuration
        // Same-sign muons
        gconf.mu_reliso_veto      = 0.4;
        gconf.mu_reliso_fo        = 0.4;
        gconf.mu_reliso_tight     = 0.03;
        gconf.mu_addlep_veto      = false;
        gconf.mu_addlep_fo        = true;
        gconf.mu_addlep_tight     = true;
        // Same-sign electrons
        gconf.el_reliso_veto      = 0.4;
        gconf.el_reliso_fo        = 0.4;
        gconf.el_reliso_tight     = 0.03;
        gconf.el_addlep_veto      = false;
        gconf.el_addlep_fo        = true;
        gconf.el_addlep_tight     = true;
        // Three-lepton muons (Shares same veto as same-sign)
        gconf.mu_reliso_3l_fo     = 0.4;
        gconf.mu_reliso_3l_tight  = 0.07;
        gconf.mu_addlep_3l_fo     = true;
        gconf.mu_addlep_3l_tight  = true;
        // Three-lepton electrons (Shares same veto as same-sign)
        gconf.el_reliso_3l_fo     = 0.4;
        gconf.el_reliso_3l_tight  = 0.05;
        gconf.el_addlep_3l_fo     = true;
        gconf.el_addlep_3l_tight  = true;
    }
    else
    {
        std::cout << "year not recognized! gconf.year = " << gconf.year << std::endl;
        FATALERROR(__FUNCTION__);
    }
}

//##############################################################################################################
void babyMaker_v2::SaveOutput()
{

    // change directory to output file
    ofile->cd();

    // This is always saved
    h_neventsinfile->Write();

    // Based on the baby mode now preselect and if it doesn't pass return
    switch (babymode)
    {
        case kWWWBaby:   SaveWWWBaby(); break;
        case kFRBaby:    SaveWWWBaby(); break;
        case kOSBaby:    SaveWWWBaby(); break;
        case kTnPBaby:   SaveTnPBaby(); break;
        case kAllBaby:   SaveWWWBaby(); break;
        case kPOGBaby:   SavePOGBaby(); break;
        case kLooseBaby: SaveWWWBaby(); break;
        default: return;
    }

}

//##############################################################################################################
void babyMaker_v2::SaveWWWBaby()
{
    ofile->cd();

    t_os      = t->CopyTree("((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]<0))");
    t_ss      = t->CopyTree("((nVlep>=3)||((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]>0)))");

    t_www     = t->CopyTree("((nVlep>=3)||((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]>0)))*(bkgtype==\"WWW\")||(bkgtype==\"WHtoWWW\")");
    t_qflip   = t->CopyTree("((nVlep>=3)||((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]>0)))*(bkgtype==\"chargeflips\")");
    t_photon  = t->CopyTree("((nVlep>=3)||((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]>0)))*(bkgtype==\"photonfakes\")");
    t_fakes   = t->CopyTree("((nVlep>=3)||((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]>0)))*(bkgtype==\"fakes\")");
    t_prompt  = t->CopyTree("((nVlep>=3)||((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]>0)))*((nVlep==2&&bkgtype==\"trueSS\")||(nVlep>=3&&bkgtype==\"trueWWW\"))");
    t_lostlep = t->CopyTree("((nVlep>=3)||((nVlep==2)&&(lep_pdgId[0]*lep_pdgId[1]>0)))*((nVlep==2&&bkgtype==\"SSLL\")||(nVlep>=3&&(bkgtype==\"3lLL\"||bkgtype==\"true3L\")))");

    t_os->SetName("t_os");
    t_ss->SetName("t_ss");
    t_www->SetName("t_www");
    t_qflip->SetName("t_qflip");
    t_photon->SetName("t_photon");
    t_fakes->SetName("t_fakes");
    t_prompt->SetName("t_prompt");
    t_lostlep->SetName("t_lostlep");

    t_os->SetTitle("Opposite Sign Events");
    t_ss->SetTitle("Same Sign Events");
    t_www->SetTitle("Signal Events");
    t_qflip->SetTitle("Charge Flip Events");
    t_photon->SetTitle("Prompt Photon Events");
    t_fakes->SetTitle("Fake Lepton Events");
    t_prompt->SetTitle("Prompt Lepton Events");
    t_lostlep->SetTitle("Lost Lepton Events or WZ in 3-lepton categories");

    t->Write();
    t_os->Write();
    t_ss->Write();
    t_www->Write();
    t_qflip->Write();
    t_photon->Write();
    t_fakes->Write();
    t_prompt->Write();
    t_lostlep->Write();

    if (isWHSUSY())
    {
        h_nevents_SMS->Write();
        h_nrawevents_SMS->Write();
    }
}

//##############################################################################################################
void babyMaker_v2::SavePOGBaby()
{
    SaveWWWBaby();
}

//##############################################################################################################
void babyMaker_v2::SaveTnPBaby()
{
    t->Write();
    ofile->cd();
}


//##############################################################################################################
void babyMaker_v2::Process()
{

    // Process leptons via CoreUtil
    ProcessLeptons();

    if (!isPass())
        return;

    // Process other non-lepton objects via CoreUtil
    ProcessNonLeptonObjects();

    // Fill the output ttree
    FillBaby();

}

//##############################################################################################################
void babyMaker_v2::ProcessLeptons()
{
    // Loop over electrons
    ProcessElectrons();

    // Loop over muons
    ProcessMuons();

}

//##############################################################################################################
bool babyMaker_v2::isPass()
{
    // Based on the baby mode now preselect and if it doesn't pass return
    switch (babymode)
    {
        case kWWWBaby:   if (!PassWWWPreselection  ()) return false; break;
        case kFRBaby:    if (!PassFRPreselection   ()) return false; break;
        case kOSBaby:    if (!PassOSPreselection   ()) return false; break;
        case kTnPBaby:   if (!PassTnPPreselection  ()) return false; break;
        case kAllBaby:   /* no cut is applied accept all */  break;
        case kPOGBaby:   if (!PassPOGPreselection  ()) return false; break;
        case kLooseBaby: if (!PassLoosePreselection()) return false; break;
        default: return false;
    }

    // If it passed then accept
    return true;
}

//##############################################################################################################
void babyMaker_v2::ProcessNonLeptonObjects()
{
    // Triggers
    ProcessTriggers();

    // Loop over gen particles
    ProcessGenParticles();

    // Loop over Jets
    ProcessJets();

    // Loop over fatJets
    ProcessFatJets();

    // Process MET (recalculate etc.)
    ProcessMET();

    // Loop over charged particle candidates
    ProcessTracks();
}

//##############################################################################################################
void babyMaker_v2::FillBaby()
{
    // Now process the rest of the stuff
    switch (babymode)
    {
        case kWWWBaby:   FillWWWBaby(); break;
        case kFRBaby:    FillWWWBaby(); break;
        case kOSBaby:    FillWWWBaby(); break;
        case kTnPBaby:   FillTnPBaby(); break;
        case kAllBaby:   FillWWWBaby(); break;
        case kPOGBaby:   FillPOGBaby(); break;
        case kLooseBaby: FillWWWBaby(); break;
    }
}

//##############################################################################################################
void babyMaker_v2::FillWWWBaby()
{
    if (coreElectron.index.size() + coreMuon.index.size() < 1)
        return;

    // Fill baby ntuple branches with event information (evt, lumi etc.)
    // These are always created and filled
    FillEventInfo();

    // Truth level variables studies for WWW
    FillTruthLevelStudyVariables();

    // Fill baby ntuple branches corresponding to lepton indices
    FillElectrons();

    // Fill baby ntuple branches corresponding to lepton indices
    FillMuons();

    // Organize leptons by sorting by their pt
    SortLeptonBranches();

    // Fill baby ntuple branches for jets
    FillJets();

    // Fill baby ntuple branches for fat-jets
    FillFatJets();

    // Organize jets by sorting by their pt
    SortJetBranches();

    // Fill baby ntuple branches for emu channel for easy lepton SF application
    FillEMuLeptons();

    // Fill baby ntuple branches for MET
    FillMET();

    // Fill baby ntuple branches for track related variables (iso track veto)
    FillTracks();

    // Fill generatore level particles
    FillGenParticles();

    // Fill vertex info
    FillVertexInfo();

    // Fill MET filter info
    FillMETFilter();

    // Fill summary variables
    FillSummaryVariables();

    // Fill trigger bits (This comes after summary variables and leptons
    FillTrigger();

    // Fill Weights
    FillWeights(); // TODO weight has ID dependency // TODO How to get POG root file? turn them into functions?

    // Fill baby ntuple branches for lbnt lepton for easy fake rate application
    FillLbntLeptons();

    // Fill year dependent info
    FillYearInfo();

    // Fill TTree (NOTE: also clears internal variables)
    FillTTree();
}

//##############################################################################################################
void babyMaker_v2::FillPOGBaby()
{
    // Fill TTree (NOTE: also clears internal variables)
    FillWWWBaby();
}

//##############################################################################################################
void babyMaker_v2::FillTnPBaby()
{
    // Fill Muons
    FillTnPMuons();

    // Fill Electrons
    FillTnPElectrons();
}

//##############################################################################################################
void babyMaker_v2::ProcessTriggers() { coreTrigger.process(); }

//##############################################################################################################
void babyMaker_v2::ProcessGenParticles() { coreGenPart.process(); }

//##############################################################################################################
void babyMaker_v2::ProcessElectrons()
{
    // Based on the baby mode Process different set of electrons
    switch (babymode)
    {
        case kWWWBaby:   ProcessNominalElectrons(); return; break;
        case kFRBaby:    ProcessNominalElectrons(); return; break;
        case kOSBaby:    ProcessNominalElectrons(); return; break;
        case kTnPBaby:   ProcessTnPElectrons();     return; break;
        case kAllBaby:   ProcessNominalElectrons(); return; break;
        case kPOGBaby:   ProcessPOGElectrons();     return; break;
        case kLooseBaby: ProcessNominalElectrons(); return; break;
        default: return;
    }
}

//##############################################################################################################
void babyMaker_v2::ProcessNominalElectrons() { coreElectron.process(isVetoElectron); }

//##############################################################################################################
void babyMaker_v2::ProcessPOGElectrons() { coreElectron.process(is2017POGVetoElectron); }

//##############################################################################################################
void babyMaker_v2::ProcessTnPElectrons() { coreElectron.process(isProbeElectron, isTagElectron); }

//##############################################################################################################
void babyMaker_v2::ProcessMuons()
{
    // Based on the baby mode now preselect and if it doesn't pass return
    switch (babymode)
    {
        case kWWWBaby:   ProcessNominalMuons(); return; break;
        case kFRBaby:    ProcessNominalMuons(); return; break;
        case kOSBaby:    ProcessNominalMuons(); return; break;
        case kTnPBaby:   ProcessTnPMuons();     return; break;
        case kAllBaby:   ProcessNominalMuons(); return; break;
        case kPOGBaby:   ProcessPOGMuons();     return; break;
        case kLooseBaby: ProcessNominalMuons(); return; break;
        default: return;
    }
}

//##############################################################################################################
void babyMaker_v2::ProcessNominalMuons() { coreMuon.process(isVetoMuon); }

//##############################################################################################################
void babyMaker_v2::ProcessPOGMuons() { coreMuon.process(is2017POGVetoMuon); }

//##############################################################################################################
void babyMaker_v2::ProcessTnPMuons() { coreMuon.process(isProbeMuon, isTagMuon); }

//##############################################################################################################
void babyMaker_v2::ProcessJets() { coreJet.process(coreJec, coreJer); }

//##############################################################################################################
void babyMaker_v2::ProcessFatJets() { coreFatJet.process(coreJec); }

//##############################################################################################################
void babyMaker_v2::ProcessMET() { coreMET.process(coreJec, coreJer); }

//##############################################################################################################
void babyMaker_v2::ProcessTracks() { coreTrack.process(); }

//##############################################################################################################
bool babyMaker_v2::PassWWWPreselection()
{
    return PassPresel();
}

//##############################################################################################################
bool babyMaker_v2::PassPresel()
{
    return PassPresel_v3();
}

//##############################################################################################################
bool babyMaker_v2::PassPresel_v1()
{
    // Select 2 SS lepton events or 3 or more lepton events
    vector<int> el_idx = coreElectron.index;
    vector<int> mu_idx = coreMuon.index;

    // If 3 or more veto lepton events. then require 3 and only 3 loose lepton events
    // Such requirement ensures simplicity in various things like SFOS counter and fake rate application
    // Also, require that it is total charge == 1
    if (el_idx.size() + mu_idx.size() > 2)
    {
        int nloose = 0;
        int chargesum = 0;
        for (auto& iel : coreElectron.index)
        {
            if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_FO_3L))
            {
                nloose++;
                chargesum += cms3.els_charge()[iel];
            }
        }
        for (auto& imu : coreMuon.index)
        {
            if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_FO_3L))
            {
                nloose++;
                chargesum += cms3.mus_charge()[imu];
            }
        }
        if (nloose == 3)
            return abs(chargesum) == 1;
        else
            return false;
    }
    // If less than 2 leptons skip
    if (el_idx.size() + mu_idx.size() < 2) return false;
    // If equal to 2 leptons then must be a same-sign and at must have two loose lepton
    int nloose = 0;
    for (auto& iel : coreElectron.index)
    {
        if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_FO_SS))
            nloose++;
    }
    for (auto& imu : coreMuon.index)
    {
        if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_FO_SS))
            nloose++;
    }
    if (nloose != 2)
        return false;
    return true;
    if (mu_idx.size() == 2)
    {
        if (cms3.mus_charge()[mu_idx[0]] * cms3.mus_charge()[mu_idx[1]] > 0)
            return true;
        else
            return false;
    }
    else if (el_idx.size() == 2)
    {
        if (cms3.els_charge()[el_idx[0]] * cms3.els_charge()[el_idx[1]] > 0)
            return true;
        else
            return false;
    }
    else if (mu_idx.size() == 1 && el_idx.size() == 1)
    {
        if (cms3.mus_charge()[mu_idx[0]] * cms3.els_charge()[el_idx[0]] > 0)
            return true;
        else
            return false;
    }
    else
    {
        FATALERROR(__FUNCTION__);
        return false;
    }
}

//##############################################################################################################
bool babyMaker_v2::PassPresel_v2()
{
    // Select 2 SS lepton events or 3 or more lepton events
    vector<int> el_idx = coreElectron.index;
    vector<int> mu_idx = coreMuon.index;

    int nveto = el_idx.size() + mu_idx.size();

    // If more than 3 veto lepton events, then dump.
    if (nveto > 3)
        return false;

    // If 3 veto lepton events, then require at least 2 tight 3lid leptons (n.b. double fakes are not possible.)
    // Also, require that it is total charge == 1
    if (nveto == 3)
    {
        int ntight = 0;
        int chargesum = 0;
        for (auto& iel : coreElectron.index)
        {
            if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_TIGHT_3L))
                ntight++;
            chargesum += cms3.els_charge()[iel];
        }
        for (auto& imu : coreMuon.index)
        {
            if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_TIGHT_3L))
                ntight++;
            chargesum += cms3.mus_charge()[imu];
        }
        if (ntight >= 2)
            return abs(chargesum) == 1;
        else
            return false;
    }

    // If 2 veto lepton events, then require at least one tight ss lepton (n.b. double fakes are not possible.)
    if (nveto == 2)
    {
        int ntight = 0;
        int chargesum = 0;
        for (auto& iel : coreElectron.index)
        {
            if (cms3.els_p4()[iel].pt() > 25. && passElectronSelection_VVV(iel, VVV_TIGHT_SS))
                ntight++;
            chargesum += cms3.els_charge()[iel];
        }
        for (auto& imu : coreMuon.index)
        {
            if (cms3.mus_p4()[imu].pt() > 25. && passMuonSelection_VVV(imu, VVV_TIGHT_SS))
                ntight++;
            chargesum += cms3.mus_charge()[imu];
        }
        if (ntight >= 1)
        {
            if (abs(chargesum) == 2) // If same-sign accept
                return true;
            else if (ntight == 2) // If two-tight accept
                return true;
            else
                return false; // If not same-sign and not two-tight discard
        }
        else
        {
            return false;
        }
    }

    // If 1 veto lepton events, then write out the whole thing
    if (nveto == 1)
    {
        // Check if data,
        bool isdata = isData();
        if (isdata)
        {
            // If data then check the triggers
            // These triggers are checked in coreutil, but to optimize the code performance I hand check them if data
            // I don't wish to run coreTrigger.process() for all events
            if (el_idx.size() == 1)
            {
                int HLT_SingleIsoEl17 = passHLTTriggerPattern("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
                if (HLT_SingleIsoEl17 > 0) return true;
                int HLT_SingleIsoEl23 = passHLTTriggerPattern("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
                if (HLT_SingleIsoEl23 > 0) return true;
                int HLT_SingleIsoEl8 = passHLTTriggerPattern("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
                if (HLT_SingleIsoEl8 > 0) return true;
            }
            else if (mu_idx.size() == 1)
            {
                int HLT_SingleIsoMu17 = passHLTTriggerPattern("HLT_Mu17_TrkIsoVVL_v");
                if (HLT_SingleIsoMu17 > 0) return true;
                int HLT_SingleIsoMu8 = passHLTTriggerPattern("HLT_Mu8_TrkIsoVVL_v");
                if (HLT_SingleIsoMu8 > 0) return true;
            }
            // If it reaches this point, then it means that none of the trigger passed
            return false;
        }
        return true;
    }

    // Anything else dump
    return false;
}

//##############################################################################################################
bool babyMaker_v2::PassPresel_v3()
{

    if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
    {
        std::cout << std::endl;
        std::cout << "This event is in the txt" << std::endl;
        std::cout <<  " coreElectron.index.size(): " << coreElectron.index.size() <<  " coreMuon.index.size(): " << coreMuon.index.size() <<  std::endl;
        std::cout <<  " coreElectron.index.size()+coreMuon.index.size(): " << coreElectron.index.size()+coreMuon.index.size() <<  std::endl;
        std::cout <<  " cms3.evt_run(): " << cms3.evt_run() <<  " cms3.evt_lumiBlock(): " << cms3.evt_lumiBlock() <<  " cms3.evt_event(): " << cms3.evt_event() <<  std::endl;
        for (auto& iel : coreElectron.index)
        {
            std::cout <<  " passElectronSelection_VVV(iel,VVV_FO_3L): " << passElectronSelection_VVV(iel,VVV_FO_3L) <<  " passElectronSelection_VVV(iel,VVV_TIGHT_3L): " << passElectronSelection_VVV(iel,VVV_TIGHT_3L) <<  std::endl;
            std::cout <<  " passElectronSelection_VVV(iel,VVV_cutbased_3l_fo_v4): " << passElectronSelection_VVV(iel,VVV_cutbased_3l_fo_v4) <<  " passElectronSelection_VVV(iel,VVV_cutbased_3l_tight_v4): " << passElectronSelection_VVV(iel,VVV_cutbased_3l_tight_v4) <<  std::endl;
            std::cout <<  " fabs(cms3.els_etaSC().at(iel)): " << fabs(cms3.els_etaSC().at(iel)) <<  std::endl;
            std::cout <<  " fabs(cms3.els_dxyPV().at(iel)): " << fabs(cms3.els_dxyPV().at(iel)) <<  std::endl;
            std::cout <<  " fabs(cms3.els_dzPV().at(iel)): " << fabs(cms3.els_dzPV().at(iel)) <<  std::endl;
            std::cout <<  " getMVAoutput(iel): " << getMVAoutput(iel) <<  std::endl;
            std::cout <<  " eleRelIso03EA(iel,2,false): " << eleRelIso03EA(iel,2,false) <<  std::endl;
            std::cout <<  " eleRelIso03EA(iel,2,true): " << eleRelIso03EA(iel,2,true) <<  std::endl;
            std::cout <<  " fabs(cms3.els_ip3d()[iel]): " << fabs(cms3.els_ip3d()[iel]) <<  std::endl;
            std::cout <<  " threeChargeAgree(iel): " << threeChargeAgree(iel) <<  std::endl;
            std::cout <<  " isTriggerSafenoIso_v1(iel): " << isTriggerSafenoIso_v1(iel) <<  std::endl;
            std::cout <<  " isTriggerSafe_v1(iel): " << isTriggerSafe_v1(iel) <<  std::endl;
        }
        for (auto& imu : coreMuon.index)
        {
            std::cout <<  " passMuonSelection_VVV(imu,VVV_FO_3L): " << passMuonSelection_VVV(imu,VVV_FO_3L) <<  " passMuonSelection_VVV(imu,VVV_TIGHT_3L): " << passMuonSelection_VVV(imu,VVV_TIGHT_3L) <<  std::endl;
        }
    }

    // Select 2 SS lepton events or 3 or more lepton events
    vector<int> el_idx = coreElectron.index;
    vector<int> mu_idx = coreMuon.index;

    // If 3 or more veto lepton events. then require the event contains 3 loose and 2 or more tight leptons
    // Such requirement ensures simplicity in various things like SFOS counter and fake rate application
    // Also, require that it is total charge == 1
    if (el_idx.size() + mu_idx.size() > 2)
    {
//        if (el_idx.size() + mu_idx.size() != 3)
//            return false;
        int nloose = 0;
        int ntight = 0;
        int chargesum = 0;
        for (auto& iel : coreElectron.index)
        {
            if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_FO_3L))
            {
                nloose++;
                chargesum += cms3.els_charge()[iel];
            }
            if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_TIGHT_3L))
            {
                ntight++;
            }
        }
        for (auto& imu : coreMuon.index)
        {
            if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_FO_3L))
            {
                nloose++;
                chargesum += cms3.mus_charge()[imu];
            }
            if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_TIGHT_3L))
            {
                ntight++;
            }
        }
        if (nloose == 3 && ntight >= 2)
            return abs(chargesum) == 1;
        else
            return false;
    }
    // If less than 2 leptons skip
    if (el_idx.size() + mu_idx.size() < 2) return false;
    // If equal to 2 leptons then must be a same-sign and at must have two loose lepton and one or more tight lepton
    int nloose = 0;
    int ntight = 0;
    for (auto& iel : coreElectron.index)
    {
        if (cms3.els_p4()[iel].pt() > 25. && passElectronSelection_VVV(iel, VVV_FO_SS))
            nloose++;
        if (cms3.els_p4()[iel].pt() > 25. && passElectronSelection_VVV(iel, VVV_TIGHT_SS))
            ntight++;
    }
    for (auto& imu : coreMuon.index)
    {
        if (cms3.mus_p4()[imu].pt() > 25. && passMuonSelection_VVV(imu, VVV_FO_SS))
            nloose++;
        if (cms3.mus_p4()[imu].pt() > 25. && passMuonSelection_VVV(imu, VVV_TIGHT_SS))
            ntight++;
    }
    if (nloose != 2 || ntight < 1)
        return false;
    if (mu_idx.size() == 2)
    {
        if (cms3.mus_charge()[mu_idx[0]] * cms3.mus_charge()[mu_idx[1]] > 0)
            return true;
        else
            return false;
    }
    else if (el_idx.size() == 2)
    {
        if (cms3.els_charge()[el_idx[0]] * cms3.els_charge()[el_idx[1]] > 0)
            return true;
        else
            return false;
    }
    else if (mu_idx.size() == 1 && el_idx.size() == 1)
    {
        if (cms3.mus_charge()[mu_idx[0]] * cms3.els_charge()[el_idx[0]] > 0)
            return true;
        else
            return false;
    }
    else
    {
        FATALERROR(__FUNCTION__);
        return false;
    }
}

//##############################################################################################################
bool babyMaker_v2::PassFRPreselection()
{
    // Select 2 SS lepton events or 3 or more lepton events
    vector<int> el_idx = coreElectron.index;
    vector<int> mu_idx = coreMuon.index;

    if (el_idx.size() + mu_idx.size() < 1)
        return false;

//    if (el_idx.size() + mu_idx.size() == 1)
//    {
        // Check if data,
        if (isData())
        {
            if (coreSample.is2016(looper.getCurrentFileName()))
            {
                // If data then check the triggers
                // These triggers are checked in coreutil, but to optimize the code performance I hand check them if data
                // I don't wish to run coreTrigger.process() for all events
                int HLT_SingleIsoEl8;
                int HLT_SingleIsoEl17;
                int HLT_SingleIsoEl23;
                int HLT_SingleIsoMu8;
                int HLT_SingleIsoMu17;
                setHLTBranch("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v" ,  true, HLT_SingleIsoEl8 );
                setHLTBranch("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",  true, HLT_SingleIsoEl17 );
                setHLTBranch("HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v",  true, HLT_SingleIsoEl23 );
                setHLTBranch("HLT_Mu8_TrkIsoVVL_v",  true, HLT_SingleIsoMu8 );
                setHLTBranch("HLT_Mu17_TrkIsoVVL_v",  true, HLT_SingleIsoMu17 );
                if (HLT_SingleIsoEl17 > 0) return true;
                if (HLT_SingleIsoEl23 > 0) return true;
                if (HLT_SingleIsoEl8 > 0) return true;
                if (HLT_SingleIsoMu17 > 0) return true;
                if (HLT_SingleIsoMu8 > 0) return true;
                // If it reaches this point, then it means that none of the trigger passed
                return false;
            }
            else if (coreSample.is2017(looper.getCurrentFileName()))
            {
                // If data then check the triggers
                // These triggers are checked in coreutil, but to optimize the code performance I hand check them if data
                // I don't wish to run coreTrigger.process() for all events
                int HLT_SingleEl8;
                int HLT_SingleEl17;
                int HLT_SingleIsoEl8;
                int HLT_SingleIsoEl23;
                int HLT_SingleIsoMu8;
                int HLT_SingleIsoMu17;
                setHLTBranch("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v"       ,  true, HLT_SingleEl8 );
                setHLTBranch("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v"      ,  true, HLT_SingleEl17 );
                setHLTBranch("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v" ,  true, HLT_SingleIsoEl8 );
                setHLTBranch("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v",  true, HLT_SingleIsoEl23 );
                setHLTBranch("HLT_Mu8_TrkIsoVVL_v",  true, HLT_SingleIsoMu8 );
                setHLTBranch("HLT_Mu17_TrkIsoVVL_v",  true, HLT_SingleIsoMu17 );
                if (HLT_SingleEl8 > 0) return true;
                if (HLT_SingleEl17 > 0) return true;
                if (HLT_SingleIsoEl23 > 0) return true;
                if (HLT_SingleIsoEl8 > 0) return true;
                if (HLT_SingleIsoMu17 > 0) return true;
                if (HLT_SingleIsoMu8 > 0) return true;
                // If it reaches this point, then it means that none of the trigger passed
                return false;
            }
        }
        return true;
//    }
//    else
//    {
//        return false;
//    }
}

//##############################################################################################################
bool babyMaker_v2::PassOSPreselection()
{
    // Select 2 SS lepton events or 3 or more lepton events
    vector<int> el_idx = coreElectron.index;
    vector<int> mu_idx = coreMuon.index;

    // If less than 2 veto leptons
    if (el_idx.size() + mu_idx.size() < 2)
        return false;

    // The 2 veto leptons of 20 GeV each at the least must pass 3L tight options
    int ntight = 0;
    int nloose = 0;
    for (auto& iel : coreElectron.index)
    {
        if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_TIGHT_3L))
            ntight++;
        if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_FO_3L))
            nloose++;
    }
    for (auto& imu : coreMuon.index)
    {
        if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_TIGHT_3L))
            ntight++;
        if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_FO_3L))
            nloose++;
    }

    if (ntight >= 1 && nloose >= 2)
        return true;
    else
        return false;
}

//##############################################################################################################
bool babyMaker_v2::PassTnPPreselection()
{
    return true;
}

//##############################################################################################################
bool babyMaker_v2::PassPOGPreselection()
{
    // Select 2 SS lepton events or 3 or more lepton events
    vector<int> el_idx = coreElectron.index;
    vector<int> mu_idx = coreMuon.index;

    if (el_idx.size() + mu_idx.size() >= 2)
    {
//        std::cout <<  " el_idx.size(): " << el_idx.size() <<  " mu_idx.size(): " << mu_idx.size() <<  std::endl;
        int nloose = 0;
        for (auto& iel : coreElectron.index) if (is2017POGLooseElectron(iel)) nloose++;
        for (auto& imu : coreMuon.index)     if (is2017POGLooseMuon(imu)) nloose++;
//        std::cout <<  " nloose: " << nloose <<  std::endl;
        if (nloose >= 2)
            return true;
        else
            return false;
    }
    else
    {
        return false;
    }
}

//##############################################################################################################
bool babyMaker_v2::PassLoosePreselection()
{

    if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
    {
        std::cout <<  " cms3.evt_run(): " << cms3.evt_run() <<  " cms3.evt_lumiBlock(): " << cms3.evt_lumiBlock() <<  " cms3.evt_event(): " << cms3.evt_event() <<  std::endl;
        std::cout << "here" << std::endl;
        coreGenPart.process();
        coreGenPart.printAllParticles();
        abort();
    }
    else
    {
        return false;
    }

    // Select 2 SS lepton events or 3 or more lepton events
    vector<int> el_idx = coreElectron.index;
    vector<int> mu_idx = coreMuon.index;

    // If 3 or more leptons
    if (el_idx.size() + mu_idx.size() > 2)
    {
        // Cut out 4 or more leptons
        if (el_idx.size() + mu_idx.size() != 3)
            return false;
        int nloose = 0;
        int ntight = 0;
        int chargesum = 0;
        for (auto& iel : coreElectron.index)
        {
            if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_FO_3L))
            {
                nloose++;
                chargesum += cms3.els_charge()[iel];
            }
            if (cms3.els_p4()[iel].pt() > 20. && passElectronSelection_VVV(iel, VVV_TIGHT_3L))
            {
                ntight++;
            }
        }
        for (auto& imu : coreMuon.index)
        {
            if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_FO_3L))
            {
                nloose++;
                chargesum += cms3.mus_charge()[imu];
            }
            if (cms3.mus_p4()[imu].pt() > 20. && passMuonSelection_VVV(imu, VVV_TIGHT_3L))
            {
                ntight++;
            }
        }
        // require 3 loose leptons
        if (nloose == 3)
            return abs(chargesum) == 1;
        else
            return false;
    }
    // If less than 2 leptons skip
    if (el_idx.size() + mu_idx.size() < 2) return false;
    // If equal to 2 leptons then must be a same-sign and at must have two loose lepton and one or more tight lepton
    int nloose = 0;
    int ntight = 0;
    for (auto& iel : coreElectron.index)
    {
        if (cms3.els_p4()[iel].pt() > 25. && passElectronSelection_VVV(iel, VVV_FO_SS))
            nloose++;
        if (cms3.els_p4()[iel].pt() > 25. && passElectronSelection_VVV(iel, VVV_TIGHT_SS))
            ntight++;
    }
    for (auto& imu : coreMuon.index)
    {
        if (cms3.mus_p4()[imu].pt() > 25. && passMuonSelection_VVV(imu, VVV_FO_SS))
            nloose++;
        if (cms3.mus_p4()[imu].pt() > 25. && passMuonSelection_VVV(imu, VVV_TIGHT_SS))
            ntight++;
    }
    // require two loose
    if (nloose != 2)
        return false;
    // require same charge
    if (mu_idx.size() == 2)
    {
        if (cms3.mus_charge()[mu_idx[0]] * cms3.mus_charge()[mu_idx[1]] > 0)
            return true;
        else
            return false;
    }
    else if (el_idx.size() == 2)
    {
        if (cms3.els_charge()[el_idx[0]] * cms3.els_charge()[el_idx[1]] > 0)
            return true;
        else
            return false;
    }
    else if (mu_idx.size() == 1 && el_idx.size() == 1)
    {
        if (cms3.mus_charge()[mu_idx[0]] * cms3.els_charge()[el_idx[0]] > 0)
            return true;
        else
            return false;
    }
    else
    {
        FATALERROR(__FUNCTION__);
        return false;
    }
}

//##############################################################################################################
void babyMaker_v2::FillTruthLevelStudyVariables()
{
    // WWW sample truth level study
    studyWWW();

    // Higgs pt dependent variable study
    studyHiggsDecay();

    // Study doubly charged higgs sample
    studyDoublyChargedHiggs();
}

//##############################################################################################################
void babyMaker_v2::FillEventInfo()
{
    tx->setBranch<Int_t>("run", cms3.evt_run());
    tx->setBranch<Int_t>("lumi", cms3.evt_lumiBlock());
    tx->setBranch<unsigned long long>("evt", cms3.evt_event());
    tx->setBranch<int>("isData", cms3.evt_isRealData());
    if (cms3.evt_isRealData())
    {
        tx->setBranch<float>("evt_scale1fb", 1);
        tx->setBranch<float>("xsec_br", 1);
        tx->setBranch<int>("evt_passgoodrunlist", goodrun(cms3.evt_run(), cms3.evt_lumiBlock()));
    }
    else
    {
        float scale1fb = 1;
        if (isDoublyChargedHiggs())
            scale1fb = 0.01;
        else if (isWprime())
            scale1fb = 0.004;
        else if (isWHSUSY())
            scale1fb = 1;
        else
            scale1fb = coreDatasetInfo.getScale1fb();

        float xsec = 1;
        if (isDoublyChargedHiggs())
            xsec = 1;
        else if (isWprime())
            xsec = 1;
        else if (isWHSUSY())
            xsec = 1;
        else
            xsec = coreDatasetInfo.getXsec();

        //==================OUTDATED OUTDATED OUTDATED OUTDATED======================================================================================
        // CMS3 www_2l_mia has scale1fb = 1/91900 * 1000. or 1/164800 * 1000.
        // CMS4 WWW samples are set to "SampleNiceName() == "www_2l_", so below snippet will be skipped
        /* deprecated */ //if (SampleNiceName().BeginsWith("www_2l_mia_cms3")     ) scale1fb *= 0.053842752 * 1.14082 *  91900. / (91900. + 164800.);
        /* deprecated */ //if (SampleNiceName().BeginsWith("www_2l_ext1_mia_cms3")) scale1fb *= 0.053842752 * 1.1402  * 164800. / (91900. + 164800.);
        //                                                                                ^^^^^^^^^^^   ^^^^^^^   ^^^^^^^   ^^^^^^^^^^^^^^^^^^
        //                                                                                 sigma*BR     neg-wgt   tot-evt   tot-evt of both sample
        // sigma(pp -> WWW) = 216 fb (MadGraph5)
        // BR (W -> e,mu,tau+neutrino) = 0.3258
        // BR (WWW->2l, l=e,mu,tau)) = 3 * (0.3258^2 * (1-0.3258)) + 0.3258^3 = 0.249272
        // sigma (pp -> WWW) * BR(WWW->2l) = 53.842752
        //
        //
        // 216/208.6 fb as the scale1fb data store is 208.6 in CORE.
        // But do this only for 2017 production.
        // This is because the 2016 production did not have this factor applied but rather applied at the looper level.
        // So for posterity (to produces same results with the tagged code) do not apply this for 2016
        //if (SampleNiceName().BeginsWith("www_") && coreSample.is2017(looper.getCurrentFileName()))
        //    scale1fb *= 1.035475;
        //==================OUTDATED OUTDATED OUTDATED OUTDATED======================================================================================

        // Sigma of WWW process for non-WH WWW process = 0.216 pb
        // BR (WWW->2l, l=e,mu,tau)) = 3 * (0.3258^2 * (1-0.3258)) + 0.3258^3 = 0.249272
        // Sigma * BR = 0.053842752
        // /PrivateWWW/www-cms4-Private80X-v1/MINIAODSIM        CMS4_V00-00-02_2017Sep27    92200      80870   0.053842752  0.0002439922782
        // /PrivateWWW/wwwext-cms4-Private80X-v1/MINIAODSIM     CMS4_V00-00-02_2017Sep27    159500     139804  0.053842752  0.0002439922782
        // Total events 80870 + 139804 = 220674
        // Per event weight = 0.053842752 * 1000 / 220674 = 0.0002439922782



        tx->setBranch<float>("evt_scale1fb", scale1fb);
        tx->setBranch<float>("xsec_br", xsec);
        tx->setBranch<int>("evt_passgoodrunlist", true);
    }
    tx->setBranch<TString>("CMS4path", looper.getCurrentFileName());
    tx->setBranch<int>("CMS4index", looper.getCurrentEventIndex());
}

//##############################################################################################################
void babyMaker_v2::FillElectrons()
{
    for (auto& idx : coreElectron.index)
    {
        float conecorrptfactorraw = coreElectron.index.size() + coreMuon.index.size() > 2 ? eleRelIso03EA(idx, 2, true) - 0.03: eleRelIso03EA(idx, 2, true) - 0.05;
        float conecorrptfactor = max(0., (double) conecorrptfactorraw) + 1.; // To clip correcting once it passes tight isolation criteria

        if (babymode != kFRBaby)
        {
            if (coreElectron.index.size() + coreMuon.index.size() == 2)
            {
                if (!( passElectronSelection_VVV(idx, VVV_FO_SS) ))
                    continue;
            }
            else if (coreElectron.index.size() + coreMuon.index.size() > 2)
            {
                if (!( passElectronSelection_VVV(idx, VVV_FO_3L) ))
                    continue;
            }
        }

        tx->pushbackToBranch<LorentzVector> ("lep_p4"                           , cms3.els_p4()[idx]);
        tx->pushbackToBranch<float>         ("lep_pt"                           , cms3.els_p4()[idx].pt());
        tx->pushbackToBranch<float>         ("lep_eta"                          , cms3.els_p4()[idx].eta());
        tx->pushbackToBranch<float>         ("lep_phi"                          , cms3.els_p4()[idx].phi());
        tx->pushbackToBranch<float>         ("lep_coneCorrPt"                   , cms3.els_p4()[idx].pt() * conecorrptfactor);
        tx->pushbackToBranch<float>         ("lep_ip3d"                         , cms3.els_ip3d()[idx]);
        tx->pushbackToBranch<float>         ("lep_ip3derr"                      , cms3.els_ip3derr()[idx]);
        tx->pushbackToBranch<int>           ("lep_isTriggerSafe_v1"             , isTriggerSafe_v1(idx));
        tx->pushbackToBranch<int>           ("lep_lostHits"                     , cms3.els_lostHits()[idx]);
        tx->pushbackToBranch<int>           ("lep_convVeto"                     , !cms3.els_conv_vtx_flag()[idx]);
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_3l_fo"      , cms3.els_p4()[idx].pt() > 20. && passElectronSelection_VVV(idx, VVV_FO_3L));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_3l_tight"   , cms3.els_p4()[idx].pt() > 20. && passElectronSelection_VVV(idx, VVV_TIGHT_3L));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_fo"         , cms3.els_p4()[idx].pt() > 25. && passElectronSelection_VVV(idx, VVV_FO_SS));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_tight"      , cms3.els_p4()[idx].pt() > 25. && passElectronSelection_VVV(idx, VVV_TIGHT_SS));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_veto"       , cms3.els_p4()[idx].pt() > 10. && passElectronSelection_VVV(idx, VVV_VETO));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_fo_noiso"   , cms3.els_p4()[idx].pt() > 20. && passElectronSelection_VVV(idx, VVV_FO_NOISO));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_tight_noiso", cms3.els_p4()[idx].pt() > 25. && passElectronSelection_VVV(idx, VVV_TIGHT_NOISO));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_veto_noiso" , cms3.els_p4()[idx].pt() > 10. && passElectronSelection_VVV(idx, VVV_VETO_NOISO));
        tx->pushbackToBranch<int>           ("lep_pass_POG_veto"                , cms3.els_p4()[idx].pt() > 10. && is2017POGVetoElectron(idx));
        tx->pushbackToBranch<int>           ("lep_pass_POG_loose"               , cms3.els_p4()[idx].pt() > 20. && is2017POGLooseElectron(idx));
        tx->pushbackToBranch<int>           ("lep_pass_POG_medium"              , cms3.els_p4()[idx].pt() > 20. && is2017POGMediumElectron(idx));
        tx->pushbackToBranch<int>           ("lep_pass_POG_tight"               , cms3.els_p4()[idx].pt() > 20. && is2017POGTightElectron(idx));
        tx->pushbackToBranch<int>           ("lep_pdgId"                        , cms3.els_charge()[idx]*(-11));
        tx->pushbackToBranch<float>         ("lep_dxy"                          , cms3.els_dxyPV()[idx]);
        tx->pushbackToBranch<float>         ("lep_dz"                           , cms3.els_dzPV()[idx]);
        tx->pushbackToBranch<float>         ("lep_pterr"                        , cms3.els_ptErr()[idx]);
        tx->pushbackToBranch<float>         ("lep_relIso04DB"                   , -999);
        tx->pushbackToBranch<float>         ("lep_relIso03EA"                   , eleRelIso03EA(idx, 4));
        tx->pushbackToBranch<float>         ("lep_relIso03EALep"                , eleRelIso03EA(idx, 4, true));
        tx->pushbackToBranch<float>         ("lep_relIso03EAv2"                 , eleRelIso03EA(idx, 2));
        tx->pushbackToBranch<float>         ("lep_relIso03EAv2Lep"              , eleRelIso03EA(idx, 2, true));
        tx->pushbackToBranch<int>           ("lep_tightCharge"                  , tightChargeEle(idx));
        tx->pushbackToBranch<float>         ("lep_trk_pt"                       , cms3.els_trk_p4()[idx].pt());
        tx->pushbackToBranch<int>           ("lep_charge"                       , cms3.els_charge()[idx]);
        tx->pushbackToBranch<float>         ("lep_etaSC"                        , cms3.els_etaSC()[idx]);
        tx->pushbackToBranch<float>         ("lep_MVA"                          , getMVAoutput(idx, true));
        tx->pushbackToBranch<int>           ("lep_isMediumPOG"                  , 1);
        tx->pushbackToBranch<int>           ("lep_isTightPOG"                   , 1);
        tx->pushbackToBranch<float>         ("lep_r9"                           , cms3.els_r9()[idx]);
        tx->pushbackToBranch<int>           ("lep_nlayers"                      , -1);
        if (!cms3.evt_isRealData())
        {
            pair<int, int> motherId_genIdx = lepMotherID_v2(Lep(cms3.els_charge()[idx] * (-11), idx)); //don't forget the sign
            tx->pushbackToBranch<int>       ("lep_motherIdSS"                   , motherId_genIdx.first);
            tx->pushbackToBranch<int>       ("lep_genPart_index"                , motherId_genIdx.second);
            tx->pushbackToBranch<int>       ("lep_isFromW"                      , isFromW(11, idx));
            tx->pushbackToBranch<int>       ("lep_isFromZ"                      , isFromZ(11, idx));
            tx->pushbackToBranch<int>       ("lep_isFromB"                      , isFromB(11, idx));
            tx->pushbackToBranch<int>       ("lep_isFromC"                      , isFromC(11, idx));
            tx->pushbackToBranch<int>       ("lep_isFromL"                      , isFromLight(11, idx));
            tx->pushbackToBranch<int>       ("lep_isFromLF"                     , isFromLightFake(11, idx));
        }
        else
        {
            tx->pushbackToBranch<int>       ("lep_motherIdSS"                   , 99);
            tx->pushbackToBranch<int>       ("lep_genPart_index"                , -1);
            tx->pushbackToBranch<int>       ("lep_isFromW"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromZ"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromB"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromC"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromL"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromLF"                     , -1);
        }

        // Due to CMS4 not having pf candidates...
        std::cout.setstate(std::ios_base::failbit); // To suppress warning about CMS4 not having PF candidates
        if (cms3.evt_CMS3tag()[0].Contains("CMS3"))
            tx->pushbackToBranch<float>     ("lep_relIso04EAv2"                 , elRelIsoCustomCone(idx, 0.4, false, 0.0, /*useDBCorr=*/false, /*useEACorr=*/true, /*mindr=*/ -1, /*eaversion=*/2));
        else
            tx->pushbackToBranch<float>     ("lep_relIso04EAv2"                 , eleRelIso03EA(idx, 2));
        std::cout.clear();
    }
}

//##############################################################################################################
void babyMaker_v2::FillMuons()
{
    for (auto& idx : coreMuon.index)
    {
        float conecorrptfactorraw = coreElectron.index.size() + coreMuon.index.size() > 2 ? muRelIso03EA(idx, 2, true) - 0.03: muRelIso03EA(idx, 2, true) - 0.07;
        float conecorrptfactor = max(0., (double) conecorrptfactorraw) + 1.; // To clip correcting once it passes tight isolation criteria

        if (babymode != kFRBaby)
        {
            if (coreElectron.index.size() + coreMuon.index.size() == 2)
            {
                if (!( passMuonSelection_VVV(idx, VVV_FO_SS) ))
                    continue;
            }
            else if (coreElectron.index.size() + coreMuon.index.size() > 2)
            {
                if (!( passMuonSelection_VVV(idx, VVV_FO_3L) ))
                    continue;
            }
        }

        tx->pushbackToBranch<LorentzVector> ("lep_p4"                           , cms3.mus_p4()[idx]);
        tx->pushbackToBranch<float>         ("lep_pt"                           , cms3.mus_p4()[idx].pt());
        tx->pushbackToBranch<float>         ("lep_eta"                          , cms3.mus_p4()[idx].eta());
        tx->pushbackToBranch<float>         ("lep_phi"                          , cms3.mus_p4()[idx].phi());
        tx->pushbackToBranch<float>         ("lep_coneCorrPt"                   , cms3.mus_p4()[idx].pt() * conecorrptfactor);
        tx->pushbackToBranch<float>         ("lep_ip3d"                         , cms3.mus_ip3d()[idx]);
        tx->pushbackToBranch<float>         ("lep_ip3derr"                      , cms3.mus_ip3derr()[idx]);
        tx->pushbackToBranch<int>           ("lep_isTriggerSafe_v1"             , true); // Electron specific branch. So muons always pass.
        tx->pushbackToBranch<int>           ("lep_lostHits"                     , cms3.mus_lostHits()[idx]);
        tx->pushbackToBranch<int>           ("lep_convVeto"                     , 1);
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_3l_fo"      , cms3.mus_p4()[idx].pt() > 20. && passMuonSelection_VVV(idx, VVV_FO_3L));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_3l_tight"   , cms3.mus_p4()[idx].pt() > 20. && passMuonSelection_VVV(idx, VVV_TIGHT_3L));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_fo"         , cms3.mus_p4()[idx].pt() > 25. && passMuonSelection_VVV(idx, VVV_FO_SS));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_tight"      , cms3.mus_p4()[idx].pt() > 25. && passMuonSelection_VVV(idx, VVV_TIGHT_SS));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_veto"       , cms3.mus_p4()[idx].pt() > 10. && passMuonSelection_VVV(idx, VVV_VETO));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_fo_noiso"   , cms3.mus_p4()[idx].pt() > 20. && passMuonSelection_VVV(idx, VVV_FO_NOISO));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_tight_noiso", cms3.mus_p4()[idx].pt() > 25. && passMuonSelection_VVV(idx, VVV_TIGHT_NOISO));
        tx->pushbackToBranch<int>           ("lep_pass_VVV_cutbased_veto_noiso" , cms3.mus_p4()[idx].pt() > 10. && passMuonSelection_VVV(idx, VVV_VETO_NOISO));
        tx->pushbackToBranch<int>           ("lep_pass_POG_veto"                , cms3.mus_p4()[idx].pt() > 10. && is2017POGVetoMuon(idx));
        tx->pushbackToBranch<int>           ("lep_pass_POG_loose"               , cms3.mus_p4()[idx].pt() > 20. && is2017POGLooseMuon(idx));
        tx->pushbackToBranch<int>           ("lep_pass_POG_medium"              , cms3.mus_p4()[idx].pt() > 20. && is2017POGMediumMuon(idx));
        tx->pushbackToBranch<int>           ("lep_pass_POG_tight"               , cms3.mus_p4()[idx].pt() > 20. && is2017POGTightMuon(idx));
        tx->pushbackToBranch<int>           ("lep_pdgId"                        , cms3.mus_charge()[idx]*(-13));
        tx->pushbackToBranch<float>         ("lep_dxy"                          , cms3.mus_dxyPV()[idx]);
        tx->pushbackToBranch<float>         ("lep_dz"                           , cms3.mus_dzPV()[idx]);
        tx->pushbackToBranch<float>         ("lep_pterr"                        , cms3.mus_ptErr()[idx]);
        tx->pushbackToBranch<float>         ("lep_relIso04DB"                   , muRelIso04DB(idx));
        tx->pushbackToBranch<float>         ("lep_relIso03EA"                   , muRelIso03EA(idx, 4));
        tx->pushbackToBranch<float>         ("lep_relIso03EALep"                , muRelIso03EA(idx, 4, true));
        tx->pushbackToBranch<float>         ("lep_relIso03EAv2"                 , muRelIso03EA(idx, 2));
        tx->pushbackToBranch<float>         ("lep_relIso03EAv2Lep"              , muRelIso03EA(idx, 2, true));
        tx->pushbackToBranch<int>           ("lep_tightCharge"                  , tightChargeMuon(idx));
        tx->pushbackToBranch<float>         ("lep_trk_pt"                       , cms3.mus_trk_p4()[idx].pt());
        tx->pushbackToBranch<int>           ("lep_charge"                       , cms3.mus_charge()[idx]);
        tx->pushbackToBranch<float>         ("lep_etaSC"                        , cms3.mus_p4()[idx].eta()); // Electron specific branch. Just take muon's regular eta.
        tx->pushbackToBranch<float>         ("lep_MVA"                          , -99);
        tx->pushbackToBranch<int>           ("lep_isMediumPOG"                  , isMediumMuonPOG(idx));
        tx->pushbackToBranch<int>           ("lep_isTightPOG"                   , isTightMuonPOG(idx));
        tx->pushbackToBranch<float>         ("lep_r9"                           , 0);
        tx->pushbackToBranch<int>           ("lep_nlayers"                      , cms3.mus_nlayers()[idx]);
        if (!cms3.evt_isRealData())
        {
            pair<int, int> motherId_genIdx = lepMotherID_v2(Lep(cms3.mus_charge()[idx] * (-13), idx)); //don't forget the sign
            tx->pushbackToBranch<int>       ("lep_motherIdSS"                   , motherId_genIdx.first);
            tx->pushbackToBranch<int>       ("lep_genPart_index"                , motherId_genIdx.second);
            tx->pushbackToBranch<int>       ("lep_isFromW"                      , isFromW(13, idx));
            tx->pushbackToBranch<int>       ("lep_isFromZ"                      , isFromZ(13, idx));
            tx->pushbackToBranch<int>       ("lep_isFromB"                      , isFromB(13, idx));
            tx->pushbackToBranch<int>       ("lep_isFromC"                      , isFromC(13, idx));
            tx->pushbackToBranch<int>       ("lep_isFromL"                      , isFromLight(13, idx));
            tx->pushbackToBranch<int>       ("lep_isFromLF"                     , isFromLightFake(13, idx));
        }
        else
        {
            tx->pushbackToBranch<int>       ("lep_motherIdSS"                   , 99);
            tx->pushbackToBranch<int>       ("lep_genPart_index"                , -1);
            tx->pushbackToBranch<int>       ("lep_isFromW"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromZ"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromB"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromC"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromL"                      , -1);
            tx->pushbackToBranch<int>       ("lep_isFromLF"                     , -1);
        }

        // Due to CMS4 not having pf candidates...
        std::cout.setstate(std::ios_base::failbit); // To suppress warning about CMS4 not having PF candidates
        if (cms3.evt_CMS3tag()[0].Contains("CMS3"))
            tx->pushbackToBranch<float>     ("lep_relIso04EAv2"                 , muRelIsoCustomCone(idx, 0.4, /*useVetoCones=*/false, 0.5, false, true, -1, 2));
        else
            tx->pushbackToBranch<float>     ("lep_relIso04EAv2"                 , muRelIso03EA(idx, 2));
        std::cout.clear();
    }
}

//##############################################################################################################
void babyMaker_v2::FillGenParticles()
{
    if (cms3.evt_isRealData())
        return;
    tx->setBranch<float>("gen_ht", coreGenPart.gen_ht);
    tx->setBranch<vector<LorentzVector>>("genPart_p4", coreGenPart.genPart_p4);
    tx->setBranch<vector<int>>("genPart_motherId", coreGenPart.genPart_motherId);
    tx->setBranch<vector<int>>("genPart_pdgId", coreGenPart.genPart_pdgId);
    tx->setBranch<vector<int>>("genPart_charge", coreGenPart.genPart_charge);
    tx->setBranch<vector<int>>("genPart_status", coreGenPart.genPart_status);
    tx->setBranch<int>("ngenLep", coreGenPart.ngenLep);
    tx->setBranch<int>("ngenLepFromTau", coreGenPart.ngenLepFromTau);
    tx->setBranch<int>("ngenLepFromBoson", coreGenPart.ngenLepFromBoson);
}

//##############################################################################################################
void babyMaker_v2::SortLeptonBranches()
{
    tx->sortVecBranchesByPt("lep_p4",
            {
            "lep_pt",
            "lep_eta",
            "lep_phi",
            "lep_dxy",
            "lep_dz",
            "lep_ip3d",
            "lep_ip3derr",
            "lep_pterr",
            "lep_relIso04DB",
            "lep_relIso03EA",
            "lep_relIso03EALep",
            "lep_relIso03EAv2",
            "lep_relIso03EAv2Lep",
            "lep_relIso04EAv2",
            "lep_trk_pt",
            "lep_etaSC",
            "lep_MVA",
            "lep_coneCorrPt",
            "lep_r9",
            },
            {
            "lep_pass_VVV_cutbased_3l_fo",
            "lep_pass_VVV_cutbased_3l_tight",
            "lep_pass_VVV_cutbased_fo",
            "lep_pass_VVV_cutbased_tight",
            "lep_pass_VVV_cutbased_veto",
            "lep_pass_VVV_cutbased_fo_noiso",
            "lep_pass_VVV_cutbased_tight_noiso",
            "lep_pass_VVV_cutbased_veto_noiso",
            "lep_pass_POG_veto",
            "lep_pass_POG_loose",
            "lep_pass_POG_medium",
            "lep_pass_POG_tight",
            "lep_isTriggerSafe_v1",
            "lep_lostHits",
            "lep_convVeto",
            "lep_motherIdSS",
            "lep_genPart_index",
            "lep_pdgId",
            "lep_tightCharge",
            "lep_charge",
            "lep_isMediumPOG",
            "lep_isTightPOG",
            "lep_isFromW",
            "lep_isFromZ",
            "lep_isFromB",
            "lep_isFromC",
            "lep_isFromL",
            "lep_isFromLF",
            "lep_nlayers",
            },
            {});
}

//##############################################################################################################
void babyMaker_v2::FillJets()
{
    coreBtagSF.clearSF();
    coreBtagSFFastSim.clearSF();
    for (unsigned ijet = 0; ijet < coreJet.index.size(); ++ijet)
    {
        int idx = coreJet.index[ijet];
        float corr = coreJet.corrs[ijet];
        float shift = coreJet.shifts[ijet];
        float current_csv_val = cms3.getbtagvalue("pfCombinedInclusiveSecondaryVertexV2BJetTags", idx);

        // Check whether this jet overlaps with any of the leptons
        if (isLeptonOverlappingWithJet(ijet))
            continue;

        LorentzVector jet = cms3.pfjets_p4()[idx] * cms3.pfjets_undoJEC()[idx] * corr;

        // Smear value for JER (only MC for data set it to nominal)
        float smear = cms3.evt_isRealData() ? 1 : coreJet.smears[ijet];
        float smearup = cms3.evt_isRealData() ? 1 : coreJet.smears_up[ijet];
        float smeardn = cms3.evt_isRealData() ? 1 : coreJet.smears_dn[ijet];

        if (jet.pt() > 20)
        {
            tx->pushbackToBranch<LorentzVector>("jets_p4", jet);
            tx->pushbackToBranch<float>("jets_csv", current_csv_val);
            if (!isWHSUSY())
                coreBtagSF.accumulateSF(idx, jet.pt(), jet.eta());
            else
                coreBtagSFFastSim.accumulateSF(idx, jet.pt(), jet.eta());
            if (jet.pt() > 30. and abs(jet.eta()) < 2.5)
                tx->pushbackToBranch<LorentzVector>("jets30_p4", jet);
        }

        LorentzVector jet_up = jet * (1. + shift);
        if (jet_up.pt() > 20)
        {
            tx->pushbackToBranch<LorentzVector>("jets_up_p4", jet_up);
            tx->pushbackToBranch<float>("jets_up_csv", current_csv_val);
            if (jet_up.pt() > 30. and abs(jet_up.eta()) < 2.5)
                tx->pushbackToBranch<LorentzVector>("jets30_up_p4", jet_up);
        }

        LorentzVector jet_dn = jet * (1. - shift);
        if (jet_dn.pt() > 20)
        {
            tx->pushbackToBranch<LorentzVector>("jets_dn_p4", jet_dn);
            tx->pushbackToBranch<float>("jets_dn_csv", current_csv_val);
            if (jet_dn.pt() > 30. and abs(jet_dn.eta()) < 2.5)
                tx->pushbackToBranch<LorentzVector>("jets30_dn_p4", jet_dn);
        }

        LorentzVector jet_jer = jet * smear;
        if (jet_jer.pt() > 20)
        {
            tx->pushbackToBranch<LorentzVector>("jets_jer_p4", jet_jer);
            tx->pushbackToBranch<float>("jets_jer_csv", current_csv_val);
            if (jet_jer.pt() > 30. and abs(jet_jer.eta()) < 2.5)
                tx->pushbackToBranch<LorentzVector>("jets30_jer_p4", jet_jer);
        }

        LorentzVector jet_jerup = jet * smearup;
        if (jet_jerup.pt() > 20)
        {
            tx->pushbackToBranch<LorentzVector>("jets_jerup_p4", jet_jerup);
            tx->pushbackToBranch<float>("jets_jerup_csv", current_csv_val);
            if (jet_jerup.pt() > 30. and abs(jet_jerup.eta()) < 2.5)
                tx->pushbackToBranch<LorentzVector>("jets30_jerup_p4", jet_jerup);
        }

        LorentzVector jet_jerdn = jet * smeardn;
        if (jet_jerdn.pt() > 20)
        {
            tx->pushbackToBranch<LorentzVector>("jets_jerdn_p4", jet_jerdn);
            tx->pushbackToBranch<float>("jets_jerdn_csv", current_csv_val);
            if (jet_jerdn.pt() > 30. and abs(jet_jerdn.eta()) < 2.5)
                tx->pushbackToBranch<LorentzVector>("jets30_jerdn_p4", jet_jerdn);
        }

    }

    if (!cms3.evt_isRealData())
    {
        if (!isWHSUSY())
        {
            tx->setBranch<float>("weight_btagsf"         , coreBtagSF.btagprob_data     / coreBtagSF.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_heavy_DN", coreBtagSF.btagprob_heavy_DN / coreBtagSF.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_heavy_UP", coreBtagSF.btagprob_heavy_UP / coreBtagSF.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_light_DN", coreBtagSF.btagprob_light_DN / coreBtagSF.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_light_UP", coreBtagSF.btagprob_light_UP / coreBtagSF.btagprob_mc);
        }
        else
        {
            tx->setBranch<float>("weight_btagsf"         , coreBtagSFFastSim.btagprob_data     / coreBtagSFFastSim.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_heavy_DN", coreBtagSFFastSim.btagprob_heavy_DN / coreBtagSFFastSim.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_heavy_UP", coreBtagSFFastSim.btagprob_heavy_UP / coreBtagSFFastSim.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_light_DN", coreBtagSFFastSim.btagprob_light_DN / coreBtagSFFastSim.btagprob_mc);
            tx->setBranch<float>("weight_btagsf_light_UP", coreBtagSFFastSim.btagprob_light_UP / coreBtagSFFastSim.btagprob_mc);
        }
    }
    else
    {
        tx->setBranch<float>("weight_btagsf"         , 1);
        tx->setBranch<float>("weight_btagsf_heavy_DN", 1);
        tx->setBranch<float>("weight_btagsf_heavy_UP", 1);
        tx->setBranch<float>("weight_btagsf_light_DN", 1);
        tx->setBranch<float>("weight_btagsf_light_UP", 1);
    }

}

//##############################################################################################################
void babyMaker_v2::FillFatJets()
{
    for (unsigned ijet = 0; ijet < coreFatJet.index.size(); ++ijet)
    {
        int idx = coreFatJet.index[ijet];
        float corr = coreFatJet.corrs[ijet];
        float shift = coreFatJet.shifts[ijet];
        LorentzVector fatjet = cms3.ak8jets_p4()[idx] * cms3.ak8jets_undoJEC()[idx] * corr;
        if (fatjet.pt() > 20)
        {
            tx->pushbackToBranch<LorentzVector>("ak8jets_p4", fatjet);
            // 2016 CMS4 has a different naming scheme
            if (coreSample.is2016(looper.getCurrentFileName()))
            {
                tx->pushbackToBranch<float>("ak8jets_softdropMass", cms3.ak8jets_softdropMass()[idx]);
                tx->pushbackToBranch<float>("ak8jets_prunedMass", cms3.ak8jets_prunedMass()[idx]);
                tx->pushbackToBranch<float>("ak8jets_nJettinessTau1", cms3.ak8jets_nJettinessTau1()[idx]);
                tx->pushbackToBranch<float>("ak8jets_nJettinessTau2", cms3.ak8jets_nJettinessTau2()[idx]);
                tx->pushbackToBranch<float>("ak8jets_mass", cms3.ak8jets_mass()[idx]);
                tx->pushbackToBranch<float>("ak8jets_puppi_nJettinessTau1", cms3.ak8jets_puppi_nJettinessTau1()[idx]);
                tx->pushbackToBranch<float>("ak8jets_puppi_nJettinessTau2", cms3.ak8jets_puppi_nJettinessTau2()[idx]);
                tx->pushbackToBranch<float>("ak8jets_puppi_softdropMass", cms3.ak8jets_puppi_softdropMass()[idx]);
                tx->pushbackToBranch<float>("ak8jets_puppi_eta", cms3.ak8jets_puppi_eta()[idx]);
                tx->pushbackToBranch<float>("ak8jets_puppi_phi", cms3.ak8jets_puppi_phi()[idx]);
                tx->pushbackToBranch<float>("ak8jets_puppi_pt", cms3.ak8jets_puppi_pt()[idx]);
                tx->pushbackToBranch<float>("ak8jets_puppi_mass", cms3.ak8jets_puppi_mass()[idx]);
                tx->pushbackToBranch<float>("ak8jets_softdropPuppiSubjet1", cms3.ak8jets_softdropPuppiSubjet1()[idx].M());
                tx->pushbackToBranch<float>("ak8jets_softdropPuppiSubjet2", cms3.ak8jets_softdropPuppiSubjet2()[idx].M());
            }
            // 2017 uses PUPPI for default and several of these mass variables are not even provided
            else if (coreSample.is2017(looper.getCurrentFileName()))
            {
                tx->pushbackToBranch<float>("ak8jets_softdropMass", cms3.ak8jets_puppi_softdropMass()[idx]);
                tx->pushbackToBranch<float>("ak8jets_softdropPuppiSubjet1", cms3.ak8jets_softdropPuppiSubjet1()[idx].M());
                tx->pushbackToBranch<float>("ak8jets_softdropPuppiSubjet2", cms3.ak8jets_softdropPuppiSubjet2()[idx].M());
                tx->pushbackToBranch<float>("ak8jets_nJettinessTau1", cms3.ak8jets_nJettinessTau1()[idx]);
                tx->pushbackToBranch<float>("ak8jets_nJettinessTau2", cms3.ak8jets_nJettinessTau2()[idx]);
                tx->pushbackToBranch<float>("ak8jets_deep_rawdisc_h4q", cms3.ak8jets_deep_rawdisc_h4q()[idx]);
                tx->pushbackToBranch<float>("ak8jets_deep_rawdisc_hbb", cms3.ak8jets_deep_rawdisc_hbb()[idx]);
                tx->pushbackToBranch<float>("ak8jets_deep_rawdisc_top", cms3.ak8jets_deep_rawdisc_top()[idx]);
                tx->pushbackToBranch<float>("ak8jets_deep_rawdisc_qcd", cms3.ak8jets_deep_rawdisc_qcd()[idx]);
                tx->pushbackToBranch<float>("ak8jets_deep_rawdisc_w", cms3.ak8jets_deep_rawdisc_w()[idx]);
                tx->pushbackToBranch<float>("ak8jets_deep_rawdisc_z", cms3.ak8jets_deep_rawdisc_z()[idx]);
                tx->pushbackToBranch<float>("ak8jets_deep_rawdisc_zbb", cms3.ak8jets_deep_rawdisc_zbb()[idx]);
            }


            // Some how cannot find the following in CMS4 2016... unfortunate
//            tx->pushbackToBranch<float>("ak8jets_trimmedMass", cms3.ak8jets_trimmedMass()[idx]);
        }
    }
}

//##############################################################################################################
void babyMaker_v2::FillMET()
{
    tx->setBranch<float>("met_pt", coreMET.met_pt);
    tx->setBranch<float>("met_phi", coreMET.met_phi);
    tx->setBranch<float>("met_up_pt", coreMET.met_up_pt);
    tx->setBranch<float>("met_up_phi", coreMET.met_up_phi);
    tx->setBranch<float>("met_dn_pt", coreMET.met_dn_pt);
    tx->setBranch<float>("met_dn_phi", coreMET.met_dn_phi);
    tx->setBranch<float>("met_gen_pt", cms3.gen_met());
    tx->setBranch<float>("met_gen_phi", cms3.gen_metPhi());
    tx->setBranch<float>("met_jer_pt", coreMET.met_jer_pt);
    tx->setBranch<float>("met_jer_phi", coreMET.met_jer_phi);
    tx->setBranch<float>("met_jerup_pt", coreMET.met_jerup_pt);
    tx->setBranch<float>("met_jerup_phi", coreMET.met_jerup_phi);
    tx->setBranch<float>("met_jerdn_pt", coreMET.met_jerdn_pt);
    tx->setBranch<float>("met_jerdn_phi", coreMET.met_jerdn_phi);
}

//##############################################################################################################
void babyMaker_v2::FillTracks()
{
    int nisotrack = 0;
    for (auto& idx : coreTrack.index)
    {
        // Check whether this track overlaps with any of the veto leptons
        if (isLeptonOverlappingWithTrack(idx))
            continue;
        nisotrack++;
    }

    // Set the branch
    tx->setBranch<int>("nisoTrack_mt2_cleaned_VVV_cutbased_veto", nisotrack);
}

//##############################################################################################################
void babyMaker_v2::SortJetBranches()
{
    tx->sortVecBranchesByPt("jets_p4", {"jets_csv"}, {}, {});
    tx->sortVecBranchesByPt("jets_up_p4", {"jets_up_csv"}, {}, {});
    tx->sortVecBranchesByPt("jets_dn_p4", {"jets_dn_csv"}, {}, {});
    tx->sortVecBranchesByPt("jets_jer_p4", {"jets_jer_csv"}, {}, {});
    tx->sortVecBranchesByPt("jets_jerup_p4", {"jets_jerup_csv"}, {}, {});
    tx->sortVecBranchesByPt("jets_jerdn_p4", {"jets_jerdn_csv"}, {}, {});
    tx->sortVecBranchesByPt("jets30_p4", {}, {}, {});
    tx->sortVecBranchesByPt("jets30_up_p4", {}, {}, {});
    tx->sortVecBranchesByPt("jets30_dn_p4", {}, {}, {});
    tx->sortVecBranchesByPt("jets30_jer_p4", {}, {}, {});
}

//##############################################################################################################
void babyMaker_v2::FillTrigger()
{
    if (cms3.evt_isRealData())
    {
        tx->setBranch<int>("HLT_DoubleMu", coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_DoubleMu_2017 : coreTrigger.HLT_DoubleMu);
        tx->setBranch<int>("HLT_DoubleEl", coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_DoubleEl_2017 : coreTrigger.HLT_DoubleEl);
        tx->setBranch<int>("HLT_DoubleEl_DZ", coreTrigger.HLT_DoubleEl_DZ);
        tx->setBranch<int>("HLT_DoubleEl_DZ_2", coreTrigger.HLT_DoubleEl_DZ_2);
        tx->setBranch<int>("HLT_MuEG", coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_MuEG_2017 : coreTrigger.HLT_MuEG);
        tx->setBranch<int>("HLT_MuEG_2016", coreTrigger.HLT_MuEG);
        tx->setBranch<int>("HLT_SingleEl8"  , coreTrigger.HLT_SingleEl8_2017 );
        tx->setBranch<int>("HLT_SingleEl17" , coreTrigger.HLT_SingleEl17_2017 );
        tx->setBranch<int>("HLT_SingleIsoEl8"  , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoEl8_2017  : coreTrigger.HLT_SingleIsoEl8  );
        tx->setBranch<int>("HLT_SingleIsoEl17" , coreTrigger.HLT_SingleIsoEl17 );
        tx->setBranch<int>("HLT_SingleIsoEl23" , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoEl23_2017 : coreTrigger.HLT_SingleIsoEl23 );
        tx->setBranch<int>("HLT_SingleIsoMu8"  , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoMu8_2017  : coreTrigger.HLT_SingleIsoMu8  );
        tx->setBranch<int>("HLT_SingleIsoMu17" , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoMu17_2017 : coreTrigger.HLT_SingleIsoMu17 );
        tx->setBranch<int>("HLT_PFMET140_PFMHT140_IDTight", coreTrigger.HLT_PFMET140_PFMHT140_IDTight);
    }
    else
    {
        tx->setBranch<int>("HLT_DoubleMu", 1);
        tx->setBranch<int>("HLT_DoubleEl", 1);
        tx->setBranch<int>("HLT_DoubleEl_DZ", 1);
        tx->setBranch<int>("HLT_DoubleEl_DZ_2", 1);
        tx->setBranch<int>("HLT_MuEG", 1);
        tx->setBranch<int>("HLT_MuEG_2016", 1);
        tx->setBranch<int>("HLT_SingleIsoEl8", 1);
        tx->setBranch<int>("HLT_SingleIsoEl17", 1);
        tx->setBranch<int>("HLT_SingleIsoEl23", 1);
        tx->setBranch<int>("HLT_SingleIsoMu8", 1);
        tx->setBranch<int>("HLT_SingleIsoMu17", 1);
        tx->setBranch<int>("HLT_PFMET140_PFMHT140_IDTight", 1);
    }
    if (isWHSUSY())
    {
        tx->setBranch<int>("mc_HLT_DoubleMu", 1);
        tx->setBranch<int>("mc_HLT_DoubleEl", 1);
        tx->setBranch<int>("mc_HLT_DoubleEl_DZ", 1);
        tx->setBranch<int>("mc_HLT_DoubleEl_DZ_2", 1);
        tx->setBranch<int>("mc_HLT_MuEG", 1);
        tx->setBranch<int>("mc_HLT_MuEG_2016", 1);
        tx->setBranch<int>("mc_HLT_SingleIsoEl8", 1);
        tx->setBranch<int>("mc_HLT_SingleIsoEl17", 1);
        tx->setBranch<int>("mc_HLT_SingleIsoEl23", 1);
        tx->setBranch<int>("mc_HLT_SingleIsoMu8", 1);
        tx->setBranch<int>("mc_HLT_SingleIsoMu17", 1);
        tx->setBranch<int>("mc_HLT_PFMET140_PFMHT140_IDTight", 1);
    }
    else
    {
        tx->setBranch<int>("mc_HLT_DoubleMu", coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_DoubleMu_2017 : coreTrigger.HLT_DoubleMu);
        tx->setBranch<int>("mc_HLT_DoubleEl", coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_DoubleEl_2017 : coreTrigger.HLT_DoubleEl);
        tx->setBranch<int>("mc_HLT_DoubleEl", coreTrigger.HLT_DoubleEl);
        tx->setBranch<int>("mc_HLT_DoubleEl_DZ", coreTrigger.HLT_DoubleEl_DZ);
        tx->setBranch<int>("mc_HLT_DoubleEl_DZ_2", coreTrigger.HLT_DoubleEl_DZ_2);
        tx->setBranch<int>("mc_HLT_MuEG", coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_MuEG_2017 : coreTrigger.HLT_MuEG);
        tx->setBranch<int>("mc_HLT_MuEG_2016", coreTrigger.HLT_MuEG);
        tx->setBranch<int>("mc_HLT_SingleEl8"  , coreTrigger.HLT_SingleEl8_2017 );
        tx->setBranch<int>("mc_HLT_SingleEl17" , coreTrigger.HLT_SingleEl17_2017 );
        tx->setBranch<int>("mc_HLT_SingleIsoEl8"  , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoEl8_2017  : coreTrigger.HLT_SingleIsoEl8  );
        tx->setBranch<int>("mc_HLT_SingleIsoEl17" , coreTrigger.HLT_SingleIsoEl17 );
        tx->setBranch<int>("mc_HLT_SingleIsoEl23" , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoEl23_2017 : coreTrigger.HLT_SingleIsoEl23 );
        tx->setBranch<int>("mc_HLT_SingleIsoMu8"  , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoMu8_2017  : coreTrigger.HLT_SingleIsoMu8  );
        tx->setBranch<int>("mc_HLT_SingleIsoMu17" , coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_SingleIsoMu17_2017 : coreTrigger.HLT_SingleIsoMu17 );
        tx->setBranch<int>("mc_HLT_PFMET140_PFMHT140_IDTight", coreTrigger.HLT_PFMET140_PFMHT140_IDTight);
    }
    if (cms3.evt_isRealData())
    {
        bool trig_ee = coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_DoubleEl_2017 : coreTrigger.HLT_DoubleEl || coreTrigger.HLT_DoubleEl_DZ;
        bool trig_em = coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_MuEG_2017     : coreTrigger.HLT_MuEG;
        bool trig_mm = coreSample.is2017(looper.getCurrentFileName()) ? coreTrigger.HLT_DoubleMu_2017 : coreTrigger.HLT_DoubleMu;
        bool is_pd_ee = looper.getCurrentFileName().Contains("DoubleEG");
        bool is_pd_em = looper.getCurrentFileName().Contains("MuonEG");
        bool is_pd_mm = looper.getCurrentFileName().Contains("DoubleMuon");
        bool pass_duplicate_ee_em_mm = false;
        bool pass_duplicate_mm_em_ee = false;
        if (is_pd_ee)
        {
            if (trig_ee)
                pass_duplicate_ee_em_mm = true;
            if (!trig_mm && !trig_em && trig_ee)
                pass_duplicate_mm_em_ee = true;
        }
        else if (is_pd_em)
        {
            if (!trig_ee && trig_em)
                pass_duplicate_ee_em_mm = true;
            if (!trig_mm && trig_em)
                pass_duplicate_mm_em_ee = true;
        }
        else if (is_pd_mm)
        {
            if (!trig_ee && !trig_em && trig_mm)
                pass_duplicate_ee_em_mm = true;
            if (trig_mm)
                pass_duplicate_mm_em_ee = true;
        }
        tx->setBranch<int>("pass_duplicate_ee_em_mm", pass_duplicate_ee_em_mm);
        tx->setBranch<int>("pass_duplicate_mm_em_ee", pass_duplicate_mm_em_ee);

        // Until something weird is resolved
        bool trig_em_2016 = coreTrigger.HLT_MuEG;
        bool pass_duplicate_ee_em2016_mm = false;
        bool pass_duplicate_mm_em2016_ee = false;
        if (is_pd_ee)
        {
            if (trig_ee)
                pass_duplicate_ee_em2016_mm = true;
            if (!trig_mm && !trig_em_2016 && trig_ee)
                pass_duplicate_mm_em2016_ee = true;
        }
        else if (is_pd_em)
        {
            if (!trig_ee && trig_em_2016)
                pass_duplicate_ee_em2016_mm = true;
            if (!trig_mm && trig_em_2016)
                pass_duplicate_mm_em2016_ee = true;
        }
        else if (is_pd_mm)
        {
            if (!trig_ee && !trig_em_2016 && trig_mm)
                pass_duplicate_ee_em2016_mm = true;
            if (trig_mm)
                pass_duplicate_mm_em2016_ee = true;
        }
        tx->setBranch<int>("pass_duplicate_ee_em2016_mm", pass_duplicate_ee_em2016_mm);
        tx->setBranch<int>("pass_duplicate_mm_em2016_ee", pass_duplicate_mm_em2016_ee);
    }
    else
    {
        tx->setBranch<int>("pass_duplicate_ee_em_mm", 1);
        tx->setBranch<int>("pass_duplicate_mm_em_ee", 1);
    }
    tx->setBranch<int>("passTrigger", passTrigger());
}

//##############################################################################################################
bool babyMaker_v2::passTrigger()
{
    if (tx->getBranch<int>("nLlep", true) < 2)
        return false;

    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const int mc_HLT_DoubleEl    = tx->getBranch<int>("mc_HLT_DoubleEl", true);
    const int mc_HLT_DoubleEl_DZ = tx->getBranch<int>("mc_HLT_DoubleEl_DZ", true);
    const int mc_HLT_MuEG        = tx->getBranch<int>("mc_HLT_MuEG", true);
    const int mc_HLT_DoubleMu    = tx->getBranch<int>("mc_HLT_DoubleMu", true);
    const int nVlep              = tx->getBranch<int>("nVlep", true);
    const int nLlep              = tx->getBranch<int>("nLlep", true);


    if (nVlep != 2 && nVlep != 3)
        return 0;

    if (nLlep != 2 && nLlep != 3)
        return 0;

    if (lep_pdgId.size() < 2)
        return 0;

    if (nVlep == 2 && nLlep == 2)
    {
        int lepprod = lep_pdgId.at(0)*lep_pdgId.at(1);
        if (abs(lepprod) == 121)
            return (mc_HLT_DoubleEl || mc_HLT_DoubleEl_DZ);
        else if (abs(lepprod) == 143)
            return mc_HLT_MuEG;
        else if (abs(lepprod) == 169)
            return mc_HLT_DoubleMu;
        else
            return 0;
    }
    else if (nVlep == 3 && nLlep == 3)
    {
        int lepprod01 = lep_pdgId.at(0)*lep_pdgId.at(1);
        if (abs(lepprod01) == 121 && (mc_HLT_DoubleEl || mc_HLT_DoubleEl_DZ))
            return true;
        else if (abs(lepprod01) == 143 && mc_HLT_MuEG)
            return true;
        else if (abs(lepprod01) == 169 && mc_HLT_DoubleMu)
            return true;

        int lepprod02 = lep_pdgId.at(0)*lep_pdgId.at(2);
        if (abs(lepprod02) == 121 && (mc_HLT_DoubleEl || mc_HLT_DoubleEl_DZ))
            return true;
        else if (abs(lepprod02) == 143 && mc_HLT_MuEG)
            return true;
        else if (abs(lepprod02) == 169 && mc_HLT_DoubleMu)
            return true;

        int lepprod12 = lep_pdgId.at(1)*lep_pdgId.at(2);
        if (abs(lepprod12) == 121 && (mc_HLT_DoubleEl || mc_HLT_DoubleEl_DZ))
            return true;
        else if (abs(lepprod12) == 143 && mc_HLT_MuEG)
            return true;
        else if (abs(lepprod12) == 169 && mc_HLT_DoubleMu)
            return true;

        return false;
    }
    else
    {
        return 0;
    }
}

//##############################################################################################################
void babyMaker_v2::FillVertexInfo()
{
    tx->setBranch<int>("firstgoodvertex", firstGoodVertex());
    if (!cms3.evt_isRealData())
        tx->setBranch<int>("nTrueInt", cms3.puInfo_trueNumInteractions().at(0));
    else
        tx->setBranch<int>("nTrueInt", -999);
    int nVert = 0;
    for (unsigned int ivtx = 0; ivtx < cms3.evt_nvtxs(); ivtx++)
        if (isGoodVertex(ivtx)) nVert++;
    tx->setBranch<int>("nVert", nVert);
}

//##############################################################################################################
void babyMaker_v2::FillMETFilter()
{
    // get CMS3 version number to use later
    TString cms3_version = cms3.evt_CMS3tag()[0];

    // convert last two digits of version number to int
    int small_cms3_version = TString(cms3_version(cms3_version.Length() - 2, cms3_version.Length())).Atoi();

    // Compute whether this is a recent cms3 version
    bool recent_cms3_version = (cms3_version.Contains("V08-00") && small_cms3_version <= 12) ? false : true;

    // Is it CMS4?
    bool is_cms4 = cms3_version.Contains("CMS4");

    // If 2017
    if (coreSample.is2017(looper.getCurrentFileName()))
    {

        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
        // Using CMSSW 94X case https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2018
        tx->setBranch<int>("Flag_goodVertices", cms3.filt_goodVertices());
        tx->setBranch<int>("Flag_globalTightHalo2016", cms3.filt_globalTightHalo2016());
        tx->setBranch<int>("Flag_HBHENoiseFilter", cms3.filt_hbheNoise());
        tx->setBranch<int>("Flag_HBHEIsoNoiseFilter", cms3.filt_hbheNoiseIso());
        tx->setBranch<int>("Flag_EcalDeadCellTriggerPrimitiveFilter", cms3.filt_ecalTP());
        tx->setBranch<int>("Flag_badMuonFilter", cms3.filt_BadPFMuonFilter());
        tx->setBranch<int>("Flag_badChargedCandidateFilter", cms3.filt_BadChargedCandidateFilter());
        tx->setBranch<int>("Flag_eeBadScFilter", cms3.filt_eeBadSc());
        tx->setBranch<int>("Flag_ecalBadCalibFilter", cms3.filt_ecalBadCalibFilter());
        if (cms3.evt_isRealData())
        {
            tx->setBranch<int>("Flag_eeBadScFilter", 1); // Not suggested for FullSim
            if (isFastSim())
            {
                tx->setBranch<int>("Flag_globalTightHalo2016", 1);
            }
        }

        // Pass all event filters if it is a data event or just pass it if it is an MC
        tx->setBranch<int>("Flag_AllEventFilters",
                tx->getBranch<int>("Flag_goodVertices")
                && tx->getBranch<int>("Flag_globalTightHalo2016")
                && tx->getBranch<int>("Flag_HBHENoiseFilter")
                && tx->getBranch<int>("Flag_HBHEIsoNoiseFilter")
                && tx->getBranch<int>("Flag_EcalDeadCellTriggerPrimitiveFilter")
                && tx->getBranch<int>("Flag_badMuonFilter")
                && tx->getBranch<int>("Flag_badChargedCandidateFilter")
                && tx->getBranch<int>("Flag_eeBadScFilter")
                && tx->getBranch<int>("Flag_ecalBadCalibFilter")
                );
    }
    else // 2016
    {

        if (cms3.evt_isRealData())
        {
            tx->setBranch<int>("Flag_ecalLaserCorrFilter", cms3.filt_ecalLaser());
            tx->setBranch<int>("Flag_hcalLaserEventFilter", cms3.filt_hcalLaser());
            tx->setBranch<int>("Flag_trackingFailureFilter", cms3.filt_trackingFailure());
            tx->setBranch<int>("Flag_CSCTightHaloFilter", cms3.filt_cscBeamHalo());
        }
        else
        {
            tx->setBranch<int>("Flag_ecalLaserCorrFilter", 1);
            tx->setBranch<int>("Flag_hcalLaserEventFilter", 1);
            tx->setBranch<int>("Flag_trackingFailureFilter", 1);
            tx->setBranch<int>("Flag_CSCTightHaloFilter", 1);
        }

        // in data and MC
        if (!isWHSUSY())
        {
            tx->setBranch<int>("Flag_HBHENoiseFilter", cms3.filt_hbheNoise());
            tx->setBranch<int>("Flag_HBHEIsoNoiseFilter", cms3.filt_hbheNoiseIso());
            tx->setBranch<int>("Flag_CSCTightHalo2015Filter", cms3.filt_cscBeamHalo2015());
            tx->setBranch<int>("Flag_EcalDeadCellTriggerPrimitiveFilter", cms3.filt_ecalTP());
            tx->setBranch<int>("Flag_goodVertices", cms3.filt_goodVertices());
            tx->setBranch<int>("Flag_eeBadScFilter", cms3.filt_eeBadSc());
        }
        else
        {
            tx->setBranch<int>("Flag_HBHENoiseFilter", 1);
            tx->setBranch<int>("Flag_HBHEIsoNoiseFilter", 1);
            tx->setBranch<int>("Flag_CSCTightHalo2015Filter", 1);
            tx->setBranch<int>("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);
            tx->setBranch<int>("Flag_goodVertices", 1);
            tx->setBranch<int>("Flag_eeBadScFilter", 1);
        }

        if (!is_cms4)
            tx->setBranch<int>("Flag_badChargedCandidateFilter", badChargedCandidateFilter());
        else
            tx->setBranch<int>("Flag_badChargedCandidateFilter", 1);

        // inputs for badMuonFilters in latest cms3 tags
        if (recent_cms3_version)
        {
            if (!isWHSUSY())
            {
                tx->setBranch<int>("Flag_globalTightHalo2016", cms3.filt_globalTightHalo2016());
            }
            else
            {
                tx->setBranch<int>("Flag_globalTightHalo2016", 1);
            }
            if (!is_cms4)
                tx->setBranch<int>("Flag_badMuonFilter", badMuonFilter());
            else
                tx->setBranch<int>("Flag_badMuonFilter", 1);
            tx->setBranch<int>("Flag_badMuonFilterv2", badMuonFilterV2());
            tx->setBranch<int>("Flag_badChargedCandidateFilterv2", badChargedCandidateFilterV2());
            if (small_cms3_version >= 18 && !is_cms4)
            {
                tx->setBranch<int>("Flag_badMuons", cms3.filt_badMuons());
                tx->setBranch<int>("Flag_duplicateMuons", cms3.filt_duplicateMuons());
                tx->setBranch<int>("Flag_noBadMuons", cms3.filt_noBadMuons());
            }
            else
            {
                tx->setBranch<int>("Flag_badMuons", 1);
                tx->setBranch<int>("Flag_duplicateMuons", 1);
                tx->setBranch<int>("Flag_noBadMuons", 1);
            }
        }

        // Pass all event filters if it is a data event or just pass it if it is an MC
        tx->setBranch<int>("Flag_AllEventFilters",
                (
                 cms3.evt_isRealData() &&
                 tx->getBranch<int>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
                 tx->getBranch<int>("Flag_HBHEIsoNoiseFilter") &&
                 tx->getBranch<int>("Flag_HBHENoiseFilter") &&
                 tx->getBranch<int>("Flag_badChargedCandidateFilter") &&
                 tx->getBranch<int>("Flag_badMuonFilter") &&
                 tx->getBranch<int>("Flag_eeBadScFilter") &&
                 tx->getBranch<int>("Flag_globalTightHalo2016") &&
                 tx->getBranch<int>("Flag_goodVertices")
                )
                ||
                (
                 !cms3.evt_isRealData()
                )
                );

    }

    // fastsim filter for SMS sample
    if (isWHSUSY())
    {
        // FastSim filter, Nominal Jets
        bool fastsimfilt = false;
        for (unsigned ijet = 0; ijet < coreJet.index.size(); ++ijet)
        {
            int idx = coreJet.index[ijet];
            float corr = coreJet.corrs[ijet];
            LV p4 = cms3.pfjets_p4()[idx] * cms3.pfjets_undoJEC()[idx] * corr;
            if (p4.Pt() < 30) continue;
            if (fabs(p4.Eta()) > 5.0) continue;
            bool isgenmatch = false;
            for (unsigned int gix = 0; gix < cms3.genjets_p4NoMuNoNu().size(); ++gix)
            {
                if (ROOT::Math::VectorUtil::DeltaR(cms3.genjets_p4NoMuNoNu()[gix], p4) < 0.3)
                {
                    isgenmatch = true;
                    break;
                }
            }
            if (isgenmatch)
                continue;
            if (cms3.pfjets_chargedHadronE().at(idx)/ (cms3.pfjets_undoJEC().at(idx) * p4.energy()) > 0.1)
                continue;
            fastsimfilt = true;
            break;
        }
        tx->setBranch<int>("fastsimfilt", fastsimfilt);
    }
}

//##############################################################################################################
void babyMaker_v2::FillEMuLeptons()
{
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4");
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId");
    const vector<float>& lep_relIso04DB = tx->getBranch<vector<float>>("lep_relIso04DB");
    const vector<float>& lep_relIso03EA = tx->getBranch<vector<float>>("lep_relIso03EA");
    const vector<float>& lep_relIso03EALep = tx->getBranch<vector<float>>("lep_relIso03EALep");
    const vector<float>& lep_ip3d = tx->getBranch<vector<float>>("lep_ip3d");

    // TODO check for two leptons only ?

    for (unsigned i = 0; i < lep_p4.size(); ++i)
    {
        if (abs(lep_pdgId[i]) == 11)
        {
            tx->setBranch<float>("el_pt", lep_p4[i].pt());
            tx->setBranch<float>("el_eta", lep_p4[i].eta());
            tx->setBranch<float>("el_phi", lep_p4[i].phi());
            tx->setBranch<float>("el_relIso03EA", lep_relIso03EA[i]);
            tx->setBranch<float>("el_relIso03EALep", lep_relIso03EALep[i]);
            tx->setBranch<float>("el_ip3d", lep_ip3d[i]);
        }
        else if (abs(lep_pdgId[i]) == 13)
        {
            tx->setBranch<float>("mu_pt", lep_p4[i].pt());
            tx->setBranch<float>("mu_eta", lep_p4[i].eta());
            tx->setBranch<float>("mu_phi", lep_p4[i].phi());
            tx->setBranch<float>("mu_relIso04DB", lep_relIso04DB[i]);
            tx->setBranch<float>("mu_relIso03EA", lep_relIso03EA[i]);
            tx->setBranch<float>("mu_relIso03EALep", lep_relIso03EALep[i]);
            tx->setBranch<float>("mu_ip3d", lep_ip3d[i]);
        }
        else
        {
            FATALERROR(__FUNCTION__);
        }
    }
}

//##############################################################################################################
void babyMaker_v2::FillLbntLeptons()
{
    // Retrieve relevant variables
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId");
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4");
    const vector<float>& lep_coneCorrPt = tx->getBranch<vector<float>>("lep_coneCorrPt");

    // If the lepton counts reveal that it is not even an event in the AR then return 0
    const int& nVlep = tx->getBranch<int>("nVlep");
    const int& nLlep = tx->getBranch<int>("nLlep");
    const int& nTlep = tx->getBranch<int>("nTlep");
    if ((nLlep - 1) != nTlep)
        return;

    // Retrieve the flags for the lepton ID
    const vector<int>& istight = nVlep == 2 ? tx->getBranch<vector<int>>(gconf.wwwcfg["tightid"]) : tx->getBranch<vector<int>>(gconf.wwwcfg["3ltightid"]);
    const vector<int>& isloose = nVlep == 2 ? tx->getBranch<vector<int>>(gconf.wwwcfg["looseid"]) : tx->getBranch<vector<int>>(gconf.wwwcfg["3llooseid"]);

    // Figure out the index by whichever one is the loose but not tight
    int index = -1;
    for (int ilep = 0; ilep < nLlep; ++ilep)
    {
        if (!istight[ilep] && isloose[ilep])
            index = ilep;
    }

    // Sanity check
    if (index < 0)
    {
        std::cout <<  " passCount(tx->getBranch<vector<int>>(gconf.wwwcfg['3llooseid'])): " << passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["3llooseid"])) <<  std::endl;
        std::cout <<  " nVlep: " << nVlep <<  std::endl;
        vector<int> el_idx = coreElectron.index;
        vector<int> mu_idx = coreMuon.index;
        std::cout <<  " el_idx.size(): " << el_idx.size() <<  " mu_idx.size(): " << mu_idx.size() <<  std::endl;
        std::cout <<  " isloose.size(): " << isloose.size() <<  std::endl;
        for (int ilep = 0; ilep < nVlep; ++ilep)
            std::cout <<  " istight[ilep]: " << istight[ilep] <<  " isloose[ilep]: " << isloose[ilep] <<  std::endl;
        FATALERROR(__FUNCTION__);
    }

    tx->setBranch<float>("lbnt_coneCorrPt", tx->getBranch<vector<float>>("lep_coneCorrPt")[index]);
    tx->setBranch<float>("lbnt_pt", tx->getBranch<vector<float>>("lep_pt")[index]);
    tx->setBranch<float>("lbnt_abseta", fabs(tx->getBranch<vector<float>>("lep_eta")[index]));
    tx->setBranch<float>("lbnt_pdgId", tx->getBranch<vector<int>>("lep_pdgId")[index]);
    if (abs(tx->getBranch<vector<int>>("lep_pdgId")[index]) == 11)
    {
        tx->setBranch<float>("lbnt_el_coneCorrPt", tx->getBranch<vector<float>>("lep_coneCorrPt")[index]);
        tx->setBranch<float>("lbnt_el_pt", tx->getBranch<vector<float>>("lep_pt")[index]);
        tx->setBranch<float>("lbnt_el_abseta", fabs(tx->getBranch<vector<float>>("lep_eta")[index]));
        tx->setBranch<float>("lbnt_mu_coneCorrPt", -999);
        tx->setBranch<float>("lbnt_mu_pt", -999);
        tx->setBranch<float>("lbnt_mu_abseta", -999);
    }
    else if (abs(tx->getBranch<vector<int>>("lep_pdgId")[index]) == 13)
    {
        tx->setBranch<float>("lbnt_el_coneCorrPt", -999);
        tx->setBranch<float>("lbnt_el_pt", -999);
        tx->setBranch<float>("lbnt_el_abseta", -999);
        tx->setBranch<float>("lbnt_mu_coneCorrPt", tx->getBranch<vector<float>>("lep_coneCorrPt")[index]);
        tx->setBranch<float>("lbnt_mu_pt", tx->getBranch<vector<float>>("lep_pt")[index]);
        tx->setBranch<float>("lbnt_mu_abseta", fabs(tx->getBranch<vector<float>>("lep_eta")[index]));
    }
    else
    {
        std::cout <<  " index: " << index <<  std::endl;
        std::cout <<  " tx->getBranch<vector<int>>('lep_pdgId').size(): " << tx->getBranch<vector<int>>("lep_pdgId").size() <<  std::endl;
        std::cout <<  " abs(tx->getBranch<vector<int>>('lep_pdgId')[index]): " << abs(tx->getBranch<vector<int>>("lep_pdgId")[index]) <<  std::endl;
        std::cout <<  " nVlep: " << nVlep <<  std::endl;
        FATALERROR(__FUNCTION__);
    }

    return;
}

//##############################################################################################################
void babyMaker_v2::FillYearInfo()
{
    tx->setBranch<int>("is2016", (coreSample.is2016(looper.getCurrentFileName())));
    tx->setBranch<int>("is2017", (coreSample.is2017(looper.getCurrentFileName())));
}

//##############################################################################################################
void babyMaker_v2::FillSummaryVariables()
{
    FillJetVariables(0);
    FillJetVariables(1);
    FillJetVariables(-1);
    FillJetVariables(2); // JER
    FillJetVariables(3); // JER up
    FillJetVariables(-2); // JER down
    FillLeptonVariables();
    FillLepJetVariables(0);
    FillLepJetVariables(1);
    FillLepJetVariables(-1);
    FillLepJetVariables(2); // JER
    FillLepJetVariables(3); // JER up
    FillLepJetVariables(-3); // JER down
    FillEventTags();
}

//##############################################################################################################
void babyMaker_v2::FillTnPMuons()
{
    for (unsigned int i = 0; i < coreMuon.index.size(); ++i)
    {
        int idx = coreMuon.index[i];
        int tag_idx = coreMuon.index.size() == coreMuon.tagindex.size() ? coreMuon.tagindex[i] : -1;

        if (tag_idx < 0)
            continue;

        // Get the tag p4
        LV dilep_p4 = cms3.mus_p4()[idx] + cms3.mus_p4()[tag_idx];
        float dilep_mass = dilep_p4.mass();

        // Skip if off the Z peak
        if (!( dilep_mass > 60 && dilep_mass < 120 ))
            continue;

        // Fill baby ntuple branches with event information (evt, lumi etc.)
        // These are always created and filled
        FillEventInfo();
        FillMuonTrigger(idx, tag_idx);
        FillMuonIDVariables(idx, tag_idx);
        FillTTree();
    }
}

//##############################################################################################################
void babyMaker_v2::FillMuonTrigger(int idx, int tag_idx)
{
    // Tag single muon trigger
    tx->setBranch<int>("tag_HLT_IsoMu24", passHLTTriggerPattern("HLT_IsoMu24_v"));
    tx->setBranch<int>("tag_HLT_IsoMu27", passHLTTriggerPattern("HLT_IsoMu27_v"));
    tx->setBranch<int>("tag_HLT_IsoTkMu24", passHLTTriggerPattern("HLT_IsoTkMu24_v"));
}

//
//##############################################################################################################
void babyMaker_v2::FillMuonIDVariables(int idx, int tag_idx)
{
    LV dilep_p4 = tag_idx >= 0 ? cms3.mus_p4()[idx] + cms3.mus_p4()[tag_idx] : LV();
    float dilep_mass = dilep_p4.mass();

    tx->setBranch<LV>("dilep_p4", dilep_p4);
    tx->setBranch<float>("dilep_mass", dilep_mass);
    tx->setBranch<LV>("p4", cms3.mus_p4()[idx]);
    tx->setBranch<LV>("tag_p4", tag_idx >= 0 ? cms3.mus_p4()[tag_idx] : LV());
    tx->setBranch<int>("tag_charge", tag_idx >= 0 ? cms3.mus_charge()[tag_idx] : 0);
    tx->setBranch<bool>("conv_vtx_flag", false);
    tx->setBranch<bool>("evt_isRealData", cms3.evt_isRealData());
    tx->setBranch<bool>("isPF", false);

    tx->setBranch<int>("id", -13.0 * cms3.mus_charge()[idx]);

    tx->setBranch<float>("RelIso03EA", muRelIso03EA(idx));
    tx->setBranch<float>("RelIso03EAv2", muRelIso03EA(idx, 2));
    tx->setBranch<float>("RelIso03EAv2wlep", muRelIso03EA(idx, 2, true));

    tx->setBranch<float>("dxyPV", cms3.mus_dxyPV()[idx]);
    tx->setBranch<float>("dZ", cms3.mus_dzPV()[idx]);
    tx->setBranch<float>("ip3d", cms3.mus_ip3d()[idx]);
    tx->setBranch<float>("ip3derr", cms3.mus_ip3derr()[idx]);

    tx->setBranch<bool>("passes_VVV_tight", muonID(idx, VVV_tight_v5));
    tx->setBranch<bool>("passes_VVV_tight_noiso", muonID(idx, VVV_tight_noiso_v5));
    tx->setBranch<bool>("passes_VVV_3l_tight", muonID(idx, VVV_3l_tight_v5));
    tx->setBranch<bool>("passes_VVV_3l_tight_noiso", muonID(idx, VVV_3l_tight_noiso_v5));
    tx->setBranch<bool>("passes_VVV_fo", muonID(idx, VVV_fo_v5));
    tx->setBranch<bool>("passes_VVV_fo_noiso", muonID(idx, VVV_fo_noiso_v5));
    tx->setBranch<bool>("passes_VVV_3l_fo", muonID(idx, VVV_3l_fo_v5));
    tx->setBranch<bool>("passes_VVV_3l_fo_noiso", muonID(idx, VVV_3l_fo_noiso_v5));
    tx->setBranch<bool>("passes_VVV_veto", muonID(idx, VVV_veto_v5));
    tx->setBranch<bool>("passes_VVV_veto_noiso", muonID(idx, VVV_veto_noiso_v5));
    tx->setBranch<bool>("passes_POG_mediumID", isMediumMuonPOG(idx));
    tx->setBranch<bool>("passes_POG_looseID", isLooseMuonPOG(idx));

    if (!cms3.evt_isRealData())
    {
        tx->setBranch<LV>("mc_motherp4", cms3.mus_mc_motherp4()[idx]);
        tx->setBranch<LV>("mc_p4", cms3.mus_mc_p4()[idx]);
        tx->setBranch<int>("mc_motherid", cms3.mus_mc_motherid()[idx]);
        tx->setBranch<int>("motherID", lepMotherID_v2(Lep(-13 * cms3.mus_charge()[idx], idx)).first);
    }

    tx->setBranch<float>("ptErr", cms3.mus_ptErr()[idx]);
    tx->setBranch<float>("trk_pt", cms3.mus_trk_p4()[idx].pt());

    tx->setBranch<int>("evt_event", tx->getBranch<unsigned long long>("evt"));
    tx->setBranch<int>("evt_lumiBlock", tx->getBranch<int>("lumi"));
    tx->setBranch<int>("evt_run", tx->getBranch<int>("run"));

    // Count number of good vertices.
    int nvtx = 0;
    for (unsigned int ivtx = 0; ivtx < cms3.evt_nvtxs(); ivtx++)
    {
        if (!isGoodVertex(ivtx)) { continue; }
        nvtx++;
    }
    tx->setBranch<int>("nvtx", nvtx);
}

//##############################################################################################################
void babyMaker_v2::FillTnPElectrons()
{
    for (unsigned int i = 0; i < coreElectron.index.size(); ++i)
    {
        int idx = coreElectron.index[i];
        int tag_idx = coreElectron.index.size() == coreElectron.tagindex.size() ? coreElectron.tagindex[i] : -1;

        if (tag_idx < 0)
            continue;

        // Get the tag p4
        LV dilep_p4 = cms3.els_p4()[idx] + cms3.els_p4()[tag_idx];
        float dilep_mass = dilep_p4.mass();

        // Skip if off the Z peak
        if (!( dilep_mass > 60 && dilep_mass < 120 ))
            continue;

        // Fill baby ntuple branches with event information (evt, lumi etc.)
        // These are always created and filled
        FillEventInfo();
        FillElectronTrigger(idx, tag_idx);
        FillElectronIDVariables(idx, tag_idx);
        FillTTree();
    }
}

//##############################################################################################################
void babyMaker_v2::FillElectronTrigger(int idx, int tag_idx)
{
    // Set the branch of the tag electrons single lepton trigger bit
    tx->setBranch<int>("tag_HLT_Ele27_eta2p1_WPTight_Gsf", passHLTTriggerPattern("HLT_Ele27_eta2p1_WPTight_Gsf_v"));
    tx->setBranch<int>("tag_HLT_Ele32_eta2p1_WPTight_Gsf", passHLTTriggerPattern("HLT_Ele32_eta2p1_WPTight_Gsf_v"));
    tx->setBranch<int>("tag_HLT_Ele32_WPTight_Gsf", passHLTTriggerPattern("HLT_Ele32_WPTight_Gsf_v"));
    tx->setBranch<int>("tag_HLT_Ele35_WPTight_Gsf", passHLTTriggerPattern("HLT_Ele35_WPTight_Gsf_v"));
}

//##############################################################################################################
void babyMaker_v2::FillElectronIDVariables(int idx, int tag_idx)
{
    LV dilep_p4 = tag_idx >= 0 ? cms3.els_p4()[idx] + cms3.els_p4()[tag_idx] : LV();
    float dilep_mass = dilep_p4.mass();

    tx->setBranch<LV>("dilep_p4", dilep_p4);
    tx->setBranch<float>("dilep_mass", dilep_mass);
    tx->setBranch<LV>("p4", cms3.els_p4()[idx]);
    tx->setBranch<LV>("tag_p4", tag_idx >= 0 ? cms3.els_p4()[tag_idx] : LV());
    tx->setBranch<int>("tag_charge", tag_idx >= 0 ? cms3.els_charge()[tag_idx] : 0);
    tx->setBranch<bool>("conv_vtx_flag", false);
    tx->setBranch<bool>("evt_isRealData", cms3.evt_isRealData());
    tx->setBranch<bool>("isPF", false);

    tx->setBranch<int>("id", -11.0 * cms3.els_charge()[idx]);

    tx->setBranch<float>("RelIso03EA", eleRelIso03EA(idx));
    tx->setBranch<float>("RelIso03EAv2", eleRelIso03EA(idx, 2));
    tx->setBranch<float>("RelIso03EAv2wlep", eleRelIso03EA(idx, 2, true));

    tx->setBranch<float>("dxyPV", cms3.els_dxyPV()[idx]);
    tx->setBranch<float>("dZ", cms3.els_dzPV()[idx]);
    tx->setBranch<float>("ip3d", cms3.els_ip3d()[idx]);
    tx->setBranch<float>("ip3derr", cms3.els_ip3derr()[idx]);

    tx->setBranch<bool>("passes_VVV_tight", electronID(idx, VVV_tight_v5));
    tx->setBranch<bool>("passes_VVV_tight_noiso", electronID(idx, VVV_tight_noiso_v5));
    tx->setBranch<bool>("passes_VVV_3l_tight", electronID(idx, VVV_3l_tight_v5));
    tx->setBranch<bool>("passes_VVV_3l_tight_noiso", electronID(idx, VVV_3l_tight_noiso_v5));
    tx->setBranch<bool>("passes_VVV_fo", electronID(idx, VVV_fo_v5));
    tx->setBranch<bool>("passes_VVV_fo_noiso", electronID(idx, VVV_fo_noiso_v5));
    tx->setBranch<bool>("passes_VVV_3l_fo", electronID(idx, VVV_3l_fo_v5));
    tx->setBranch<bool>("passes_VVV_3l_fo_noiso", electronID(idx, VVV_3l_fo_noiso_v5));
    tx->setBranch<bool>("passes_VVV_veto", electronID(idx, VVV_veto_v5));
    tx->setBranch<bool>("passes_VVV_veto_noiso", electronID(idx, VVV_veto_noiso_v5));
    tx->setBranch<bool>("passes_POG_mediumID", isMediumElectronPOG(idx));
    tx->setBranch<bool>("passes_POG_looseID", isLooseElectronPOG(idx));
    tx->setBranch<bool>("passes_POG_MVA_HZZ", isMVAHZZNoIsofall17(idx));
    tx->setBranch<bool>("passes_POG_MVA_wp80", isMVAwp80NoIsofall17(idx));
    tx->setBranch<bool>("passes_POG_MVA_wp90", isMVAwp90NoIsofall17(idx));

    if (!cms3.evt_isRealData())
    {
        tx->setBranch<LV>("mc_motherp4", cms3.els_mc_motherp4()[idx]);
        tx->setBranch<LV>("mc_p4", cms3.els_mc_p4()[idx]);
        tx->setBranch<int>("mc_motherid", cms3.els_mc_motherid()[idx]);
        tx->setBranch<int>("motherID", lepMotherID_v2(Lep(-11 * cms3.els_charge()[idx], idx)).first);
    }

    tx->setBranch<bool>("threeChargeAgree", threeChargeAgree(idx));
    tx->setBranch<float>("dEtaIn", cms3.els_dEtaIn()[idx]);
    tx->setBranch<float>("dPhiIn", cms3.els_dPhiIn()[idx]);
    tx->setBranch<float>("eOverPIn", cms3.els_eOverPIn()[idx]);
    tx->setBranch<float>("ecalEnergy", cms3.els_ecalEnergy()[idx]);
    tx->setBranch<float>("etaSC", cms3.els_etaSC()[idx]);
    tx->setBranch<float>("hOverE", cms3.els_hOverE()[idx]);
    tx->setBranch<float>("sigmaIEtaIEta_full5x5", cms3.els_sigmaIEtaIEta_full5x5()[idx]);
    tx->setBranch<float>("mva", getMVAoutput(idx));
    tx->setBranch<float>("mva_25ns", getMVAoutput(idx));
    tx->setBranch<float>("mva_2017", getMVAoutput(idx, true));

    tx->setBranch<int>("evt_event", tx->getBranch<unsigned long long>("evt"));
    tx->setBranch<int>("evt_lumiBlock", tx->getBranch<int>("lumi"));
    tx->setBranch<int>("evt_run", tx->getBranch<int>("run"));

    // Count number of good vertices.
    int nvtx = 0;
    for (unsigned int ivtx = 0; ivtx < cms3.evt_nvtxs(); ivtx++)
    {
        if (!isGoodVertex(ivtx)) { continue; }
        nvtx++;
    }
    tx->setBranch<int>("nvtx", nvtx);
}

//##############################################################################################################
void babyMaker_v2::FillTTree()
{
    tx->fill();
    tx->clear();
}

//##############################################################################################################
bool babyMaker_v2::isLeptonOverlappingWithJet(int ijet)
{
    bool is_overlapping = false;

    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);
    const vector<int>& lep_id = tx->getBranch<vector<int>>(gconf.wwwcfg["3llooseid"], true);

    int idx = coreJet.index[ijet];

    for (unsigned ilep = 0; ilep < lep_id.size(); ++ilep)
    {
        const LV& p4 = lep_p4[ilep];
        const int& id = lep_id[ilep];
        if (id == 0) continue;
        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfjets_p4()[idx], p4) < 0.4)
        {
            is_overlapping = true;
            break;
        }
    }
//
//    int idx = coreJet.index[ijet];
//
//    for (auto& imu : coreMuon.index)
//    {
//        if (!(isLooseMuon(imu)))
//            continue;
//
//        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfjets_p4()[idx], cms3.mus_p4()[imu]) < 0.4)
//        {
//            is_overlapping = true;
//            break;
//        }
//    }
//
//    if (is_overlapping)
//        return true;
//
//    for (auto& iel : coreElectron.index)
//    {
//        if (!(isLooseElectron(iel)))
//            continue;
//
//        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfjets_p4()[idx], cms3.els_p4()[iel]) < 0.4)
//        {
//            is_overlapping = true;
//            break;
//        }
//    }

    if (is_overlapping)
        return true;

    return false;
}

//##############################################################################################################
bool babyMaker_v2::isLeptonOverlappingWithTrack(int idx)
{
    bool is_overlapping = false;

    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);
    const vector<int>& lep_id = tx->getBranch<vector<int>>(gconf.wwwcfg["vetoid"], true);

    for (unsigned ilep = 0; ilep < lep_id.size(); ++ilep)
    {
        const LV& p4 = lep_p4[ilep];
        const int& id = lep_id[ilep];
        if (id == 0) continue;
        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfcands_p4()[idx], p4) < 0.4)
        {
            is_overlapping = true;
            break;
        }
    }

//    for (auto& imu : coreMuon.index)
//    {
//        if (!(isVetoMuon(imu)))
//            continue;
//
//        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfcands_p4()[idx], cms3.mus_p4()[imu]) < 0.01)
//        {
//            is_overlapping = true;
//            break;
//        }
//    }
//
//    if (is_overlapping)
//        return true;
//
//    for (auto& iel : coreElectron.index)
//    {
//        if (!(isVetoElectron(iel)))
//            continue;
//
//        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfcands_p4()[idx], cms3.els_p4()[iel]) < 0.01)
//        {
//            is_overlapping = true;
//            break;
//        }
//    }

    if (is_overlapping)
        return true;

    return false;
}

//##############################################################################################################
// Used to overlap remova against jets
bool babyMaker_v2::isLooseMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 20.         )) return false;
    if (!( passMuonSelection_VVV(idx, VVV_FO_3L) )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against jets
bool babyMaker_v2::isLooseElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 20.             )) return false;
    if (!( passElectronSelection_VVV(idx, VVV_FO_3L) )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::isVetoMuon(int idx)
{
    if (!( passMuonSelection_VVV(idx, VVV_VETO) )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::isVetoElectron(int idx)
{
    if (!( passElectronSelection_VVV(idx, VVV_VETO) )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::isPt10Electron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 10.        )) return false;
    if (!( fabs(cms3.els_p4()[idx].eta()) < 2.5 )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::isPt20Electron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 20.        )) return false;
    if (!( fabs(cms3.els_p4()[idx].eta()) < 2.5 )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGVetoElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 10.        )) return false;
    if (!( fabs(cms3.els_p4()[idx].eta()) < 2.5 )) return false;
    if (!( isVetoElectronPOGfall17_v2(idx)      )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGLooseElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 20.        )) return false;
    if (!( fabs(cms3.els_p4()[idx].eta()) < 2.5 )) return false;
    if (!( isLooseElectronPOGfall17_v2(idx)     )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGMediumElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 20.        )) return false;
    if (!( fabs(cms3.els_p4()[idx].eta()) < 2.5 )) return false;
    if (!( isMediumElectronPOGfall17_v2(idx)    )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGTightElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 20.        )) return false;
    if (!( fabs(cms3.els_p4()[idx].eta()) < 2.5 )) return false;
    if (!( isTightElectronPOGfall17_v2(idx)     )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGVetoMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 10.        )) return false;
    if (!( fabs(cms3.mus_p4()[idx].eta()) < 2.4 )) return false;
    if (!( isLooseMuonPOG(idx)                  )) return false;
    if (!( muRelIso04DB(idx) < 0.25             )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGLooseMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 20.        )) return false;
    if (!( fabs(cms3.mus_p4()[idx].eta()) < 2.4 )) return false;
    if (!( isLooseMuonPOG(idx)                  )) return false;
    if (!( muRelIso04DB(idx) < 0.25             )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGMediumMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 20.        )) return false;
    if (!( fabs(cms3.mus_p4()[idx].eta()) < 2.4 )) return false;
    if (!( isMediumMuonPOG(idx)                 )) return false;
    if (!( muRelIso04DB(idx) < 0.25             )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker_v2::is2017POGTightMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 20.        )) return false;
    if (!( fabs(cms3.mus_p4()[idx].eta()) < 2.4 )) return false;
    if (!( isTightMuonPOG(idx)                  )) return false;
    if (!( muRelIso04DB(idx) < 0.05             )) return false;
    return true;
}

//##############################################################################################################
// Used for probe in Tag-and-Probe method (i.e. no cuts whatsoever other than pt)
bool babyMaker_v2::isProbeMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 10. )) return false;
    return true;
}

//##############################################################################################################
// Used for probe in Tag-and-Probe method (i.e. no cuts whatsoever other than pt)
bool babyMaker_v2::isProbeElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 10. )) return false;
    return true;
}

//##############################################################################################################
// Used for tag in Tag-and-Probe method (i.e. some tight requirement)
bool babyMaker_v2::isTagMuon(int i, int j)
{
    // Tag muon selection
    if (!( cms3.mus_p4()[j].pt()                                    >= 20.0  )) return false;
    if (!( fabs(cms3.mus_p4()[j].eta())                             <=  2.4  )) return false;
    if (!( fabs(cms3.mus_dxyPV()[j])                                <=  0.02 )) return false;
    if (!( fabs(cms3.mus_dzPV()[j])                                 <=  0.05 )) return false;
    if (!( fabs(cms3.mus_ip3d()[j] / cms3.mus_ip3derr()[j])         <=  4    )) return false;
    if (!( isTightMuonPOG(j)                                                 )) return false;
    if (!( muRelIso03EA(j)                                          <=  0.2  )) return false;
    return true;
}

//##############################################################################################################
// Used for tag in Tag-and-Probe method (i.e. some tight requirement)
bool babyMaker_v2::isTagElectron(int i, int j)
{
    if (!( cms3.els_p4()[j].pt()                                    >= 20.0  )) return false;
    if (!( fabs(cms3.els_etaSC()[j])                                <=  2.5  )) return false;
    if (!( isMediumElectronPOG(j)                                            )) return false;
    if (!( fabs(cms3.els_ip3d()[j] / cms3.els_ip3derr()[j])         <=  4    )) return false;
    return true;
}

//##############################################################################################################
int babyMaker_v2::passCount(const vector<int>& v)
{
    return std::count_if(v.begin(), v.end(), [](int i){return i > 0;});
}

//##############################################################################################################
void babyMaker_v2::sortVecP4(vector<LV>& v)
{
    std::sort(v.begin(), v.end(), [](const LV& a, const LV& b){return a.pt() > b.pt();});
}

//##############################################################################################################
int babyMaker_v2::nSFOS()
{
    int nSFOS = 0;
    for (auto& pdgid0 : tx->getBranch<vector<int>>("lep_pdgId"))
    {
        for (auto& pdgid1 : tx->getBranch<vector<int>>("lep_pdgId"))
        {
            if (pdgid0 == -pdgid1)
                nSFOS++;
        }
    }
    if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
    {
        std::cout <<  " nSFOS/2: " << nSFOS/2 <<  std::endl;
    }
    return nSFOS/2; // Because the nested for loop double counts
}

//##############################################################################################################
void babyMaker_v2::FillJetVariables(int variation)
{
    //
    // Assumes FillJets and SortJetBranches are already called
    //

    TString jets_csv     = variation == 0 ? "jets_csv"          : variation == 1 ? "jets_up_csv"          : variation ==-1 ?  "jets_dn_csv"         : variation == 2 ? "jets_jer_csv"         : variation == 3 ? "jets_jerup_csv"        : /*variation ==-3 ?*/ "jets_jerdn_csv"        ;
    TString jets_p4      = variation == 0 ? "jets_p4"           : variation == 1 ? "jets_up_p4"           : variation ==-1 ?  "jets_dn_p4"          : variation == 2 ? "jets_jer_p4"          : variation == 3 ? "jets_jerup_p4"         : /*variation ==-3 ?*/ "jets_jerdn_p4"         ;
    TString nj_bn        = variation == 0 ? "nj"                : variation == 1 ? "nj_up"                : variation ==-1 ?  "nj_dn"               : variation == 2 ? "nj_jer"               : variation == 3 ? "nj_jerup"              : /*variation ==-3 ?*/ "nj_jerdn"              ;
    TString nj30_bn      = variation == 0 ? "nj30"              : variation == 1 ? "nj30_up"              : variation ==-1 ?  "nj30_dn"             : variation == 2 ? "nj30_jer"             : variation == 3 ? "nj30_jerup"            : /*variation ==-3 ?*/ "nj30_jerdn"            ;
    TString nb_bn        = variation == 0 ? "nb"                : variation == 1 ? "nb_up"                : variation ==-1 ?  "nb_dn"               : variation == 2 ? "nb_jer"               : variation == 3 ? "nb_jerup"              : /*variation ==-3 ?*/ "nb_jerdn"              ;
    TString Mjj_bn       = variation == 0 ? "Mjj"               : variation == 1 ? "Mjj_up"               : variation ==-1 ?  "Mjj_dn"              : variation == 2 ? "Mjj_jer"              : variation == 3 ? "Mjj_jerup"             : /*variation ==-3 ?*/ "Mjj_jerdn"             ;
    TString MjjL_bn      = variation == 0 ? "MjjL"              : variation == 1 ? "MjjL_up"              : variation ==-1 ?  "MjjL_dn"             : variation == 2 ? "MjjL_jer"             : variation == 3 ? "MjjL_jerup"            : /*variation ==-3 ?*/ "MjjL_jerdn"            ;
    TString DetajjL_bn   = variation == 0 ? "DetajjL"           : variation == 1 ? "DetajjL_up"           : variation ==-1 ?  "DetajjL_dn"          : variation == 2 ? "DetajjL_jer"          : variation == 3 ? "DetajjL_jerup"         : /*variation ==-3 ?*/ "DetajjL_jerdn"         ;
    TString MjjVBF_bn    = variation == 0 ? "MjjVBF"            : variation == 1 ? "MjjVBF_up"            : variation ==-1 ?  "MjjVBF_dn"           : variation == 2 ? "MjjVBF_jer"           : variation == 3 ? "MjjVBF_jerup"          : /*variation ==-3 ?*/ "MjjVBF_jerdn"          ;
    TString DetajjVBF_bn = variation == 0 ? "DetajjVBF"         : variation == 1 ? "DetajjVBF_up"         : variation ==-1 ?  "DetajjVBF_dn"        : variation == 2 ? "DetajjVBF_jer"        : variation == 3 ? "DetajjVBF_jerup"       : /*variation ==-3 ?*/ "DetajjVBF_jerdn"       ;
    TString MjjDR1_bn    = variation == 0 ? "MjjDR1"            : variation == 1 ? "MjjDR1_up"            : variation ==-1 ?  "MjjDR1_dn"           : variation == 2 ? "MjjDR1_jer"           : variation == 3 ? "MjjDR1_jerup"          : /*variation ==-3 ?*/ "MjjDR1_jerdn"          ;
    TString DRjj_bn      = variation == 0 ? "DRjj"              : variation == 1 ? "DRjj_up"              : variation ==-1 ?  "DRjj_dn"             : variation == 2 ? "DRjj_jer"             : variation == 3 ? "DRjj_jerup"            : /*variation ==-3 ?*/ "DRjj_jerdn"            ;
    TString DRjjDR1_bn   = variation == 0 ? "DRjjDR1"           : variation == 1 ? "DRjjDR1_up"           : variation ==-1 ?  "DRjjDR1_dn"          : variation == 2 ? "DRjjDR1_jer"          : variation == 3 ? "DRjjDR1_jerup"         : /*variation ==-3 ?*/ "DRjjDR1_jerdn"         ;
    TString j0_p4_bn     = variation == 0 ? "jet0_wtag_p4"      : variation == 1 ? "jet0_wtag_p4_up"      : variation ==-1 ?  "jet0_wtag_p4_dn"     : variation == 2 ? "jet0_wtag_p4_jer"     : variation == 3 ? "jet0_wtag_p4_jerup"    : /*variation ==-3 ?*/ "jet0_wtag_p4_jerdn"    ;
    TString j1_p4_bn     = variation == 0 ? "jet1_wtag_p4"      : variation == 1 ? "jet1_wtag_p4_up"      : variation ==-1 ?  "jet1_wtag_p4_dn"     : variation == 2 ? "jet1_wtag_p4_jer"     : variation == 3 ? "jet1_wtag_p4_jerup"    : /*variation ==-3 ?*/ "jet1_wtag_p4_jerdn"    ;
    TString j0_p4_DR1_bn = variation == 0 ? "jet0_wtag_p4_DR1"  : variation == 1 ? "jet0_wtag_p4_DR1_up"  : variation ==-1 ?  "jet0_wtag_p4_DR1_dn" : variation == 2 ? "jet0_wtag_p4_DR1_jer" : variation == 3 ? "jet0_wtag_p4_DR1_jerup": /*variation ==-3 ?*/ "jet0_wtag_p4_DR1_jerdn";
    TString j1_p4_DR1_bn = variation == 0 ? "jet1_wtag_p4_DR1"  : variation == 1 ? "jet1_wtag_p4_DR1_up"  : variation ==-1 ?  "jet1_wtag_p4_DR1_dn" : variation == 2 ? "jet1_wtag_p4_DR1_jer" : variation == 3 ? "jet1_wtag_p4_DR1_jerup": /*variation ==-3 ?*/ "jet1_wtag_p4_DR1_jerdn";

    int nj = 0;
    int nj30 = 0;
    int nb = 0;
    float Mjj = 0;
    float MjjL = -999;
    float DetajjL = -999;
    float MjjVBF = -999;
    float MjjDR1 = -999;
    float DetajjVBF = -999;
    float tmpDR = 9999;
    float tmpDR_DR1 = 9999;
    float DRjj = -999;
    float DRjjDR1 = -999;
    LV j0_p4;
    LV j1_p4;
    LV j0_p4_DR1;
    LV j1_p4_DR1;
    for (unsigned int i = 0; i < tx->getBranch<vector<float>>(jets_csv, false).size(); ++i)
    {
        const LV& p4 = tx->getBranch<vector<LV>>(jets_p4, false)[i];
        const float& csv= tx->getBranch<vector<float>>(jets_csv, false)[i];

        // nb jets
        if (!( p4.pt() > 20. ))
            continue;
        if (fabs(p4.eta()) < 2.4 && csv >= 0.5426)
            nb++;

        // njets across all eta
        if (!( p4.pt() > 30. ))
            continue;
        if (fabs(p4.eta()) < 5.0)
            nj++;

        // njets for central for SS
        if (fabs(p4.eta()) < 2.5)
            nj30++;

        // Compute Mjj using the closest two jets
        for (unsigned int j = i + 1; j < tx->getBranch<vector<float>>(jets_csv, false).size(); ++j)
        {
            const LV& p4_2 = tx->getBranch<vector<LV>>(jets_p4, false)[j];

            if (!( p4_2.pt() > 30. ))
                continue;

            // central eta
            if (fabs(p4.eta()) < 2.5 && fabs(p4_2.eta()) < 2.5)
            {
                // Choose the closest two jets
                float this_dR = ROOT::Math::VectorUtil::DeltaR(p4, p4_2);
                if (this_dR < tmpDR)
                {
                    tmpDR = this_dR;
                    Mjj = (p4 + p4_2).M();
                    DRjj = tmpDR;
                    j0_p4 = p4.pt() > p4_2.pt() ? p4 : p4_2;
                    j1_p4 = p4.pt() > p4_2.pt() ? p4_2 : p4;
                }

                // Choose the jets with angle closest to dR = 1
                if (fabs(this_dR - 1.) < fabs(tmpDR_DR1 - 1.))
                {
                    tmpDR_DR1 = this_dR;
                    MjjDR1 = (p4 + p4_2).M();
                    DRjjDR1 = tmpDR_DR1;
                    j0_p4_DR1 = p4.pt() > p4_2.pt() ? p4 : p4_2;
                    j1_p4_DR1 = p4.pt() > p4_2.pt() ? p4_2 : p4;
                }

                // If they were not set then set (this makes it choose two leading ones)
                if (MjjL < 0)
                {
                    MjjL = (p4 + p4_2).M();
                    DetajjL = fabs(p4.eta() - p4_2.eta());
                }
            }

            // all eta
            if (fabs(p4.eta()) < 5.0 && fabs(p4_2.eta()) < 5.0)
            {
                // If they were not set then set (this makes it choose two leading ones)
                if (MjjVBF < 0)
                {
                    MjjVBF = (p4 + p4_2).M();
                    DetajjVBF = fabs(p4.eta() - p4_2.eta());
                }
            }
        }
    }
    tx->setBranch<int>(nj_bn, nj);
    tx->setBranch<int>(nj30_bn, nj30);
    tx->setBranch<int>(nb_bn, nb);
    tx->setBranch<float>(Mjj_bn, Mjj);
    tx->setBranch<float>(MjjL_bn, MjjL);
    tx->setBranch<float>(DetajjL_bn, DetajjL);
    tx->setBranch<float>(MjjVBF_bn, MjjVBF);
    tx->setBranch<float>(MjjDR1_bn, MjjDR1);
    tx->setBranch<float>(DetajjVBF_bn, DetajjVBF);
    tx->setBranch<float>(DRjj_bn, DRjj);
    tx->setBranch<float>(DRjjDR1_bn, DRjjDR1);
    tx->setBranch<LV>(j0_p4_bn, j0_p4);
    tx->setBranch<LV>(j1_p4_bn, j1_p4);
    tx->setBranch<LV>(j0_p4_DR1_bn, j0_p4_DR1);
    tx->setBranch<LV>(j1_p4_DR1_bn, j1_p4_DR1);
}

//##############################################################################################################
void babyMaker_v2::FillLeptonVariables()
{

    //
    // Assumes FillElectrons and FillLeptons as well as SortLeptonBranches are called
    //

    // Calculate tight and loose lepton based on number of veto leptons
    int nVlep = coreElectron.index.size() + coreMuon.index.size();
    int nTlep = nVlep == 2 ? passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["tightid"])) : passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["3ltightid"]));
    int nLlep = nVlep == 2 ? passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["looseid"])) : passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["3llooseid"]));
    int nTlepSS = passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["tightid"]));
    int nLlep3L = passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["3llooseid"]));
    int nTlep3L = passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["3ltightid"]));

    // Set the number of lepton counter
    tx->setBranch<int>("nVlep", nVlep);
    tx->setBranch<int>("nTlep", nTlep);
    tx->setBranch<int>("nTlepSS", nTlepSS);
    tx->setBranch<int>("nLlep", nLlep);
    tx->setBranch<int>("nLlep3L", nLlep3L);
    tx->setBranch<int>("nTlep3L", nTlep3L);

    // nVlep == 2
    //  nLlep == 2 isSS
    // nVlep >= 3
    //  nLlep == 3 is3l
    //    nLlepSS >= 2 (is the pair a same-sign?) need idx isSS3lCR

    //bool isSS = (nVlep == 2 && nLlep == 2);
    //bool is3L = (nVlep >= 3 && nLlep == 3);
    bool isSS = (nLlep == 2);
    bool is3L = (nLlep == 3);

    // is SS events then set the SS related variables
    if (isSS)
    {
        FillSSLeptonVariables(0, 1);
    }
    // is 3L events then set the 3L related variables
    else if (is3L)
    {
        Fill3LLeptonVariables();

        // If it additionally passes isSSCR then set the "Same-sign" variables for the 3L events
        bool issscr;
        int index0;
        int index1;
        tie(issscr, index0, index1) = isSSCR();
        if (issscr)
        {
            FillSSLeptonVariables(index0, index1);
        }
    }
}

//##############################################################################################################
void babyMaker_v2::FillLepJetVariables(int variation)
{

    // check if two leptons and two jets
    TString nj30  = variation == 0 ? "nj30"  : variation == 1 ? "nj30_up"  : variation ==-1 ? "nj30_dn" : variation== 2 ? "nj30_jer" : variation == 3 ? "nj30_jerup" : /*variation ==-3 ?*/ "nj30_jerdn";
    if (tx->getBranch<int>("nVlep") == 2 and tx->getBranch<int>(nj30) >= 2)
    {

        // Get leptons
        const LV& lep_p4_0 = tx->getBranch<vector<LV>>("lep_p4")[0];
        const LV& lep_p4_1 = tx->getBranch<vector<LV>>("lep_p4")[1];

        //
        // Assumes FillJets and SortJetBranches are already called as well as leptons
        //

        TString Ml0j0_bn   = variation == 0 ? "Ml0j0"      : variation == 1 ? "Ml0j0_up"      : variation ==-1 ? "Ml0j0_dn"       : variation == 2 ? "Ml0j0_jer"       : variation == 3 ? "Ml0j0_jerup"       : /*variation ==-3 ?*/ "Ml0j0_jerdn"      ;
        TString Ml0j1_bn   = variation == 0 ? "Ml0j1"      : variation == 1 ? "Ml0j1_up"      : variation ==-1 ? "Ml0j1_dn"       : variation == 2 ? "Ml0j1_jer"       : variation == 3 ? "Ml0j1_jerup"       : /*variation ==-3 ?*/ "Ml0j1_jerdn"      ;
        TString Ml1j0_bn   = variation == 0 ? "Ml1j0"      : variation == 1 ? "Ml1j0_up"      : variation ==-1 ? "Ml1j0_dn"       : variation == 2 ? "Ml1j0_jer"       : variation == 3 ? "Ml1j0_jerup"       : /*variation ==-3 ?*/ "Ml1j0_jerdn"      ;
        TString Ml1j1_bn   = variation == 0 ? "Ml1j1"      : variation == 1 ? "Ml1j1_up"      : variation ==-1 ? "Ml1j1_dn"       : variation == 2 ? "Ml1j1_jer"       : variation == 3 ? "Ml1j1_jerup"       : /*variation ==-3 ?*/ "Ml1j1_jerdn"      ;
        TString MinMlj_bn  = variation == 0 ? "MinMlj"     : variation == 1 ? "MinMlj_up"     : variation ==-1 ? "MinMlj_dn"      : variation == 2 ? "MinMlj_jer"      : variation == 3 ? "MinMlj_jerup"      : /*variation ==-3 ?*/ "MinMlj_jerdn"     ;
        TString MaxMlj_bn  = variation == 0 ? "MaxMlj"     : variation == 1 ? "MaxMlj_up"     : variation ==-1 ? "MaxMlj_dn"      : variation == 2 ? "MaxMlj_jer"      : variation == 3 ? "MaxMlj_jerup"      : /*variation ==-3 ?*/ "MaxMlj_jerdn"     ;
        TString SumMlj_bn  = variation == 0 ? "SumMlj"     : variation == 1 ? "SumMlj_up"     : variation ==-1 ? "SumMlj_dn"      : variation == 2 ? "SumMlj_jer"      : variation == 3 ? "SumMlj_jerup"      : /*variation ==-3 ?*/ "SumMlj_jerdn"     ;
        TString Ml0jj_bn   = variation == 0 ? "Ml0jj"      : variation == 1 ? "Ml0jj_up"      : variation ==-1 ? "Ml0jj_dn"       : variation == 2 ? "Ml0jj_jer"       : variation == 3 ? "Ml0jj_jerup"       : /*variation ==-3 ?*/ "Ml0jj_jerdn"      ;
        TString Ml1jj_bn   = variation == 0 ? "Ml1jj"      : variation == 1 ? "Ml1jj_up"      : variation ==-1 ? "Ml1jj_dn"       : variation == 2 ? "Ml1jj_jer"       : variation == 3 ? "Ml1jj_jerup"       : /*variation ==-3 ?*/ "Ml1jj_jerdn"      ;
        TString MinMljj_bn = variation == 0 ? "MinMljj"    : variation == 1 ? "MinMljj_up"    : variation ==-1 ? "MinMljj_dn"     : variation == 2 ? "MinMljj_jer"     : variation == 3 ? "MinMljj_jerup"     : /*variation ==-3 ?*/ "MinMljj_jerdn"    ;
        TString MaxMljj_bn = variation == 0 ? "MaxMljj"    : variation == 1 ? "MaxMljj_up"    : variation ==-1 ? "MaxMljj_dn"     : variation == 2 ? "MaxMljj_jer"     : variation == 3 ? "MaxMljj_jerup"     : /*variation ==-3 ?*/ "MaxMljj_jerdn"    ;
        TString SumMljj_bn = variation == 0 ? "SumMljj"    : variation == 1 ? "SumMljj_up"    : variation ==-1 ? "SumMljj_dn"     : variation == 2 ? "SumMljj_jer"     : variation == 3 ? "SumMljj_jerup"     : /*variation ==-3 ?*/ "SumMljj_jerdn"    ;
        TString SumMinMlj01= variation == 0 ? "SumMinMlj01": variation == 1 ? "SumMinMlj01_up": variation ==-1 ? "SumMinMlj01_dn" : variation == 2 ? "SumMinMlj01_jer" : variation == 3 ? "SumMinMlj01_jerup" : /*variation ==-3 ?*/ "SumMinMlj01_jerdn";
        TString jets30_p4  = variation == 0 ? "jets30_p4"  : variation == 1 ? "jets30_up_p4"  : variation ==-1 ? "jets30_dn_p4"   : variation == 2 ? "jets30_jer_p4"   : variation == 3 ? "jets30_jerup_p4"   : /*variation ==-3 ?*/ "jets30_jerdn_p4"  ;

        const LV& jet_p4_0 = tx->getBranch<vector<LV>>(jets30_p4)[0];
        const LV& jet_p4_1 = tx->getBranch<vector<LV>>(jets30_p4)[1];

        float Ml0j0 = (lep_p4_0 + jet_p4_0).mass();
        float Ml0j1 = (lep_p4_0 + jet_p4_1).mass();
        float Ml1j0 = (lep_p4_1 + jet_p4_0).mass();
        float Ml1j1 = (lep_p4_1 + jet_p4_1).mass();

        tx->setBranch<float>(Ml0j0_bn, Ml0j0);
        tx->setBranch<float>(Ml0j1_bn, Ml0j1);
        tx->setBranch<float>(Ml1j0_bn, Ml1j0);
        tx->setBranch<float>(Ml1j1_bn, Ml1j1);

        tx->setBranch<float>(MinMlj_bn, std::min( std::min( std::min(Ml0j0, Ml0j1), Ml1j0), Ml1j1));
        tx->setBranch<float>(MaxMlj_bn, std::max( std::max( std::max(Ml0j0, Ml0j1), Ml1j0), Ml1j1));
        tx->setBranch<float>(SumMlj_bn, Ml0j0 + Ml0j1 + Ml1j0 + Ml1j1);

        float Ml0jj = (lep_p4_0 + jet_p4_0 + jet_p4_1).mass();
        float Ml1jj = (lep_p4_1 + jet_p4_0 + jet_p4_1).mass();

        tx->setBranch<float>(Ml0jj_bn, Ml0jj);
        tx->setBranch<float>(Ml1jj_bn, Ml1jj);

        tx->setBranch<float>(MinMljj_bn, std::min(Ml0jj, Ml1jj));
        tx->setBranch<float>(MaxMljj_bn, std::min(Ml0jj, Ml1jj));
        tx->setBranch<float>(SumMljj_bn, Ml0jj + Ml1jj);

    }

}

//##############################################################################################################
void babyMaker_v2::FillSSLeptonVariables(int idx0, int idx1)
{
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);

    // Fill the same-sign pairs lep id tag
    tx->setBranch<int>("passSSee", (lep_pdgId[idx0] * lep_pdgId[idx1]) == 121);
    tx->setBranch<int>("passSSem", (lep_pdgId[idx0] * lep_pdgId[idx1]) == 143);
    tx->setBranch<int>("passSSmm", (lep_pdgId[idx0] * lep_pdgId[idx1]) == 169);

    // Save the lepton index for later purpose
    tx->setBranch<int>("lep_idx0_SS", idx0);
    tx->setBranch<int>("lep_idx1_SS", idx1);

    tx->setBranch<float>("MllSS", (lep_p4[idx0] + lep_p4[idx1]).mass());
    if (abs(lep_pdgId[idx0]) == 11 && abs(lep_pdgId[idx1]) == 11)
        tx->setBranch<float>("MeeSS", (lep_p4[idx0] + lep_p4[idx1]).mass());

    // Set up MET variables
    const float& met_pt = tx->getBranch<float>("met_pt");
    const float& met_phi = tx->getBranch<float>("met_phi");
    const float& met_up_pt = tx->getBranch<float>("met_up_pt");
    const float& met_up_phi = tx->getBranch<float>("met_up_phi");
    const float& met_dn_pt = tx->getBranch<float>("met_dn_pt");
    const float& met_dn_phi = tx->getBranch<float>("met_dn_phi");
    const float& met_gen_pt = tx->getBranch<float>("met_gen_pt");
    const float& met_gen_phi = tx->getBranch<float>("met_gen_phi");
    const float& met_jer_pt = tx->getBranch<float>("met_jer_pt");
    const float& met_jer_phi = tx->getBranch<float>("met_jer_phi");
    const float& met_jerup_pt = tx->getBranch<float>("met_jerup_pt");
    const float& met_jerup_phi = tx->getBranch<float>("met_jerup_phi");
    const float& met_jerdn_pt = tx->getBranch<float>("met_jerdn_pt");
    const float& met_jerdn_phi = tx->getBranch<float>("met_jerdn_phi");
    LV MET;
    LV MET_up;
    LV MET_dn;
    LV MET_gen;
    LV MET_jer;
    LV MET_jerup;
    LV MET_jerdn;
    MET.SetPxPyPzE(met_pt * TMath::Cos(met_phi), met_pt * TMath::Sin(met_phi), 0, met_pt);
    MET_up.SetPxPyPzE(met_up_pt * TMath::Cos(met_up_phi), met_up_pt * TMath::Sin(met_up_phi), 0, met_up_pt);
    MET_dn.SetPxPyPzE(met_dn_pt * TMath::Cos(met_dn_phi), met_dn_pt * TMath::Sin(met_dn_phi), 0, met_dn_pt);
    MET_gen.SetPxPyPzE(met_gen_pt * TMath::Cos(met_gen_phi), met_gen_pt * TMath::Sin(met_gen_phi), 0, met_gen_pt);
    MET_jer.SetPxPyPzE(met_jer_pt * TMath::Cos(met_jer_phi), met_jer_pt * TMath::Sin(met_jer_phi), 0, met_jer_pt);
    MET_jerup.SetPxPyPzE(met_jerup_pt * TMath::Cos(met_jerup_phi), met_jerup_pt * TMath::Sin(met_jerup_phi), 0, met_jerup_pt);
    MET_jerdn.SetPxPyPzE(met_jerdn_pt * TMath::Cos(met_jerdn_phi), met_jerdn_pt * TMath::Sin(met_jerdn_phi), 0, met_jerdn_pt);

    float MT0    = mT(lep_p4[idx0], MET);
    float MT1    = mT(lep_p4[idx1], MET);
    float MT0_up = mT(lep_p4[idx0], MET_up);
    float MT1_up = mT(lep_p4[idx1], MET_up);
    float MT0_dn = mT(lep_p4[idx0], MET_dn);
    float MT1_dn = mT(lep_p4[idx1], MET_dn);
    float MT0_gen = mT(lep_p4[idx0], MET_gen);
    float MT1_gen = mT(lep_p4[idx1], MET_gen);
    float MT0_jer = mT(lep_p4[idx0], MET_jer);
    float MT1_jer = mT(lep_p4[idx1], MET_jer);
    float MT0_jerup = mT(lep_p4[idx0], MET_jerup);
    float MT1_jerup = mT(lep_p4[idx1], MET_jerup);
    float MT0_jerdn = mT(lep_p4[idx0], MET_jerdn);
    float MT1_jerdn = mT(lep_p4[idx1], MET_jerdn);

    tx->setBranch<float>("MTmax", MT0 > MT1 ? MT0 : MT1);
    tx->setBranch<float>("MTmax_up", MT0_up > MT1_up ? MT0_up : MT1_up);
    tx->setBranch<float>("MTmax_dn", MT0_dn > MT1_dn ? MT0_dn : MT1_dn);
    tx->setBranch<float>("MTmax_gen", MT0_gen > MT1_gen ? MT0_gen : MT1_gen);
    tx->setBranch<float>("MTmax_jer", MT0_jer > MT1_jer ? MT0_jer : MT1_jer);
    tx->setBranch<float>("MTmax_jerup", MT0_jerup > MT1_jerup ? MT0_jerup : MT1_jerup);
    tx->setBranch<float>("MTmax_jerdn", MT0_jerdn > MT1_jerdn ? MT0_jerdn : MT1_jerdn);
    tx->setBranch<float>("MTmin", MT0 < MT1 ? MT0 : MT1);
    tx->setBranch<float>("MTmin_up", MT0_up < MT1_up ? MT0_up : MT1_up);
    tx->setBranch<float>("MTmin_dn", MT0_dn < MT1_dn ? MT0_dn : MT1_dn);
    tx->setBranch<float>("MTmin_gen", MT0_gen < MT1_gen ? MT0_gen : MT1_gen);
    tx->setBranch<float>("MTmin_jer", MT0_jer < MT1_jer ? MT0_jer : MT1_jer);
    tx->setBranch<float>("MTmin_jerup", MT0_jerup < MT1_jerup ? MT0_jerup : MT1_jerup);
    tx->setBranch<float>("MTmin_jerdn", MT0_jerdn < MT1_jerdn ? MT0_jerdn : MT1_jerdn);
}

//##############################################################################################################
void babyMaker_v2::Fill3LLeptonVariables()
{
    int nSFOS = babyMaker_v2::nSFOS();
    int nSFOSinZ = 0;
    if (nSFOS == 0)
    {
        tx->setBranch<float>("Mll3L", get0SFOSMll());
        tx->setBranch<float>("Mee3L", get0SFOSMee());
    }
    else if (nSFOS == 1)
    {
        tx->setBranch<float>("Mll3L", get1SFOSMll());
        if (fabs(tx->getBranch<float>("Mll3L") - 91.1876) < 20.)
            nSFOSinZ = 1;
    }
    else if (nSFOS == 2)
    {
        tx->setBranch<float>("Mll3L", get2SFOSMll0());
        tx->setBranch<float>("Mll3L1", get2SFOSMll1());
        if (fabs(tx->getBranch<float>("Mll3L") - 91.1876) < 20.)
            nSFOSinZ++;
        if (fabs(tx->getBranch<float>("Mll3L1") - 91.1876) < 20.)
            nSFOSinZ++;
    }
    tx->setBranch<int>("nSFOS", nSFOS);
    tx->setBranch<int>("nSFOSinZ", nSFOSinZ);

    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);

    // Extra variables related to the three leptons
    tx->setBranch<float>("M3l", (lep_p4[0] + lep_p4[1] + lep_p4[2]).mass());
    tx->setBranch<float>("Pt3l", (lep_p4[0] + lep_p4[1] + lep_p4[2]).pt());
    tx->setBranch<float>("M01", (lep_p4[0] + lep_p4[1]).mass());
    tx->setBranch<float>("M02", (lep_p4[0] + lep_p4[2]).mass());
    tx->setBranch<float>("M12", (lep_p4[1] + lep_p4[2]).mass());
    tx->setBranch<int>("isSFOS01", lep_pdgId[0] == -lep_pdgId[1]);
    tx->setBranch<int>("isSFOS02", lep_pdgId[0] == -lep_pdgId[2]);
    tx->setBranch<int>("isSFOS12", lep_pdgId[1] == -lep_pdgId[2]);

    // Set up MET variables
    const float& met_pt = tx->getBranch<float>("met_pt");
    const float& met_phi = tx->getBranch<float>("met_phi");
    const float& met_up_pt = tx->getBranch<float>("met_up_pt");
    const float& met_up_phi = tx->getBranch<float>("met_up_phi");
    const float& met_dn_pt = tx->getBranch<float>("met_dn_pt");
    const float& met_dn_phi = tx->getBranch<float>("met_dn_phi");
    const float& met_gen_pt = tx->getBranch<float>("met_gen_pt");
    const float& met_gen_phi = tx->getBranch<float>("met_gen_phi");
    const float& met_jer_pt = tx->getBranch<float>("met_jer_pt");
    const float& met_jer_phi = tx->getBranch<float>("met_jer_phi");
    const float& met_jerup_pt = tx->getBranch<float>("met_jerup_pt");
    const float& met_jerup_phi = tx->getBranch<float>("met_jerup_phi");
    const float& met_jerdn_pt = tx->getBranch<float>("met_jerdn_pt");
    const float& met_jerdn_phi = tx->getBranch<float>("met_jerdn_phi");
    LV MET;
    LV MET_up;
    LV MET_dn;
    LV MET_gen;
    LV MET_jer;
    LV MET_jerup;
    LV MET_jerdn;
    MET.SetPxPyPzE( met_pt * TMath::Cos(met_phi), met_pt * TMath::Sin(met_phi), 0, met_pt);
    MET_up.SetPxPyPzE( met_up_pt * TMath::Cos(met_up_phi), met_up_pt * TMath::Sin(met_up_phi), 0, met_up_pt);
    MET_dn.SetPxPyPzE( met_dn_pt * TMath::Cos(met_dn_phi), met_dn_pt * TMath::Sin(met_dn_phi), 0, met_dn_pt);
    MET_gen.SetPxPyPzE(met_gen_pt * TMath::Cos(met_gen_phi), met_gen_pt * TMath::Sin(met_gen_phi), 0, met_gen_pt);
    MET_jer.SetPxPyPzE(met_jer_pt * TMath::Cos(met_jer_phi), met_jer_pt * TMath::Sin(met_jer_phi), 0, met_jer_pt);
    MET_jerup.SetPxPyPzE(met_jerup_pt * TMath::Cos(met_jerup_phi), met_jerup_pt * TMath::Sin(met_jerup_phi), 0, met_jerup_pt);
    MET_jerdn.SetPxPyPzE(met_jerdn_pt * TMath::Cos(met_jerdn_phi), met_jerdn_pt * TMath::Sin(met_jerdn_phi), 0, met_jerdn_pt);

    // Set MET lep related variables
    tx->setBranch<float>("DPhi3lMET"   , fabs(ROOT::Math::VectorUtil::DeltaPhi((lep_p4[0] + lep_p4[1] + lep_p4[2]), MET   )));
    tx->setBranch<float>("DPhi3lMET_up", fabs(ROOT::Math::VectorUtil::DeltaPhi((lep_p4[0] + lep_p4[1] + lep_p4[2]), MET_up)));
    tx->setBranch<float>("DPhi3lMET_dn", fabs(ROOT::Math::VectorUtil::DeltaPhi((lep_p4[0] + lep_p4[1] + lep_p4[2]), MET_dn)));
    tx->setBranch<float>("DPhi3lMET_gen", fabs(ROOT::Math::VectorUtil::DeltaPhi((lep_p4[0] + lep_p4[1] + lep_p4[2]), MET_gen)));
    tx->setBranch<float>("DPhi3lMET_jer", fabs(ROOT::Math::VectorUtil::DeltaPhi((lep_p4[0] + lep_p4[1] + lep_p4[2]), MET_jer)));
    tx->setBranch<float>("DPhi3lMET_jerup", fabs(ROOT::Math::VectorUtil::DeltaPhi((lep_p4[0] + lep_p4[1] + lep_p4[2]), MET_jerup)));
    tx->setBranch<float>("DPhi3lMET_jerdn", fabs(ROOT::Math::VectorUtil::DeltaPhi((lep_p4[0] + lep_p4[1] + lep_p4[2]), MET_jerdn)));

    if (nSFOS == 1)
    {
        int idx = -1;
        if (lep_pdgId[0] == -lep_pdgId[1]) idx = 2;
        if (lep_pdgId[0] == -lep_pdgId[2]) idx = 1;
        if (lep_pdgId[1] == -lep_pdgId[2]) idx = 0;
        tx->setBranch<float>("MT3rd"   , mT(lep_p4[idx], MET   ));
        tx->setBranch<float>("MT3rd_up", mT(lep_p4[idx], MET_up));
        tx->setBranch<float>("MT3rd_dn", mT(lep_p4[idx], MET_dn));
        tx->setBranch<float>("MT3rd_gen", mT(lep_p4[idx], MET_gen));
        tx->setBranch<float>("MT3rd_jer", mT(lep_p4[idx], MET_jer));
        tx->setBranch<float>("MT3rd_jerup", mT(lep_p4[idx], MET_jerup));
        tx->setBranch<float>("MT3rd_jerdn", mT(lep_p4[idx], MET_jerdn));
    }

    float MT0    = mT(lep_p4[0], MET);
    float MT1    = mT(lep_p4[1], MET);
    float MT2    = mT(lep_p4[2], MET);
    float MT0_up = mT(lep_p4[0], MET_up);
    float MT1_up = mT(lep_p4[1], MET_up);
    float MT2_up = mT(lep_p4[2], MET_up);
    float MT0_dn = mT(lep_p4[0], MET_dn);
    float MT1_dn = mT(lep_p4[1], MET_dn);
    float MT2_dn = mT(lep_p4[2], MET_dn);
    float MT0_gen = mT(lep_p4[0], MET_gen);
    float MT1_gen = mT(lep_p4[1], MET_gen);
    float MT2_gen = mT(lep_p4[2], MET_gen);
    float MT0_jer = mT(lep_p4[0], MET_jer);
    float MT1_jer = mT(lep_p4[1], MET_jer);
    float MT2_jer = mT(lep_p4[2], MET_jer);
    float MT0_jerup = mT(lep_p4[0], MET_jerup);
    float MT1_jerup = mT(lep_p4[1], MET_jerup);
    float MT2_jerup = mT(lep_p4[2], MET_jerup);
    float MT0_jerdn = mT(lep_p4[0], MET_jerdn);
    float MT1_jerdn = mT(lep_p4[1], MET_jerdn);
    float MT2_jerdn = mT(lep_p4[2], MET_jerdn);

    tx->setBranch<float>("MTmax3L", TMath::Max(MT0, TMath::Max(MT1, MT2)));
    tx->setBranch<float>("MTmax3L_up", TMath::Max(MT0_up, TMath::Max(MT1_up, MT2_up)));
    tx->setBranch<float>("MTmax3L_dn", TMath::Max(MT0_dn, TMath::Max(MT1_dn, MT2_dn)));
    tx->setBranch<float>("MTmax3L_gen", TMath::Max(MT0_gen, TMath::Max(MT1_gen, MT2_gen)));
    tx->setBranch<float>("MTmax3L_jer", TMath::Max(MT0_jer, TMath::Max(MT1_jer, MT2_jer)));
    tx->setBranch<float>("MTmax3L_jerup", TMath::Max(MT0_jerup, TMath::Max(MT1_jerup, MT2_jerup)));
    tx->setBranch<float>("MTmax3L_jerdn", TMath::Max(MT0_jerdn, TMath::Max(MT1_jerdn, MT2_jerdn)));
}

//##############################################################################################################
void babyMaker_v2::FillEventTags()
{
    tx->setBranch<TString>("bkgtype", process());
    tx->setBranch<int>("vetophoton", vetophotonprocess());
}

//##############################################################################################################
void babyMaker_v2::FillWeights()
{

    if (tx->getBranch<int>("nLlep") < 2)
        return;

    if (cms3.evt_isRealData())
    {
        tx->setBranch<float>("purewgt", 1);
        tx->setBranch<float>("purewgt_up", 1);
        tx->setBranch<float>("purewgt_dn", 1);
    }
    else
    {
        int nTrueInt = cms3.puInfo_trueNumInteractions().at(0);
        tx->setBranch<float>("purewgt", puWeight(nTrueInt));
        tx->setBranch<float>("purewgt_up", puWeight(nTrueInt, 1));
        tx->setBranch<float>("purewgt_dn", puWeight(nTrueInt, -1));
    }

    // fakerate
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    float ffwgt;
    float ffwgterr;
    int loose_lep_idx;
    int loose_lep_idx_qcd;
    float ffwgtqcd;
    float ffwgtqcderr;
    std::tie(ffwgt, ffwgterr, loose_lep_idx) = getlepFakeRateandErrorandLooseLepIdx(true);
    std::tie(ffwgtqcd, ffwgtqcderr, loose_lep_idx_qcd) = getlepFakeRateandErrorandLooseLepIdx(false);
    // The closure is taken from the statistical uncertainty of the closure test
    // The closure central values are better than 27% and 28%
    // electron closure fSumw[1]=1, x=0.5, error=0.279466
    // muon closure     fSumw[1]=1, x=0.5, error=0.268032
    const double el_closure = 0.279466;
    const double mu_closure = 0.268032;
    tx->setBranch<float>("ffwgt", ffwgt);
    tx->setBranch<float>("ffwgt_up", ffwgt + ffwgterr);
    tx->setBranch<float>("ffwgt_dn", ffwgt - ffwgterr);
    tx->setBranch<float>("ffwgt_el_up", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt + ffwgterr : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_el_dn", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt - ffwgterr : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_mu_up", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 13 ? ffwgt + ffwgterr : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_mu_dn", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 13 ? ffwgt - ffwgterr : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_closure_el_up", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt + (ffwgt * el_closure) : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_closure_el_dn", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt - (ffwgt * el_closure) : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_closure_mu_up", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 13 ? ffwgt + (ffwgt * mu_closure) : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_closure_mu_dn", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 13 ? ffwgt - (ffwgt * mu_closure) : ffwgt) : 0);
    tx->setBranch<float>("ffwgt_closure_up", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt + (ffwgt * el_closure) : ffwgt + (ffwgt * mu_closure)) : 0);
    tx->setBranch<float>("ffwgt_closure_dn", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt - (ffwgt * el_closure) : ffwgt - (ffwgt * mu_closure)) : 0);
    tx->setBranch<float>("ffwgt_full_up", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt + sqrt(pow(ffwgt * el_closure, 2) + pow(ffwgterr, 2)) : ffwgt + sqrt(pow(ffwgt * mu_closure, 2) + pow(ffwgterr, 2))) : 0);
    tx->setBranch<float>("ffwgt_full_dn", loose_lep_idx >= 0 ? (abs(lep_pdgId[loose_lep_idx]) == 11 ? ffwgt - sqrt(pow(ffwgt * el_closure, 2) + pow(ffwgterr, 2)) : ffwgt - sqrt(pow(ffwgt * mu_closure, 2) + pow(ffwgterr, 2))) : 0);
    tx->setBranch<float>("ffwgtqcd", ffwgtqcd);
    tx->setBranch<float>("ffwgtqcd_up", ffwgtqcd + ffwgtqcderr);
    tx->setBranch<float>("ffwgtqcd_dn", ffwgtqcd - ffwgtqcderr);

    // lepsf
    if (!cms3.evt_isRealData())
    {
        float lepsf;
        float lepsferr;
        std::tie(lepsf, lepsferr) = getlepSFWeightandError();
        tx->setBranch<float>("lepsf", lepsf);
        tx->setBranch<float>("lepsf_up", lepsf + lepsferr);
        tx->setBranch<float>("lepsf_dn", lepsf - lepsferr);
    }
    else
    {
        tx->setBranch<float>("lepsf", 1);
        tx->setBranch<float>("lepsf_up", 1);
        tx->setBranch<float>("lepsf_dn", 1);
    }

    // trigger eff
    if (!cms3.evt_isRealData())
    {
//        float trigeff;
//        float trigefferr;
//        std::tie(trigeff, trigefferr) = getTrigEffandError();
//        tx->setBranch<float>("trigeff", trigeff);
//        tx->setBranch<float>("trigeff_up", trigeff + trigefferr);
//        tx->setBranch<float>("trigeff_dn", trigeff - trigefferr);
        tx->setBranch<float>("trigeff", 1);
        tx->setBranch<float>("trigeff_up", 1);
        tx->setBranch<float>("trigeff_dn", 1);
        float trigsf;
        float trigsferr;
        std::tie(trigsf, trigsferr) = getTrigSFandError();
        tx->setBranch<float>("trigsf", trigsf);
        tx->setBranch<float>("trigsf_up", trigsf + trigsferr);
        tx->setBranch<float>("trigsf_dn", trigsf - trigsferr);
    }
    else
    {
        tx->setBranch<float>("trigeff", 1);
        tx->setBranch<float>("trigeff_up", 1);
        tx->setBranch<float>("trigeff_dn", 1);
        tx->setBranch<float>("trigsf", 1);
        tx->setBranch<float>("trigsf_up", 1);
        tx->setBranch<float>("trigsf_dn", 1);
    }

//    // Smearing weights
//    if (cms3.evt_isRealData())
//    {
//        double dtSF = 1;
//        for (unsigned int i = 0; i < tx->getBranch<vector<int>>("lep_pdgId").size(); ++i)
//        {
//            int Q = tx->getBranch<vector<int>>("lep_charge")[i];
//            double pt = tx->getBranch<vector<float>>("lep_pt")[i];
//            double eta = tx->getBranch<vector<float>>("lep_eta")[i];
//            double phi = tx->getBranch<vector<float>>("lep_phi")[i];
//            dtSF *= rc.kScaleDT(Q, pt, eta, phi);
//        }
//        tx->setBranch<float>("musmear_sf", dtSF);
//    }
//    else
//    {
//        double mcSF = 1;
//        for (unsigned int i = 0; i < tx->getBranch<vector<int>>("lep_pdgId").size(); ++i)
//        {
//            int Q = tx->getBranch<vector<int>>("lep_charge")[i];
//            double pt = tx->getBranch<vector<float>>("lep_pt")[i];
//            double eta = tx->getBranch<vector<float>>("lep_eta")[i];
//            double phi = tx->getBranch<vector<float>>("lep_phi")[i];
//            int nl = tx->getBranch<vector<int>>("lep_nlayers")[i];
//            rnd.SetSeed(cms3.evt_event());
//            double u1 = TMath::Min(rnd.Rndm(), 1.0);
//            double u2 = TMath::Min(rnd.Rndm(), 1.0);
//            mcSF *= rc.kScaleAndSmearMC(Q, pt, eta, phi, nl, u1, u2);
//        }
//        tx->setBranch<float>("musmear_sf", mcSF);
//    }
}

//##############################################################################################################
void babyMaker_v2::FillGenWeights()
{
    if (cms3.evt_isRealData())
        return;

    // scale pdf variation
    if (isWHSUSY())
        FillGenWeightsForWHSUSY();
    else
        FillGenWeightsNominal();
}

//##############################################################################################################
void babyMaker_v2::FillGenWeightsForWHSUSY()
{
    using namespace tas;
    setWHSMSMass();
    const float& mass_chargino = tx->getBranch<float>("chimass");
    const float& mass_lsp = tx->getBranch<float>("lspmass");
    if (mass_chargino < 0 || mass_lsp < 0)
        FATALERROR(__FUNCTION__);

    if (genps_weight() > 0) h_nrawevents_SMS -> Fill(mass_chargino, mass_lsp, 1);
    else if (genps_weight() < 0) h_nrawevents_SMS -> Fill(mass_chargino, mass_lsp, -1);
    float xsec = hxsec->GetBinContent(hxsec->FindBin(mass_chargino));
    float xsec_uncert = hxsec->GetBinError(hxsec->FindBin(mass_chargino));

    float SMSpdf_weight_up = 1;
    float SMSpdf_weight_down = 1;
    float SMSsum_of_weights = 0;
    float SMSaverage_of_weights = 0;
    //error on pdf replicas
    //fastsim has first genweights bin being ==1
    if (genweights().size() > 111) { //fix segfault
        for (int ipdf = 10; ipdf < 110; ipdf++) {
            SMSaverage_of_weights += cms3.genweights().at(ipdf);
        } // average of weights
        SMSaverage_of_weights = SMSaverage_of_weights / 100.;
        for (int ipdf = 10; ipdf < 110; ipdf++) {
            SMSsum_of_weights += pow(cms3.genweights().at(ipdf) - SMSaverage_of_weights, 2);
        } //std of weights.
        SMSpdf_weight_up = (SMSaverage_of_weights + sqrt(SMSsum_of_weights / 99.));
        SMSpdf_weight_down = (SMSaverage_of_weights - sqrt(SMSsum_of_weights / 99.));
        // NOTE: The SMS samples have a shifted index (instead of 0 being fr_r1_f1 it is 1)
        tx->setBranch<float>("weight_fr_r1_f1", cms3.genweights()[1]);
        tx->setBranch<float>("weight_fr_r1_f2", cms3.genweights()[2]);
        tx->setBranch<float>("weight_fr_r1_f0p5", cms3.genweights()[3]);
        tx->setBranch<float>("weight_fr_r2_f1", cms3.genweights()[4]);
        tx->setBranch<float>("weight_fr_r2_f2", cms3.genweights()[5]);
        tx->setBranch<float>("weight_fr_r2_f0p5", cms3.genweights()[6]);
        tx->setBranch<float>("weight_fr_r0p5_f1", cms3.genweights()[7]);
        tx->setBranch<float>("weight_fr_r0p5_f2", cms3.genweights()[8]);
        tx->setBranch<float>("weight_fr_r0p5_f0p5", cms3.genweights()[9]);
        tx->setBranch<float>("weight_pdf_up", SMSpdf_weight_up);
        tx->setBranch<float>("weight_pdf_down", SMSpdf_weight_up);
        tx->setBranch<float>("weight_alphas_down", cms3.genweights()[110]);
        tx->setBranch<float>("weight_alphas_up", cms3.genweights()[111]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 1, genweights()[1]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 2, genweights()[2]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 3, genweights()[3]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 4, genweights()[4]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 5, genweights()[5]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 6, genweights()[6]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 7, genweights()[7]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 8, genweights()[8]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 9, genweights()[9]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 10, SMSpdf_weight_up);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 11, SMSpdf_weight_down);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 12, genweights()[110]); // α_s variation.
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 13, genweights()[111]); // α_s variation.

        // ISR weights
        LV hardsystem;
        for (unsigned int genx = 0; genx < genps_p4().size(); genx++)
        {
            if ((abs(genps_id().at(genx)) == 1000024 || abs(genps_id().at(genx)) == 1000023) && genps_isLastCopy().at(genx))
                hardsystem += genps_p4().at(genx);
        }

        float weight_ISR = 0;
        float weight_ISRup = 0;
        float weight_ISRdown = 0;
        if (hardsystem.pt() > 600.)
            weight_ISR = .783;
        else if (hardsystem.pt() > 400.)
            weight_ISR = .912;
        else if (hardsystem.pt() > 300.)
            weight_ISR = 1.;
        else if (hardsystem.pt() > 200.)
            weight_ISR = 1.057;
        else if (hardsystem.pt() > 150.)
            weight_ISR = 1.150;
        else if (hardsystem.pt() > 100.)
            weight_ISR = 1.179;
        else if (hardsystem.pt() > 50.)
            weight_ISR = 1.052;
        else
            weight_ISR = 1.;
        float isrweight_unc = fabs(weight_ISR - 1);
        weight_ISRup = 1. + isrweight_unc;
        weight_ISRdown = 1. - isrweight_unc;
        tx->setBranch<float>("weight_isr", weight_ISR);
        tx->setBranch<float>("weight_isr_up", weight_ISRup);
        tx->setBranch<float>("weight_isr_down", weight_ISRdown);
        h_nevents_SMS->Fill(mass_chargino, mass_lsp, 19, weight_ISR);
        h_nevents_SMS->Fill(mass_chargino, mass_lsp, 20, weight_ISRup);
        h_nevents_SMS->Fill(mass_chargino, mass_lsp, 21, weight_ISRdown);
    }
}

//##############################################################################################################
void babyMaker_v2::FillGenWeightsNominal()
{
    float sum_of_pdf_weights = 0;
    float average_of_pdf_weights = 0;
    //error on pdf replicas
    if (cms3.genweights().size() > 110)
    {
        // average of weights
        for (int ipdf = 9; ipdf < 109; ipdf++)
        {
            average_of_pdf_weights += cms3.genweights().at(ipdf);
        }
        average_of_pdf_weights =  average_of_pdf_weights / 100;
        //std of weights.
        for (int ipdf = 9; ipdf < 109; ipdf++)
        {
            sum_of_pdf_weights += (cms3.genweights().at(ipdf) - average_of_pdf_weights) * (cms3.genweights().at(ipdf) - average_of_pdf_weights);
        }
        tx->setBranch<float>("weight_fr_r1_f1", cms3.genweights()[0]);
        tx->setBranch<float>("weight_fr_r1_f2", cms3.genweights()[1]);
        tx->setBranch<float>("weight_fr_r1_f0p5", cms3.genweights()[2]);
        tx->setBranch<float>("weight_fr_r2_f1", cms3.genweights()[3]);
        tx->setBranch<float>("weight_fr_r2_f2", cms3.genweights()[4]);
        tx->setBranch<float>("weight_fr_r2_f0p5", cms3.genweights()[5]);
        tx->setBranch<float>("weight_fr_r0p5_f1", cms3.genweights()[6]);
        tx->setBranch<float>("weight_fr_r0p5_f2", cms3.genweights()[7]);
        tx->setBranch<float>("weight_fr_r0p5_f0p5", cms3.genweights()[8]);
        tx->setBranch<float>("weight_pdf_up", (average_of_pdf_weights + sqrt(sum_of_pdf_weights / 99)));
        tx->setBranch<float>("weight_pdf_down", (average_of_pdf_weights - sqrt(sum_of_pdf_weights / 99)));
        tx->setBranch<float>("weight_alphas_down", cms3.genweights()[109]);
        tx->setBranch<float>("weight_alphas_up", cms3.genweights()[110]);
        h_neventsinfile->Fill(1, tx->getBranch<float>("weight_fr_r1_f1"));
        h_neventsinfile->Fill(2, tx->getBranch<float>("weight_fr_r1_f2"));
        h_neventsinfile->Fill(3, tx->getBranch<float>("weight_fr_r1_f0p5"));
        h_neventsinfile->Fill(4, tx->getBranch<float>("weight_fr_r2_f1"));
        h_neventsinfile->Fill(5, tx->getBranch<float>("weight_fr_r2_f2"));
        h_neventsinfile->Fill(6, tx->getBranch<float>("weight_fr_r2_f0p5"));
        h_neventsinfile->Fill(7, tx->getBranch<float>("weight_fr_r0p5_f1"));
        h_neventsinfile->Fill(8, tx->getBranch<float>("weight_fr_r0p5_f2"));
        h_neventsinfile->Fill(9, tx->getBranch<float>("weight_fr_r0p5_f0p5"));
        h_neventsinfile->Fill(10, tx->getBranch<float>("weight_pdf_up"));
        h_neventsinfile->Fill(11, tx->getBranch<float>("weight_pdf_down"));
        h_neventsinfile->Fill(12, tx->getBranch<float>("weight_alphas_down"));
        h_neventsinfile->Fill(13, tx->getBranch<float>("weight_alphas_up"));
        h_neventsinfile->Fill(14, 1);
        h_neventsinfile->Fill(15, (cms3.genps_weight() > 0) - (cms3.genps_weight() < 0));
    }
}

//##############################################################################################################
float babyMaker_v2::get0SFOSMll()
{
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);

    int pdgid0 = lep_pdgId[0];
    int pdgid1 = lep_pdgId[1];
    int pdgid2 = lep_pdgId[2];

    if (pdgid0 == pdgid1)
        return (lep_p4[0] + lep_p4[1]).mass();
    else if (pdgid0 == pdgid2)
        return (lep_p4[0] + lep_p4[2]).mass();
    else if (pdgid1 == pdgid2)
        return (lep_p4[1] + lep_p4[2]).mass();

    cout << "Warning: get0SFOSMll: Shouldn't be here if function call are at the right places." << endl;

    return -999;
}

//##############################################################################################################
float babyMaker_v2::get0SFOSMee()
{
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);

    int pdgid0 = lep_pdgId[0];
    int pdgid1 = lep_pdgId[1];
    int pdgid2 = lep_pdgId[2];

    if (pdgid0 == pdgid1 && abs(pdgid0) == 11)
        return (lep_p4[0] + lep_p4[1]).mass();
    else if (pdgid0 == pdgid2 && abs(pdgid0) == 11)
        return (lep_p4[0] + lep_p4[2]).mass();
    else if (pdgid1 == pdgid2 && abs(pdgid1) == 11)
        return (lep_p4[1] + lep_p4[2]).mass();
    return -999;
}

//##############################################################################################################
float babyMaker_v2::get1SFOSMll()
{
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);

    int pdgid0 = lep_pdgId[0];
    int pdgid1 = lep_pdgId[1];
    int pdgid2 = lep_pdgId[2];

    if (pdgid0 == -pdgid1)
        return (lep_p4[0] + lep_p4[1]).mass();
    else if (pdgid0 == -pdgid2)
        return (lep_p4[0] + lep_p4[2]).mass();
    else if (pdgid1 == -pdgid2)
        return (lep_p4[1] + lep_p4[2]).mass();

    cout << "Warning: get1SFOSMll() Shouldn't be here if function call are at the right places." << endl;
    std::cout <<  " pdgid0: " << pdgid0 <<  " pdgid1: " << pdgid1 <<  " pdgid2: " << pdgid2 <<  std::endl;
    std::cout <<  " nSFOS(): " << nSFOS() <<  std::endl;

    for (auto& pdgid0 : tx->getBranch<vector<int>>("lep_pdgId"))
    {
        std::cout <<  " pdgid0: " << pdgid0 <<  std::endl;
    }

    int nSFOS = 0;
    for (auto& pdgid0 : tx->getBranch<vector<int>>("lep_pdgId"))
    {
        for (auto& pdgid1 : tx->getBranch<vector<int>>("lep_pdgId"))
        {
            std::cout <<  " pdgid0: " << pdgid0 <<  " pdgid1: " << pdgid1 <<  std::endl;
            if (pdgid0 == -pdgid1)
                nSFOS++;
        }
    }
    std::cout <<  " nSFOS: " << nSFOS <<  " nSFOS/2: " << nSFOS/2 <<  std::endl;

    return -999;
}

//##############################################################################################################
float babyMaker_v2::get2SFOSMll0()
{
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);

    int pdgid0 = lep_pdgId[0];
    int pdgid1 = lep_pdgId[1];
    int pdgid2 = lep_pdgId[2];

    if (pdgid0 == -pdgid1)
        return (lep_p4[0] + lep_p4[1]).mass();
    else if (pdgid0 == -pdgid2)
        return (lep_p4[0] + lep_p4[2]).mass();
    else if (pdgid1 == -pdgid2)
        return (lep_p4[1] + lep_p4[2]).mass();

    cout << "Warning: get2SFOSMll0() Shouldn't be here if function call are at the right places." << endl;

    return -999;
}

//##############################################################################################################
float babyMaker_v2::get2SFOSMll1()
{
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);

    int pdgid0 = lep_pdgId[0];
    int pdgid1 = lep_pdgId[1];
    int pdgid2 = lep_pdgId[2];

    if (pdgid2 == -pdgid1)
        return (lep_p4[1] + lep_p4[2]).mass();
    else if (pdgid0 == -pdgid2)
        return (lep_p4[0] + lep_p4[2]).mass();
    else if (pdgid1 == -pdgid0)
        return (lep_p4[0] + lep_p4[1]).mass();

    cout << "Warning: get2SFOSMll1() Shouldn't be here if function call are at the right places." << endl;

    return -999;
}

//##############################################################################################################
float babyMaker_v2::mT(LV p4, LV met)
{
    float phi1 = p4.Phi();
    float phi2 = met.Phi();
    float Et1  = p4.Et();
    float Et2  = met.Et();
    return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

//##############################################################################################################
// Returns true + indices of the SSID passing leptons that makes a SS pair.
tuple<bool, int, int> babyMaker_v2::isSSCR()
{
    // Assumes that the FillElectron/FillMuon is done
    // Retrieve some info we need from the lepton containers
    const vector<int>& lep_pdgId  = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<int>& lep_charge = tx->getBranch<vector<int>>("lep_charge", true);
    const vector<int>& lep_tight  = tx->getBranch<vector<int>>(gconf.wwwcfg["tightid"]);

    // Count the number of tight "SSID" leptons
    int ntight = passCount(lep_tight);

    // If there are no two "SSID" passing leptons then return false, and -1 indices
    if (ntight < 2)
    {
        return make_tuple(false, -1, -1);
    }
    else if (ntight >= 2)
    {
        // Determine whether I have two positive or two negative
        // Assumes abs charge == 1 ensured
        int totalcharge = lep_charge[0] + lep_charge[1] + lep_charge[2];

        // Save the indices here
        vector<int> lep_index;

        // Loop over leptons and pick the two indices that passes tight ID and are "same-sign"
        for (unsigned int i = 0; i < lep_pdgId.size(); ++i)
        {
            if (totalcharge == 1 ? lep_pdgId[i] < 0 && lep_tight[i] : lep_pdgId[i] > 0 && lep_tight[i])
                lep_index.push_back(i);
        }

        // If the lep_index does not have two leptons
        // this event did not pass the SSCR requirement from the 3L region
        if (lep_index.size() != 2)
            return make_tuple(false, -1, -1);
        else
            return make_tuple(true, lep_index[0], lep_index[1]);
    }
    else
    {
        FATALERROR(__FUNCTION__);
        return make_tuple(false, -1, -1);
    }
}

//##############################################################################################################
TString babyMaker_v2::process()
{
    if (coreSample.is2016(looper.getCurrentFileName()))
    {
        if (cms3.evt_isRealData()) return "Data";
        if (splitVH()) return "WHtoWWW";
        if (isSMWWW()) return "WWW";
        if (isSMWWW() && coreSample.is2017(looper.getCurrentFileName())) return "WWW";
        if (isSMWWW()) return "WWWv2";
        if (isData()) return "Data";
        if (tx->getBranch<int>("nVlep") == 2 && tx->getBranch<int>("nLlep") == 2)
        {
            if (tx->getBranch<int>("nLlep") < 2) return "not2l";
            int gentype = gentype_v2();
            if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
            {
                std::cout <<  " gentype: " << gentype <<  std::endl;
            }
            if      (gentype == 0) return "trueSS";
            else if (gentype == 2) return "chargeflips";
            else if (gentype == 3) return "SSLL";
            else if (gentype == 4) return "fakes";
            else if (gentype == 5) return "photonfakes";
            else                   return "others";
        }
        //this is 3l
        if (tx->getBranch<int>("nLlep") < 3) return "not3l";
        int gentype = gentype_v2();
        if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
        {
            std::cout <<  " gentype: " << gentype <<  std::endl;
        }
        if      (gentype == 0) return "trueWWW";
        else if (gentype == 1) return "true3L";
        else if (gentype == 2) return "chargeflips";
        else if (gentype == 3) return "3lLL";
        else if (gentype == 4) return "fakes";
        else if (gentype == 5) return "photonfakes";
        else                   return "others";
    }
    else if (coreSample.is2017(looper.getCurrentFileName()))
    {
        if (cms3.evt_isRealData())
            return "Data";

        if (splitVH())
            return "WHtoWWW";

        if (isSMWWW())
            return "WWW";

        if (isSMWWW() && coreSample.is2017(looper.getCurrentFileName()))
            return "WWW";

        if (isData())
            return "Data";

        return gentype_v3();
    }
}

//##############################################################################################################
bool babyMaker_v2::splitVH()
{
    if (!isVH()) return false; //file is certainly no WHtoWWW
    bool isHtoWW = false;
    bool isWnotFromH = false;
    bool isZthere = false;
    const vector<int>& genPart_pdgId = coreGenPart.genPart_pdgId;
    const vector<int>& genPart_status = coreGenPart.genPart_status;
    const vector<int>& genPart_motherId = coreGenPart.genPart_motherId;
    for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
    {

//        if (genPart_status[i] == 22)
//            std::cout <<  " genPart_pdgId[i]: " << genPart_pdgId[i] <<  " genPart_status[i]: " << genPart_status[i] <<  " genPart_motherId[i]: " << genPart_motherId[i] <<  std::endl;

        if (abs(genPart_pdgId[i]) == 23 && genPart_status[i] == 22) isZthere = true;

        if (abs(genPart_pdgId[i]) != 24) continue;
        if (genPart_status[i] != 22) continue;
        if (abs(genPart_motherId[i]) == 25)  isHtoWW     = true;
        if (abs(genPart_motherId[i]) != 2)   isWnotFromH = true;
        if (isHtoWW && isWnotFromH && !isZthere) break;
    }
    return isHtoWW && isWnotFromH && !isZthere;
}

//##############################################################################################################
bool babyMaker_v2::studyDoublyChargedHiggs()
{
    ProcessGenParticles();
    if (!isDoublyChargedHiggs()) return false; //file is certainly not doubly charged higgs
    const vector<int>& genPart_pdgId = coreGenPart.genPart_pdgId;
    const vector<int>& genPart_status = coreGenPart.genPart_status;
    const vector<int>& genPart_motherId = coreGenPart.genPart_motherId;
    const vector<int>& genPart_grandmaId = coreGenPart.genPart_grandmaId;
    const vector<LV>& genPart_p4 = coreGenPart.genPart_p4;
    bool iswhpmpm = false;
    int mother_of_w = 0;
    int wid = 0;
    for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
    {
        int pdgId = genPart_pdgId[i];
        int status = genPart_status[i];
        int motherId = genPart_motherId[i];
        int grandmaId = genPart_grandmaId[i];
        if (status < 21 || status > 23)
            continue;
        if (grandmaId != 2212)
            continue;
        if (abs(pdgId) != 24)
            continue;
        iswhpmpm = true;
        mother_of_w = motherId;
        wid = pdgId;
    }

    if (!iswhpmpm)
    {
        tx->clear();
        tx->setBranch<int>("iswhchannel", 0);
        return false;
    }

    vector<LV> qs;
    vector<LV> h;
    vector<LV> w;
    vector<LV> ws;
    vector<LV> l;
    vector<LV> v;

    std::cout.setstate(std::ios_base::failbit);
    std::cout << std::endl;
    int hid = 0;
    for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
    {
        int pdgId = genPart_pdgId[i];
        int status = genPart_status[i];
        int motherId = genPart_motherId[i];
        int grandmaId = genPart_grandmaId[i];
        LV p4 = genPart_p4[i];
        // The quarks from the associated W decay
        if (status == 23 && grandmaId == wid && motherId == wid && (abs(pdgId) >= 1 && abs(pdgId) <= 4))
        {
            qs.push_back(p4);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        else if (status == 22 && grandmaId == 2212 && abs(pdgId) == 255)
        {
            h.push_back(p4);
            hid = pdgId;
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        else if (status == 22 && grandmaId == 2212 && abs(pdgId) == 24)
        {
            w.push_back(p4);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        // Assume that in the particle record doubly charged higgs will always show up before product
        else if (abs(pdgId) == 24 && motherId == hid)
        {
            ws.push_back(p4);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        // Assume that in the particle record doubly charged higgs will always show up before product
        else if (abs(motherId) == 24 && grandmaId == hid && (abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15))
        {
            l.push_back(p4);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        else if (abs(motherId) == 24 && grandmaId == hid && (abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16))
        {
            v.push_back(p4);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
    }
    std::cout.clear();

    // sanity check
    if (qs.size() != 2)
        std::cout << "Found different from 2 quarks " << qs.size() << std::endl;
    if (h.size() != 1)
        std::cout << "Found different from 1 doubly charged higgs " << h.size() << std::endl;
    if (w.size() != 1)
        std::cout << "Found different from 1 associated w" << w.size() << std::endl;
    if (ws.size() != 2)
        std::cout << "Found different from 2 w decays" << ws.size() << std::endl;
    if (l.size() != 2)
        std::cout << "Found different from 2 leptons" << l.size() << std::endl;
    if (v.size() != 2)
        std::cout << "Found different from 2 neutrinos" << v.size() << std::endl;

    // Fill in tree variable
    tx->setBranch<int>("iswhchannel", true);
    tx->setBranch<LV>("w_p4", w[0]);
    tx->setBranch<LV>("q0_p4", qs[0].pt() > qs[1].pt() ? qs[0] : qs[1]);
    tx->setBranch<LV>("q1_p4", qs[0].pt() > qs[1].pt() ? qs[1] : qs[0]);
    tx->setBranch<LV>("h_p4", h[0]);
    tx->setBranch<LV>("w0_p4", ws[0].pt() > ws[1].pt() ? ws[0] : ws[1]);
    tx->setBranch<LV>("w1_p4", ws[0].pt() > ws[1].pt() ? ws[1] : ws[0]);
    tx->setBranch<LV>("l0_p4", l[0].pt() > l[1].pt() ? l[0] : l[1]);
    tx->setBranch<LV>("l1_p4", l[0].pt() > l[1].pt() ? l[1] : l[0]);
    tx->setBranch<LV>("v0_p4", v[0].pt() > v[1].pt() ? v[0] : v[1]);
    tx->setBranch<LV>("v1_p4", v[0].pt() > v[1].pt() ? v[1] : v[0]);
    return true;
}

//##############################################################################################################
bool babyMaker_v2::studyHiggsDecay()
{
    ProcessGenParticles();
    if (!isWprime() && !isVH()) return false; //file is certainly not doubly charged higgs
    const vector<int>& genPart_idx = coreGenPart.genPart_idx;
    const vector<int>& genPart_pdgId = coreGenPart.genPart_pdgId;
    const vector<int>& genPart_status = coreGenPart.genPart_status;
    const vector<int>& genPart_motherId = coreGenPart.genPart_motherId;
    const vector<int>& genPart_grandmaId = coreGenPart.genPart_grandmaId;
    const vector<LV>& genPart_p4 = coreGenPart.genPart_p4;

    // Find Higgs and the 2 decay 4 vector
    vector<int> higgs_idx;
    vector<LV> higgs_p4;
    vector<int> higgsdecay_idx;
    vector<LV> higgsdecay_p4;
    vector<int> higgsdecay_id;
    vector<int> higgsdecay_isstar;
    // WW decay
    bool isWW = false;
    vector<int> lepton_idx;
    vector<int> lepton_id;
    vector<int> lepton_mid;
    vector<int> lepton_isstar;
    vector<LV> lepton_p4;
    vector<int> neutrino_idx;
    vector<int> neutrino_id;
    vector<LV> neutrino_p4;
    vector<int> quark_idx;
    vector<int> quark_id;
    vector<int> quark_mid;
    vector<int> quark_isstar;
    vector<LV> quark_p4;
    vector<int> lepton_from_tau_idx;
    vector<int> lepton_from_tau_id;
    vector<LV> lepton_from_tau_p4;
    for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
    {
        int idx = genPart_idx[i];
        int pdgId = genPart_pdgId[i];
        int status = genPart_status[i];
        int motherId = genPart_motherId[i];
        int grandmaId = genPart_grandmaId[i];
        LV p4 = genPart_p4[i];
        if (abs(pdgId) == 25 && ((abs(motherId) == 9000002 && isWprime())||(status == 22 && isVH())))
        {
            higgs_idx.push_back(idx);
            higgs_p4.push_back(p4);
        }
        if (motherId == 25 && pdgId != 25)
        {
            higgsdecay_idx.push_back(idx);
            higgsdecay_p4.push_back(p4);
            higgsdecay_id.push_back(pdgId);
            higgsdecay_isstar.push_back(p4.mass() < 63.);
            if (abs(pdgId) == 24) isWW = true;
        }
        if ((abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15) && ((status == 23 || status == 1) || (abs(pdgId) == 15 && status == 2)) && abs(motherId) == 24 && abs(grandmaId) == 25)
        {
            lepton_idx.push_back(idx);
            lepton_id.push_back(pdgId);
            lepton_mid.push_back(motherId);
            lepton_p4.push_back(p4);
        }
        if ((abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16) && (status == 23 || status == 1) && abs(motherId) == 24 && abs(grandmaId) == 25)
        {
            neutrino_idx.push_back(idx);
            neutrino_id.push_back(pdgId);
            neutrino_p4.push_back(p4);
        }
        if ((abs(pdgId) >= 1 && abs(pdgId) <= 4))
        {
//            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        if ((abs(pdgId) >= 11 && abs(pdgId) <= 16))
        {
//            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        if ((abs(pdgId) >= 1 && abs(pdgId) <= 5) && (status == 23 || status == 1) && abs(motherId) == 24) // up to pdgid = 5 because W "can" decay to bc
        {
            if ( (abs(grandmaId) == 25 && isVH()) || (isWprime()) )
            {
                quark_idx.push_back(idx);
                quark_id.push_back(pdgId);
                quark_mid.push_back(motherId);
                quark_p4.push_back(p4);
//                std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
            }
        }
    }

//    std::cout <<  " higgs_idx.size(): " << higgs_idx.size() <<  " higgsdecay_idx.size(): " << higgsdecay_idx.size() <<  std::endl;
//    for (auto& hd_id : higgsdecay_id)
//        std::cout <<  " hd_id: " << hd_id <<  std::endl;
//    for (auto& h_idx : higgs_idx)
//        std::cout <<  " h_idx: " << h_idx <<  std::endl;
//
//    LV sumdecay = higgsdecay_p4[0] + higgsdecay_p4[1];
//
//    std::cout <<  " sumdecay.pt(): " << sumdecay.pt() <<  " sumdecay.eta(): " << sumdecay.eta() <<  " sumdecay.phi(): " << sumdecay.phi() <<  " sumdecay.energy(): " << sumdecay.energy() <<  " sumdecay.mass(): " << sumdecay.mass() <<  std::endl;
//    std::cout <<  " higgs_p4[0].pt(): " << higgs_p4[0].pt() <<  " higgs_p4[0].eta(): " << higgs_p4[0].eta() <<  " higgs_p4[0].phi(): " << higgs_p4[0].phi() <<  " higgs_p4[0].energy(): " << higgs_p4[0].energy() <<  " higgs_p4[0].mass(): " << higgs_p4[0].mass() <<  std::endl;

    if (isWW)
    {
        if (lepton_idx.size() + neutrino_idx.size() + quark_idx.size() != 4)
        {
            std::cout <<  " looper.getCurrentEventIndex(): " << looper.getCurrentEventIndex() <<  std::endl;
            std::cout <<  " lepton_idx.size(): " << lepton_idx.size() <<  " quark_idx.size(): " << quark_idx.size() <<  std::endl;
            for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
            {
                int idx = genPart_idx[i];
                int pdgId = genPart_pdgId[i];
                int status = genPart_status[i];
                int motherId = genPart_motherId[i];
                int grandmaId = genPart_grandmaId[i];
                LV p4 = genPart_p4[i];
//                std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
                if (abs(motherId) == 15 && abs(grandmaId) == 24)
                {
//                    std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
                }
            }
        }
        bool hastau = false;
        for (auto& id : lepton_id)
        {
            if (abs(id) == 15)
            {
                hastau = true;
            }
        }
        if (hastau)
        {
            for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
            {
                int idx = genPart_idx[i];
                int pdgId = genPart_pdgId[i];
                int status = genPart_status[i];
                int motherId = genPart_motherId[i];
                int grandmaId = genPart_grandmaId[i];
                LV p4 = genPart_p4[i];
                //                std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
                if (abs(motherId) == 15)
                {
//                    std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
                }
            }
        }
    }

    // Save information for H
    tx->setBranch<LV>("higgs_p4", higgs_p4[0]);

    // Save information for H->XX
    tx->setBranch<int>("higgsdecay", abs(higgsdecay_id[0]));
    tx->setBranch<vector<LV>>("decay_p4", higgsdecay_p4);
    tx->setBranch<vector<int>>("decay_id", higgsdecay_id);
    tx->setBranch<vector<int>>("decay_isstar", higgsdecay_isstar);
    tx->sortVecBranchesByPt("decay_p4", {}, {"decay_id", "decay_isstar"}, {});

    // Compute and fill the boosted variables as well as toe rotation
    // This also returns the reference vector
    float deta0, dphi0;
    float deta250, dphi250;
    float deta500, dphi500;
    float deta1000, dphi1000;
    float deta1500, dphi1500;
    std::tie(deta0, dphi0) = FillReBoostedVariables("decay"  , 0  , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );
    std::tie(deta250, dphi250) = FillReBoostedVariables("decay"  , 250  , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );
    std::tie(deta500, dphi500) = FillReBoostedVariables("decay"  , 500  , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );
    std::tie(deta1000, dphi1000) = FillReBoostedVariables("decay"  , 1000 , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );
    std::tie(deta1500, dphi1500) = FillReBoostedVariables("decay"  , 1500 , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );

    // If not WW stop here and move on
    if (abs(higgsdecay_id[0]) != 24) return true;

    // First save the information about the leptonic/hadronic channel information
    tx->setBranch<int>("nlep", lepton_id.size());
    tx->setBranch<int>("nquark", quark_id.size());

    // Determining which lepton or quarks had the mother of off-shell W
    int wstar_id = 0;
    for (unsigned i = 0; i < higgsdecay_id.size(); ++i)
    {
        if (higgsdecay_p4[i].mass() < 63.) // Perhaps you "could" have 62.5 and 62.5 ..? probably very unlikely
        {
            wstar_id = higgsdecay_id[i];
            break;
        }
    }

    // Finding off-shell failed. Sanity check.
    if (wstar_id == 0)
    {
        std::cout << "Failed to find off-shell W" << std::endl;
        FATALERROR(__FUNCTION__);
    }

    // Compute whether specific quark or lepton had a off-shell mother or not
    for (auto& id : lepton_mid) lepton_isstar.push_back(id == wstar_id);
    for (auto& id : quark_mid ) quark_isstar .push_back(id == wstar_id);

    // Save lepton and quark information
    tx->setBranch<vector<LV>>("lepton_p4", lepton_p4);
    tx->setBranch<vector<int>>("lepton_id", lepton_id);
    tx->setBranch<vector<int>>("lepton_isstar", lepton_isstar);
    tx->sortVecBranchesByPt("lepton_p4", {}, {"lepton_id", "lepton_isstar"}, {});
    tx->setBranch<vector<LV>>("quark_p4", quark_p4);
    tx->setBranch<vector<int>>("quark_id", quark_id);
    tx->setBranch<vector<int>>("quark_isstar", quark_isstar);
    tx->sortVecBranchesByPt("quark_p4", {}, {"quark_id", "quark_isstar"}, {});

    // Now start boosting things to normalize to certain set higgs pt
    FillReBoostedVariables("lepton" , 0    , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta0   , dphi0   );
    FillReBoostedVariables("lepton" , 250  , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta250 , dphi250 );
    FillReBoostedVariables("lepton" , 500  , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta500 , dphi500 );
    FillReBoostedVariables("lepton" , 1000 , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta1000, dphi1000);
    FillReBoostedVariables("lepton" , 1500 , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta1500, dphi1500);
    FillReBoostedVariables("quark"  , 0    , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta0   , dphi0    );
    FillReBoostedVariables("quark"  , 250  , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta250 , dphi250  );
    FillReBoostedVariables("quark"  , 500  , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta500 , dphi500  );
    FillReBoostedVariables("quark"  , 1000 , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta1000, dphi1000 );
    FillReBoostedVariables("quark"  , 1500 , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta1500, dphi1500 );

    FillDecayProductsDRVariables(250);
    FillDecayProductsDRVariables(500);
    FillDecayProductsDRVariables(1000);
    FillDecayProductsDRVariables(1500);


    return true;
}

//##############################################################################################################
std::tuple<float, float> babyMaker_v2::FillReBoostedVariables(TString bname, int ptBoost, const LV& higgs_p4, const vector<LV>& decay_p4, const vector<int>& decay_id, const vector<int>& decay_isstar, float ref_deta, float ref_dphi)
{
    vector<LV> decay_boosted_p4;
    vector<float> decay_h_dr;
    vector<float> decay_h_deta;
    vector<float> decay_h_dphi;
    vector<float> decay_h_deta_rotated;
    vector<float> decay_h_dphi_rotated;
    float ref_vec_deta;
    float ref_vec_dphi;

    // Get all the various re-boosted p4, dr, deta, dphi values
    std::tie(decay_boosted_p4, decay_h_dr, decay_h_deta, decay_h_dphi, decay_h_deta_rotated, decay_h_dphi_rotated, ref_vec_deta, ref_vec_dphi) = getReBoostedDRDEtaDPhi(ptBoost, higgs_p4, decay_p4, ref_deta, ref_dphi);

    // Set branches
    tx->setBranch<vector<LV>>   (TString::Format("boosted%d_%s_p4"    , ptBoost, bname.Data()), decay_boosted_p4);
    tx->setBranch<vector<int>>  (TString::Format("boosted%d_%s_id"    , ptBoost, bname.Data()), decay_id);
    tx->setBranch<vector<int>>  (TString::Format("boosted%d_%s_isstar", ptBoost, bname.Data()), decay_isstar);
    tx->setBranch<vector<float>>(TString::Format("boosted%d_%s_h_dr"  , ptBoost, bname.Data()), decay_h_dr);
    tx->setBranch<vector<float>>(TString::Format("boosted%d_%s_h_deta", ptBoost, bname.Data()), decay_h_deta);
    tx->setBranch<vector<float>>(TString::Format("boosted%d_%s_h_dphi", ptBoost, bname.Data()), decay_h_dphi);
    tx->setBranch<vector<float>>(TString::Format("boosted%d_%s_h_deta_rotated", ptBoost, bname.Data()), decay_h_deta_rotated);
    tx->setBranch<vector<float>>(TString::Format("boosted%d_%s_h_dphi_rotated", ptBoost, bname.Data()), decay_h_dphi_rotated);

    tx->sortVecBranchesByPt(
            TString::Format("boosted%d_%s_p4", ptBoost, bname.Data()),
            {
                TString::Format("boosted%d_%s_h_dr", ptBoost, bname.Data()),
                TString::Format("boosted%d_%s_h_deta", ptBoost, bname.Data()),
                TString::Format("boosted%d_%s_h_dphi", ptBoost, bname.Data()),
                TString::Format("boosted%d_%s_h_deta_rotated", ptBoost, bname.Data()),
                TString::Format("boosted%d_%s_h_dphi_rotated", ptBoost, bname.Data())
            },
            {
                TString::Format("boosted%d_%s_id", ptBoost, bname.Data()),
                TString::Format("boosted%d_%s_isstar", ptBoost, bname.Data())
            },
            {});

    return make_tuple(ref_vec_deta, ref_vec_dphi);
}

//##############################################################################################################
std::tuple<std::vector<LV>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, float, float> babyMaker_v2::getReBoostedDRDEtaDPhi(int ptBoost, const LV& higgs_p4, const vector<LV>& higgsdecay_p4, float ref_deta, float ref_dphi)
{

    // Return the values
    if (higgsdecay_p4.size() == 0)
        return make_tuple(std::vector<LV>(), std::vector<float>(), std::vector<float>(), std::vector<float>(), std::vector<float>(), std::vector<float>(), -999, -999);

    // Compute the boost vectors
    TLorentzVector target_tlv = RooUtil::Calc::getTLV(higgs_p4);
    if (ptBoost == 0)
        target_tlv.SetPtEtaPhiM(target_tlv.Pt(), target_tlv.Eta(), target_tlv.Phi(), target_tlv.M());
    else
        target_tlv.SetPtEtaPhiM(ptBoost, target_tlv.Eta(), target_tlv.Phi(), target_tlv.M());
    TVector3 to_target = target_tlv.BoostVector();
    TVector3 to_higgs = RooUtil::Calc::boostVector(higgs_p4);

    // Re-compute all the decay product with the provided boost vector
    vector<LV> boosted_decay_p4;
    for (unsigned int i = 0; i < higgsdecay_p4.size(); ++i)
    {
        LV decay = RooUtil::Calc::getBoosted(higgsdecay_p4[i], -to_higgs);
        RooUtil::Calc::boost(decay, to_target);
        boosted_decay_p4.push_back(decay);
    }

    // Compute the higgs and re-boosted objects' dr, deta, dphi
    vector<float> boosted_decay_h_dr;
    vector<float> boosted_decay_h_deta;
    vector<float> boosted_decay_h_dphi;

    //for (auto& p4 : boosted_decay_p4)
    for (unsigned int i = 0; i < boosted_decay_p4.size(); ++i)
    {
        LV p4 = boosted_decay_p4[i];
        boosted_decay_h_dr.push_back(RooUtil::Calc::DeltaR(p4, higgs_p4));
        boosted_decay_h_deta.push_back(RooUtil::Calc::DeltaEta(p4, higgs_p4));
        boosted_decay_h_dphi.push_back(RooUtil::Calc::DeltaPhi(p4, higgs_p4));
    }

    // Compute the rotation angle based off of the leading decay particle
    int lead_idx = boosted_decay_p4[0].pt() > boosted_decay_p4[1].pt() ? 0 : 1;
    TVector2 ref_vec(boosted_decay_h_deta[lead_idx], boosted_decay_h_dphi[lead_idx]);
    if (ref_dphi != -999)
        ref_vec.Set(ref_deta, ref_dphi);

    vector<float> boosted_decay_h_deta_rotated;
    vector<float> boosted_decay_h_dphi_rotated;
    for (unsigned i = 0; i < boosted_decay_h_deta.size(); ++i)
    {
        TVector2 ang_vec;
        ang_vec.Set(boosted_decay_h_deta[i], boosted_decay_h_dphi[i]);
        TVector2 new_vec = ang_vec.Rotate(-ref_vec.Phi() + TMath::Pi() / 2.);
        boosted_decay_h_deta_rotated.push_back(new_vec.Px());
        boosted_decay_h_dphi_rotated.push_back(new_vec.Py());
    }

    // Return the values
    return make_tuple(boosted_decay_p4, boosted_decay_h_dr, boosted_decay_h_deta, boosted_decay_h_dphi, boosted_decay_h_deta_rotated, boosted_decay_h_dphi_rotated, (float) ref_vec.Px(), (float) ref_vec.Py());
}

//##############################################################################################################
void babyMaker_v2::FillDecayProductsDRVariables(int pt)
{
    // if not lvqq decay then no point of calculating this max dr
    if (tx->getBranch<int>("nlep", true) != 1)
        return;
    const vector<LV>& lepton_p4 = tx->getBranch<vector<LV>>(TString::Format("boosted%d_lepton_p4", pt), true);
    const vector<LV>& quark_p4 = tx->getBranch<vector<LV>>(TString::Format("boosted%d_quark_p4", pt), true);
    const LV& l_p4 = lepton_p4[0];
    const LV& q0_p4 = quark_p4[0];
    const LV& q1_p4 = quark_p4[1];
    float lq0_dr = RooUtil::Calc::DeltaR(l_p4, q0_p4);
    float lq1_dr = RooUtil::Calc::DeltaR(l_p4, q1_p4);
    float q0q1_dr = RooUtil::Calc::DeltaR(q0_p4, q1_p4);
    tx->setBranch<float>(TString::Format("boosted%d_lqq_max_dr", pt), max(lq0_dr, max(lq1_dr, q0q1_dr)));
    tx->setBranch<float>(TString::Format("boosted%d_lq0_dr", pt), lq0_dr);
    tx->setBranch<float>(TString::Format("boosted%d_lq1_dr", pt), lq1_dr);
    tx->setBranch<float>(TString::Format("boosted%d_qq_dr", pt), q0q1_dr);
}

//##############################################################################################################
bool babyMaker_v2::studyWWW()
{
    ProcessGenParticles();
    if (!isWWW() && !isWHSUSY()) return false;
    if (isWHWWW()) tx->setBranch<int>("iswhwww", 1);
    if (isSMWWW()) tx->setBranch<int>("iswhwww", 0);
    const vector<int>& genPart_pdgId = coreGenPart.genPart_pdgId;
    const vector<int>& genPart_status = coreGenPart.genPart_status;
    const vector<int>& genPart_motherId = coreGenPart.genPart_motherId;
    const vector<int>& genPart_grandmaId = coreGenPart.genPart_grandmaId;
    const vector<int>& genPart_mother_idx = coreGenPart.genPart_mother_idx;
    const vector<LV>& genPart_p4 = coreGenPart.genPart_p4;

    std::cout.setstate(std::ios_base::failbit);

    std::cout << std::endl;
    vector<LV> h;
    vector<LV> w;
    vector<int> w_m_id;
    vector<LV> l;
    vector<int> l_id;
    vector<int> l_m_idx;
    vector<LV> v;
    vector<int> v_m_idx;
    vector<LV> q;
    vector<int> q_m_idx;
    for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
    {
        int pdgId = genPart_pdgId[i];
        int status = genPart_status[i];
        int motherId = genPart_motherId[i];
        int grandmaId = genPart_grandmaId[i];
        int mother_idx = genPart_mother_idx[i];
        LV p4 = genPart_p4[i];

        // The W bosons (status 22)
        if (status == 22 && abs(pdgId) == 24)
        {
            w.push_back(p4);
            w_m_id.push_back(motherId);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        // W decays
        else if ((abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15) && ((status == 23 || status == 1) || (abs(pdgId) == 15 && status == 2)) && abs(motherId) == 24)
        {
            l.push_back(p4);
            l_id.push_back(pdgId);
            l_m_idx.push_back(mother_idx);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        else if ((abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16) && ((status == 23 || status == 1) || (abs(pdgId) == 16 && status == 2)) && abs(motherId) == 24)
        {
            v.push_back(p4);
            v_m_idx.push_back(mother_idx);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        else if ((abs(pdgId) >= 1 && abs(pdgId) <= 4) && (status == 23 || status == 1) && abs(motherId) == 24)
        {
            q.push_back(p4);
            q_m_idx.push_back(mother_idx);
            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
    }

    std::cout <<  " w.size(): " << w.size() <<  std::endl;
    std::cout <<  " q.size(): " << q.size() <<  std::endl;
    std::cout <<  " l.size(): " << l.size() <<  std::endl;
    std::cout <<  " v.size(): " << v.size() <<  std::endl;

    std::cout.clear();

    if (isWHSUSY() && w.size() != 3)
        return false;

//    // Sanity check
//    if (w.size() != 3)
//    {
//        std::cout << "Found different from 3 w = " << w.size() << std::endl;
//    }

//    // Sanity check
//    if (w.size() != 3)
//    {
//        std::cout << "Found different from 3 w = " << w.size() << std::endl;
//        for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
//        {
//            int pdgId = genPart_pdgId[i];
//            int status = genPart_status[i];
//            int motherId = genPart_motherId[i];
//            int grandmaId = genPart_grandmaId[i];
//            int mother_idx = genPart_mother_idx[i];
//            LV p4 = genPart_p4[i];
//            //if (abs(pdgId) == 24)
//            //    std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
//            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
//        }
//
//        bool isHtoWW = false;
//        bool isZthere = false;
//        for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
//        {
//            if (abs(genPart_pdgId[i]) == 23) isZthere = true;
//            if (abs(genPart_pdgId[i]) == 24 && genPart_status[i] == 22 && abs(genPart_motherId[i]) == 25) isHtoWW = true;
//            if (isHtoWW && isZthere) break;
//        }
//        std::cout <<  " isHtoWW: " << isHtoWW <<  " isZthere: " << isZthere <<  std::endl;
//
//    }

    // Set has_tau
    bool has_tau = false;
    for (unsigned i = 0; i < l_id.size(); ++i) if (abs(l_id[i]) == 15) has_tau = true;
    tx->setBranch<int>("has_tau", has_tau);

    // Set www_channel
    if      (l_id.size() == 2 && l_id[0] * l_id[1] > 0) tx->setBranch<int>("www_channel", 2);
    else if (l_id.size() == 2 && l_id[0] * l_id[1] < 0) tx->setBranch<int>("www_channel", 4);
    else if (l_id.size() == 3) tx->setBranch<int>("www_channel", 3);
    else if (l_id.size() == 1) tx->setBranch<int>("www_channel", 1);
    else if (l_id.size() == 0) tx->setBranch<int>("www_channel", 0);
    else FATALERROR(__FUNCTION__);

    // Set w boson variables
    for (unsigned i = 0; i < w.size(); ++i)
    {
        tx->pushbackToBranch<LV>("w_p4", w[i]);
        tx->pushbackToBranch<int>("w_isH", abs(w_m_id[i]) == 25);
        tx->pushbackToBranch<int>("w_isstar", w[i].mass() < 65.);
    }

    // Set leptons
    for (unsigned i = 0; i < l.size(); ++i)
    {
        tx->pushbackToBranch<LV>("l_p4", l[i]);
        tx->pushbackToBranch<float>("l_w_pt", genPart_p4[l_m_idx[i]].pt());
        tx->pushbackToBranch<float>("l_w_eta", genPart_p4[l_m_idx[i]].eta());
        tx->pushbackToBranch<float>("l_w_phi", genPart_p4[l_m_idx[i]].phi());
        tx->pushbackToBranch<float>("l_w_mass", genPart_p4[l_m_idx[i]].mass());
        tx->pushbackToBranch<int>("l_w_id", genPart_pdgId[l_m_idx[i]]);
        tx->pushbackToBranch<int>("l_isH", abs(genPart_motherId[l_m_idx[i]]) == 25);
        tx->pushbackToBranch<int>("l_isstar", genPart_p4[l_m_idx[i]].mass() < 65.);
        tx->pushbackToBranch<int>("l_istau", abs(l_id[i]) == 15);
    }

    // Set quarks
    for (unsigned i = 0; i < q.size(); ++i)
    {
        tx->pushbackToBranch<LV>("q_p4", q[i]);
        tx->pushbackToBranch<float>("q_w_pt", genPart_p4[q_m_idx[i]].pt());
        tx->pushbackToBranch<float>("q_w_eta", genPart_p4[q_m_idx[i]].eta());
        tx->pushbackToBranch<float>("q_w_phi", genPart_p4[q_m_idx[i]].phi());
        tx->pushbackToBranch<float>("q_w_mass", genPart_p4[q_m_idx[i]].mass());
        tx->pushbackToBranch<int>("q_w_id", genPart_pdgId[q_m_idx[i]]);
        tx->pushbackToBranch<int>("q_isH", abs(genPart_motherId[q_m_idx[i]]) == 25);
        tx->pushbackToBranch<int>("q_isstar", genPart_p4[q_m_idx[i]].mass() < 65.);
    }

    tx->sortVecBranchesByPt("w_p4", {}, {"w_isH", "w_isstar"}, {});
    tx->sortVecBranchesByPt("l_p4", {"l_w_pt", "l_w_eta", "l_w_phi", "l_w_mass"}, {"l_isH", "l_isstar", "l_istau"}, {});
    tx->sortVecBranchesByPt("q_p4", {"q_w_pt", "q_w_eta", "q_w_phi", "q_w_mass"}, {"q_isH", "q_isstar"}, {});

    if (l.size() == 2)
    {
        tx->setBranch<float>("dRllSS", ROOT::Math::VectorUtil::DeltaR(l[0], l[1]));
    }
    if (q.size() == 2)
    {
        tx->setBranch<float>("dRqqSS", ROOT::Math::VectorUtil::DeltaR(q[0], q[1]));
    }

    if (l.size() == 3 && isWHWWW())
    {
        LV lep0;
        LV lep1;
        for (unsigned i = 0; i < l.size(); ++i)
        {
            if (abs(genPart_motherId[l_m_idx[i]]) == 25)
            {
                if (lep0.pt() == 0)
                {
                    lep0 = l[i];
                }
                else
                {
                    lep1 = l[i];
                }
            }
        }

        LV neu0;
        LV neu1;
        for (unsigned i = 0; i < v.size(); ++i)
        {
            if (abs(genPart_motherId[v_m_idx[i]]) == 25)
            {
                if (neu0.pt() == 0)
                {
                    neu0 = v[i];
                }
                else
                {
                    neu1 = v[i];
                }
            }
        }

        tx->setBranch<float>("Mll_higgs", (lep0 + lep1).mass());
        tx->setBranch<float>("DPhill_higgs", fabs(ROOT::Math::VectorUtil::DeltaPhi(lep0, lep1)));
        tx->setBranch<float>("MT_higgs", mT(lep0 + lep1, neu0 + neu1));
    }


//    tx->pushbackToBranch<LV>("v_p4");
//    tx->pushbackToBranch<LV>("v_w_p4");
//    tx->pushbackToBranch<int>("v_isstar");
//    tx->pushbackToBranch<int>("v_isH");
//
//    tx->pushbackToBranch<LV>("q_p4");
//    tx->pushbackToBranch<LV>("q_w_p4");
//    tx->pushbackToBranch<int>("q_isstar");
//    tx->pushbackToBranch<int>("q_isH");

//    tx->setBranch<LV>("w0_p4", w[0]);
//    tx->setBranch<LV>("w1_p4", w[1]);
//    tx->setBranch<LV>("w2_p4", w[2]);
//    tx->setBranch<int>("w0_islep");
//    tx->setBranch<int>("w1_islep");
//    tx->setBranch<int>("w2_islep");
//    tx->setBranch<int>("w0_isstar");
//    tx->setBranch<int>("w1_isstar");
//    tx->setBranch<int>("w2_isstar");
//    tx->setBranch<int>("w0_isH");
//    tx->setBranch<int>("w1_isH");
//    tx->setBranch<int>("w2_isH");
//    tx->setBranch<LV>("l0_p4");
//    tx->setBranch<LV>("l1_p4");
//    tx->setBranch<LV>("l2_p4");
//    tx->setBranch<LV>("v0_p4");
//    tx->setBranch<LV>("v1_p4");
//    tx->setBranch<LV>("v2_p4");
//    tx->setBranch<LV>("q0_p4");
//    tx->setBranch<LV>("q1_p4");
//    tx->setBranch<LV>("q2_p4");
//    tx->setBranch<LV>("q3_p4");
//    tx->setBranch<LV>("q4_p4");
//    tx->setBranch<LV>("q5_p4");
//    tx->setBranch<int>("l0_windex");
//    tx->setBranch<int>("l1_windex");
//    tx->setBranch<int>("l2_windex");
//    tx->setBranch<int>("v0_windex");
//    tx->setBranch<int>("v1_windex");
//    tx->setBranch<int>("v2_windex");
//    tx->setBranch<int>("q0_windex");
//    tx->setBranch<int>("q1_windex");
//    tx->setBranch<int>("q2_windex");
//    tx->setBranch<int>("q3_windex");
//    tx->setBranch<int>("q4_windex");
//    tx->setBranch<int>("q5_windex");
//    tx->setBranch<LV>("l0_w_p4");
//    tx->setBranch<LV>("l1_w_p4");
//    tx->setBranch<LV>("l2_w_p4");
//    tx->setBranch<LV>("v0_w_p4");
//    tx->setBranch<LV>("v1_w_p4");
//    tx->setBranch<LV>("v2_w_p4");
//    tx->setBranch<LV>("q0_w_p4");
//    tx->setBranch<LV>("q1_w_p4");
//    tx->setBranch<LV>("q2_w_p4");
//    tx->setBranch<LV>("q3_w_p4");
//    tx->setBranch<LV>("q4_w_p4");
//    tx->setBranch<LV>("q5_w_p4");
//    tx->setBranch<int>("l0_isstar");
//    tx->setBranch<int>("l1_isstar");
//    tx->setBranch<int>("l2_isstar");
//    tx->setBranch<int>("v0_isstar");
//    tx->setBranch<int>("v1_isstar");
//    tx->setBranch<int>("v2_isstar");
//    tx->setBranch<int>("q0_isstar");
//    tx->setBranch<int>("q1_isstar");
//    tx->setBranch<int>("q2_isstar");
//    tx->setBranch<int>("q3_isstar");
//    tx->setBranch<int>("q4_isstar");
//    tx->setBranch<int>("q5_isstar");
//    tx->setBranch<int>("l0_isH");
//    tx->setBranch<int>("l1_isH");
//    tx->setBranch<int>("l2_isH");
//    tx->setBranch<int>("v0_isH");
//    tx->setBranch<int>("v1_isH");
//    tx->setBranch<int>("v2_isH");
//    tx->setBranch<int>("q0_isH");
//    tx->setBranch<int>("q1_isH");
//    tx->setBranch<int>("q2_isH");
//    tx->setBranch<int>("q3_isH");
//    tx->setBranch<int>("q4_isH");
//    tx->setBranch<int>("q5_isH");
//    tx->setBranch<int>("l0_istau");
//    tx->setBranch<int>("l1_istau");
//    tx->setBranch<int>("l2_istau");

}

//##############################################################################################################
bool babyMaker_v2::studyWHWWW()
{
    ProcessGenParticles();

    if (!isVH()) return false;

    if (isWHWWW()) tx->setBranch<int>("iswhwww", 1);

    const vector<int>& genPart_idx = coreGenPart.genPart_idx;
    const vector<int>& genPart_pdgId = coreGenPart.genPart_pdgId;
    const vector<int>& genPart_status = coreGenPart.genPart_status;
    const vector<int>& genPart_motherId = coreGenPart.genPart_motherId;
    const vector<int>& genPart_grandmaId = coreGenPart.genPart_grandmaId;
    const vector<LV>& genPart_p4 = coreGenPart.genPart_p4;

    // Find Higgs and the 2 decay 4 vector
    vector<int> higgs_idx;
    vector<LV> higgs_p4;
    vector<int> higgsdecay_idx;
    vector<LV> higgsdecay_p4;
    vector<int> higgsdecay_id;
    vector<int> higgsdecay_isstar;
    // WW decay
    bool isWW = false;
    vector<int> lepton_idx;
    vector<int> lepton_id;
    vector<int> lepton_mid;
    vector<int> lepton_isstar;
    vector<LV> lepton_p4;
    vector<int> neutrino_idx;
    vector<int> neutrino_id;
    vector<LV> neutrino_p4;
    vector<int> quark_idx;
    vector<int> quark_id;
    vector<int> quark_mid;
    vector<int> quark_isstar;
    vector<LV> quark_p4;
    vector<int> lepton_from_tau_idx;
    vector<int> lepton_from_tau_id;
    vector<LV> lepton_from_tau_p4;
    for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
    {
        int idx = genPart_idx[i];
        int pdgId = genPart_pdgId[i];
        int status = genPart_status[i];
        int motherId = genPart_motherId[i];
        int grandmaId = genPart_grandmaId[i];
        LV p4 = genPart_p4[i];
        if (abs(pdgId) == 25)
        {
            higgs_idx.push_back(idx);
            higgs_p4.push_back(p4);
        }
        if (motherId == 25 && pdgId != 25)
        {
            higgsdecay_idx.push_back(idx);
            higgsdecay_p4.push_back(p4);
            higgsdecay_id.push_back(pdgId);
            higgsdecay_isstar.push_back(p4.mass() < 63.);
            if (abs(pdgId) == 24) isWW = true;
        }
        if ((abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15) && ((status == 23 || status == 1) || (abs(pdgId) == 15 && status == 2)) && abs(motherId) == 24 && abs(grandmaId) == 25)
        {
            lepton_idx.push_back(idx);
            lepton_id.push_back(pdgId);
            lepton_mid.push_back(motherId);
            lepton_p4.push_back(p4);
        }
        if ((abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16) && (status == 23 || status == 1) && abs(motherId) == 24 && abs(grandmaId) == 25)
        {
            neutrino_idx.push_back(idx);
            neutrino_id.push_back(pdgId);
            neutrino_p4.push_back(p4);
        }
        if ((abs(pdgId) >= 1 && abs(pdgId) <= 4))
        {
//            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        if ((abs(pdgId) >= 11 && abs(pdgId) <= 16))
        {
//            std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
        }
        if ((abs(pdgId) >= 1 && abs(pdgId) <= 5) && (status == 23 || status == 1) && abs(motherId) == 24) // up to pdgid = 5 because W "can" decay to bc
        {
            quark_idx.push_back(idx);
            quark_id.push_back(pdgId);
            quark_mid.push_back(motherId);
            quark_p4.push_back(p4);
        }
    }


    if (isWW)
    {
//        if (lepton_idx.size() + neutrino_idx.size() + quark_idx.size() != 4)
//        {
//            std::cout <<  " looper.getCurrentEventIndex(): " << looper.getCurrentEventIndex() <<  std::endl;
//            std::cout <<  " lepton_idx.size(): " << lepton_idx.size() <<  " quark_idx.size(): " << quark_idx.size() <<  std::endl;
//            for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
//            {
//                int idx = genPart_idx[i];
//                int pdgId = genPart_pdgId[i];
//                int status = genPart_status[i];
//                int motherId = genPart_motherId[i];
//                int grandmaId = genPart_grandmaId[i];
//                LV p4 = genPart_p4[i];
////                std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
//                if (abs(motherId) == 15 && abs(grandmaId) == 24)
//                {
////                    std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
//                }
//            }
//        }
        bool hastau = false;
        for (auto& id : lepton_id)
        {
            if (abs(id) == 15)
            {
                hastau = true;
            }
        }
        if (hastau)
        {
            for (unsigned int i = 0; i < genPart_pdgId.size(); ++i)
            {
                int idx = genPart_idx[i];
                int pdgId = genPart_pdgId[i];
                int status = genPart_status[i];
                int motherId = genPart_motherId[i];
                int grandmaId = genPart_grandmaId[i];
                LV p4 = genPart_p4[i];
                //                std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
                if (abs(motherId) == 15)
                {
//                    std::cout <<  " pdgId: " << pdgId <<  " motherId: " << motherId <<  " grandmaId: " << grandmaId <<  " status: " << status <<  " p4.pt(): " << p4.pt() <<  " p4.eta(): " << p4.eta() <<  " p4.phi(): " << p4.phi() <<  std::endl;
                }
            }
        }
    }

    // Save information for H
    tx->setBranch<LV>("higgs_p4", higgs_p4[0]);

    // Save information for H->XX
    tx->setBranch<int>("higgsdecay", abs(higgsdecay_id[0]));
    tx->setBranch<vector<LV>>("decay_p4", higgsdecay_p4);
    tx->setBranch<vector<int>>("decay_id", higgsdecay_id);
    tx->setBranch<vector<int>>("decay_isstar", higgsdecay_isstar);
    tx->sortVecBranchesByPt("decay_p4", {}, {"decay_id", "decay_isstar"}, {});

    // Compute and fill the boosted variables as well as toe rotation
    // This also returns the reference vector
    float deta0, dphi0;
    float deta500, dphi500;
    float deta1000, dphi1000;
    float deta1500, dphi1500;
    std::tie(deta0, dphi0) = FillReBoostedVariables("decay"  , 0  , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );
    std::tie(deta500, dphi500) = FillReBoostedVariables("decay"  , 500  , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );
    std::tie(deta1000, dphi1000) = FillReBoostedVariables("decay"  , 1000 , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );
    std::tie(deta1500, dphi1500) = FillReBoostedVariables("decay"  , 1500 , higgs_p4[0] , higgsdecay_p4 , higgsdecay_id , higgsdecay_isstar );

    // If not WW stop here and move on
    if (abs(higgsdecay_id[0]) != 24) return true;

    // First save the information about the leptonic/hadronic channel information
    tx->setBranch<int>("nlep", lepton_id.size());
    tx->setBranch<int>("nquark", quark_id.size());

    // Determining which lepton or quarks had the mother of off-shell W
    int wstar_id = 0;
    for (unsigned i = 0; i < higgsdecay_id.size(); ++i)
    {
        if (higgsdecay_p4[i].mass() < 63.) // Perhaps you "could" have 62.5 and 62.5 ..? probably very unlikely
        {
            wstar_id = higgsdecay_id[i];
            break;
        }
    }

    // Finding off-shell failed. Sanity check.
    if (wstar_id == 0)
    {
        std::cout << "Failed to find off-shell W" << std::endl;
        FATALERROR(__FUNCTION__);
    }

    // Compute whether specific quark or lepton had a off-shell mother or not
    for (auto& id : lepton_mid) lepton_isstar.push_back(id == wstar_id);
    for (auto& id : quark_mid ) quark_isstar .push_back(id == wstar_id);

    // Save lepton and quark information
    tx->setBranch<vector<LV>>("lepton_p4", lepton_p4);
    tx->setBranch<vector<int>>("lepton_id", lepton_id);
    tx->setBranch<vector<int>>("lepton_isstar", lepton_isstar);
    tx->sortVecBranchesByPt("lepton_p4", {}, {"lepton_id", "lepton_isstar"}, {});
    tx->setBranch<vector<LV>>("quark_p4", quark_p4);
    tx->setBranch<vector<int>>("quark_id", quark_id);
    tx->setBranch<vector<int>>("quark_isstar", quark_isstar);
    tx->sortVecBranchesByPt("quark_p4", {}, {"quark_id", "quark_isstar"}, {});

    // Now start boosting things to normalize to certain set higgs pt
    FillReBoostedVariables("lepton" , 0    , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta0   , dphi0   );
    FillReBoostedVariables("lepton" , 500  , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta500 , dphi500 );
    FillReBoostedVariables("lepton" , 1000 , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta1000, dphi1000);
    FillReBoostedVariables("lepton" , 1500 , higgs_p4[0] , lepton_p4     , lepton_id     , lepton_isstar    , deta1500, dphi1500);
    FillReBoostedVariables("quark"  , 0    , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta0   , dphi0    );
    FillReBoostedVariables("quark"  , 500  , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta500 , dphi500  );
    FillReBoostedVariables("quark"  , 1000 , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta1000, dphi1000 );
    FillReBoostedVariables("quark"  , 1500 , higgs_p4[0] , quark_p4      , quark_id      , quark_isstar     , deta1500, dphi1500 );


}

//##############################################################################################################
void babyMaker_v2::setWHSMSMassAndWeights()
{
    using namespace tas;
    //get susy particle masses from sparms
    float mass_chargino;
    float mass_lsp;
    float mass_gluino;
    if (!SampleNiceName().BeginsWith("whsusy-2l"))
    {
        for (unsigned int nsparm = 0; nsparm < sparm_names().size(); ++nsparm) {
            if (sparm_names().at(nsparm).Contains("mCh")) mass_chargino = sparm_values().at(nsparm);
            if (sparm_names().at(nsparm).Contains("mLSP")) mass_lsp = sparm_values().at(nsparm);
            if (sparm_names().at(nsparm).Contains("mGl")) mass_gluino = sparm_values().at(nsparm);
        }
    }
    else
    {
        mass_chargino = 200;
        mass_lsp = 50;
    }
    tx->clear();
    tx->setBranch<float>("chimass", mass_chargino);
    tx->setBranch<float>("lspmass", mass_lsp);
//    std::cout <<  " h_nrawevents_SMS: " << h_nrawevents_SMS <<  std::endl;
    if (genps_weight() > 0) h_nrawevents_SMS -> Fill(mass_chargino, mass_lsp, 1);
    else if (genps_weight() < 0) h_nrawevents_SMS -> Fill(mass_chargino, mass_lsp, -1);
    float xsec = hxsec->GetBinContent(hxsec->FindBin(mass_chargino));
    float xsec_uncert = hxsec->GetBinError(hxsec->FindBin(mass_chargino));

    float SMSpdf_weight_up = 1;
    float SMSpdf_weight_down = 1;
    float SMSsum_of_weights = 0;
    float SMSaverage_of_weights = 0;
    //error on pdf replicas
    //fastsim has first genweights bin being ==1
    if (genweights().size() > 111) { //fix segfault
        for (int ipdf = 10; ipdf < 110; ipdf++) {
            SMSaverage_of_weights += cms3.genweights().at(ipdf);
        } // average of weights
        SMSaverage_of_weights = SMSaverage_of_weights / 100.;
        for (int ipdf = 10; ipdf < 110; ipdf++) {
            SMSsum_of_weights += pow(cms3.genweights().at(ipdf) - SMSaverage_of_weights, 2);
        } //std of weights.
        SMSpdf_weight_up = (SMSaverage_of_weights + sqrt(SMSsum_of_weights / 99.));
        SMSpdf_weight_down = (SMSaverage_of_weights - sqrt(SMSsum_of_weights / 99.));
//        std::cout <<  " h_nevents_SMS: " << h_nevents_SMS <<  std::endl;
//        std::cout <<  " genweights().size(): " << genweights().size() <<  std::endl;
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 1, genweights()[1]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 2, genweights()[2]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 3, genweights()[3]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 4, genweights()[4]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 5, genweights()[5]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 6, genweights()[6]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 7, genweights()[7]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 8, genweights()[8]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 9, genweights()[9]);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 10, SMSpdf_weight_up);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 11, SMSpdf_weight_down);
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 12, genweights()[110]); // α_s variation.
        h_nevents_SMS -> Fill(mass_chargino, mass_lsp, 13, genweights()[111]); // α_s variation.
    }
}

//##############################################################################################################
void babyMaker_v2::setWHSMSMass()
{
    using namespace tas;
    //get susy particle masses from sparms
    float mass_chargino;
    float mass_lsp;
    float mass_gluino;
    if (!SampleNiceName().BeginsWith("whsusy-2l"))
    {
        for (unsigned int nsparm = 0; nsparm < sparm_names().size(); ++nsparm) {
            if (sparm_names().at(nsparm).Contains("mCh")) mass_chargino = sparm_values().at(nsparm);
            if (sparm_names().at(nsparm).Contains("mLSP")) mass_lsp = sparm_values().at(nsparm);
            if (sparm_names().at(nsparm).Contains("mGl")) mass_gluino = sparm_values().at(nsparm);
        }
    }
    else
    {
        mass_chargino = 200;
        mass_lsp = 50;
    }
    tx->clear();
    tx->setBranch<float>("chimass", mass_chargino);
    tx->setBranch<float>("lspmass", mass_lsp);
}

//##############################################################################################################
bool babyMaker_v2::filterWHMass(float chimass, float lspmass)
{
    using namespace tas;
    //get susy particle masses from sparms
    float mass_chargino;
    float mass_lsp;
    for (unsigned int nsparm = 0; nsparm < sparm_names().size(); ++nsparm) {
        if (chimass >= 0)
            if (sparm_names().at(nsparm).Contains("mCh"))
                if (sparm_values().at(nsparm) != chimass)
                    return false;
        if (lspmass >= 0)
            if (sparm_names().at(nsparm).Contains("mLSP"))
                if (sparm_values().at(nsparm) != lspmass)
                    return false;
    }
    return true;
}

//##############################################################################################################
int babyMaker_v2::gentype_v2()
{
    if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
    {
        int ngenLepFromTau = tx->getBranch<int>("ngenLepFromTau");
        int ngenLepFromBoson = tx->getBranch<int>("ngenLepFromBoson");
        int ngenLep = tx->getBranch<int>("ngenLep");
        int nVlep = tx->getBranch<int>("nVlep");
        int nLlep = tx->getBranch<int>("nLlep");
        vector<LV> lep_p4 = tx->getBranch<vector<LV>>("lep_p4");
        vector<int> lep_motherIdSS = tx->getBranch<vector<int>>("lep_motherIdSS");
        vector<int> lep_isFromW = tx->getBranch<vector<int>>("lep_isFromW");
        vector<int> lep_isFromZ = tx->getBranch<vector<int>>("lep_isFromZ");
        vector<int> genPart_charge = tx->getBranch<vector<int>>("genPart_charge");
        vector<int> lep_charge = tx->getBranch<vector<int>>("lep_charge");
        vector<int> lep_genPart_index = tx->getBranch<vector<int>>("lep_genPart_index");
        std::cout <<  " ngenLepFromTau: " << ngenLepFromTau <<  std::endl;
        std::cout <<  " ngenLepFromBoson: " << ngenLepFromBoson <<  std::endl;
        std::cout <<  " ngenLep: " << ngenLep <<  std::endl;
        std::cout <<  " nVlep: " << nVlep <<  std::endl;
        std::cout <<  " nLlep: " << nLlep <<  std::endl;
        std::cout <<  " lep_p4[0].pt(): " << lep_p4[0].pt() <<  " lep_p4[1].pt(): " << lep_p4[1].pt() <<  " lep_p4[2].pt(): " << lep_p4[2].pt() <<  std::endl;
        std::cout <<  " lep_motherIdSS[0]: " << lep_motherIdSS[0] <<  " lep_motherIdSS[1]: " << lep_motherIdSS[1] <<  " lep_motherIdSS[2]: " << lep_motherIdSS[2] <<  std::endl;
        std::cout <<  " lep_isFromW[0]: " << lep_isFromW[0] <<  " lep_isFromW[1]: " << lep_isFromW[1] <<  " lep_isFromW[2]: " << lep_isFromW[2] <<  std::endl;
        std::cout <<  " lep_isFromZ[0]: " << lep_isFromZ[0] <<  " lep_isFromZ[1]: " << lep_isFromZ[1] <<  " lep_isFromZ[2]: " << lep_isFromZ[2] <<  std::endl;
//        std::cout <<  " genPart_charge[lep_genPart_index[0]]: " << genPart_charge[lep_genPart_index[0]] <<  " genPart_charge[lep_genPart_index[1]]: " << genPart_charge[lep_genPart_index[1]] <<  " genPart_charge[lep_genPart_index[2]]: " << genPart_charge[lep_genPart_index[2]] <<  std::endl;
        std::cout <<  " lep_charge[0]: " << lep_charge[0] <<  " lep_charge[1]: " << lep_charge[1] <<  " lep_charge[2]: " << lep_charge[2] <<  std::endl;
        std::cout <<  " lep_genPart_index[0]: " << lep_genPart_index[0] <<  " lep_genPart_index[1]: " << lep_genPart_index[1] <<  " lep_genPart_index[2]: " << lep_genPart_index[2] <<  std::endl;
        std::cout <<  " (lep_p4[0]+lep_p4[1]).mass(): " << (lep_p4[0]+lep_p4[1]).mass() <<  std::endl;

        for (unsigned int iGen = 0; iGen < cms3.genps_p4().size(); iGen++)
        {
            std::cout <<  " iGen: " << iGen <<  std::endl;
            std::cout <<  " cms3.genps_p4()[iGen]: " << cms3.genps_p4()[iGen] <<  std::endl;
            std::cout <<  " cms3.genps_p4()[iGen].pt(): " << cms3.genps_p4()[iGen].pt() <<  std::endl;
            std::cout <<  " cms3.genps_p4()[iGen].eta(): " << cms3.genps_p4()[iGen].eta() <<  std::endl;
            std::cout <<  " cms3.genps_p4()[iGen].phi(): " << cms3.genps_p4()[iGen].phi() <<  std::endl;
            std::cout <<  " cms3.genps_p4()[iGen].mass(): " << cms3.genps_p4()[iGen].mass() <<  std::endl;
            std::cout <<  " cms3.genps_id()[iGen]: " << cms3.genps_id()[iGen] <<  std::endl;
            std::cout <<  " cms3.genps_status()[iGen]: " << cms3.genps_status()[iGen] <<  std::endl;
            std::cout <<  " int(cms3.genps_charge()[iGen]): " << int(cms3.genps_charge()[iGen]) <<  std::endl;
            std::cout <<  " cms3.genps_id_simplemother()[iGen]: " << cms3.genps_id_simplemother()[iGen] <<  std::endl;
            std::cout <<  " cms3.genps_id_simplegrandma()[iGen]: " << cms3.genps_id_simplegrandma()[iGen] <<  std::endl;
            std::cout <<  " cms3.genps_idx_simplemother()[iGen]: " << cms3.genps_idx_simplemother()[iGen] <<  std::endl;
            std::cout <<  " cms3.genps_isMostlyLikePythia6Status3()[iGen]: " << cms3.genps_isMostlyLikePythia6Status3()[iGen] <<  std::endl;
        }
    }

    bool gammafake = false;
    bool jetfake   = false;
    unsigned int ngenlep = tx->getBranch<int>("ngenLepFromTau") + tx->getBranch<int>("ngenLepFromBoson");
    unsigned int nW(0), nZ(0);
    bool lep1_real = tx->getBranch<vector<int>>("lep_motherIdSS")[0] > 0;
    bool lep2_real = tx->getBranch<vector<int>>("lep_motherIdSS")[1] > 0;
    bool lep3_real = false;
    if (tx->getBranch<int>("nLlep") >= 3) lep3_real = tx->getBranch<vector<int>>("lep_motherIdSS")[2] > 0;
    vector<int> reallepindex;
    for (unsigned int lepindex = 0; lepindex < tx->getBranch<vector<LV>>("lep_p4").size(); ++lepindex)
    {
        if (tx->getBranch<vector<int>>("lep_motherIdSS")[lepindex] > 0)
            reallepindex.push_back(lepindex);
        else if (tx->getBranch<vector<int>>("lep_motherIdSS")[lepindex] == -3)
            gammafake = true;
        else
            jetfake = true;
        if (tx->getBranch<vector<int>>("lep_isFromW")[lepindex]) nW++;
        if (tx->getBranch<vector<int>>("lep_isFromZ")[lepindex]) nZ++;
    }
    if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
    {
        std::cout <<  " gammafake: " << gammafake <<  std::endl;
        std::cout <<  " jetfake: " << jetfake <<  std::endl;
        std::cout <<  " ngenlep: " << ngenlep <<  std::endl;
        std::cout <<  " nW: " << nW <<  " nZ: " << nZ <<  std::endl;
        std::cout <<  " lep1_real: " << lep1_real <<  std::endl;
        std::cout <<  " lep2_real: " << lep2_real <<  std::endl;
        std::cout <<  " lep3_real: " << lep3_real <<  std::endl;
    }
    //found two real leptons
    if (tx->getBranch<int>("nLlep") == 2)
    {
        bool ischargeflip = false;
        bool isSS = false;
        if (lep1_real && lep2_real)
        {
            int ilep1 = tx->getBranch<vector<int>>("lep_genPart_index").at(0);
            int ilep2 = tx->getBranch<vector<int>>("lep_genPart_index").at(1);
            bool lep1_chargeflip  = tx->getBranch<vector<int>>("genPart_charge").at(ilep1) != tx->getBranch<vector<int>>("lep_charge").at(0);
            bool lep2_chargeflip  = tx->getBranch<vector<int>>("genPart_charge").at(ilep2) != tx->getBranch<vector<int>>("lep_charge").at(1);
            if (!lep1_chargeflip && !lep2_chargeflip && nW == 2) return 0; // true SS
            else if (!lep1_chargeflip && !lep2_chargeflip)        isSS = true; // true SS - but could be still lost lepton WZ
            if (lep1_chargeflip || lep2_chargeflip)        ischargeflip = true;      // charge flip
        }
        if (ngenlep > 2 || reallepindex.size() > 2 || (nW > 0 && nZ > 0)) return 3; // lostlep
        if ((ngenlep < 2 || !lep1_real || !lep2_real) &&       jetfake) return 4; // jetfake - if double fake with one jet fake and one gamma fake call it jet fake
        if ((ngenlep < 2 || !lep1_real || !lep2_real) &&     gammafake) return 5; // gammafake
        if ((ngenlep < 2 || !lep1_real || !lep2_real) &&    !gammafake) return 4; // call all without gamma fake jetfake - safety cut
        if (isSS)         return 0;
        if (ischargeflip) return 2;
        cout << "This event was not classified - 2 lepton event - v2" << endl;
        return 1;
    }
    else
    {
        //found three real leptons
        bool ischargeflip = false;
        bool isthreelep = false;
        if (lep1_real && lep2_real && lep3_real)
        {
            int ilep1 = tx->getBranch<vector<int>>("lep_genPart_index").at(0);
            int ilep2 = tx->getBranch<vector<int>>("lep_genPart_index").at(1);
            int ilep3 = tx->getBranch<vector<int>>("lep_genPart_index").at(2);
            bool lep1_chargeflip  = tx->getBranch<vector<int>>("genPart_charge").at(ilep1) != tx->getBranch<vector<int>>("lep_charge").at(0);
            bool lep2_chargeflip  = tx->getBranch<vector<int>>("genPart_charge").at(ilep2) != tx->getBranch<vector<int>>("lep_charge").at(1);
            bool lep3_chargeflip  = tx->getBranch<vector<int>>("genPart_charge").at(ilep3) != tx->getBranch<vector<int>>("lep_charge").at(2);
            if (!lep1_chargeflip && !lep2_chargeflip && !lep3_chargeflip && nW == 3) return 0; // true WWW
            else if (!lep1_chargeflip && !lep2_chargeflip && !lep3_chargeflip)        isthreelep = true; // true 3l, but could be lost lepton ZZ
            if (lep1_chargeflip || lep2_chargeflip || lep3_chargeflip)        ischargeflip = true;     // charge flip
        }
        if (ngenlep > 3 || reallepindex.size() > 3 || (nW >= 2 && nZ >= 1) || (nZ >= 3)) return 3; // lostlep (2 lep from W and 2 from Z, or 4 from Z)
        //there is the case of having WZZ with two lost leptons --> ngenlep>3 - correctly put has lostlep
        if ((ngenlep < 3 || !lep1_real || !lep2_real || !lep3_real) &&        jetfake) return 4; // jetfake
        if ((ngenlep < 3 || !lep1_real || !lep2_real || !lep3_real) &&      gammafake) return 5; // gammafake
        if ((ngenlep < 3 || !lep1_real || !lep2_real || !lep3_real) &&     !gammafake) return 4; // jetfake
        if (isthreelep)   return 1;
        if (ischargeflip) return 2;
        cout << "This event was not classified - 3 lepton event - v2" << endl;
        return 0;
    }
}

//##############################################################################################################
TString babyMaker_v2::gentype_v3()
{

    // Mainly used for 2017 analysis (maybe apply for 2016 for clealiness?)
    // Goal: classify events by the following categories

    // lep_motherIdSS convention (defined in CORE/SSSelections.cc)
    //  2 = charge flip    good lepton
    //  1 = charge correct good lepton (perfect)
    //  0 = unmatched
    // -1 = fake from bottom
    // -2 = fake from charm
    // -3 = fake from photon (I think?)
    // -4 = fkae from light

//    if (eventlist.has(cms3.evt_run(), cms3.evt_lumiBlock(), cms3.evt_event()))
//    {
//        coreGenPart.printParticleOfInterest();
//        std::cout <<  " tx->getBranch<vector<LV>>('lep_p4')[0].pt(): " << tx->getBranch<vector<LV>>("lep_p4")[0].pt() <<  " tx->getBranch<vector<LV>>('lep_p4')[1].pt(): " << tx->getBranch<vector<LV>>("lep_p4")[1].pt() <<  std::endl;
//        std::cout <<  " tx->getBranch<vector<LV>>('lep_p4')[0].phi(): " << tx->getBranch<vector<LV>>("lep_p4")[0].phi() <<  " tx->getBranch<vector<LV>>('lep_p4')[1].phi(): " << tx->getBranch<vector<LV>>("lep_p4")[1].phi() <<  std::endl;
//        std::cout <<  " tx->getBranch<vector<int>>('lep_motherIdSS')[0]: " << tx->getBranch<vector<int>>("lep_motherIdSS")[0] <<  " tx->getBranch<vector<int>>('lep_motherIdSS')[1]: " << tx->getBranch<vector<int>>("lep_motherIdSS")[1] <<  std::endl;
//        std::cout <<  " tx->getBranch<vector<int>>('lep_isFromW')[0]: " << tx->getBranch<vector<int>>("lep_isFromW")[0] <<  " tx->getBranch<vector<int>>('lep_isFromW')[1]: " << tx->getBranch<vector<int>>("lep_isFromW")[1] <<  std::endl;
//        std::cout <<  " tx->getBranch<vector<int>>('lep_isFromZ')[0]: " << tx->getBranch<vector<int>>("lep_isFromZ")[0] <<  " tx->getBranch<vector<int>>('lep_isFromZ')[1]: " << tx->getBranch<vector<int>>("lep_isFromZ")[1] <<  std::endl;
//        if (tx->getBranch<int>("nLlep") >= 3)
//        {
//            std::cout <<  " tx->getBranch<vector<LV>>('lep_p4')[2].pt(): " << tx->getBranch<vector<LV>>("lep_p4")[2].pt() <<  std::endl;
//            std::cout <<  " tx->getBranch<vector<LV>>('lep_p4')[2].phi(): " << tx->getBranch<vector<LV>>("lep_p4")[2].phi() <<  std::endl;
//            std::cout <<  " tx->getBranch<vector<int>>('lep_motherIdSS')[2]: " << tx->getBranch<vector<int>>("lep_motherIdSS")[2] <<  std::endl;
//            std::cout <<  " tx->getBranch<vector<int>>('lep_isFromW')[2]: " << tx->getBranch<vector<int>>("lep_isFromW")[2] <<  std::endl;
//            std::cout <<  " tx->getBranch<vector<int>>('lep_isFromZ')[2]: " << tx->getBranch<vector<int>>("lep_isFromZ")[2] <<  std::endl;
//            std::cout <<  " (tx->getBranch<vector<LV>>('lep_p4')[0]+tx->getBranch<vector<LV>>('lep_p4')[1]).mass(): " << (tx->getBranch<vector<LV>>("lep_p4")[0]+tx->getBranch<vector<LV>>("lep_p4")[1]).mass() <<  std::endl;
//        }
//    }

    int lep1_mid = tx->getBranch<vector<int>>("lep_motherIdSS")[0];
    int lep2_mid = tx->getBranch<vector<int>>("lep_motherIdSS")[1];

    // For SS channel (i.e. nLlep == 2 nVlep == 2)
    if (tx->getBranch<int>("nLlep") == 2 && tx->getBranch<int>("nVlep") == 2)
    {
        // trueSS       = i.e. Irreducible
        // chargeflips  = i.e. charge mis-id
        // SSLL         = i.e. lost-lepton
        // fakes        = i.e. non-prompt
        // photonfakes  = i.e. gamma->ell
        // others       = i.e. unmatched
        if (lep1_mid == 1 && lep2_mid == 1)
        {
            if (sampleIsZX())
                return "SSLL";
            else
                return "trueSS";
        }
        else if (lep1_mid == 2 || lep2_mid == 2)
        {
            return "chargeflips";
        }
        else if (lep1_mid ==-3 || lep2_mid ==-3)
        {
            return "photonfakes";
        }
        else if (lep1_mid <= 0 || lep2_mid <= 0)
        {
            return "fakes";
        }
        else
        {
            return "others";
        }
    }
    // For 3L channel (i.e. nLlep == 3 nVlep == 3)
    else if (tx->getBranch<int>("nLlep") == 3 && tx->getBranch<int>("nVlep") == 3)
    {
        // trueWWW      = i.e. Irreducible
        // 3lLL         = e.g. ZZ-> lost lepton + 3 leptons (has at least one Z)
        // chargeflips  = i.e. charge mis-id
        // fakes        = i.e. non-prompt
        // photonfakes  = i.e. gamma->ell
        // others       = i.e. unmatched
        int lep3_mid = tx->getBranch<int>("nLlep") >= 3 ? tx->getBranch<vector<int>>("lep_motherIdSS")[2] : -999;
        if (lep1_mid == 1 && lep2_mid == 1 && lep3_mid == 1)
        {
            if (sampleIsZX())
                return "3lLL";
            else
                return "trueWWW";
        }
        else if (lep1_mid == 2 || lep2_mid == 2 || lep3_mid == 2)
        {
            return "chargeflips";
        }
        else if (lep1_mid ==-3 || lep2_mid ==-3 || lep3_mid ==-3)
        {
            return "photonfakes";
        }
        else if (lep1_mid <= 0 || lep2_mid <= 0 || lep3_mid <= 0)
        {
            return "fakes";
        }
        else
        {
            return "others";
        }
    }
    else
    {
        return "others";
    }
}


//##############################################################################################################
bool babyMaker_v2::sampleIsZX()
{
    TString fname = looper.getCurrentFileName();
    //if (fname.Contains("VHToNonbb")) return true;
    if (fname.Contains("WZ")) return true;
    if (fname.Contains("ZZ")) return true;
    if (fname.Contains("TZ")) return true;
    if (fname.Contains("tZq")) return true;
    return false;
}

//##############################################################################################################
bool babyMaker_v2::vetophotonprocess()
{
    bool process = tx->getBranch<TString>("bkgtype").EqualTo("photonfakes");
    if (
        (looper.getCurrentFileName().Contains("/WJetsToLNu")
       ||looper.getCurrentFileName().Contains("/DYJets")
       ||looper.getCurrentFileName().Contains("/TT_")
       ||looper.getCurrentFileName().Contains("/TTJets_")
       ||looper.getCurrentFileName().Contains("/WW_")
       ||looper.getCurrentFileName().Contains("/WWTo")
       ||looper.getCurrentFileName().Contains("/WZ_")
       ||looper.getCurrentFileName().Contains("/WZTo")
       )
        &&(process)
       ) return true;
    if (
        (looper.getCurrentFileName().Contains("/WGTo")
       ||looper.getCurrentFileName().Contains("/WGstar")
       ||looper.getCurrentFileName().Contains("/ZGTo")
       ||looper.getCurrentFileName().Contains("/TTGamma")
       ||looper.getCurrentFileName().Contains("/WWG")
       ||looper.getCurrentFileName().Contains("/WZG")
       )
        &&(!process)
       )
        return true;
    return false;
}

//##############################################################################################################
std::tuple<float, float, int> babyMaker_v2::getlepFakeRateandErrorandLooseLepIdx(bool data, int lepton_id_version)
{
    // Retrieve relevant variables
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);
    const vector<float>& lep_coneCorrPt = tx->getBranch<vector<float>>("lep_coneCorrPt", true);

    // If the lepton counts reveal that it is not even an event in the AR then return 0
    const int& nVlep = tx->getBranch<int>("nVlep");
    const int& nLlep = tx->getBranch<int>("nLlep");
    const int& nTlep = tx->getBranch<int>("nTlep");
    if ((nLlep - 1) != nTlep)
        return make_tuple(0., 0., -1);

    // Create an int with -1 or 1 to indicate if 3l or ss event
    int version_control = 1;
    if (nLlep >= 3)
        version_control = -1;

    // This is going to take either -1 or 1 (newID) or 0 (oldID).
    int FR_version = version_control * lepton_id_version;

    // Retrieve the flags for the lepton ID
    const vector<int>& istight = nVlep == 2 ? tx->getBranch<vector<int>>(gconf.wwwcfg["tightid"]) : tx->getBranch<vector<int>>(gconf.wwwcfg["3ltightid"]);
    const vector<int>& isloose = nVlep == 2 ? tx->getBranch<vector<int>>(gconf.wwwcfg["looseid"]) : tx->getBranch<vector<int>>(gconf.wwwcfg["3llooseid"]);

    // Figure out the index by whichever one is the loose but not tight
    int index = -1;
    for (int ilep = 0; ilep < nVlep; ++ilep)
    {
        if (!istight[ilep] && isloose[ilep])
            index = ilep;
    }

    // Sanity check
    if (index < 0)
    {
        std::cout <<  " passCount(tx->getBranch<vector<int>>(gconf.wwwcfg['3llooseid'])): " << passCount(tx->getBranch<vector<int>>(gconf.wwwcfg["3llooseid"])) <<  std::endl;
        std::cout <<  " nVlep: " << nVlep <<  std::endl;
        vector<int> el_idx = coreElectron.index;
        vector<int> mu_idx = coreMuon.index;
        std::cout <<  " el_idx.size(): " << el_idx.size() <<  " mu_idx.size(): " << mu_idx.size() <<  std::endl;
        std::cout <<  " isloose.size(): " << isloose.size() <<  std::endl;
        for (int ilep = 0; ilep < nVlep; ++ilep)
            std::cout <<  " istight[ilep]: " << istight[ilep] <<  " isloose[ilep]: " << isloose[ilep] <<  std::endl;
        FATALERROR(__FUNCTION__);
    }

    // Variables to read out fake rate
    float error = 0;
    float faker = 0;
    float conept = lep_coneCorrPt[index];

    // Get the actual fakerates
    if (data)
    {
        if (abs(lep_pdgId[index]) == 11)
        {
            error = fakerate_el_data_unc(conept, lep_p4[index].Eta(), FR_version);
            faker = fakerate_el_data    (conept, lep_p4[index].Eta(), FR_version);
        }
        else
        {
            error = fakerate_mu_data_unc(conept, lep_p4[index].Eta(), FR_version);
            faker = fakerate_mu_data    (conept, lep_p4[index].Eta(), FR_version);
        }
    }
    else
    {
        if (abs(lep_pdgId[index]) == 11)
        {
            error = fakerate_el_qcd_unc(conept, lep_p4[index].Eta(), FR_version);
            faker = fakerate_el_qcd    (conept, lep_p4[index].Eta(), FR_version);
        }
        else
        {
            error = fakerate_mu_qcd_unc(conept, lep_p4[index].Eta(), FR_version);
            faker = fakerate_mu_qcd    (conept, lep_p4[index].Eta(), FR_version);
        }
    }

    // return
    return make_tuple(faker, error, index);
}

////##############################################################################################################
//std::tuple<float, float> babyMaker_v2::getlepSFandError(int index, int lepton_id_version) // deprecated
//{
//    // Retrieve relevant variables
//    const vector<int>&   lep_pdgId = tx->getBranch<vector<int>>  ("lep_pdgId");
//    const vector<LV>&    lep_p4    = tx->getBranch<vector<LV>>   ("lep_p4");
//    const vector<float>& lep_etaSC = tx->getBranch<vector<float>>("lep_etaSC");
//
//    // If the provided index is out of bound return null
//    if (index < 0) return make_tuple(1., 0.);
//    if (index >= (int) lep_pdgId.size()) return make_tuple(1., 0.);
//
//    // If electron
//    if (abs(lep_pdgId[index]) == 11)
//    {
//        // Parametrized kinematic quantity for scale factors
//        float pt = lep_p4[index].Pt();
//        float eta = lep_etaSC[index];
//        // Scale Factors
//        float sfreco   = lepsf_EGammaReco               (pt, eta);
//        float sftight  = lepsf_EGammaTightID            (pt, eta);
//        float sfWWW    = lepsf_EGammaTightPOG_EGammaVVV (pt, eta, lepton_id_version);
//        float sfWWWIso = lepsf_EGammaVVV_Isolation      (pt, eta, lepton_id_version);
//        float sf       = sfreco * sftight * sfWWW * sfWWWIso;
//        // Errors
//        float errreco   = lepsf_EGammaReco_unc           (pt, eta);
//        float errtight  = lepsf_EGammaTightID_unc        (pt, eta);
//        float errWWW    = fabs(lepsf_EGammaTightPOG_EGammaVVV(pt, eta, lepton_id_version, 1) - sfWWW   ); // difference against syst = +1 is the error
//        float errWWWIso = fabs(lepsf_EGammaVVV_Isolation     (pt, eta, lepton_id_version, 1) - sfWWWIso); // difference against syst = +1 is the error
//        // Following the formula from https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas, cause I always forget
//        float fracerrreco   = sfreco   > 0 ? errreco   / sfreco   : 0;
//        float fracerrtight  = sftight  > 0 ? errtight  / sftight  : 0;
//        float fracerrWWW    = sfWWW    > 0 ? errWWW    / sfWWW    : 0;
//        float fracerrWWWIso = sfWWWIso > 0 ? errWWWIso / sfWWWIso : 0;
//        float fracerrtotal  = sqrt(pow(fracerrreco, 2) + pow(fracerrtight, 2) + pow(fracerrWWW, 2) + pow(fracerrWWWIso, 2));
//        // Assume the "tight" and "WWW" error (which are MVA based for lepton_id_version >= 2) are 100% correlated
//        if (lepton_id_version >= 2)
//            fracerrtotal = sqrt(pow(fracerrreco, 2) + pow(fracerrtight + fracerrWWW, 2) + pow(fracerrWWWIso, 2));
//        float error         = fracerrtotal * sf;
//        // Return the value
//        return make_tuple(sf, error);
//    }
//    else  //muon
//    {
//        // Parametrized kinematic quantity for scale factors
//        float pt = lep_p4[index].Pt();
//        float eta = lep_p4[index].Eta();
//        // Scale Factors
//        float sfreco   = lepsf_MuReco                (pt, eta);
//        float sftight1 = lepsf_MuMediumID_BtoF       (pt, eta);
//        float sftight2 = lepsf_MuMediumID_GH         (pt, eta);
//        float sfWWW    = lepsf_MuMediumPOG_MuTightVVV(pt, eta, lepton_id_version);
//        float sftight  = sftight1 * 0.549833 + sftight2 * 0.450167; //luminosity weights B-F vs. G+H
//        float sf       = sftight * sfWWW;
//        // Errors
//        // Recommendation from https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults and linked twiki
//        // From the recommendations what we conclude is to add is a 1% additional uncertainty.
//        float errreco   = lepsf_MuReco_unc         (pt, eta);
//        float errtight1 = lepsf_MuMediumID_BtoF_unc(pt, eta); // Additional 1%
//        float errtight2 = lepsf_MuMediumID_GH_unc  (pt, eta); // Additional 1%
//        float errWWW    = fabs(lepsf_MuMediumPOG_MuTightVVV(pt, eta, lepton_id_version, 1) - sfWWW); // difference against syst = +1 is the error
//        float errtight  = sqrt(pow(errtight1 * 0.549833, 2) + pow(errtight2 * 0.450167, 2));
//        // Following the formula from https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas, cause I always forget
//        float fracerrreco  = sfreco  > 0 ? errreco  / sfreco : 0;
//        float fracerrtight = sftight > 0 ? errtight / sftight : 0;
//        float fracerrWWW   = sfWWW   > 0 ? errWWW   / sfWWW   : 0;
//        float fracerrtotal = sqrt(pow(fracerrreco, 2) + pow(fracerrtight, 2) + pow(fracerrWWW, 2) + pow(0.01, 2)); // flat additional 1% uncertainty
//        float error        = fracerrtotal * sf;
//        // Return the value
//        return make_tuple(sf, error);
//    }
//}

//##############################################################################################################
std::tuple<float, float> babyMaker_v2::getlepSFandError(int index, int lepton_id_version)
{
    // Retrieve relevant variables
    const vector<int>&   lep_pdgId = tx->getBranch<vector<int>>  ("lep_pdgId", true);
    const vector<LV>&    lep_p4    = tx->getBranch<vector<LV>>   ("lep_p4", true);
    const vector<float>& lep_etaSC = tx->getBranch<vector<float>>("lep_etaSC", true);

    // If the provided index is out of bound return null
    if (index < 0) return make_tuple(1., 0.);
    if (index >= (int) lep_pdgId.size()) return make_tuple(1., 0.);

    // Warning message if the leptonid version is not 2
    if (abs(lepton_id_version) != 2)
        RooUtil::error("getlepSFandError: lepton_id_version != 2");

    // If electron
    if (abs(lep_pdgId[index]) == 11)
    {
        // Parametrized kinematic quantity for scale factors
        float pt = lep_p4[index].Pt();
        float eta = lep_etaSC[index];
        // Scale Factors
        float sfreco   = lepsf_elec_reco(pt, eta);
        float sfpogid  = lepton_id_version > 0 ? lepsf_elec_mva80 (pt, eta) : lepsf_elec_mva90 (pt, eta);
        float sfWWW    = lepton_id_version > 0 ? lepsf_elec_ss_id (pt, eta) : lepsf_elec_3l_id (pt, eta);
        float sfWWWIso = lepton_id_version > 0 ? lepsf_elec_ss_iso(pt, eta) : lepsf_elec_3l_iso(pt, eta);
        float sf       = sfreco * sfpogid * sfWWW * sfWWWIso;
        // Errors
        float errreco   = lepsf_elec_reco(pt, eta, 1) - sfreco;
        float errpogid  = (lepton_id_version > 0 ? lepsf_elec_mva80 (pt, eta, 1) : lepsf_elec_mva90 (pt, eta, 1)) - sfpogid;
        float errWWW    = (lepton_id_version > 0 ? lepsf_elec_ss_id (pt, eta, 1) : lepsf_elec_3l_id (pt, eta, 1)) - sfWWW; // difference against syst = +1 is the error
        float errWWWIso = (lepton_id_version > 0 ? lepsf_elec_ss_iso(pt, eta, 1) : lepsf_elec_3l_iso(pt, eta, 1)) - sfWWWIso; // difference against syst = +1 is the error
        // Following the formula from https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas, cause I always forget
        float fracerrreco   = sfreco   > 0 ? errreco   / sfreco   : 0;
        float fracerrpogid  = sfpogid  > 0 ? errpogid  / sfpogid  : 0;
        float fracerrWWW    = sfWWW    > 0 ? errWWW    / sfWWW    : 0;
        float fracerrWWWIso = sfWWWIso > 0 ? errWWWIso / sfWWWIso : 0;
        float fracerrtotal  = sqrt(pow(fracerrreco, 2) + pow(fracerrpogid, 2) + pow(fracerrWWW, 2) + pow(fracerrWWWIso, 2));
        // Assume the "tight" and "WWW" error (which are MVA based for lepton_id_version >= 2) are 100% correlated
        if (lepton_id_version >= 2)
            fracerrtotal = sqrt(pow(fracerrreco, 2) + pow(fracerrpogid + fracerrWWW, 2) + pow(fracerrWWWIso, 2));
        float error         = fracerrtotal * sf;
        // Return the value
        return make_tuple(sf, error);
    }
    else  //muon
    {
        // Parametrized kinematic quantity for scale factors
        float pt = lep_p4[index].Pt();
        float eta = lep_p4[index].Eta();
        // Scale Factors
        float sftrk    = lepsf_muon_trk(pt, eta);
        float sfpogid  = lepsf_muon_id (pt, eta);
        float sfWWW    = lepton_id_version > 0 ? lepsf_muon_ss(pt, eta) : lepsf_muon_3l(pt, eta);
        float sf       = sftrk * sfpogid * sfWWW;
        // Errors
        // Recommendation from https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults and linked twiki
        // From the recommendations what we conclude is to add is a 1% additional uncertainty.
        float errtrk    = lepsf_muon_trk(pt, eta, 1) - sftrk;
        float errpogid  = lepsf_muon_id (pt, eta, 1) - sfpogid;
        float errWWW    = (lepton_id_version > 0 ? lepsf_muon_ss(pt, eta, 1) : lepsf_muon_3l(pt, eta, 1)) - sfWWW;
        // Following the formula from https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas, cause I always forget
        float fracerrtrk   = sftrk   > 0 ? errtrk   / sftrk   : 0;
        float fracerrpogid = sfpogid > 0 ? errpogid / sfpogid : 0;
        float fracerrWWW   = sfWWW   > 0 ? errWWW   / sfWWW   : 0;
        float fracerrtotal = sqrt(pow(fracerrtrk, 2) + pow(fracerrpogid, 2) + pow(fracerrWWW, 2) + pow(0.01, 2)); // flat additional 1% uncertainty
        float error        = fracerrtotal * sf;
        // Return the value
        return make_tuple(sf, error);
    }
}

//##############################################################################################################
std::tuple<float, float> babyMaker_v2::getlepSFWeightandError(int lepton_id_version)
{
    // Retrieve relevant variables
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);

    // Check the bounds
    if (lep_pdgId.size() == 0)
        return make_tuple(1., 0.);

    // Return values
    float SF = 1.;
    float error = 0.;

    // To aggregate the SFs
    vector<float> myeff, myerr;

    // if 3lepton use a different scale factor
    bool is3l = lep_pdgId.size() > 2;
    int multiplier_factor = is3l ? -1 : 1;

    // Loop over leptons and aggregate lepton SFs
    for (unsigned int i = 0; i < lep_pdgId.size(); ++i)
    {
        float tempeff = 0;
        float temperr = 0;
        std::tie(tempeff, temperr) = getlepSFandError(i, multiplier_factor * lepton_id_version);
        myeff.push_back(tempeff);
        myerr.push_back(temperr);
        SF *= tempeff;
    }

    // If SF is = 0
    if (SF == 0)
        return make_tuple(0., 0.); //don't compute error

    // Combine the error by square sum
    for (unsigned int i = 0; i < lep_pdgId.size(); ++i)
        error += pow(SF / myeff[i] * myerr[i], 2);

    // Calculate the error
    error = sqrt(error);

    // Return the result
    return make_tuple(SF, error);
}

//##############################################################################################################
std::tuple<float, float> babyMaker_v2::getTrigEffandError(int lepton_id_version) // Deprecated
{
    return make_tuple(0.0, 0.0);
//    // Retrieve relevant variables
//    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId");
//    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4");
//
//    // If less than two leptons (should never happen anyways)
//    // return dummy value
//    if (lep_pdgId.size() < 2)
//        return make_tuple(0.0, 0.0);
//
//    // Return values
//    float eff = 1.0;
//    float err = 0.0;
//
//    // even if there are more than 2 leptons apply to the first two only
//    float pt0 = lep_p4[0].pt();
//    float pt1 = lep_p4[1].pt();
//    float eta0 = lep_p4[0].eta();
//    float eta1 = lep_p4[1].eta();
//    float leadeta  = std::min(fabs(eta0), (float)2.39);
//    float traileta = std::min(fabs(eta1), (float)2.39);
//    float leadpt   = std::min(pt0, (float)499.);
//    float trailpt  = std::min(pt1, (float)499.);
//
//    float bigeta   = leadeta > traileta ? leadeta  : traileta;
//    float smalleta = leadeta > traileta ? traileta : leadeta;
//
//    // is ee events
//    if (abs(lep_pdgId[0]) == 11 && abs(lep_pdgId[1]) == 11)
//    {
//        // related to lepton legs
//        float e_l0 = trigeff_el_lead(leadeta, leadpt);
//        float e_t1 = trigeff_el_trail(traileta, trailpt);
//        float e_l1 = trigeff_el_lead(traileta, trailpt);
//        float e_t0 = trigeff_el_trail(leadeta, leadpt);
//        float d_l0 = trigeff_el_lead(leadeta, leadpt, 1) - trigeff_el_lead(leadeta, leadpt);
//        float d_t1 = trigeff_el_trail(traileta, trailpt, 1) - trigeff_el_trail(traileta, trailpt);
//        float d_l1 = trigeff_el_lead(traileta, trailpt, 1) - trigeff_el_lead(traileta, trailpt);
//        float d_t0 = trigeff_el_trail(leadeta, leadpt, 1) - trigeff_el_trail(leadeta, leadpt);
//        float tempeff = 1.0;
//        float temperr = 0.0;
//        std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, e_l1, e_t0, e_t1, d_l0, d_l1, d_t0, d_t1);
//        // dz
//        float dzeff = trigeff_diel_dz(smalleta, bigeta);
//        float dzerr = trigeff_diel_dz(smalleta, bigeta, 1) - trigeff_diel_dz(smalleta, bigeta, 1);
//        eff = tempeff * dzeff;
//        err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
//    }
//
//    // emu trigger's DZ filter was near 100% given statistics error also same-sign analysis observes the same.
//    // So apply only a flat err of 2%
//
//    // is em events
//    if (abs(lep_pdgId[0]) == 11 && abs(lep_pdgId[1]) == 13)
//    {
//        // related to lepton legs
//        float e_l0 = trigeff_el_lead(leadeta, leadpt);
//        float e_t1 = trigeff_mu_trail(traileta, trailpt);
//        float d_l0 = trigeff_el_lead(leadeta, leadpt, 1) - trigeff_el_lead(leadeta, leadpt);
//        float d_t1 = trigeff_mu_trail(traileta, trailpt, 1) - trigeff_mu_trail(traileta, trailpt);
//        float tempeff = 1.0;
//        float temperr = 0.0;
//        std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
//        // dz
//        float dzeff = 1.0;
//        float dzerr = 0.02;
//        eff = tempeff * dzeff;
//        err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
//    }
//
//    // is me events
//    if (abs(lep_pdgId[0]) == 13 && abs(lep_pdgId[1]) == 11)
//    {
//        // nominal
//        float e_l0 = trigeff_mu_lead(leadeta, leadpt);
//        float e_t1 = trigeff_el_trail(traileta, trailpt);
//        float d_l0 = trigeff_mu_lead(leadeta, leadpt, 1) - trigeff_mu_lead(leadeta, leadpt);
//        float d_t1 = trigeff_el_trail(traileta, trailpt, 1) - trigeff_el_trail(traileta, trailpt);
//        float tempeff = 1.0;
//        float temperr = 0.0;
//        std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
//        // dz
//        float dzeff = 1.0;
//        float dzerr = 0.02;
//        eff = tempeff * dzeff;
//        err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
//    }
//
//    // is mm events
//    if (abs(lep_pdgId[0]) == 13 && abs(lep_pdgId[1]) == 13)
//    {
//        // related to lepton legs
//        float e_l0 = trigeff_mu_lead(leadeta, leadpt);
//        float e_t1 = trigeff_mu_trail(traileta, trailpt);
//        float e_l1 = trigeff_mu_lead(traileta, trailpt);
//        float e_t0 = trigeff_mu_trail(leadeta, leadpt);
//        float d_l0 = trigeff_mu_lead(leadeta, leadpt, 1) - trigeff_mu_lead(leadeta, leadpt);
//        float d_t1 = trigeff_mu_trail(traileta, trailpt, 1) - trigeff_mu_trail(traileta, trailpt);
//        float d_l1 = trigeff_mu_lead(traileta, trailpt, 1) - trigeff_mu_lead(traileta, trailpt);
//        float d_t0 = trigeff_mu_trail(leadeta, leadpt, 1) - trigeff_mu_trail(leadeta, leadpt);
//        float tempeff = 1.0;
//        float temperr = 0.0;
//        std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, e_l1, e_t0, e_t1, d_l0, d_l1, d_t0, d_t1);
//        // dz
//        float dzeff = 0.241 * trigeff_dimu_dz(smalleta, bigeta) + (1 - 0.241) * 1; // Because DZ filter only affects Period H
//        float dzerr = 0.241 * (trigeff_dimu_dz(smalleta, bigeta, 1) - trigeff_dimu_dz(smalleta, bigeta));
//        eff = tempeff * dzeff;
//        err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
//        // And the fractino of period H is calculated from here: http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/Run2_Data2016
//        // 8.636 + 0.221 / 36.814 = 0.241
//    }
//
//    // Return result
//    return make_tuple(eff, err);
}

//##############################################################################################################
std::tuple<float, float> babyMaker_v2::getTrigSFandError(int lepton_id_version)
{
    // Retrieve relevant variables
    const vector<int>& lep_pdgId = tx->getBranch<vector<int>>("lep_pdgId", true);
    const vector<LV>& lep_p4 = tx->getBranch<vector<LV>>("lep_p4", true);

    // If less than two leptons (should never happen anyways)
    // return dummy value
    if (lep_pdgId.size() < 2)
        return make_tuple(0.0, 0.0);

    // Return values
    float eff = 1.0;
    float err = 0.0;

    // even if there are more than 2 leptons apply to the first two only
    float pt0 = lep_p4[0].pt();
    float pt1 = lep_p4[1].pt();
    float eta0 = lep_p4[0].eta();
    float eta1 = lep_p4[1].eta();
    float leadeta  = std::min(fabs(eta0), (float)2.39);
    float traileta = std::min(fabs(eta1), (float)2.39);
    float leadpt   = std::min(pt0, (float)199.);
    float trailpt  = std::min(pt1, (float)199.);

    float bigeta   = leadeta > traileta ? leadeta  : traileta;
    float smalleta = leadeta > traileta ? traileta : leadeta;

    // NOTE: "eff" variable names are really "efficiency scale factors"

    // For Fast Sim the trigger emulator is not used so apply the efficiency
    if (isWHSUSY())
    {
        // is ee events
        if (abs(lep_pdgId[0]) == 11 && abs(lep_pdgId[1]) == 11)
        {
            // related to lepton legs
            float e_l0 = trigdata_el_lead(leadpt, leadeta);
            float e_t1 = trigdata_el_trail(trailpt, traileta);
            float d_l0 = trigdata_el_lead(leadpt, leadeta, 1) - trigdata_el_lead(leadpt, leadeta);
            float d_t1 = trigdata_el_trail(trailpt, traileta, 1) - trigdata_el_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = trigdata_diel_dz(smalleta, bigeta);
            float dzerr = trigdata_diel_dz(smalleta, bigeta, 1) - trigdata_diel_dz(smalleta, bigeta);
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
        }

        // emu trigger's DZ filter was near 100% given statistics error also same-sign analysis observes the same.
        // So apply only a flat err of 2%

        // is em events
        if (abs(lep_pdgId[0]) == 11 && abs(lep_pdgId[1]) == 13)
        {
            // related to lepton legs
            float e_l0 = trigdata_el_lead(leadpt, leadeta);
            float e_t1 = trigdata_mu_trail(trailpt, traileta);
            float d_l0 = trigdata_el_lead(leadpt, leadeta, 1) - trigdata_el_lead(leadpt, leadeta);
            float d_t1 = trigdata_mu_trail(trailpt, traileta, 1) - trigdata_mu_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = 1.0;
            float dzerr = 0.02;
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
        }

        // is me events
        if (abs(lep_pdgId[0]) == 13 && abs(lep_pdgId[1]) == 11)
        {
            // nominal
            float e_l0 = trigdata_mu_lead(leadpt, leadeta);
            float e_t1 = trigdata_el_trail(trailpt, traileta);
            float d_l0 = trigdata_mu_lead(leadpt, leadeta, 1) - trigdata_mu_lead(leadpt, leadeta);
            float d_t1 = trigdata_el_trail(trailpt, traileta, 1) - trigdata_el_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = 1.0;
            float dzerr = 0.02;
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
        }

        // is mm events
        if (abs(lep_pdgId[0]) == 13 && abs(lep_pdgId[1]) == 13)
        {
            // related to lepton legs
            float e_l0 = trigdata_mu_lead(leadpt, leadeta);
            float e_t1 = trigdata_mu_trail(trailpt, traileta);
            float d_l0 = trigdata_mu_lead(leadpt, leadeta, 1) - trigdata_mu_lead(leadpt, leadeta);
            float d_t1 = trigdata_mu_trail(trailpt, traileta, 1) - trigdata_mu_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = 0.241 * trigdata_dimu_dz(smalleta, bigeta) + (1 - 0.241) * 1; // Because DZ filter only affects Period H
            float dzerr = 0.241 * (trigdata_dimu_dz(smalleta, bigeta, 1) - trigdata_dimu_dz(smalleta, bigeta));
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
            // And the fractino of period H is calculated from here: http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/Run2_Data2016
            // 8.636 + 0.221 / 36.814 = 0.241
        }
    }
    else
    {
        // is ee events
        if (abs(lep_pdgId[0]) == 11 && abs(lep_pdgId[1]) == 11)
        {
            // related to lepton legs
            float e_l0 = trigsf_el_lead(leadpt, leadeta);
            float e_t1 = trigsf_el_trail(trailpt, traileta);
            float d_l0 = trigsf_el_lead(leadpt, leadeta, 1) - trigsf_el_lead(leadpt, leadeta);
            float d_t1 = trigsf_el_trail(trailpt, traileta, 1) - trigsf_el_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = trigsf_diel_dz(smalleta, bigeta);
            float dzerr = trigsf_diel_dz(smalleta, bigeta, 1) - trigsf_diel_dz(smalleta, bigeta);
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
        }

        // emu trigger's DZ filter was near 100% given statistics error also same-sign analysis observes the same.
        // So apply only a flat err of 2%

        // is em events
        if (abs(lep_pdgId[0]) == 11 && abs(lep_pdgId[1]) == 13)
        {
            // related to lepton legs
            float e_l0 = trigsf_el_lead(leadpt, leadeta);
            float e_t1 = trigsf_mu_trail(trailpt, traileta);
            float d_l0 = trigsf_el_lead(leadpt, leadeta, 1) - trigsf_el_lead(leadpt, leadeta);
            float d_t1 = trigsf_mu_trail(trailpt, traileta, 1) - trigsf_mu_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = 1.0;
            float dzerr = 0.02;
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
        }

        // is me events
        if (abs(lep_pdgId[0]) == 13 && abs(lep_pdgId[1]) == 11)
        {
            // nominal
            float e_l0 = trigsf_mu_lead(leadpt, leadeta);
            float e_t1 = trigsf_el_trail(trailpt, traileta);
            float d_l0 = trigsf_mu_lead(leadpt, leadeta, 1) - trigsf_mu_lead(leadpt, leadeta);
            float d_t1 = trigsf_el_trail(trailpt, traileta, 1) - trigsf_el_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = 1.0;
            float dzerr = 0.02;
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
        }

        // is mm events
        if (abs(lep_pdgId[0]) == 13 && abs(lep_pdgId[1]) == 13)
        {
            // related to lepton legs
            float e_l0 = trigsf_mu_lead(leadpt, leadeta);
            float e_t1 = trigsf_mu_trail(trailpt, traileta);
            float d_l0 = trigsf_mu_lead(leadpt, leadeta, 1) - trigsf_mu_lead(leadpt, leadeta);
            float d_t1 = trigsf_mu_trail(trailpt, traileta, 1) - trigsf_mu_trail(trailpt, traileta);
            float tempeff = 1.0;
            float temperr = 0.0;
            std::tie(tempeff, temperr) = getCombinedTrigEffandError(e_l0, 0., 0., e_t1, d_l0, 0., 0., d_t1);
            // dz
            float dzeff = 0.241 * trigsf_dimu_dz(smalleta, bigeta) + (1 - 0.241) * 1; // Because DZ filter only affects Period H
            float dzerr = 0.241 * (trigsf_dimu_dz(smalleta, bigeta, 1) - trigsf_dimu_dz(smalleta, bigeta));
            eff = tempeff * dzeff;
            err = eff * sqrt(pow(temperr / tempeff, 2) + pow(dzerr / dzeff, 2));
            // And the fractino of period H is calculated from here: http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/Run2_Data2016
            // 8.636 + 0.221 / 36.814 = 0.241
        }
    }

    // Return result
    return make_tuple(eff, err);
}

std::tuple<float, float> babyMaker_v2::getCombinedTrigEffandError(
        float e_l0,
        float e_l1,
        float e_t0,
        float e_t1,
        float d_l0,
        float d_l1,
        float d_t0,
        float d_t1)
{
    bool scheme_that_i_think_is_correct = true;
    if (scheme_that_i_think_is_correct)
    {
        float e_lt = e_l0 * e_t1;
        float f_l0 = e_l0 > 0 ? d_l0 / e_l0 : 0;
        float f_t1 = e_t1 > 0 ? d_t1 / e_t1 : 0;
        float f_lt = sqrt(pow(f_l0, 2) + pow(f_t1, 2));
        float d_lt = e_lt * f_lt;
        float eff = e_lt;
        float err = d_lt;
        return make_tuple(eff, err);
    }
    else
    {
        float e_lt = e_l0 * e_t1;
        float e_tl = e_t0 * e_l1;
        float f_l0 = e_l0 > 0 ? d_l0 / e_l0 : 0;
        float f_l1 = e_l1 > 0 ? d_l1 / e_l1 : 0;
        float f_t0 = e_t0 > 0 ? d_t0 / e_t0 : 0;
        float f_t1 = e_t1 > 0 ? d_t1 / e_t1 : 0;
        float f_lt = sqrt(pow(f_l0, 2) + pow(f_t1, 2));
        float f_tl = sqrt(pow(f_t0, 2) + pow(f_l1, 2));
        float d_lt = e_lt * f_lt;
        float d_tl = e_tl * f_tl;
        float eff = e_lt + (1.0 - e_lt) * e_tl;
        float err2 = pow(d_lt, 2) + pow(d_tl, 2) + pow((e_lt * e_tl) * sqrt(pow(f_lt, 2) + pow(f_tl, 2)), 2);
        float err = sqrt(err2);
        return make_tuple(eff, err);
    }
}

//---------------==================-----------------==================-----------------=================----------
//---------------==================-----------------==================-----------------=================----------
//---------------==================-----------------==================-----------------=================----------
//---------------==================-----------------==================-----------------=================----------
//---------------==================-----------------==================-----------------=================----------




//
//
// Old version overlap removal settings and the lepton definition used to perform the overlap removal
// This is kept only for documentation purpose only
//
//


//##############################################################################################################
bool babyMaker_v2::isLeptonOverlappingWithJet_OldVersion(int ijet)
{
    bool is_overlapping = false;

    int idx = coreJet.index[ijet];
    float corr = coreJet.corrs[ijet];
    float pt = (cms3.pfjets_p4()[idx] * corr).pt();
    float eta = (cms3.pfjets_p4()[idx] * corr).eta();

    if (!(pt > JET_PT_MIN || (pt > BJET_PT_MIN && cms3.getbtagvalue("pfCombinedInclusiveSecondaryVertexV2BJetTags", idx) >= BJET_CSV_MED))) { return false; }
    if (fabs(eta) > JET_ETA_MAX) { return false; }
    if (!(isLoosePFJet_Summer16_v1(idx))) { return false; }

    for (auto& imu : coreMuon.index)
    {
        if (!(cms3.mus_p4()[imu].pt() > 20))
            continue;

        if (!(isVetoMuonNoIso_OldVersion(imu)))
            continue;

        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfjets_p4()[idx], cms3.mus_p4()[imu]) < 0.4)
        {
            is_overlapping = true;
            break;
        }
    }

    if (is_overlapping)
        return true;

    for (auto& iel : coreElectron.index)
    {
        if (!(cms3.els_p4()[iel].pt() > 20))
            continue;

        if (!(isVetoElectronNoIso_OldVersion(iel)))
            continue;

        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfjets_p4()[idx], cms3.els_p4()[iel]) < 0.4)
        {
            is_overlapping = true;
            break;
        }
    }

    if (is_overlapping)
        return true;

    return false;
}

//##############################################################################################################
// Used to "pre"select leptons from CMS3
bool babyMaker_v2::isVetoMuonNoIso_OldVersion(int idx)
{
    if (!( cms3.mus_p4()[idx].pt()            > 10.    )) return false;
    if (!( isLooseMuonPOG(idx)                         )) return false;
    return true;
}

//##############################################################################################################
// Used to "pre"select leptons from CMS3
bool babyMaker_v2::isVetoElectronNoIso_OldVersion(int idx)
{
    if (!( cms3.els_p4()[idx].pt()            >  10.   )) return false;
    if (!(isVetoElectronPOGspring16noIso_v1(idx) || passElectronSelection_VVV(idx, VVV_MVAbased_tight_noiso))) return false;
    return true;
}

//##############################################################################################################
// FATAL error code
void babyMaker_v2::FATALERROR(const char* funcname)
{
    std::cout <<  " cms3.evt_run(): " << cms3.evt_run() <<  " cms3.evt_lumiBlock(): " << cms3.evt_lumiBlock() <<  " cms3.evt_event(): " << cms3.evt_event() <<  std::endl;
    cout << "FATAL ERROR: I Should never be here! In function : " << funcname << endl;
}

//eof
