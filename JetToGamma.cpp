#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <TLorentzVector.h>

#include "AtlasUtils.C"
#include "AtlasLabels.C"
#include "AtlasStyle.C"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TError.h"


/// МК фоны + сигнал
 const char *fname[111] = {
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361045.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361046.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361047.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361048.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361049.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361050.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361051.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361052.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361053.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361054.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361055.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361056.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361045.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361046.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361047.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361048.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361049.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361050.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361051.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361052.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361053.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361054.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361055.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361056.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361045.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361046.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361047.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361048.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361049.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361050.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361051.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361052.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361053.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361054.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361055.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361056.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_ttgamma_MC16a.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_ttgamma_MC16d.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_ttgamma_MC16e.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16a_361273.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16a_361274.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16a_361275.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16d_361273.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16d_361274.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16d_361275.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16e_361273.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16e_361274.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16e_361275.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16a_364525.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16a_364530.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16a_364535.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16d_364525.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16d_364530.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16d_364535.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16e_364525.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16e_364530.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16e_364535.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364184.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364185.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364186.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364187.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364188.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364189.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364190.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364191.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364192.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364193.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364194.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364195.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364196.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364184.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364185.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364186.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364187.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364188.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364189.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364190.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364191.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364192.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364193.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364194.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364195.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364196.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364184.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364185.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364186.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364187.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364188.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364189.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364190.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364191.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364192.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364193.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364194.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364195.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364196.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16a_364504.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16a_364509.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16a_364514.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16d_364504.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16d_364509.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16d_364514.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16e_364504.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16e_364509.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16e_364514.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16a.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16d.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16e.root" };

  /// данные
   const char *fname_data = "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Data.root";

  /// Работа с кодом (что менять):
  // 1) Изоляцию
  // 2) Рабочие точки
  // 3) Отборы и Предотборы
  // 4) Пределы по изоляции
  // 5) Значение W(enu)

  // Электрослабый анализ (используются только предотброры для увеличеия статистики):
// (fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
// if(mc_ph_type < 13 || mc_ph_type > 16) continue;
   // for signal files we use seletion on mc_ph_type [13, 16]
   // for background we don't use this seletion
// (metTST_pt <= 120) continue; // EWK analysis
// (ph_pt <= 150) continue;
// (n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

// Инклюзивный анализ:
// preselections
// (fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
// (metTST_pt <= 130) continue; // EWK analysis
// (ph_pt <= 150) continue;
// (n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;
// selections
// (metTSTsignif <= 11) continue;
// if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
// if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
// if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;


  bool LoosePrime2 = true;
  bool LoosePrime3 = false;
  bool LoosePrime4 = false;
  bool LoosePrime5 = false;

  int Loose = 0x27fc01;

  bool ElectroWeak = false;
  bool Inclusive = false;
  bool InclusiveSelections = false;

  double sum_A = 0, sum_B = 0, sum_D = 0, sum_C = 0, R_sum = 0, del_R_sum = 0;
  double sum_err_A = 0, sum_err_B = 0, sum_err_C = 0, sum_err_D = 0;
  double sum = 0, ph_pt, errA_data, errB_data, errC_data, errD_data, errB, errD;

  void JetToGamma(){

    double N_A_data, N_B_data, N_C_data, N_D_data;

    TFile *file_data = new TFile(fname_data, "READ");
    cout<<fname_data<<endl;
    TTree *tree = (TTree*)file_data->Get("output_tree");

    TH1D *hist_A_data = new TH1D ("hist_A_data", "hist_A_data", 12, -3, 3);
    TH1D *hist_B_data = new TH1D ("hist_B_data", "hist_B_data", 12, -3, 3);
    TH1D *hist_C_data = new TH1D ("hist_C_data", "hist_C_data", 12, -3, 3);
    TH1D *hist_D_data = new TH1D ("hist_D_data", "hist_D_data", 12, -3, 3);

    double sum = 0, ph_pt, errA_data, errB_data, errC_data, errD_data, errB, errD;
    double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, ph_phi, ph_eta;
    UInt_t ph_isem, n_ph, n_mu, n_e_medium, n_jet;
    Int_t mc_ph_type;
    double jet_lead_phi, jet_sublead_phi, metTST_phi, metTSTsignif, weight, soft_term_pt;
    double jet_lead_eta, jet_lead_pt, jet_lead_E,jet_sublead_pt, jet_sublead_eta, jet_sublead_E;
    TLorentzVector met, ph, jet, jet2;

    tree->SetBranchAddress("weight",&weight);
    tree->SetBranchAddress("ph_pt",&ph_pt);
    tree->SetBranchAddress("ph_phi",&ph_phi);
    tree->SetBranchAddress("ph_eta",&ph_eta);

    tree->SetBranchAddress("jet_lead_pt", &jet_lead_pt);  //leading jet p_x
    tree->SetBranchAddress("jet_lead_eta", &jet_lead_eta);  //p_y
    tree->SetBranchAddress("jet_lead_phi", &jet_lead_phi);  //p_z
    tree->SetBranchAddress("jet_lead_E", &jet_lead_E);    //E

    tree->SetBranchAddress("jet_sublead_pt", &jet_sublead_pt);  //leading jet p_x
    tree->SetBranchAddress("jet_sublead_eta", &jet_sublead_eta);  //p_y
    tree->SetBranchAddress("jet_sublead_phi", &jet_sublead_phi);  //p_z
    tree->SetBranchAddress("jet_sublead_E", &jet_sublead_E);    //E

    tree->SetBranchAddress("metTST_pt", &metTST_pt);  //MET p_x
    tree->SetBranchAddress("metTST_phi", &metTST_phi);  //p_y

    tree->SetBranchAddress("ph_iso_et40", &ph_iso_et40);
    tree->SetBranchAddress("ph_iso_et20", &ph_iso_et20);
    tree->SetBranchAddress("ph_iso_pt", &ph_iso_pt);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("n_ph", &n_ph);
    tree->SetBranchAddress("n_jet", &n_jet);
    tree->SetBranchAddress("n_mu", &n_mu);
    tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
    tree->SetBranchAddress("ph_isem", &ph_isem);
    tree->SetBranchAddress("ph_z_point", &ph_z_point);
    tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
    tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
    tree->SetBranchAddress("soft_term_pt", &soft_term_pt);

    int N_data = (int)tree->GetEntries();

    double MinCut = 2.45, MediumCut = 4.45, MaxCut = 29.45;
    double EtoGam_A = 111, EtoGam_A_err = 0.295466, EtoGam_B = 4.33, EtoGam_B_err = 0.118058;
    double EtoGam_C = 1.38, EtoGam_C_err = 0.02, EtoGam_D = 0.095, EtoGam_D_err = 0.01;


    /// начало цикла по дереву данных

    for(int i = 0; i < N_data; i++){

     tree->GetEntry(i);
     jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
     jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
     met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
     ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

     if(ph_iso_pt/ph_pt >= 0.05 || fabs(ph_z_point)>=250) continue;
     if(ph_pt <= 150) continue;
     if(metTST_pt <= 120) continue;
     if(n_jet < 2) continue;
     if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

     if(metTSTsignif <= 12) continue;
     if(fabs(met.DeltaPhi(jet)) <= 0.3 ) continue;
     if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
     if(fabs(met.DeltaPhi(ph)) <= 0.4) continue;
     if(soft_term_pt >= 16) continue;


     if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 ) hist_A_data->Fill(ph_eta, 1.0);
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B_data->Fill(ph_eta, 1.0);
     else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C_data->Fill(ph_eta, 1.0);
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D_data->Fill(ph_eta, 1.0);
   }

   N_A_data =  hist_A_data->IntegralAndError(1, hist_A_data->GetNbinsX(), errA_data, "");
   N_B_data =  hist_B_data->IntegralAndError(1, hist_B_data->GetNbinsX(), errB_data, "");
   N_C_data =  hist_C_data->IntegralAndError(1, hist_C_data->GetNbinsX(), errC_data, "");
   N_D_data =  hist_D_data->IntegralAndError(1, hist_D_data->GetNbinsX(), errD_data, "");


   cout<<"N_A_data: "<<N_A_data<<" +- "<<errA_data<<endl;
   cout<<"N_B_data: "<<N_B_data<<" +- "<<errB_data<<endl;
   cout<<"N_B_data: "<<N_C_data<<" +- "<<errC_data<<endl;
   cout<<"N_F_data: "<<N_D_data<<" +- "<<errD_data<<endl;

   file_data->Close();





  //------------------------------------------------------------
  //------------------------------------------------------------
  //------------------------------------------------------------
  //--------------------------MK--------------------------------
  //------------------------------------------------------------
  //------------------------------------------------------------
  //------------------------------------------------------------






  double sum_A = 0, sum_B = 0, sum_C = 0, sum_D = 0, R_sum = 0, del_R_sum = 0;
  double sum_err_A = 0, sum_err_B = 0, sum_err_C = 0, sum_err_D = 0;

  double c_A, c_B, c_C, c_D;
  double err_c_A, err_c_B, err_c_C, err_c_D;

  for(int i = 0; i<108; i++){

    char ftempname[104]{};
    sprintf( ftempname, "%s", fname[i] );
    TFile *file = new TFile(ftempname, "READ");
    cout<<ftempname<<endl;

  double sum_of_weights_bk_xAOD, sumw_MC16a = 0, weight, sum = 0, ph_pt, sum_koef = 0, koef;
  double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, ph_phi, ph_eta;
  UInt_t ph_isem, n_ph, n_mu, n_e_medium, n_jet;
  Int_t mc_ph_type;
  double jet_lead_phi, jet_sublead_phi, metTST_phi, metTSTsignif, soft_term_pt;
  double jet_lead_eta, jet_lead_pt, jet_lead_E,jet_sublead_pt, jet_sublead_eta, jet_sublead_E;
  TLorentzVector met, ph, jet, jet2;
  TTree *tree_MC_sw = (TTree*)file->Get("output_tree_sw");
  TTree *tree = (TTree*)file->Get("output_tree");
  TTree *tree_norm = (TTree*)file->Get("norm_tree");
  tree_MC_sw->SetBranchAddress("sum_of_weights_bk_xAOD",&sum_of_weights_bk_xAOD);
  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("ph_pt",&ph_pt);

  tree->SetBranchAddress("ph_phi",&ph_phi);
  tree->SetBranchAddress("ph_eta",&ph_eta);

  tree->SetBranchAddress("jet_lead_pt", &jet_lead_pt);  //leading jet p_x
  tree->SetBranchAddress("jet_lead_eta", &jet_lead_eta);  //p_y
  tree->SetBranchAddress("jet_lead_phi", &jet_lead_phi);  //p_z
  tree->SetBranchAddress("jet_lead_E", &jet_lead_E);    //E

  tree->SetBranchAddress("jet_sublead_pt", &jet_sublead_pt);  //leading jet p_x
  tree->SetBranchAddress("jet_sublead_eta", &jet_sublead_eta);  //p_y
  tree->SetBranchAddress("jet_sublead_phi", &jet_sublead_phi);  //p_z
  tree->SetBranchAddress("jet_sublead_E", &jet_sublead_E);    //E

  tree->SetBranchAddress("metTST_pt", &metTST_pt);  //MET p_x
  tree->SetBranchAddress("metTST_phi", &metTST_phi);  //p_y

  tree->SetBranchAddress("ph_iso_et40", &ph_iso_et40);
  tree->SetBranchAddress("ph_iso_et20", &ph_iso_et20);
  tree->SetBranchAddress("ph_iso_pt", &ph_iso_pt);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("n_ph", &n_ph);
  tree->SetBranchAddress("n_jet", &n_jet);
  tree->SetBranchAddress("n_mu", &n_mu);
  tree->SetBranchAddress("metTST_pt", &metTST_pt);
  tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
  tree->SetBranchAddress("ph_isem", &ph_isem);
  tree->SetBranchAddress("ph_z_point", &ph_z_point);
  tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
  tree_norm->SetBranchAddress("koef",&koef);
  tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
  tree->SetBranchAddress("soft_term_pt", &soft_term_pt);

  int entry = (int)tree_MC_sw->GetEntries();
  int N = (int)tree->GetEntries();
  int N_koef = (int)tree_norm->GetEntries();
  for (int i=0; i<entry; i++) {
   tree_MC_sw->GetEntry(i);
   sumw_MC16a += sum_of_weights_bk_xAOD;
  }
   for(int i = 0; i < 1; i++){
     tree_norm->GetEntry(i);
     sum_koef = koef;
   }


  TH1D *hist_A = new TH1D ("hist_A", "hist_A", 12, -3, 3);
  TH1D *hist_B = new TH1D ("hist_B", "hist_B", 12, -3, 3);
  TH1D *hist_C = new TH1D ("hist_C", "hist_C", 12, -3, 3);
  TH1D *hist_D = new TH1D ("hist_D", "hist_D", 12, -3, 3);


  //LoosePrime2 = ph_isem & 0x27fc01;
  //LoosePrime3 = ph_isem & 0x25fc01;
  //LoosePrime4 = ph_isem & 0x5fc01;
  //LoosePrime5 = ph_isem & 0x1fc01100

  Double_t lumi_mc16a = 36214.96;
  Double_t lumi_mc16d = 44307.4;
  Double_t lumi_mc16e = 58450.1;


  for(int i = 0; i < N; i++){

    tree->GetEntry(i);

    jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
    jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
    met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
    ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

    if(ph_iso_pt/ph_pt >= 0.05 || fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
    //if(mc_ph_type < 13 || mc_ph_type > 16) continue;
    if(ph_pt <= 150) continue;
    if(metTST_pt <= 120) continue;
    if(n_jet < 2) continue;
    if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

    if(metTSTsignif <= 12) continue;
    if(fabs(met.DeltaPhi(jet)) <= 0.3 ) continue;
    if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
    if(fabs(met.DeltaPhi(ph)) <= 0.4) continue;
    if(soft_term_pt >= 16) continue;

    TString new_ftempname = TString(ftempname);
    if(new_ftempname.Contains("MC16a")){
      if((ph_iso_et40 - 0.022*ph_pt) < 2.45 && ph_isem == 0 ) hist_A->Fill(ph_eta, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > 4.45 && (ph_iso_et40 - 0.022*ph_pt) < 29.45 && ph_isem == 0) hist_B->Fill(ph_eta, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) < 2.45 && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C->Fill(ph_eta, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > 4.45 && (ph_iso_et40 - 0.022*ph_pt) < 29.45 && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D->Fill(ph_eta, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
    }
   else if(new_ftempname.Contains("MC16d")){
      if((ph_iso_et40 - 0.022*ph_pt) < 2.45 && ph_isem == 0 ) hist_A->Fill(ph_eta, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > 4.45 && (ph_iso_et40 - 0.022*ph_pt) < 29.45 && ph_isem == 0) hist_B->Fill(ph_eta, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) < 2.45 && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C->Fill(ph_eta, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > 4.45 && (ph_iso_et40 - 0.022*ph_pt) < 29.45 && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D->Fill(ph_eta, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
   }
   else if(new_ftempname.Contains("MC16e")){
      if((ph_iso_et40 - 0.022*ph_pt) < 2.45 && ph_isem == 0 ) hist_A->Fill(ph_eta, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > 4.45 && (ph_iso_et40 - 0.022*ph_pt) < 29.45 && ph_isem == 0) hist_B->Fill(ph_eta, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) < 2.45 && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C->Fill(ph_eta, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > 4.45 && (ph_iso_et40 - 0.022*ph_pt) < 29.45 && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D->Fill(ph_eta, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
   }

  }

  Double_t errA, errB, errC, errD;

  double N_A = hist_A->IntegralAndError(1, hist_A->GetNbinsX(), errA, "");
  double N_B = hist_B->IntegralAndError(1, hist_B->GetNbinsX(), errB, "");
  double N_C = hist_C->IntegralAndError(1, hist_C->GetNbinsX(), errC, "");
  double N_D = hist_D->IntegralAndError(1, hist_D->GetNbinsX(), errD, "");


  cout<<"N_A = "<<N_A<<" +- "<<errA<<endl;
  cout<<"N_B = "<<N_B<<" +- "<<errB<<endl;
  cout<<"N_C = "<<N_C<<" +- "<<errC<<endl;
  cout<<"N_D = "<<N_D<<" +- "<<errD<<endl;

  /// couting sum of events with weights
  sum_A += N_A;
  sum_B += N_B;
  sum_C += N_C;
  sum_D += N_D;

  sum_err_A += errA*errA;
  sum_err_B += errB*errB;
  sum_err_C += errC*errC;
  sum_err_D += errD*errD;

  cout<<"loose'2:"<<endl;
  cout<<"Sum in region A = "<<sum_A<<" +- "<<sqrt(sum_err_A)<<endl;
  cout<<"Sum in region B = "<<sum_B<<" +- "<<sqrt(sum_err_B)<<endl;
  cout<<"Sum in region C = "<<sum_C<<" +- "<<sqrt(sum_err_C)<<endl;
  cout<<"Sum in region D = "<<sum_D<<" +- "<<sqrt(sum_err_D)<<endl;

  c_B = 0.0313752;
  c_C = 0.00847338;
  c_D = 0.000334873;

  err_c_B = sqrt(pow(sqrt(sum_err_A)*sum_B/(sum_A*sum_A) , 2) + pow(sqrt(sum_err_B)/(sum_A), 2));
  err_c_C = sqrt(pow(sqrt(sum_err_A)*sum_C/(sum_A*sum_A) , 2) + pow(sqrt(sum_err_C)/(sum_A), 2));
  err_c_D = sqrt(pow(sqrt(sum_err_A)*sum_D/(sum_A*sum_A) , 2) + pow(sqrt(sum_err_D)/(sum_A), 2));

  cout<<"Leakage parameters:"<<endl;
  cout<<"c_B = "<<c_B<<" +- "<<err_c_B<<endl;
  cout<<"c_C = "<<c_C<<" +- "<<err_c_C<<endl;
  cout<<"c_D = "<<c_D<<" +- "<<err_c_D<<endl;


  if(i == 107){

    sum_A = -sum_A + N_A_data - EtoGam_A;
    sum_err_A = sum_err_A + errA_data*errA_data + EtoGam_A_err*EtoGam_A_err;
    sum_B = -sum_B + N_B_data - EtoGam_B;
    sum_err_B = sum_err_B + errB_data*errB_data + EtoGam_B_err*EtoGam_B_err;
    sum_C = -sum_C + N_C_data - EtoGam_C;
    sum_err_C = sum_err_C + errC_data*errC_data + EtoGam_C_err*EtoGam_C_err;
    sum_D = -sum_D + N_D_data - EtoGam_D;
    sum_err_D = sum_err_D + errD_data*errD_data + EtoGam_D_err*EtoGam_D_err;
  }

  double N_jet_to_gam;
  double b, c, a;
  a = c_D - c_B*c_C;
  b = sum_D + (c_D*sum_A) - (c_B*sum_C + c_C*sum_B);
  c = (sum_D*sum_A) - (sum_C*sum_B);
  N_jet_to_gam = (sum_A) - ((b - sqrt(b*b - (4*a*c)))/(2*a));
  cout<<"Jet to Gam = "<<N_jet_to_gam<<endl;
  cout<<"a = "<<a<<endl;
  cout<<"b = "<<b<<endl;
  cout<<"c = "<<c<<endl;

  cout<<"----Results----"<<endl;
  cout<<"Sum in region A = "<<fabs(sum_A)<<" +- "<<sqrt(sum_err_A)<<endl;
  cout<<"Sum in region B = "<<fabs(sum_B)<<" +- "<<sqrt(sum_err_B)<<endl;
  cout<<"Sum in region C = "<<fabs(sum_C)<<" +- "<<sqrt(sum_err_C)<<endl;
  cout<<"Sum in region D = "<<fabs(sum_D)<<" +- "<<sqrt(sum_err_D)<<endl;

  file->Close();
  }


}
