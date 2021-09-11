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
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16a.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16d.root",
   "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16e.root",
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
// (mc_ph_type >= 13 && mc_ph_type <= 15) continue;
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

  //TString PhotonIsolationName = "FixedCutTight_Tight";
  //TString PhotonIsolationName = "FixedCutTightCaloOnly";
  TString  PhotonIsolationName = "FixedCutLoose";

  bool LoosePrime2 = true;
  bool LoosePrime3 = false;
  bool LoosePrime4 = false;
  bool LoosePrime5 = false;

  int Loose = 0x5fc01;


  void R_factor_Data_Driven(){

    double N_B_data, N_D_data, N_E_data, N_F_data;

    TFile *file_data = new TFile(fname_data, "READ");
    cout<<fname_data<<endl;
    TTree *tree = (TTree*)file_data->Get("output_tree");

    TH1F *hist_BE_data = new TH1F ("hist_BE_data", "hist_BE_data", 50, -1, 49);
    TH1F *hist_DF_data = new TH1F ("hist_DF_data", "hist_DF_data", 50, -1, 49);
    TH1F *hist_E_data = new TH1F ("hist_E_data", "hist_E_data", 50, -1, 49);
    TH1F *hist_F_data = new TH1F ("hist_F_data", "hist_F_data", 50, -1, 49);

    double sum = 0, ph_pt, errB_data, errD_data, errE_data, errF_data, errB, errD;
    double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, ph_phi, ph_eta;
    UInt_t ph_isem, n_ph, n_mu, n_e_medium;
    Int_t mc_ph_type;
    double jet_lead_phi, jet_sublead_phi, metTST_phi, metTSTsignif, weight;
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
    tree->SetBranchAddress("n_mu", &n_mu);
    tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
    tree->SetBranchAddress("ph_isem", &ph_isem);
    tree->SetBranchAddress("ph_z_point", &ph_z_point);
    tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
    tree->SetBranchAddress("metTSTsignif", &metTSTsignif);

    int N_data = (int)tree->GetEntries();

    double MinCut, MediumCut, MaxCut;
    double EtoGam_BE, EtoGam_BE_err, EtoGam_E, EtoGam_E_err;
    double EtoGam_DF, EtoGam_DF_err, EtoGam_F, EtoGam_F_err;
    cout<<"Enter low cut: ";
    cin>>MinCut;
    cout<<"Enter medium cut: ";
    cin>>MediumCut;
    cout<<"Enter hight cut: ";
    cin>>MaxCut;
    cout<<"You entered: ["<<MinCut<<","<<MediumCut<<"] and ["<<MediumCut<<","<<MaxCut<<"];"<<endl;
    cout<<"Enter the values of W(enu):"<<endl;
    cout<<"B-E: ";
    cin>>EtoGam_BE;
    cout<<"B-E error: ";
    cin>>EtoGam_BE_err;
    cout<<"E: ";
    cin>>EtoGam_E;
    cout<<"E error: ";
    cin>>EtoGam_E_err;

    cout<<"D-F: ";
    cin>>EtoGam_DF;
    cout<<"D-F error: ";
    cin>>EtoGam_DF_err;
    cout<<"F: ";
    cin>>EtoGam_F;
    cout<<"F error: ";
    cin>>EtoGam_F_err;


    /// начало цикла по дереву данных
  if(PhotonIsolationName.Contains("FixedCutTight_Tight")){

    for(int i = 0; i < N_data; i++){

     tree->GetEntry(i);
     jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
     jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
     met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
     ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

     if(ph_iso_pt/ph_pt >= 0.05) continue;

     //pre-selections
     if(fabs(ph_z_point)>=250) continue;
     if(metTST_pt <= 120) continue; // EWK analysis
     //if(metTST_pt <= 130) continue; // Inclusive analysis
     if(ph_pt <= 150) continue;
     if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

     //selections
     //if(metTSTsignif <= 11) continue;
     //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
     // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
     // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

     if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0 ) hist_BE_data->Fill(n_ph, 1.0);
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_E_data->Fill(n_ph, 1.0);
      else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF_data->Fill(n_ph, 1.0);
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F_data->Fill(n_ph, 1.0);
   }

   N_B_data =  hist_BE_data->IntegralAndError(-1, 49, errB_data, "");
   N_D_data =  hist_DF_data->IntegralAndError(-1, 49, errD_data, "");
   N_E_data =  hist_E_data->IntegralAndError(-1, 49, errE_data, "");
   N_F_data =  hist_F_data->IntegralAndError(-1, 49, errF_data, "");

   cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
   if(LoosePrime2) cout<<"Loose'2: "<<endl;
   else if(LoosePrime3) cout<<"Loose'3: "<<endl;
   else if(LoosePrime4) cout<<"Loose'4: "<<endl;
   else if(LoosePrime5) cout<<"Loose'5: "<<endl;
   cout<<"N_B-E_data: "<<N_B_data<<" +- "<<errB_data<<endl;
   cout<<"N_D-F_data: "<<N_D_data<<" +- "<<errD_data<<endl;
   cout<<"N_E_data: "<<N_E_data<<" +- "<<errE_data<<endl;
   cout<<"N_F_data: "<<N_F_data<<" +- "<<errF_data<<endl;

   file_data->Close();

 } else if(PhotonIsolationName.Contains("FixedCutTightCaloOnly")){

   for(int i = 0; i < N_data; i++){

    tree->GetEntry(i);
    jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
    jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
    met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
    ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

    //pre-selections
    if(fabs(ph_z_point)>=250) continue;
    if(metTST_pt <= 120) continue; // EWK analysis
    //if(metTST_pt <= 130) continue; // Inclusive analysis
    if(ph_pt <= 150) continue;
    if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

    //selections
    //if(metTSTsignif <= 11) continue;
    //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
    // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
    // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

    if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0 ) hist_BE_data->Fill(n_ph, 1.0);
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_E_data->Fill(n_ph, 1.0);
     else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF_data->Fill(n_ph, 1.0);
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F_data->Fill(n_ph, 1.0);

  }

  N_B_data =  hist_BE_data->IntegralAndError(-1, 49, errB_data, "");
  N_D_data =  hist_DF_data->IntegralAndError(-1, 49, errD_data, "");
  N_E_data =  hist_E_data->IntegralAndError(-1, 49, errE_data, "");
  N_F_data =  hist_F_data->IntegralAndError(-1, 49, errF_data, "");

  cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
  if(LoosePrime2) cout<<"Loose'2: "<<endl;
  else if(LoosePrime3) cout<<"Loose'3: "<<endl;
  else if(LoosePrime4) cout<<"Loose'4: "<<endl;
  else if(LoosePrime5) cout<<"Loose'5: "<<endl;
  cout<<"N_B-E_data: "<<N_B_data<<" +- "<<errB_data<<endl;
  cout<<"N_D-F_data: "<<N_D_data<<" +- "<<errD_data<<endl;
  cout<<"N_E_data: "<<N_E_data<<" +- "<<errE_data<<endl;
  cout<<"N_F_data: "<<N_F_data<<" +- "<<errF_data<<endl;

  file_data->Close();

} else if(PhotonIsolationName.Contains("FixedCutLoose")){

  for(int i = 0; i < N_data; i++){

   tree->GetEntry(i);
   jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
   jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
   met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
   ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

   if(ph_iso_pt/ph_pt >= 0.05) continue;

   //pre-selections
   if(fabs(ph_z_point)>=250) continue;
   if(metTST_pt <= 120) continue; // EWK analysis
   //if(metTST_pt <= 130) continue; // Inclusive analysis
   if(ph_pt <= 150) continue;
   if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

   //selections
   //if(metTSTsignif <= 11) continue;
   //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
   // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
   // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

   if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0 ) hist_BE_data->Fill(n_ph, 1.0);
   else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0) hist_E_data->Fill(n_ph, 1.0);
   else if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF_data->Fill(n_ph, 1.0);
   else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F_data->Fill(n_ph, 1.0);

  }

  N_B_data =  hist_BE_data->IntegralAndError(-1, 49, errB_data, "");
  N_D_data =  hist_DF_data->IntegralAndError(-1, 49, errD_data, "");
  N_E_data =  hist_E_data->IntegralAndError(-1, 49, errE_data, "");
  N_F_data =  hist_F_data->IntegralAndError(-1, 49, errF_data, "");

  cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
  if(LoosePrime2) cout<<"Loose'2: "<<endl;
  else if(LoosePrime3) cout<<"Loose'3: "<<endl;
  else if(LoosePrime4) cout<<"Loose'4: "<<endl;
  else if(LoosePrime5) cout<<"Loose'5: "<<endl;
  cout<<"N_B-E_data: "<<N_B_data<<" +- "<<errB_data<<endl;
  cout<<"N_D-F_data: "<<N_D_data<<" +- "<<errD_data<<endl;
  cout<<"N_E_data: "<<N_E_data<<" +- "<<errE_data<<endl;
  cout<<"N_F_data: "<<N_F_data<<" +- "<<errF_data<<endl;
}
   file_data->Close();






  //------------------------------------------------------------
  //------------------------------------------------------------
  //------------------------------------------------------------
  //--------------------------MK--------------------------------
  //------------------------------------------------------------
  //------------------------------------------------------------
  //------------------------------------------------------------


  float sum_BE = 0, sum_DF = 0, sum_E = 0, sum_F = 0, R_sum = 0, del_R_sum = 0;
  float sum_err_BE = 0, sum_err_DF = 0, sum_err_E = 0, sum_err_F = 0;
  float sum_B = 0, sum_D = 0, sum_err_B = 0, sum_err_D = 0;


  for(int i = 0; i<111; i++){

     char ftempname[111]{};
     sprintf( ftempname, "%s", fname[i] );
     TFile *file = new TFile(ftempname, "READ");
     cout<<ftempname<<endl;

   double sum_of_weights_bk_xAOD, sumw_MC16a = 0, weight, sum = 0, ph_pt, sum_koef = 0, koef;
   double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, ph_phi;
   UInt_t ph_isem, n_ph, n_mu, n_e_medium, n_jets;
   double jet_lead_phi, jet_sublead_phi, metTST_phi, metTSTsignif;
   Int_t mc_ph_type;
   double jet_lead_eta, jet_lead_pt, jet_lead_E,jet_sublead_pt, jet_sublead_eta, jet_sublead_E;
   TLorentzVector met, ph, jet, jet2;

   TTree *tree_MC_sw = (TTree*)file->Get("output_tree_sw");
   TTree *tree = (TTree*)file->Get("output_tree");
   TTree *tree_norm = (TTree*)file->Get("norm_tree");
   tree_MC_sw->SetBranchAddress("sum_of_weights_bk_xAOD", &sum_of_weights_bk_xAOD);
   tree->SetBranchAddress("weight",&weight);
   tree->SetBranchAddress("ph_pt",&ph_pt);
   tree->SetBranchAddress("ph_phi",&ph_phi);

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

   tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
   tree->SetBranchAddress("ph_iso_et40", &ph_iso_et40);
   tree->SetBranchAddress("ph_iso_et20", &ph_iso_et20);
   tree->SetBranchAddress("ph_iso_pt", &ph_iso_pt);
   tree->SetBranchAddress("weight", &weight);
   tree->SetBranchAddress("n_ph", &n_ph);
   tree->SetBranchAddress("n_mu", &n_mu);
   tree->SetBranchAddress("n_jet", &n_jets);
   tree->SetBranchAddress("metTST_pt", &metTST_pt);
   tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
   tree->SetBranchAddress("ph_isem", &ph_isem);
   tree->SetBranchAddress("ph_z_point", &ph_z_point);
   tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
   tree_norm->SetBranchAddress("koef",&koef);

   int entry = (int)tree_MC_sw->GetEntries();
   int N = (int)tree->GetEntries();
   int N_koef = (int)tree_norm->GetEntries();
   for (int i=0; i<entry; i++) {
    tree_MC_sw->GetEntry(i);
    sumw_MC16a += sum_of_weights_bk_xAOD;
   }
   for(int i = 0; i < N; i++){
     tree->GetEntry(i);
     sum += weight;
   }

    for(int i = 0; i < 1; i++){
      tree_norm->GetEntry(i);
      sum_koef = koef;
    }


   TH1F *hist_BE = new TH1F ("hist_BE", "hist_BE", 100, -1, 50);
   TH1F *hist_DF = new TH1F ("hist_DF", "hist_DF", 100, -1, 50);
   TH1F *hist_E = new TH1F ("hist_E", "hist_E", 100, -1, 50);
   TH1F *hist_F = new TH1F ("hist_F", "hist_F", 100, -1, 50);

   Double_t lumi_mc16a = 36214.96;
   Double_t lumi_mc16d = 44307.4;
   Double_t lumi_mc16e = 58450.1;

   if(PhotonIsolationName.Contains("FixedCutTight_Tight")){

     //topoetcone40 < 0.022 pT + 2.45
     //Track isolation - ptcone20/pT < 0.05

   for(int i = 0; i < N; i++){

      tree->GetEntry(i);
      jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
      jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
      met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
      ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

      if(ph_iso_pt/ph_pt >= 0.05) continue;

      //pre-selections
      if(fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
      if(metTST_pt <= 120) continue; // EWK analysis
      //if(metTST_pt <= 130) continue; // Inclusive analysis
      if(ph_pt <= 150) continue;
      if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

      //selections
      //if(metTSTsignif <= 11) continue;
      //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
      // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
      // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;


      TString new_ftempname = TString(ftempname);
     if(new_ftempname.Contains("MC16a")){
       if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0  ) hist_BE->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_DF->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem == 0) ) hist_E->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_F->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
     }
    else if(new_ftempname.Contains("MC16d")){
       if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0   ) hist_BE->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_DF->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem == 0) ) hist_E->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_F->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
    }
    else if(new_ftempname.Contains("MC16e")){
       if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0   ) hist_BE->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem == 0)) hist_E->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
    }

     }
   } else if(PhotonIsolationName.Contains("FixedCutTightCaloOnly")){

     //topoetcone40 < 0.022 pT + 2.45
     //Track isolation - none

     for(int i = 0; i < N; i++){

        tree->GetEntry(i);

        jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
        jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
        met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
        ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

        //pre-selections
        if(fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
        if(metTST_pt <= 120) continue; // EWK analysis
        //if(metTST_pt <= 130) continue; // Inclusive analysis
        if(ph_pt <= 150) continue;
        if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

        //selections
        //if(metTSTsignif <= 11) continue;
        //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
        // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
        // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

        TString new_ftempname = TString(ftempname);
       if(new_ftempname.Contains("MC16a")){
         if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0  ) hist_BE->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_DF->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem == 0) ) hist_E->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_F->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       }
      else if(new_ftempname.Contains("MC16d")){
         if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0   ) hist_BE->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_DF->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem == 0) ) hist_E->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) ) hist_F->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      }
      else if(new_ftempname.Contains("MC16e")){
         if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0   ) hist_BE->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem == 0)) hist_E->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      }

       }

   }else if(PhotonIsolationName.Contains("FixedCutLoose")){

     //topoetcone20 < 0.065 pT
     //ptcone20/pT < 0.05

     for(int i = 0; i < N; i++){

        tree->GetEntry(i);
        jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
        jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
        met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
        ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

        if(ph_iso_pt/ph_pt >= 0.05) continue;

        //pre-selections
        if(fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
        if(metTST_pt <= 120) continue; // EWK analysis
        //if(metTST_pt <= 130) continue; // Inclusive analysis
        if(ph_pt <= 150) continue;
        if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

        //selections
        //if(metTSTsignif <= 11) continue;
        //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
        // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
        // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

        TString new_ftempname = TString(ftempname);
        if(new_ftempname.Contains("MC16a")){
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0 ) hist_BE->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem == 0)) hist_E->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        }
       else if(new_ftempname.Contains("MC16d")){
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0 ) hist_BE->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem == 0)) hist_E->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       }
       else if(new_ftempname.Contains("MC16e")){
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0) hist_BE->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem == 0)) hist_E->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
          if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
       }

       }


    }

    double errBE, errDF, errE, errF, errD;

    double N_BE = hist_BE->IntegralAndError(-1, 50, errBE, "");
    double N_DF = hist_DF->IntegralAndError(-1, 50, errDF, "");
    double N_E = hist_E->IntegralAndError(-1, 50, errE, "");
    double N_F = hist_F->IntegralAndError(-1, 50, errF, "");

    double R;
    R = N_BE*N_F/(N_E*N_DF);

    sum_BE += N_BE;
    sum_DF += N_DF;
    sum_E += N_E;
    sum_F += N_F;

    sum_err_BE += errBE*errBE;
    sum_err_DF += errDF*errDF;
    sum_err_E += errE*errE;
    sum_err_F += errF*errF;

    if(i == 110){
      cout<<"Sum in region B-E = "<<sum_BE<<" +- "<<sqrt(sum_err_BE)<<endl;
      cout<<"Sum in region E = "<<sum_E<<" +- "<<sqrt(sum_err_E)<<endl;
      cout<<"Sum in region D-F = "<<sum_DF<<" +- "<<sqrt(sum_err_DF)<<endl;
      cout<<"Sum in region F = "<<sum_F<<" +- "<<sqrt(sum_err_F)<<endl;
      cout<<"Summing for R factor = "<<R_sum<<" +- "<<del_R_sum<<endl;
    }


    if(i == 110){

      sum_BE = sum_BE - N_B_data + EtoGam_BE;
      sum_err_BE = sum_err_BE + errB_data*errB_data + EtoGam_BE_err*EtoGam_BE_err;
      sum_E = sum_E - N_E_data + EtoGam_E;
      sum_err_E = sum_err_E + errE_data*errE_data + EtoGam_E_err*EtoGam_E_err;
      sum_DF = sum_DF - N_D_data + EtoGam_DF;
      sum_err_DF = sum_err_DF + errD_data*errD_data + EtoGam_DF_err*EtoGam_DF_err;
      sum_F = sum_F - N_F_data + EtoGam_F;
      sum_err_F = sum_err_F + errF_data*errF_data + EtoGam_F_err*EtoGam_F_err;
    }

    file->Close();
  }

  R_sum = sum_BE*sum_F/(sum_E*sum_DF);
  del_R_sum = sqrt(pow(sqrt(sum_err_BE)*sum_F/(sum_DF*sum_E) , 2) + pow(sqrt(sum_err_F)*sum_BE/(sum_DF*sum_E), 2)
    + pow(sqrt(sum_err_DF)*sum_F*sum_BE/(sum_DF*sum_E*sum_DF), 2) + pow(sqrt(sum_err_E)*sum_F*sum_BE/(sum_DF*sum_E*sum_E) , 2));

    cout<<"----Results----"<<endl;
    cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
    if(LoosePrime2) cout<<"loose'2: "<<endl;
    if(LoosePrime3) cout<<"loose'3: "<<endl;
    else if(LoosePrime4) cout<<"loose'4: "<<endl;
    else if(LoosePrime5) cout<<"loose'5: "<<endl;
    cout<<"Sum in region B-E = "<<fabs(sum_BE)<<" +- "<<sqrt(sum_err_BE)<<endl;
    cout<<"Sum in region E = "<<fabs(sum_E)<<" +- "<<sqrt(sum_err_E)<<endl;
    cout<<"Sum in region D-F = "<<fabs(sum_DF)<<" +- "<<sqrt(sum_err_DF)<<endl;
    cout<<"Sum in region F = "<<fabs(sum_F)<<" +- "<<sqrt(sum_err_F)<<endl;
    cout<<"Summing for R factor = "<<R_sum<<" +- "<<del_R_sum<<endl;

}
