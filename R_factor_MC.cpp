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


  const char *fname[104] = {
    /*"/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet32.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet33.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet34.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet35.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet36.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet37.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet38.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet39.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet40.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet41.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet42.root",*/
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj364222.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj364223.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366011.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366012.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366013.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366014.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366015.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366016.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366017.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366020.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366021.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366022.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366023.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366024.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366025.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366026.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366029.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366030.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366031.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366032.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366033.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366034.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366035.root",
    /*"/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet31.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet32.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet33.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet34.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet35.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet36.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet37.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet38.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet39.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet40.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet41.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet42.root",*/
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj364222.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj364223.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366011.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366012.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366013.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366014.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366015.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366016.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366017.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366020.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366021.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366022.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366023.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366024.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366025.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366026.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366029.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366030.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366031.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366032.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366033.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366034.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366035.root",
    /*"/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet31.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet32.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet33.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet34.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet35.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet36.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet37.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet38.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet39.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet40.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet41.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet42.root",*/
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj364222.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj364223.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366011.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366012.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366013.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366014.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366015.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366016.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366017.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366020.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366021.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366022.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366023.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366024.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366025.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366026.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366029.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366030.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366031.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366032.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366033.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366034.root",
    "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366035.root"};


 void R_factor_MC(){
   SetAtlasStyle();

   //в программе нужно менять:
   // 1) Рабочие точки
   // 2) Изоляцию
   // 3) Верхний предел по энергии
   // 4) Предотборы и отборы

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


   double sum_A = 0, sum_B = 0, sum_C = 0, sum_D = 0, R_sum = 0, del_R_sum = 0;
   double sum_err_A = 0, sum_err_B = 0, sum_err_C = 0, sum_err_D = 0;
   TString PhotonIsolationName = "FixedCutTight_Tight";
   //TString PhotonIsolationName = "FixedCutTightCaloOnly";
   //TString  PhotonIsolationName = "FixedCutLoose";

   bool LoosePrime2 = true;
   bool LoosePrime3 = false;
   bool LoosePrime4 = false;
   bool LoosePrime5 = false;


   for(int i = 0; i<69; i++){

     char ftempname[104]{};
     sprintf( ftempname, "%s", fname[i] );
     TFile *file = new TFile(ftempname, "READ");
     cout<<ftempname<<endl;

   double sum_of_weights_bk_xAOD, sumw_MC16a = 0, weight, sum = 0, ph_pt, sum_koef = 0, koef;
   double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, metTSTsignif;
   double ph_phi, jet_lead_phi, jet_sublead_phi, metTST_phi, ph_eta;
   UInt_t ph_isem, n_ph, n_mu, n_e_medium;
   Int_t mc_ph_type;
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
   tree->SetBranchAddress("n_mu", &n_mu);
   tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
   tree->SetBranchAddress("ph_isem", &ph_isem);
   tree->SetBranchAddress("ph_z_point", &ph_z_point);
   tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
   tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
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
   for (int i=0; i<N_koef; i++) {
    tree_norm->GetEntry(i);
    sum_koef += koef;
   }


   //cout<<"sumw_MC16a = "<<sumw_MC16a<<endl;
   //cout<<"weight = "<<sum_koef<<endl;

   TH1F *hist_A = new TH1F ("hist_A", "hist_A", 100, -1, 50);
   TH1F *hist_B = new TH1F ("hist_B", "hist_B", 100, -1, 50);
   TH1F *hist_C = new TH1F ("hist_C", "hist_C", 100, -1, 50);
   TH1F *hist_D = new TH1F ("hist_D", "hist_D", 100, -1, 50);


   //LoosePrime2 = ph_isem & 0x27fc01;
   //LoosePrime3 = ph_isem & 0x25fc01;
   //LoosePrime4 = ph_isem & 0x5fc01;
   //LoosePrime5 = ph_isem & 0x1fc01;

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
     if(mc_ph_type >= 13 && mc_ph_type <= 15) continue;
     if(metTST_pt <= 120) continue; // EWK analysis
     //if(metTST_pt <= 130) continue; // Inclusive analysis
     if(ph_pt <= 150) continue;
     if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

     //selections
     //if(metTSTsignif <= 11) continue;
     //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
     // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
     // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

     double MinCut = 2.45, MediumCut = 2.45 + 2.0, MaxCut = 29.45;

     TString new_ftempname = TString(ftempname);
     if(new_ftempname.Contains("MC16a")){
       if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
     }
    else if(new_ftempname.Contains("MC16d")){
       if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
    }
    else if(new_ftempname.Contains("MC16e")){
       if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
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
       if(mc_ph_type >= 13 && mc_ph_type <= 15) continue;
       if(metTST_pt <= 120) continue; // EWK analysis
       //if(metTST_pt <= 130) continue; // Inclusive analysis
       if(ph_pt <= 150) continue;
       if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

       //selections
       //if(metTSTsignif <= 11) continue;
       //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
       // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
       // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

       double MinCut = 2.45, MediumCut = 2.45 + 2.0, MaxCut = 29.45;

       TString new_ftempname = TString(ftempname);
       if(new_ftempname.Contains("MC16a")){
         if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       }
      else if(new_ftempname.Contains("MC16d")){
         if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      }
      else if(new_ftempname.Contains("MC16e")){
         if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut&& (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x1fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
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
       if(mc_ph_type >= 13 && mc_ph_type <= 15) continue;
       if(metTST_pt <= 120) continue; // EWK analysis
       //if(metTST_pt <= 130) continue; // Inclusive analysis
       if(ph_pt <= 150) continue;
       if(n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

       //selections
       //if(metTSTsignif <= 11) continue;
       //if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
       // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
       // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

       double MinCut = 0.0, MediumCut = 0.0 + 2.0, MaxCut = 10000000000.0;

       TString new_ftempname = TString(ftempname);
       if(new_ftempname.Contains("MC16a")){
         if((ph_iso_et20 - 0.065*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et20 - 0.065*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x25fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
         else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x25fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       }
      else if(new_ftempname.Contains("MC16d")){
        if((ph_iso_et20 - 0.065*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et20 - 0.065*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x25fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x25fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      }
      else if(new_ftempname.Contains("MC16e")){
        if((ph_iso_et20 - 0.065*ph_pt) < MinCut && ph_isem == 0 ) hist_A->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et20 - 0.065*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & 0x25fc01) == 0 )) hist_C->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & 0x25fc01) == 0 )) hist_D->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      }

      }


  }

   Double_t errA, errB, errC, errD;

   double N_A = hist_A->IntegralAndError(-1, 50, errA, "");
   double N_B = hist_B->IntegralAndError(-1, 50, errB, "");
   double N_C = hist_C->IntegralAndError(-1, 50, errC, "");
   double N_D = hist_D->IntegralAndError(-1, 50, errD, "");



   double R;
   R = N_A*N_D/(N_C*N_B);

   cout<<"N_A = "<<N_A<<" +- "<<errA<<endl;
   cout<<"N_B = "<<N_B<<" +- "<<errB<<endl;
   cout<<"N_C = "<<N_C<<" +- "<<errC<<endl;
   cout<<"N_D = "<<N_D<<" +- "<<errD<<endl;
   double deltaR;
   deltaR = sqrt(pow(errA*N_D/(N_B*N_C) , 2) + pow(errD*N_A/(N_B*N_C), 2) + pow(errB*N_D*N_A/(N_B*N_C*N_B), 2) + pow(errC*N_D*N_A/(N_B*N_C*N_C) , 2));
   cout<<"R factor = "<<R<<" +- "<<deltaR<<endl;


   /// couting sum of events with weights
   sum_A += N_A;
   sum_B += N_B;
   sum_C += N_C;
   sum_D += N_D;

   sum_err_A += errA*errA;
   sum_err_B += errB*errB;
   sum_err_C += errC*errC;
   sum_err_D += errD*errD;

   cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
   if(LoosePrime2) cout<<"Loose'2: "<<endl;
   else if(LoosePrime3) cout<<"Loose'3: "<<endl;
   else if(LoosePrime4) cout<<"Loose'4: "<<endl;
   else if(LoosePrime5) cout<<"Loose'5: "<<endl;
   cout<<"Events in region A = "<<sum_A<<" +- "<<sqrt(sum_err_A)<<endl;
   cout<<"Events in region B = "<<sum_B<<" +- "<<sqrt(sum_err_B)<<endl;
   cout<<"Events in region C = "<<sum_C<<" +- "<<sqrt(sum_err_C)<<endl;
   cout<<"Events in region D = "<<sum_D<<" +- "<<sqrt(sum_err_D)<<endl;
   R_sum = sum_A*sum_D/(sum_C*sum_B);
   del_R_sum = sqrt(pow(sqrt(sum_err_A)*sum_D/(sum_B*sum_C) , 2) + pow(sqrt(sum_err_D)*sum_A/(sum_B*sum_C), 2) + pow(sqrt(sum_err_B)*sum_D*sum_A/(sum_B*sum_C*sum_B), 2) + pow(sqrt(sum_err_C)*sum_D*sum_A/(sum_B*sum_C*sum_C) , 2));

   cout<<"R factor = "<<R_sum<<" +- "<<del_R_sum<<endl;
   file->Close();
  }
 }
