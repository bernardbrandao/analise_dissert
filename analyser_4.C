{
#include<TROOT.h>
#include<fstream>
#include<TTree.h>
#include<TFile.h>
#include<math.h>
#include<TMath.h>
#include<TLorentzVector.h>
///////////////////////
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH1.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
///////////////////////
#include "Math/Vector4D.h"
#include "THStack.h"
/////////////////////// Utilizando esta classe para selecionar o global muon
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"

TFile *f; 
f = new TFile("ggtree_mc_ALL.root");
TDirectoryFile *pasta = (TDirectoryFile*)f->Get("ggNtuplizer");
TFile *f1= new TFile("tentativa.root","RECREATE");

TTree* tree = (TTree*) pasta->Get("EventTree");
TTree* tree1 = new TTree("tree1","tentativa");	

int nMu, nEle, nMC;
Double_t zMass = 91.2; //em GeV
Double_t muMassMC,jPsiMassMC,ZMassMC;
Double_t muMass,muEt,mu1Pt;
vector <float>  *muPt = 0 ,*muEta = 0, *muEn = 0, *muPhi = 0,*muCharge = 0;
vector <float>  *mcMomMass = 0;
vector<int>     *muType; //https://github.com/cms-cvs-history/DataFormats-MuonReco/blob/master/interface/MuonSelectors.h
vector<int>     *mcMomPID;
TLorentzVector muon;                                  
TLorentzVector muon1, muon2, muon3, muon4, dimuon, jPsi, Z;
vector <UShort_t> *muIDbit;
ULong64_t HLTEleMuX;
///////////////

///////////////

TBranch *bmuPt = 0;
TBranch *bmuEta = 0;
TBranch *bmuEn = 0;
TBranch *bmuPhi = 0;
TBranch *bmuCharge = 0;
TBranch *bnMu = 0;
///////////////
TBranch *bmcMomMass = 0;
TBranch *bnMC = 0;
TBranch *bmuType = 0;
TBranch *bmcMomPID = 0;
TBranch *bmuIDbit = 0;
TBranch *bHLTEleMuX = 0;
///////////////
TBranch *bnEle = 0;
///////////////

///obtendo informacoes das Trees dos muons
tree->SetBranchAddress("nMu",&nMu,&bnMu);
tree->SetBranchAddress("muPt",&muPt,&bmuPt);		
tree->SetBranchAddress("muEta",&muEta,&bmuEta);		
tree->SetBranchAddress("muEn",&muEn,&bmuEn);		 
tree->SetBranchAddress("muPhi",&muPhi,&bmuPhi);	
tree->SetBranchAddress("muCharge",&muCharge,&bmuCharge);
tree->SetBranchAddress("mcMomMass",&mcMomMass,&bmcMomMass);
tree->SetBranchAddress("nMC",&nMC,&bnMC);
tree->SetBranchAddress("muType",&muType,&bmuType);
tree->SetBranchAddress("mcMomPID",&mcMomPID,&bmcMomPID);
tree->SetBranchAddress("muIDbit",&muIDbit,&bmuIDbit);
tree->SetBranchAddress("HLTEleMuX",&HLTEleMuX,&bHLTEleMuX);
///////////////
tree->SetBranchAddress("nEle",&nEle,&bnEle);
///////////////

///////////////
tree1->Branch("nMu",&nMu);
tree1->Branch("muEta",&muPt);
tree1->Branch("muEta",&muEta);
tree1->Branch("muEn",&muEn);
tree1->Branch("muPhi",&muPhi);
tree1->Branch("muCharge",&muCharge);
///////////////
tree1->Branch("nEle",&nEle);
///////////////        Verificar, pois nao esta enchendo
//tree1->Branch("muMassMC",&muMassMC);
//tree1->Branch("jPsiMassMC",&jPsiMassMC);
//tree1->Branch("ZMassMC",&ZMassMC);
///////////////
tree1->Branch("mcMomMass",&mcMomMass);
tree1->Branch("nMC",&nMC);
tree1->Branch("mcMomPID",&mcMomPID);
tree1->Branch("muType",&muType);
tree1->Branch("muIDbit",&muIDbit);
tree1->Branch("HLTEleMuX",&HLTEleMuX);
///////////////

///////////////
TH1F *Muon_Number = new TH1F("Muon_Number","Numero de Muons",24,0,12);// Número de Múons X Eventos
Muon_Number->GetXaxis()->SetTitle("N_{#mu}");
Muon_Number->GetYaxis()->SetTitle("Eventos");
TH1F *Muon_Number_Trigger = new TH1F("Muon_Number_Trigger","Numero de Muons",24,0,12);// Número de Múons X Eventos
Muon_Number_Trigger->GetXaxis()->SetTitle("N_{#mu}");
Muon_Number_Trigger->GetYaxis()->SetTitle("Eventos");
TH1F *Muon_Number_Global = new TH1F("Muon_Number_Global","Numero de Muons",24,0,12);// Número de Múons X Eventos
Muon_Number_Global->GetXaxis()->SetTitle("N_{#mu}");
Muon_Number_Global->GetYaxis()->SetTitle("Eventos");
TH1F *Muon_Number_Soft = new TH1F("Muon_Number_Soft","Numero de Muons",24,0,12);// Número de Múons X Eventos
Muon_Number_Soft->GetXaxis()->SetTitle("N_{#mu}");
Muon_Number_Soft->GetYaxis()->SetTitle("Eventos");
TH1F *Muon_Number_Tight = new TH1F("Muon_Number_Tight","Numero de Muons",24,0,12);// Número de Múons X Eventos
Muon_Number_Tight->GetXaxis()->SetTitle("N_{#mu}");
Muon_Number_Tight->GetYaxis()->SetTitle("Eventos");
///////////////
TH1F *Mass_mu_noCut = new TH1F("Mass_mu_noCut","Distribuic#tilde{a}o de Massa dos Muons ",100,0.1,0.11);// Massa do muons sem Combinações
TH1F *Pt_mu_noCut = new TH1F("Pt_mu_noCut","Distribuic#tilde{a}o de p_{T} dos Muons ",200,0,100);// Pt dos muons sem Combinações
Pt_mu_noCut->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_mu_noCut->GetYaxis()->SetTitle("Eventos");
TH1F *Eta_mu_noCut = new TH1F("Eta_mu_noCut","Distribuic#tilde{a}o de #eta dos Muons ",100,-4,4);// Eta do muons sem Combinações
Eta_mu_noCut->GetXaxis()->SetTitle("#eta");
Eta_mu_noCut->GetYaxis()->SetTitle("Eventos");
TH1F *Phi_mu_noCut = new TH1F("Phi_mu_noCut","Distribuic#tilde{a}o de #varphi dos Muons ",100,-4,4);// Phi dos muons sem Combinações
Phi_mu_noCut->GetXaxis()->SetTitle("#varphi");
Phi_mu_noCut->GetYaxis()->SetTitle("Eventos");
///////////////
//TH1F *Mass_mu_cut_Eta = new TH1F("Mass_mu_cut_Eta","Distribuic#tilde{a}o de Massa dos Muons ",100,0.1,0.11);// Massa do muons com corte de Eta
//TH1F *Pt_mu_cut_Eta = new TH1F("Pt_mu_cut_Eta","Distribuic#tilde{a}o de p_{T} dos Muons ",200,0,100);// Pt dos muons com corte de Eta
//TH1F *Eta_mu_cut_Eta = new TH1F("Eta_mu_cut_Eta","Distribuic#tilde{a}o de #eta dos Muons ",100,-4,4);// Eta do muons com corte de Eta
//TH1F *Phi_mu_cut_Eta = new TH1F("Phi_mu_cut_Eta","Distribuic#tilde{a}o de #varphi dos Muons ",100,-4,4);// Phi dos muons com corte de Eta
///////////////
TH1F *Mass_mu_cut_Eta_Pt = new TH1F("Mass_mu_cut_Eta_Pt","Distribuic#tilde{a}o de Massa dos Muons ",100,0.1,0.11);// Massa do muons com corte de Eta e Pt
TH1F *Pt_mu_cut_Eta_Pt = new TH1F("Pt_mu_cut_Eta_Pt","Distribuic#tilde{a}o de p_{T} dos Muons ",200,0,100);// Pt dos muons com corte de Eta e Pt
Pt_mu_cut_Eta_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_mu_cut_Eta_Pt->GetYaxis()->SetTitle("Eventos");
TH1F *Eta_mu_cut_Eta_Pt = new TH1F("Eta_mu_cut_Eta_Pt","Distribuic#tilde{a}o de #eta dos Muons ",100,-4,4);// Eta do muons com corte de Eta e Pt
Eta_mu_cut_Eta_Pt->GetXaxis()->SetTitle("#eta");
Eta_mu_cut_Eta_Pt->GetYaxis()->SetTitle("Eventos");
TH1F *Phi_mu_cut_Eta_Pt = new TH1F("Phi_mu_cut_Eta_Pt","Distribuic#tilde{a}o de #varphi dos Muons ",100,-4,4);// Phi dos muons com corte de Eta e Pt
Phi_mu_cut_Eta_Pt->GetXaxis()->SetTitle("#varphi");
Phi_mu_cut_Eta_Pt->GetYaxis()->SetTitle("Eventos");
///////////////
TH1F *Pt_mu_Z = new TH1F("Pt_mu_Z","Distribuic#tilde{a}o de p_{T} dos Muons ",200,0,100);// Pt dos muons com corte de Eta e Pt
Pt_mu_Z->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_mu_Z->GetYaxis()->SetTitle("Eventos");
TH1F *Eta_mu_Z = new TH1F("Eta_mu_Z","Distribuic#tilde{a}o de #eta dos Muons ",100,-4,4);// Eta do muons com corte de Eta e Pt
Eta_mu_Z->GetXaxis()->SetTitle("#eta");
Eta_mu_Z->GetYaxis()->SetTitle("Eventos");
TH1F *Phi_mu_Z = new TH1F("Phi_mu_Z","Distribuic#tilde{a}o de #varphi dos Muons ",100,-4,4);// Phi dos muons com corte de Eta e Pt
Phi_mu_Z->GetXaxis()->SetTitle("#varphi");
Phi_mu_Z->GetYaxis()->SetTitle("Eventos");
///////////////
TH1F *Mass_jPsi_DiferentCharge = new TH1F("Mass_jPsi_DiferentCharge","Distribuic#tilde{a}o de Massa dos J/#psi ",150,0,100);// Massa do J/#psi muons de cargas opostas
Mass_jPsi_DiferentCharge->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Mass_jPsi_DiferentCharge->GetYaxis()->SetTitle("Eventos");
TH1F *Mass_jPsi_DiferentCharge_window = new TH1F("Mass_jPsi_DiferentCharge_window","Distribuic#tilde{a}o de Massa dos J/#psi ",150,2.6,3.6);// Massa do J/#psi muons de cargas opostas
Mass_jPsi_DiferentCharge_window->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Mass_jPsi_DiferentCharge_window->GetYaxis()->SetTitle("Eventos");
TH1F *Pt_mu_cut_diferentCharge = new TH1F("Pt_mu_cut_diferentCharge","Distribuic#tilde{a}o de p_{T} dos #mu ",200,0,100);// Pt dos muons com com cargas opostas
Pt_mu_cut_diferentCharge->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_mu_cut_diferentCharge->GetYaxis()->SetTitle("Eventos");
TH1F *Eta_mu_cut_diferentCharge = new TH1F("Eta_mu_cut_diferentCharge","Distribuic#tilde{a}o de #eta dos #mu ",100,-4,4);// Eta do muons com com cargas opostas
Eta_mu_cut_diferentCharge->GetXaxis()->SetTitle("#eta");
Eta_mu_cut_diferentCharge->GetYaxis()->SetTitle("Eventos");
TH1F *Phi_mu_cut_diferentCharge = new TH1F("Phi_mu_cut_diferentCharge","#varphi Ditrisribution of #mu ",100,-4,4);// Phi dos muons com com cargas opostas
Phi_mu_cut_diferentCharge->GetXaxis()->SetTitle("#varphi");
Phi_mu_cut_diferentCharge->GetYaxis()->SetTitle("Eventos");
///////////////
///////////////////////////////////////////////
TH1F *Mass_jPsi_Global = new TH1F("Mass_jPsi_Global","Distribuic#tilde{a}o de Massa dos J/#psi ",150,0,100);// Massa do J/#psi muons de cargas opostas
Mass_jPsi_Global->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Mass_jPsi_Global->GetYaxis()->SetTitle("Eventos");
TH1F *Pt_jPsi_Global = new TH1F("Pt_jPsi_Global","Distribuic#tilde{a}o de p_{T} dos J/#psi ",200,0,100);// Pt dos muons com com cargas opostas
Pt_jPsi_Global->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_jPsi_Global->GetYaxis()->SetTitle("Eventos");
///////////////////////////////////////////////
TH1F *Mass_jPsi_Soft = new TH1F("Mass_jPsi_Soft","Distribuic#tilde{a}o de Massa dos J/#psi ",100,0,76);// Massa do J/#psi muons de cargas opostas
Mass_jPsi_Soft->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Mass_jPsi_Soft->GetYaxis()->SetTitle("Eventos");
TH1F *Pt_jPsi_Soft = new TH1F("Pt_jPsi_Soft","Distribuic#tilde{a}o de p_{T} dos J/#psi ",200,0,60);// Pt dos muons com com cargas opostas
Pt_jPsi_Soft->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_jPsi_Soft->GetYaxis()->SetTitle("Eventos");
///////////////////////////////////////////////
///////////////
TH1F *Mass_jPsi_cut1 = new TH1F("Mass_jPsi_cut1","Distribuic#tilde{a}o de Massa dos J/#psi ",150,2.6,3.5);// Massa do J/#psi muons de cargas opostas
Mass_jPsi_cut1->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Mass_jPsi_cut1->GetYaxis()->SetTitle("Eventos");
TH1F *Pt_jPsi_cut_1 = new TH1F("Pt_jPsi_cut_1","Distribuic#tilde{a}o de p_{T} dos J/#psi ",200,0,100);// Pt dos muons com com cargas opostas
Pt_jPsi_cut_1->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_jPsi_cut_1->GetYaxis()->SetTitle("Eventos");
TH1F *Mass_jPsi_cut3 = new TH1F("Mass_jPsi_cut3","Distribuic#tilde{a}o de Massa dos J/#psi ",150,2.6,3.5);// Massa do J/#psi muons de cargas opostas
Mass_jPsi_cut3->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Mass_jPsi_cut3->GetYaxis()->SetTitle("Eventos");
TH1F *Pt_jPsi_cut_3 = new TH1F("Pt_jPsi_cut_3","Distribuic#tilde{a}o de p_{T} dos J/#psi ",200,0,100);// Pt dos muons com com cargas opostas
Pt_jPsi_cut_3->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_jPsi_cut_3->GetYaxis()->SetTitle("Eventos");
TH1F *Mass_jPsi_cut4 = new TH1F("Mass_jPsi_cut4","Distribuic#tilde{a}o de Massa dos J/#psi ",150,2.6,3.5);// Massa do J/#psi muons de cargas opostas
Mass_jPsi_cut4->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Mass_jPsi_cut4->GetYaxis()->SetTitle("Eventos");
TH1F *Pt_jPsi_cut_4 = new TH1F("Pt_jPsi_cut_4","Distribuic#tilde{a}o de p_{T} dos J/#psi ",200,0,100);// Pt dos muons com com cargas opostas
Pt_jPsi_cut_4->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_jPsi_cut_4->GetYaxis()->SetTitle("Eventos");
///////////////////////////////////////////////
TH1F *Z_mass_mc = new TH1F("Z_mass_mc","Distribuic#tilde{a}o de Massa dos Z ",100,0,140);// Massa do Z
Z_mass_mc->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mc->GetYaxis()->SetTitle("Eventos");
TH1F *Z_mass_mcA = new TH1F("Z_mass_mcA","Distribuic#tilde{a}o de Massa dos Z ",100,70,130);// Massa do Z na janela de massa de 70 a 110
Z_mass_mcA->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mcA->GetYaxis()->SetTitle("Eventos");
TH1F *Z_mass_mcB = new TH1F("Z_mass_mcB","Distribuic#tilde{a}o de Massa dos Z ",100,70,130);// Massa do Z na janela de massa de 70 a 110
Z_mass_mcB->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mcB->GetYaxis()->SetTitle("Eventos");
TH1F *Z_mass_mcC = new TH1F("Z_mass_mcC","Distribuic#tilde{a}o de Massa dos Z ",100,70,130);// Massa do Z na janela de massa de 70 a 110
Z_mass_mcC->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mcC->GetYaxis()->SetTitle("Eventos");
TH1F *Z_mass_mcD = new TH1F("Z_mass_mcD","Distribuic#tilde{a}o de Massa dos Z ",100,70,130);// Massa do Z na janela de massa de 70 a 110
Z_mass_mcD->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mcD->GetYaxis()->SetTitle("Eventos");
TH1F *Pt6 = new TH1F("Pt6","Distribuic#tilde{a}o de p_{T} dos Muons ",200,0,100);// 
TH1F *Eta6 = new TH1F("Eta6","Distribuic#tilde{a}o de #eta dos Muons ",100,-4,4);// 
TH1F *Phi6 = new TH1F("Phi6","Distribuic#tilde{a}o de #varphi dos Muons ",100,-4,4);// 
///////////////
/////////       Histogramas para as informacoes dos geradores      ///////////////
TH1F *Gen_mass = new TH1F("Gen_mass","Distribuic#tilde{a}o de Massa do Gerador ",100,0,130);// Massa do gerador
TH1F *Gen_mass_Z = new TH1F("Gen_mass_Z","Distribuic#tilde{a}o de Massa dos Z ",100,70,130);// Massa do Z
Gen_mass_Z->GetXaxis()->SetTitle("M_Z [GeV]");
Gen_mass_Z->GetYaxis()->SetTitle("Eventos");
///////////////
TH1F *Gen_mass_JPsi = new TH1F("Gen_mass_JPsi","Distribuic#tilde{a}o de Massa dos J/#psi ",100,2.6,3.6);// Massa do JPsi
Gen_mass_JPsi->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Gen_mass_JPsi->GetYaxis()->SetTitle("Eventos");
///////////////
THStack *hs = new THStack("hs","Histograma");
///////////////

Long64_t n = tree->GetEntries();

  for(Int_t j=0;j<n;j++){

    Int_t count1 = 0;
    Int_t count2 = 0;
    Int_t count3 = 0;
    Int_t count4 = 0;

    Long64_t tentry = tree->LoadTree(j); 
    bmuPt->GetEntry(tentry);
    bmuEta->GetEntry(tentry);
    bmuPhi->GetEntry(tentry);
    bmuEn->GetEntry(tentry);
    bnMu->GetEntry(tentry);
    bmuCharge->GetEntry(tentry);

    bnMC->GetEntry(tentry);
    bmcMomMass->GetEntry(tentry);
    bmuType->GetEntry(tentry);
    bmcMomPID->GetEntry(tentry);
    bmuIDbit->GetEntry(tentry);
    bHLTEleMuX->GetEntry(tentry);

    bool goodTriggerEvt1 = true;
    bool goodTriggerEvt2 = true;
    bool goodTriggerEvt3 = true;

	goodTriggerEvt1 = (((HLTEleMuX >> 17) & 1) == 1) ? true : false; // HLT_Mu30_TKMu11_v*
	goodTriggerEvt2 = (((HLTEleMuX >> 19) & 1) == 1) ? true : false; // HLT_IsoMu24_v*
	goodTriggerEvt3 = (((HLTEleMuX >> 22) & 1) == 1) ? true : false; // HLT_TripleMu_12_10_5_v*

    Muon_Number->Fill(nMu);

//  std::cout<<"Evento["<<j<<"]-------------------- Numero de muons "<< nMu <<endl;
//  std::cout<<"Evento["<<j<<"]-------------------- nMonte Carlo "<< nMC <<endl;

    for (Int_t i = 0; i < nMC ; i++) {
      tree1->Fill();
      mcMomMass->push_back(i);
      Gen_mass->Fill(mcMomMass->at(i));
      if((mcMomPID->at(i)) == 23){
        Gen_mass_Z->Fill(mcMomMass->at(i)); 
	  }
      if((mcMomPID->at(i)) == 443){
        Gen_mass_JPsi->Fill(mcMomMass->at(i)); 
      }
    }

///////////////////////////////////////
    int nGoodMuons = 0; //number of good muons

	int leadingMuonIndex = -99; //leading muon index
	bool leadingMuonIsTight = false; //leading muon id (tight)

	int trailingMuonIndex = -99; //leading muon index	
	bool trailingMuonIsTight = false; //leading muon id (tight)

    int teste_index = -99;
    bool teste = false;
    
    Muon_Number->Fill(nMu);
    if(!(goodTriggerEvt1 || goodTriggerEvt2 || goodTriggerEvt3))continue;
//      std::cout << "Passaram " << nMu << " por um dos tres Triggers" << endl;
      Muon_Number_Trigger->Fill(nMu);
      for(Int_t i = 0; i < nMu ; i++){
        if(((muType->at(i) >> 1) & 1) == 1){ //selecionando global muons
          count1 = count1 + 1;
          if(muIDbit->at(i) == 3){//selecionando global muons
            count2 = count2 + 1;
          }
//          if(muIDbit->at(i) == 2){
//            count3 = count3 + 1;
//            std::cout << "Temos muons Tight" << endl;
//          }
          if(((muIDbit->at(i) >> 2) & 1) == 1){//selecionando tight muons
            count3 = count3 + 1;
          }
        }
      }//for i
//      std::cout << "Numero de Global muons: " << count1 << " e Soft muons" << count2 << endl;
      Muon_Number_Global->Fill(count1);
      Muon_Number_Soft->Fill(count2);
      Muon_Number_Tight->Fill(count3);
      if(count2 >= 2){
        for(Int_t i = 0; i < nMu; i++){
          for(Int_t k = i+1; k <nMu; k++){
            muon1.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
            muon2.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));
            if(!((muon1.Pt() > 3.5 && (abs(muon1.Eta())< 2.4)) && (muon2.Pt() > 3.5 && (abs(muon2.Eta())< 2.4))))continue;
            if((muCharge->at(i) != muCharge->at(k))){
              jPsi = muon1 + muon2;
              Mass_jPsi_DiferentCharge->Fill(jPsi.M());
              if(!((jPsi.M()) - 3.0969 < 0.5))continue;
                Mass_jPsi_cut1->Fill(jPsi.M());
                if(!((jPsi.Pt()) > 8.5))continue;
                  Mass_jPsi_cut3->Fill(jPsi.M());
                  if(!((jPsi.M()) > 2.6 && (jPsi.M()) < 3.5))continue;
                    Mass_jPsi_cut4->Fill(jPsi.M());
            }//charge
          }//for k
        }//for i
        if(count3 >= 2){
          for(Int_t i = 0; i < nMu; i++){
            for(Int_t k = i+1; k < nMu; k++){
              if(!(((muIDbit->at(i) >> 2) & 1) == 1) && (((muIDbit->at(k) >> 2) & 1) == 1))continue;
              muon3.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
              muon4.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));
              dimuon = muon3 + muon4;
              Z = jPsi + dimuon;
              if((muon3 != muon1) && (muon4 != muon1) && (muon3 != muon2) && (muon4 != muon2)){
                if(!(muCharge->at(i) != muCharge->at(k)))continue;
                  if(!( ((abs(muon3.Eta())) < 2.4) && ((abs(muon4.Eta())) < 2.4) ))continue;
                    if(!(((muon3.Pt() > 30) && (muon4.Pt() > 15))||((muon4.Pt() > 30) && (muon3.Pt() > 15))))continue;
                      std::cout<<"-----------------Os valores de pt são: " << muon3.Pt()<< " e "  <<  muon4.Pt() <<endl;
                      std::cout<<"--------- " << muon1.Pt()<<" ----------- " << muon2.Pt() << " A massa do dimuon e "<< dimuon.M() <<endl ;
                      Z_mass_mc->Fill(Z.M());
                      if(!(dimuon.M() < 80))continue;
                        Z_mass_mcA->Fill(Z.M());
                        if(!((abs(Z.M()-91.1876)) < 25))continue;
                          Z_mass_mcB->Fill(Z.M());                        
              }                         
            }//for k
          } //for i
        } // if count 3       
      }//if count2
	}//for j
}
