#define TreeAnalyser_cxx
#include "TreeAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TreeAnalyser::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TreeAnalyser.C
//      root> TreeAnalyser t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

//////////////////// Criando os Histogramas - Posteriormente estarei criando uma funcao para não precisar declarar todos estes histogramas, diminuindo reduzindo o codigo
///////////////
TH1F *Mu_Number = new TH1F("Mu_Number","Numero de Muons",24,0,12);// Número de Múons X Eventos
Mu_Number->GetXaxis()->SetTitle("N_{#mu}");
Mu_Number->GetYaxis()->SetTitle("Eventos");
TH1F *Mu_Number_Trigger = new TH1F("Mu_Number_Trigger","Numero de Muons",24,0,12);// Número de Múons X Eventos
Mu_Number_Trigger->GetXaxis()->SetTitle("N_{#mu}");
Mu_Number_Trigger->GetYaxis()->SetTitle("Eventos");
TH1F *Mu_Number_Global = new TH1F("Mu_Number_Global","Numero de Muons",24,0,12);// Número de Múons X Eventos
Mu_Number_Global->GetXaxis()->SetTitle("N_{#mu}");
Mu_Number_Global->GetYaxis()->SetTitle("Eventos");
TH1F *Mu_Number_Soft = new TH1F("Mu_Number_Soft","Numero de Muons",24,0,12);// Número de Múons X Eventos
Mu_Number_Soft->GetXaxis()->SetTitle("N_{#mu}");
Mu_Number_Soft->GetYaxis()->SetTitle("Eventos");
TH1F *Mu_Number_Tight = new TH1F("Mu_Number_Tight","Numero de Muons",24,0,12);// Número de Múons X Eventos
Mu_Number_Tight->GetXaxis()->SetTitle("N_{#mu}");
Mu_Number_Tight->GetYaxis()->SetTitle("Eventos");

///////////////

TH1F *Global_Mu_Pt = new TH1F("Global_Mu_Pt","Distribucao de p_{T} dos Muons ",200,0,100);
Global_Mu_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
Global_Mu_Pt->GetYaxis()->SetTitle("Eventos");
TH1F *Global_Mu_Eta = new TH1F("Global_Mu_Eta","Distribucao de #eta dos Muons ",100,-4,4);
Global_Mu_Eta->GetXaxis()->SetTitle("#eta");
Global_Mu_Eta->GetYaxis()->SetTitle("Eventos");
TH1F *Global_Mu_Phi = new TH1F("Global_Mu_Phi","Distribucao de #varphi dos Muons ",100,-4,4);
Global_Mu_Phi->GetXaxis()->SetTitle("#varphi");
Global_Mu_Phi->GetYaxis()->SetTitle("Eventos");

TH1F *Soft_Mu_Pt = new TH1F("Soft_Mu_Pt","Distribucao de p_{T} dos Muons ",200,0,100);
Soft_Mu_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
Soft_Mu_Pt->GetYaxis()->SetTitle("Eventos");
TH1F *Soft_Mu_Eta = new TH1F("Soft_Mu_Eta","Distribucao de #eta dos Muons ",100,-4,4);
Soft_Mu_Eta->GetXaxis()->SetTitle("#eta");
Soft_Mu_Eta->GetYaxis()->SetTitle("Eventos");
TH1F *Soft_Mu_Phi = new TH1F("Soft_Mu_Phi","Distribucao de #varphi dos Muons ",100,-4,4);
Soft_Mu_Phi->GetXaxis()->SetTitle("#varphi");
Soft_Mu_Phi->GetYaxis()->SetTitle("Eventos");

TH1F *Tight_Mu_Pt = new TH1F("Tight_Mu_Pt","Distribucao de p_{T} dos Muons ",200,0,100);
Tight_Mu_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
Tight_Mu_Pt->GetYaxis()->SetTitle("Eventos");
TH1F *Tight_Mu_Eta = new TH1F("Tight_Mu_Eta","Distribucao de #eta dos Muons ",100,-4,4);
Tight_Mu_Eta->GetXaxis()->SetTitle("#eta");
Tight_Mu_Eta->GetYaxis()->SetTitle("Eventos");
TH1F *Tight_Mu_Phi = new TH1F("Tight_Mu_Phi","Distribucao de #varphi dos Muons ",100,-4,4);
Tight_Mu_Phi->GetXaxis()->SetTitle("#varphi");
Tight_Mu_Phi->GetYaxis()->SetTitle("Eventos");
          
TH1F *Pt_cut_Pt_Eta = new TH1F("Pt_cut_Pt_Eta","Distribucao de p_{T} dos Muons ",200,0,100);
Pt_cut_Pt_Eta->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt_cut_Pt_Eta->GetYaxis()->SetTitle("Eventos");
TH1F *Eta_cut_Pt_Eta = new TH1F("Eta_cut_Pt_Eta","Distribucao de #eta dos Muons ",100,-4,4);
Eta_cut_Pt_Eta->GetXaxis()->SetTitle("#eta");
Eta_cut_Pt_Eta->GetYaxis()->SetTitle("Eventos");
TH1F *Phi_cut_Pt_Eta = new TH1F("Phi_cut_Pt_Eta","Distribucao de #varphi dos Muons ",100,-4,4);
Phi_cut_Pt_Eta->GetXaxis()->SetTitle("#varphi");
Phi_cut_Pt_Eta->GetYaxis()->SetTitle("Eventos");

/////////////////////////////////////////////////
TH1F *jPsi_Mass = new TH1F("jPsi_Mass","Distribucao de Massa do J/#psi ",150,0,100);
jPsi_Mass->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
jPsi_Mass->GetYaxis()->SetTitle("Eventos");

TH1F *jPsi_Pt = new TH1F("jPsi_Pt","Distribucao de p_{T} do J/#psi ",200,0,100);
jPsi_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
jPsi_Pt->GetYaxis()->SetTitle("Eventos");

TH1F *jPsi_Mass_cut = new TH1F("jPsi_Mass_cut","Distribucao de Massa do J/#psi ",100,2.6,3.6);
jPsi_Mass_cut->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
jPsi_Mass_cut->GetYaxis()->SetTitle("Eventos");

TH1F *jPsi_Pt_cut = new TH1F("jPsi_Pt_cut","Distribucao de p_{T} do J/#psi ",200,0,100);
jPsi_Pt_cut->GetXaxis()->SetTitle("p_{T} [GeV]");
jPsi_Pt_cut->GetYaxis()->SetTitle("Eventos");

TH1F *jPsi_Mass_cut_Pt = new TH1F("jPsi_Mass_cut_Pt","Distribucao de Massa do J/#psi ",100,2.6,3.6);
jPsi_Mass_cut_Pt->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
jPsi_Mass_cut_Pt->GetYaxis()->SetTitle("Eventos");

TH1F *jPsi_Pt_cut_Pt = new TH1F("jPsi_Pt_cut_Pt","Distribucao de p_{T} do J/#psi ",200,0,100);
jPsi_Pt_cut_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
jPsi_Pt_cut_Pt->GetYaxis()->SetTitle("Eventos");

TH1F *jPsi_Mass_window = new TH1F("jPsi_Mass_window","Distribucao de Massa do J/#psi ",100,2.6,3.6);
jPsi_Mass_window->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
jPsi_Mass_window->GetYaxis()->SetTitle("Eventos");

TH1F *jPsi_Pt_window = new TH1F("jPsi_Pt_window","Distribucao de p_{T} do J/#psi ",200,0,100);
jPsi_Pt_window->GetXaxis()->SetTitle("p_{T} [GeV]");
jPsi_Pt_window->GetYaxis()->SetTitle("Eventos");

///////////////////////////////////////////////  Z

TH1F *Z_mass_test = new TH1F("Z_mass_test","Distribucao de Massa do Z",100,0,140);// Massa do Z
Z_mass_test->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_test->GetYaxis()->SetTitle("Eventos");

TH1F *Z_Pt_test = new TH1F("Z_Pt_test","Distribucao de p_{T} do Z ",200,0,100);
Z_Pt_test->GetXaxis()->SetTitle("p_{T} [GeV]");
Z_Pt_test->GetYaxis()->SetTitle("Eventos");

TH1F *Z_mass_mc = new TH1F("Z_mass_mc","Distribucao de Massa do Z ",100,0,140);// Massa do Z
Z_mass_mc->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mc->GetYaxis()->SetTitle("Eventos");

TH1F *Z_Pt_mc = new TH1F("Z_Pt_mc","Distribucao de p_{T} do Z ",200,0,100);
Z_Pt_mc->GetXaxis()->SetTitle("p_{T} [GeV]");
Z_Pt_mc->GetYaxis()->SetTitle("Eventos");

TH1F *Z_mass_mc_Eta = new TH1F("Z_mass_mc_Eta","Distribucao de Massa do Z ",100,0,140);// Massa do Z
Z_mass_mc_Eta->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mc_Eta->GetYaxis()->SetTitle("Eventos");

TH1F *Z_Pt_mc_Eta = new TH1F("Z_Pt_mc_Eta","Distribucao de p_{T} do Z ",200,0,100);
Z_Pt_mc_Eta->GetXaxis()->SetTitle("p_{T} [GeV]");
Z_Pt_mc_Eta->GetYaxis()->SetTitle("Eventos");

TH1F *Z_mass_mc_Eta_Pt = new TH1F("Z_mass_mc_Eta_Pt","Distribucao de Massa do Z ",100,0,140);// Massa do Z
Z_mass_mc_Eta_Pt->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_mc_Eta_Pt->GetYaxis()->SetTitle("Eventos");

TH1F *Z_Pt_mc_Eta_Pt = new TH1F("Z_Pt_mc_Eta_Pt","Distribucao de p_{T} do Z ",200,0,100);
Z_Pt_mc_Eta_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
Z_Pt_mc_Eta_Pt->GetYaxis()->SetTitle("Eventos");

TH1F *Z_mass_cut_dimuon = new TH1F("Z_mass_cut_dimuon","Distribucao de Massa do Z ",100,0,140);// Massa do Z
Z_mass_cut_dimuon->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_cut_dimuon->GetYaxis()->SetTitle("Eventos");

TH1F *Z_Pt_cut_dimuon = new TH1F("Z_Pt_cut_dimuon","Distribucao de p_{T} do Z",200,0,100);
Z_Pt_cut_dimuon->GetXaxis()->SetTitle("p_{T} [GeV]");
Z_Pt_cut_dimuon->GetYaxis()->SetTitle("Eventos");

TH1F *Z_mass_final = new TH1F("Z_mass_final","Distribucao de Massa do Z ",100,0,140);// Massa do Z
Z_mass_final->GetXaxis()->SetTitle("M_{Z} [GeV]");
Z_mass_final->GetYaxis()->SetTitle("Eventos");

TH1F *Z_Pt_final = new TH1F("Z_Pt_final","Distribucao de p_{T} do Z ",200,0,100);
Z_Pt_final->GetXaxis()->SetTitle("p_{T} [GeV]");
Z_Pt_final->GetYaxis()->SetTitle("Eventos");

///////////////////////////////////////////////  dimuon

TH1F *dimuon_mass_test = new TH1F("dimuon_mass_test","Distribucao de Massa do dimuon ",100,0,140);// Massa do dimuon
dimuon_mass_test->GetXaxis()->SetTitle("M_{#mu #mu} [GeV]");
dimuon_mass_test->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_Pt_test = new TH1F("dimuon_Pt_test","Distribucao de p_{T} do dimuon ",200,0,100);
dimuon_Pt_test->GetXaxis()->SetTitle("p_{T} [GeV]");
dimuon_Pt_test->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_mass_mc = new TH1F("dimuon_mass_mc","Distribucao de Massa do dimuon ",100,0,140);// Massa do dimuon
dimuon_mass_mc->GetXaxis()->SetTitle("M_{#mu #mu} [GeV]");
dimuon_mass_mc->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_Pt_mc = new TH1F("dimuon_Pt_mc","Distribucao de p_{T} do dimuon ",200,0,100);
dimuon_Pt_mc->GetXaxis()->SetTitle("p_{T} [GeV]");
dimuon_Pt_mc->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_mass_mc_Eta = new TH1F("dimuon_mass_mc_Eta","Distribucao de Massa do dimuon",100,0,140);// Massa do dimuon
dimuon_mass_mc_Eta->GetXaxis()->SetTitle("M_{#mu #mu} [GeV]");
dimuon_mass_mc_Eta->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_Pt_mc_Eta = new TH1F("dimuon_Pt_mc_Eta","Distribucao de p_{T} do dimuon ",200,0,100);
dimuon_Pt_mc_Eta->GetXaxis()->SetTitle("p_{T} [GeV]");
dimuon_Pt_mc_Eta->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_mass_mc_Eta_Pt = new TH1F("dimuon_mass_mc_Eta_Pt","Distribucao de Massa do dimuon ",100,0,140);// Massa do dimuon
dimuon_mass_mc_Eta_Pt->GetXaxis()->SetTitle("M_{#mu #mu} [GeV]");
dimuon_mass_mc_Eta_Pt->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_Pt_mc_Eta_Pt = new TH1F("dimuon_Pt_mc_Eta_Pt","Distribucao de p_{T} do dimuon ",200,0,100);
dimuon_Pt_mc_Eta_Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
dimuon_Pt_mc_Eta_Pt->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_mass_cut_dimuon = new TH1F("dimuon_mass_cut_dimuon","Distribucao de Massa do dimuon ",160,0,80);// Massa do dimuon
dimuon_mass_cut_dimuon->GetXaxis()->SetTitle("M_{#mu #mu} [GeV]");
dimuon_mass_cut_dimuon->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_Pt_cut_dimuon = new TH1F("dimuon_Pt_cut_dimuon","Distribucao de p_{T} do dimuon ",200,0,100);
dimuon_Pt_cut_dimuon->GetXaxis()->SetTitle("p_{T} [GeV]");
dimuon_Pt_cut_dimuon->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_mass_final = new TH1F("dimuon_mass_final","Distribucao de Massa do dimuon ",100,0,80);// Massa do dimuon
dimuon_mass_final->GetXaxis()->SetTitle("M_{#mu #mu} [GeV]");
dimuon_mass_final->GetYaxis()->SetTitle("Eventos");

TH1F *dimuon_Pt_final = new TH1F("dimuon_Pt_final","Distribucao de p_{T} do dimuon ",200,0,100);
dimuon_Pt_final->GetXaxis()->SetTitle("p_{T} [GeV]");
dimuon_Pt_final->GetYaxis()->SetTitle("Eventos");

///////////////
/////////       Histogramas para as informacoes dos geradores      ///////////////
TH1F *Gen_mass = new TH1F("Gen_mass","Distribucao de Massa do Gerador ",100,0,130);// Massa do gerador
TH1F *Gen_mass_Z = new TH1F("Gen_mass_Z","Distribucao de Massa do Z ",100,70,130);// Massa do Z
Gen_mass_Z->GetXaxis()->SetTitle("M_Z [GeV]");
Gen_mass_Z->GetYaxis()->SetTitle("Eventos");
///////////////
TH1F *Gen_mass_JPsi = new TH1F("Gen_mass_JPsi","Distribucao de Massa do J/#psi ",100,2.6,3.6);// Massa do JPsi
Gen_mass_JPsi->GetXaxis()->SetTitle("M_{ J/#psi} [GeV]");
Gen_mass_JPsi->GetYaxis()->SetTitle("Eventos");

THStack *hs = new THStack("hs","Histograma");
///////////////contador e histograma para verificar se estou coletando mais de 1 J/psi no mesmo evento
Int_t countJPsi = 0;
TH1F *jPsi_Event = new TH1F("jPsi_Event","Numero de J/#psi ",10,0,10);// 


////////////////////

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      /////////////Aqui coeça a minha analise
      bool goodTriggerEvt1 = true;
      bool goodTriggerEvt2 = true;
      bool goodTriggerEvt3 = true;

	  goodTriggerEvt1 = (((HLTEleMuX >> 17) & 1) == 1) ? true : false; // HLT_Mu30_TKMu11_v*
	  goodTriggerEvt2 = (((HLTEleMuX >> 19) & 1) == 1) ? true : false; // HLT_IsoMu24_v*
	  goodTriggerEvt3 = (((HLTEleMuX >> 22) & 1) == 1) ? true : false; // HLT_TripleMu_12_10_5_v*

      Int_t count1 = 0;
      Int_t count2 = 0;
      Int_t count3 = 0;
      Int_t count4 = 0;

      Mu_Number->Fill(nMu);
      if(!(goodTriggerEvt1 || goodTriggerEvt2 || goodTriggerEvt3))continue;
      std::cout << "Passaram " << nMu << " por um dos tres Triggers" << endl;
      Mu_Number_Trigger->Fill(nMu);
      for(Int_t i = 0; i < nMu ; i++){
        if(((muType->at(i) >> 1) & 1) == 1){ //selecionando global muons
          count1 = count1 + 1;
          Global_Mu_Pt->Fill(muPt->at(i));
          Global_Mu_Eta->Fill(muEta->at(i));
          Global_Mu_Phi->Fill(muPhi->at(i));
          if(muIDbit->at(i) == 3){//selecionando soft muons
            count2 = count2 + 1;
            Soft_Mu_Pt->Fill(muPt->at(i));
            Soft_Mu_Eta->Fill(muEta->at(i));
            Soft_Mu_Phi->Fill(muPhi->at(i));         
          }
          if(muIDbit->at(i) == 2){
            count3 = count3 + 1;
            std::cout << "Temos muons Tight" << endl;
          }
          if(((muIDbit->at(i) >> 2) & 1) == 1){//selecionando tight muons
            count3 = count3 + 1;
            Tight_Mu_Pt->Fill(muPt->at(i));
            Tight_Mu_Eta->Fill(muEta->at(i));
            Tight_Mu_Phi->Fill(muPhi->at(i));            
          }
        }
      }//for i
      std::cout << "---Numero de Global muons: " << count1 << endl;
      std::cout << "------Numero de Soft muons: " << count1 << endl;
      std::cout << "---------Numero de Tight muons: " << count1 << endl;
      Mu_Number_Global->Fill(count1);
      Mu_Number_Soft->Fill(count2);
      Mu_Number_Tight->Fill(count3);

//////////////////////////
      if(count2 >= 2){//Pelo menos 2 global muons
        for(Int_t i = 0; i < nMu; i++){
          if(!(((muIDbit->at(i) >> 1) & 1) == 1))continue;//global
          if(!((muIDbit->at(i) == 3)))continue;//soft
          if(!( muPt->at(i) >= 3.5 &&  (abs(muEta->at(i))< 2.4)))continue;
          Pt_cut_Pt_Eta->Fill(muPt->at(i));//muon
          Eta_cut_Pt_Eta->Fill(muEta->at(i));//muon
          Phi_cut_Pt_Eta->Fill(muPhi->at(i));//muon
          for(Int_t k = i+1; k < nMu; k++){
            if(!(((muIDbit->at(k) >> 1) & 1) == 1))continue;//global
            if(!((muIDbit->at(k) == 3)))continue;//soft
            muon1.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
            muon2.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));
            if(!((muon1.Pt() > 3.5 && (abs(muon1.Eta())< 2.4)) && (muon2.Pt() > 3.5 && (abs(muon2.Eta())< 2.4))))continue;
            if((muCharge->at(i) != muCharge->at(k))){
              jPsi = muon1 + muon2;
              jPsi_Mass->Fill(jPsi.M());
              jPsi_Pt->Fill(jPsi.Pt());
              if(!(((abs(jPsi.M()) - 3.0969)) < 0.5))continue;
                jPsi_Mass_cut->Fill(jPsi.M());
                jPsi_Pt_cut->Fill(jPsi.Pt());
                if(!((jPsi.Pt()) > 8.5))continue;
                  jPsi_Mass_cut_Pt->Fill(jPsi.M());
                  jPsi_Pt_cut_Pt->Fill(jPsi.Pt());
                  if(!(jPsi.M() > 2.6 && jPsi.M() < 3.6))continue;
                  jPsi_Mass_window->Fill(jPsi.M());
                  jPsi_Pt_window->Fill(jPsi.Pt());
                  countJPsi = countJPsi +1;
            }//charge
          }//for k
        }//for i
        jPsi_Event->Fill(countJPsi); 
        countJPsi = 0;
        if(count3 >= 2){
          for(Int_t i = 0; i < nMu; i++){
            for(Int_t k = i+1; k < nMu; k++){
              if(!(((muIDbit->at(i) >> 2) & 1) == 1) && (((muIDbit->at(k) >> 2) & 1) == 1))continue;//Tight muons
              muon3.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
              muon4.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));
              dimuon = muon3 + muon4;
              Z = jPsi + dimuon;
              if((muon3 != muon1) && (muon4 != muon1) && (muon3 != muon2) && (muon4 != muon2)){
                Z_mass_test->Fill(Z.M());
                Z_Pt_test->Fill(Z.Pt());
                dimuon_mass_test->Fill(dimuon.M());
                dimuon_Pt_test->Fill(dimuon.Pt());
                if(!(muCharge->at(i) != muCharge->at(k)))continue;
                  Z_mass_mc->Fill(Z.M());
                  Z_Pt_mc->Fill(Z.Pt());
                  dimuon_mass_mc->Fill(dimuon.M());
                  dimuon_Pt_mc->Fill(dimuon.Pt());
                  if(!( ((abs(muon3.Eta())) < 2.4) && ((abs(muon4.Eta())) < 2.4) ))continue;
                    Z_mass_mc_Eta->Fill(Z.M());
                    Z_Pt_mc_Eta->Fill(Z.Pt());
                    dimuon_mass_mc_Eta->Fill(dimuon.M());
                    dimuon_Pt_mc_Eta->Fill(dimuon.Pt());
                    if(!(((muon3.Pt() > 30) && (muon4.Pt() > 15))||((muon4.Pt() > 30) && (muon3.Pt() > 15))))continue;
                      std::cout<<"-----------------Os valores de pt são: " << muon3.Pt()<< " e "  <<  muon4.Pt() <<endl;
                      std::cout<<"--------- " << muon1.Pt()<<" ----------- " << muon2.Pt() << " A massa do dimuon e "<< dimuon.M() <<endl ;
                      Z_mass_mc_Eta_Pt->Fill(Z.M());
                      Z_Pt_mc_Eta_Pt->Fill(Z.Pt());
                      dimuon_mass_mc_Eta_Pt->Fill(dimuon.M());
                      dimuon_Pt_mc_Eta_Pt->Fill(dimuon.Pt());
                      if(!(dimuon.M() < 80))continue;
                        Z_mass_cut_dimuon->Fill(Z.M());
                        Z_Pt_cut_dimuon->Fill(Z.Pt());
                        dimuon_mass_cut_dimuon->Fill(dimuon.M());
                        dimuon_Pt_cut_dimuon->Fill(dimuon.Pt());
                        if(!((abs(Z.M()-91.1876)) < 25))continue;
                          Z_mass_final->Fill(Z.M());
                          Z_Pt_final->Fill(Z.Pt());
                          dimuon_mass_final->Fill(dimuon.M());
                          dimuon_Pt_final->Fill(dimuon.Pt());
              }                         
            }//for k
          } //for i
        } // if count 3       
      }//if count2
//////////////////////////////////////////////////////////



   }
}
