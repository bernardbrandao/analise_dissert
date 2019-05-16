#define TreeAnalyser_cxx
#include "TreeAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
////////
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "THStack.h"
//#include "deltaPhi.h"
#include "deltaR_deltaPhi.h"
//#include "VertexDistanceXY.h"
//#include "KalmanVertexFitter.h"


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

/////////////////////////////////

      int n1 = 0;
      int n2 = 0;
      int n3 = 0;
      int n4 = 0;

      int m1 = 0;
      int m2 = 0;
      int m3 = 0;
      int m4 = 0;
      int m5 = 0;
      int m6 = 0;

   double w = 36.8 * 7.029 * (pow(10,2)) * (pow(10,3))/ 204091 ;///em (<fb^(-1)>*<fb>)/eventos = <eventos>^(-1)

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
      Int_t indice1 = 0;
      Int_t indice2 = 0;
      Int_t countJPsi = 0;
      Int_t countZ, countZ_2 = 0;
     
      int nGoodMuons = 0; //number of good muons

      Mu_Number->Fill(nMu);
      if(!(goodTriggerEvt1 || goodTriggerEvt2 || goodTriggerEvt3))continue;
      //std::cout << "Passaram " << nMu << " por um dos tres Triggers" << endl;
      Mu_Number_Trigger->Fill(nMu);
      for(Int_t i = 0; i < nMu ; i++){
          Trigger_Mu_Pt->Fill(muPt->at(i));
          Trigger_Mu_Eta->Fill(muEta->at(i));
          Trigger_Mu_Phi->Fill(muPhi->at(i));
        if(((muType->at(i) >> 1) & 1) == 1){ //selecionando global muons
          count1 = count1 + 1;
          Global_Mu_Pt->Fill(muPt->at(i));
          Global_Mu_Eta->Fill(muEta->at(i));
          Global_Mu_Phi->Fill(muPhi->at(i));
          if(((muIDbit->at(i) >> 3) & 1) == 1){//selecionando soft muons
            count2 = count2 + 1;
            Soft_Mu_Pt->Fill(muPt->at(i));
            Soft_Mu_Eta->Fill(muEta->at(i));
            Soft_Mu_Phi->Fill(muPhi->at(i));         
          }
          if(((muIDbit->at(i) >> 2) & 1) == 1){//selecionando tight muons
            count3 = count3 + 1;
            Tight_Mu_Pt->Fill(muPt->at(i));
            Tight_Mu_Eta->Fill(muEta->at(i));
            Tight_Mu_Phi->Fill(muPhi->at(i));            
          }
        }
      }//for i
      Mu_Number_Global->Fill(count1);
      Mu_Number_Soft->Fill(count2);
      Mu_Number_Tight->Fill(count3);

//////////////////////////
      if(!((nMu >= 4)||(nMu>=2 && nEle >=2)))continue;//posso inserir um "|| nMu>=2 && nEle >=2"
      if(!(count2 >= 2))continue;//Pelo menos 2 global muons
      for(Int_t i = 0; i < nMu; i++){
        for(Int_t k = i+1; k < nMu; k++){
          if(!((((muIDbit->at(k) >> 1) & 1) == 1) & (((muIDbit->at(i) >> 1) & 1) == 1)))continue;//global
          if(!((((muIDbit->at(k) >> 3) & 1) == 1) & (((muIDbit->at(i) >> 3) & 1) == 1)))continue;//selecionando soft muons
          if(muPt->at(i) >= muPt->at(k)){
            muon3.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
            muon4.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));            
          }else{
            muon4.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
            muon3.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));
          }
          hist_D0->Fill(muD0->at(i));
          hist_Dz->Fill(muDz->at(i));
          if(!((fabs(muD0->at(i))) < 0.1 & (fabs(muDz->at(i))) < 0.1))continue;// corte de dxy e dz
            hist_D0_cut->Fill(muD0->at(i));
            hist_Dz_cut->Fill(muDz->at(i));
          if((muCharge->at(i) != muCharge->at(k))){
            Pt_muon3->Fill(muon3.Pt());
            Pt_muon4->Fill(muon4.Pt());
            if(!((deltaR(muEta->at(i), muPhi->at(i), muEta->at(k), muPhi->at(k)) < 0.3)))continue;
            if(!((muon3.Pt() > 3.5 && (abs(muon3.Eta())< 2.4)) && (muon4.Pt() > 3.5 && (abs(muon4.Eta())< 2.4))))continue;
            Pt_muon3_eta_pt->Fill(muon3.Pt());
            Pt_muon4_eta_pt->Fill(muon4.Pt());
            n1 = n1 + 1;
            jPsi = muon3 + muon4;
            jPsi_Mass->Fill(jPsi.M());
            jPsi_Pt->Fill(jPsi.Pt());
            if(!(((abs(jPsi.M()) - 3.0969)) < 0.5))continue;
            Pt_muon3_mass_ref->Fill(muon3.Pt());
            Pt_muon4_mass_ref->Fill(muon4.Pt());
            n2 = n2 + 1;
            jPsi_Mass_cut->Fill(jPsi.M());
            jPsi_Pt_cut->Fill(jPsi.Pt());
            if(!((jPsi.Pt()) > 8.5))continue;
            Pt_muon3_jpsi_pt->Fill(muon3.Pt());
            Pt_muon4_jpsi_pt->Fill(muon4.Pt());
            n3 = n3 + 1;
            jPsi_Mass_cut_Pt->Fill(jPsi.M());
            jPsi_Pt_cut_Pt->Fill(jPsi.Pt());
            if(!(jPsi.M() > 2.6 && jPsi.M() < 3.6))continue;
            n4 = n4 + 1;
            Pt_muon3_window->Fill(muon3.Pt());
            Pt_muon4_window->Fill(muon4.Pt());
            jPsi_Mass_window->Fill(jPsi.M(),1/w);
            jPsi_Pt_window->Fill(jPsi.Pt());
            countJPsi = countJPsi +1;
            indice1 = i;
            indice2 = k;
//////
        teste_Chi2->Fill(muChi2NDF->at(i));
        teste_Chi2_pos->Fill(muchi2LocalPosition->at(i));
//////
            nGoodMuons++;
          }//charge
        }//for k
      }//for i
      jPsi_Event->Fill(countJPsi);
      std::cout<< "--- n1 = " << n1  /w << " --- n2 = " << n2 /w << " --- n3 = " << n3 /w  << " --- n4 = " << n4 /w  << endl;
/////////////////////
//      (deltaR(muEta->at(i),muPhi->at(i),muEta->at(k),muPhi->at(k)))<0.3
//muonsCandCollection[dimuonCandCollection[iDimuon].first].Eta(), muonsCandCollection[dimuonCandCollection[iDimuon].first].Phi(), muonsCandCollection[dimuonCandCollection[iDimuon].second].Eta(), muonsCandCollection[dimuonCandCollection[iDimuon].second].Phi()) > 0.3
////////////////////
      if(!(nGoodMuons > 0 ))continue;      
      if(!(count3 >= 2))continue;
      cout << " \t\t nMu = "  << nMu ; 
      for(Int_t i = 0; i < nMu; i++){
        cout << " \n\t i = "  << i ; 
        for(Int_t k = i+1; k < nMu; k++){
          if (i==indice1||i==indice2||k==indice1||k==indice2) continue; 
          cout << "\t k = "  << k ; 
          if(!(((muIDbit->at(i) >> 2) & 1) == 1) && (((muIDbit->at(k) >> 2) & 1) == 1))continue;//Tight muons
          if(muPt->at(i) >= muPt->at(k)){
            muon1.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
            muon2.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));            
          }else{
            muon2.SetPtEtaPhiE(muPt->at(i),muEta->at(i),muPhi->at(i),muEn->at(i));
            muon1.SetPtEtaPhiE(muPt->at(k),muEta->at(k),muPhi->at(k),muEn->at(k));
          }
          dimuon = muon1 + muon2;
          Z = jPsi + dimuon;
          Z_mass_test->Fill(Z.M());
          Z_Pt_test->Fill(Z.Pt());
          dimuon_mass_test->Fill(dimuon.M());
          dimuon_Pt_test->Fill(dimuon.Pt());
          Pt_muon1_tight->Fill(muon1.Pt());
          Pt_muon2_tight->Fill(muon2.Pt());
          if(!(muCharge->at(i) != muCharge->at(k)))continue;
          m1 = m1 + 1;
          Z_mass_mc->Fill(Z.M());
          Z_Pt_mc->Fill(Z.Pt());
          dimuon_mass_mc->Fill(dimuon.M());
          dimuon_Pt_mc->Fill(dimuon.Pt());
          Pt_muon1_charge->Fill(muon1.Pt());
          Pt_muon2_charge->Fill(muon2.Pt());
          if(!( ((abs(muon1.Eta())) < 2.4) && ((abs(muon2.Eta())) < 2.4) ))continue;
          m2 = m2 + 1;
          Z_mass_mc_Eta->Fill(Z.M());
          Z_Pt_mc_Eta->Fill(Z.Pt());
          dimuon_mass_mc_Eta->Fill(dimuon.M());
          dimuon_Pt_mc_Eta->Fill(dimuon.Pt());
            Pt_muon1_eta->Fill(muon1.Pt());
            Pt_muon2_eta->Fill(muon2.Pt());
          if(!(((muon1.Pt() > 30) && (muon2.Pt() > 15))||((muon2.Pt() > 30) && (muon1.Pt() > 15))))continue;
          m3 = m3 + 1;
          std::cout<<"\n-----------------Os valores de pt são: " << muon1.Pt()<< " e "  <<  muon2.Pt() <<endl;
          std::cout<<"\n--------- " << muon3.Pt()<<" ----------- " << muon4.Pt() << " A massa do dimuon e "<< dimuon.M() <<endl ;
          Z_mass_mc_Eta_Pt->Fill(Z.M());
          Z_Pt_mc_Eta_Pt->Fill(Z.Pt());
          dimuon_mass_mc_Eta_Pt->Fill(dimuon.M());
          dimuon_Pt_mc_Eta_Pt->Fill(dimuon.Pt());
            Pt_muon1_pt->Fill(muon1.Pt());
            Pt_muon2_pt->Fill(muon2.Pt());
          if(!(dimuon.M() < 80))continue;
          m4 = m4 + 1;
          Z_mass_cut_dimuon->Fill(Z.M());
          Z_Pt_cut_dimuon->Fill(Z.Pt());
          dimuon_mass_cut_dimuon->Fill(dimuon.M());
          dimuon_Pt_cut_dimuon->Fill(dimuon.Pt());
            Pt_muon1_small->Fill(muon1.Pt());
            Pt_muon2_small->Fill(muon2.Pt());
          if(!((abs(Z.M()-91.1876)) < 25))continue;
          m5 = m5 + 1;
            Pt_muon1_mass_ref->Fill(muon1.Pt());
            Pt_muon2_mass_ref->Fill(muon2.Pt());
          if(!(Z.M()>66 && Z.M()<116))continue;
            m6 = m6 + 1;
            countZ = countZ + 1;
            Pt_muon3_true->Fill(muon3.Pt());////
            Pt_muon4_true->Fill(muon4.Pt());////
            En_muon3->Fill(muon3.E());
            En_muon4->Fill(muon4.E());
            En_muon1->Fill(muon1.E());
            En_muon2->Fill(muon2.E());
          if((muon3.E() < muon1.E())&(muon3.E() < muon2.E())&(muon4.E() < muon1.E())&(muon4.E() < muon2.E())){
            En_muon3_v->Fill(muon3.E());
            En_muon4_v->Fill(muon4.E());
            En_muon1_v->Fill(muon1.E());
            En_muon2_v->Fill(muon2.E());
            countZ_2 = countZ_2 + 1;            
          }
//            delta_phi1->Fill(deltaPhi(Wj.Phi(),Wlep.Phi()),weigth); //modelo para delta phi
            delta_phi1->Fill( Z.M(), deltaPhi( muon1.Phi(), jPsi.Phi()));
            delta_phi2->Fill( Z.M(), deltaPhi( muon2.Phi(), jPsi.Phi()));
            Pt_muon1_true->Fill(muon1.Pt());
            Pt_muon2_true->Fill(muon2.Pt());
          Z_mass_final->Fill(Z.M());
          JPsi_mass_final->Fill(jPsi.M(),1/w);
          Z_Pt_final->Fill(Z.Pt());
          dimuon_mass_final->Fill(dimuon.M());
          dimuon_Pt_final->Fill(dimuon.Pt());
        }//for k
        cout << " " << endl;
      } //for i
      std::cout<< "--- m1 = " << m1 / w << " --- m2 = " << m2 /w << " --- m3 = " << m3  /w << " --- m4 = " << m4  /w  << " --- m5 = " << m5 /w  << " --- m6 = " << m6  /w << endl;
//     if(countZ > 0 & countZ_2 >0){
      Z_Event->Fill(countZ);
      Z_Event_2->Fill(countZ_2);
//     }
//////////////////////////////////////////////////////////
   }
//  Pt_muon4->SetLineColor(kRed);
//  Pt_muon4->Draw();
//  Pt_muon3->SetLineColor(kBlue);
//  Pt_muon3->Draw("same");
//  Pt_muon1->SetLineColor(kGreen);
//  Pt_muon1->Draw("same");
//  Pt_muon2->SetLineColor(kOrange);
//  Pt_muon2->Draw("same");

//  auto legend = new TLegend(0.1,0.7,0.48,0.9);
//  legend->AddEntry(Pt_muon1,"Pt_muon1","l");
//  legend->AddEntry(Pt_muon2,"Pt_muon2","l");
//  legend->AddEntry(Pt_muon3,"Pt_muon3","l");
//  legend->AddEntry(Pt_muon4,"Pt_muon4","l");
//  legend->Draw();
//   histo->Draw("colz");
jPsi_Mass_window->SetOption("HIST");
Z_mass_final->SetOption("HIST");
JPsi_mass_final->SetOption("HIST");

//gStyle->SetPalette(55);
}
