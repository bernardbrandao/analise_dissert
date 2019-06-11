#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TMath.h"
#include "math.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include <boost/foreach.hpp>
#include <string>
////////////////////
#include <iostream>
#include "ggAnalysis/ggNtuplizer/interface/deltaR.h"
#include "ggAnalysis/ggNtuplizer/interface/deltaPhi.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
////////////////////
using namespace std;

// (local) variables associated with tree branches
Int_t            ndiMu_;
vector<int>      diMuIndex1_;
vector<int>      diMuIndex2_;
vector<int>      diMuVtxIsValid_;
vector<float>    diMuVtx_;
vector<float>    diMuVty_;
vector<float>    diMuVtz_;
vector<float>    diMuChi2_;
vector<int>      diMuNDF_;
vector<float>    diMuVtxProb_;
vector<float>    diMu_CosAlpha_; //the angle between the reconstructed momentum vector of the dimuon system and the vector from the PV to the dimuon vertex
vector<float>    diMu_Lxy_; // Transverse decay length (Lxy) between the dimuon vertex and the primary vertex
vector<float>    diMu_Rxy_;
vector<float>    diMu_eLxy_;
vector<float>    diMu_SLxy_; // Significance of Lxy ( Lxy divided by its uncertainty )
vector<float>    diMu_ctau_; // lifetime 
vector<float>    diMu_ctauErr_;
////////////////////////////
vector<int>	 global_;
vector<int>      soft_;
vector<float>    dxyMu1_;
vector<float>    dxyMu2_;
vector<float>    dzMu1_;
vector<float>    dzMu2_;
vector<float>	 SIP1_;
vector<float>	 SIP2_;	
vector<int>  	 RIso_;	
vector<float>	 deltaR_;
vector<int>      charge_;
vector<float>    ptMu1_;
vector<float>    ptMu2_;
vector<float>    etaMu1_;
vector<float>    etaMu2_;
vector<float>    phiMu1_;
vector<float>    phiMu2_;
///////////////////////////
TLorentzVector Mu1,Mu2, jpsiCand;
vector<float>	 ptJPsi_;
vector<float>	 massJPsi_;
vector<float>	 vtxprobJPsi_;
//////////////////////////
Int_t            ndiLead_;
vector<int>	 difIndex_;
vector<int>      globalZ_;
vector<int>      tight_;
vector<float>    dxyMu3_;
vector<float>    dxyMu4_;
vector<float>    dzMu3_;
vector<float>    dzMu4_;
vector<float>    SIP3_;
vector<float>    SIP4_;
vector<int>      RIsoZ_;
vector<float>    deltaRZ_;
vector<int>	 chargeZ_;
vector<float>    ptMu3_;
vector<float>    ptMu4_;
vector<float>    etaMu3_;
vector<float>    etaMu4_;
vector<float>    phiMu3_;
vector<float>    phiMu4_;
TLorentzVector Mu3,Mu4, LeadMu, Z;
vector<float>	 massdiMu_;
vector<float>	 massZ_;
vector<float>	 vtxprobZ_;
//////////////////////////

void ggNtuplizer::branchesMuonPairs(TTree* tree) {

  tree->Branch("ndiMu",       &ndiMu_);
  tree->Branch("diMuIndex1",   &diMuIndex1_);
  tree->Branch("diMuIndex2",   &diMuIndex2_);
  tree->Branch("diMuVtxIsValid",   &diMuVtxIsValid_);
  tree->Branch("diMuVtx",   &diMuVtx_);
  tree->Branch("diMuVty",   &diMuVty_);
  tree->Branch("diMuVtz",   &diMuVtz_);
  tree->Branch("diMuChi2",   &diMuChi2_);
  tree->Branch("diMuNDF",   &diMuNDF_);
  tree->Branch("diMuVtxProb",   &diMuVtxProb_);
  tree->Branch("diMu_CosAlpha",   &diMu_CosAlpha_);
  tree->Branch("diMu_Lxy",   &diMu_Lxy_);
  tree->Branch("diMu_Rxy",   &diMu_Rxy_);
  tree->Branch("diMu_eLxy",   &diMu_eLxy_);
  tree->Branch("diMu_SLxy",   &diMu_SLxy_);
  tree->Branch("diMu_ctau",   &diMu_ctau_);
  tree->Branch("diMu_ctauErr",   &diMu_ctauErr_);
  ///////////////////////////////////////
  tree->Branch("global",    &global_);
  tree->Branch("soft",    &soft_);
  tree->Branch("dxyMu1",    &dxyMu1_);
  tree->Branch("dxyMu2",    &dxyMu2_);
  tree->Branch("dzMu1",    &dzMu1_);
  tree->Branch("dzMu2",    &dzMu2_);
  tree->Branch("SIP1", 	   &SIP1_);
  tree->Branch("SIP2",     &SIP2_);
  tree->Branch("RIso",    &RIso_);
  tree->Branch("deltaR",    &deltaR_);
  tree->Branch("charge",    &charge_);
  tree->Branch("ptMu1",    &ptMu1_);
  tree->Branch("ptMu2",    &ptMu2_);
  tree->Branch("etaMu1",    &etaMu1_);
  tree->Branch("etaMu2",    &etaMu2_);
  tree->Branch("phiMu1",    &phiMu1_);
  tree->Branch("phiMu2",    &phiMu2_);
  tree->Branch("ptJPsi",    &ptJPsi_);
  tree->Branch("massJPsi",    &massJPsi_);
  tree->Branch("vtxprobJPsi",    &vtxprobJPsi_);
  ///////////////////////////////////////////////
  tree->Branch("ndiLead",       &ndiLead_);
  tree->Branch("difIndex",  &difIndex_);
  tree->Branch("globalZ",  &globalZ_);
  tree->Branch("tight",  &tight_);
  tree->Branch("dxyMu3",  &dxyMu3_);
  tree->Branch("dxyMu4",  &dxyMu4_);
  tree->Branch("dzMu3",  &dzMu3_);
  tree->Branch("dzMu4",  &dzMu4_);
  tree->Branch("SIP3",  &SIP3_);
  tree->Branch("SIP4",  &SIP4_);
  tree->Branch("RIsoZ",  &RIsoZ_);
  tree->Branch("deltaRZ",  &deltaRZ_);
  tree->Branch("chargeZ",  &chargeZ_);
  tree->Branch("ptMu3",  &ptMu3_);
  tree->Branch("ptMu4",  &ptMu4_);
  tree->Branch("etaMu3",  &etaMu3_);
  tree->Branch("etaMu4",  &etaMu4_);
  tree->Branch("phiMu3",  &phiMu3_);
  tree->Branch("phiMu4",  &phiMu4_);
  tree->Branch("massdiMu",  &massdiMu_);
  tree->Branch("massZ",  &massZ_);
  tree->Branch("vtxprobZ",  &vtxprobZ_);
}

void ggNtuplizer::fillMuonsPairs(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  diMuIndex1_.clear();
  diMuIndex2_.clear();
  diMuVtxIsValid_.clear();
  diMuVtx_.clear();
  diMuVty_.clear();
  diMuVtz_.clear();
  diMuChi2_.clear();
  diMuNDF_.clear();
  diMuVtxProb_.clear();
  diMu_CosAlpha_.clear();
  diMu_Lxy_.clear();
  diMu_Rxy_.clear();
  diMu_eLxy_.clear();
  diMu_SLxy_.clear();
  diMu_ctau_.clear();
  diMu_ctauErr_.clear();
  ////////////////////////////////////////
  global_.clear();
  soft_.clear();
  dxyMu1_.clear();
  dxyMu2_.clear();
  dzMu1_.clear();
  dzMu2_.clear();
  SIP1_.clear();
  SIP2_.clear();
  RIso_.clear();
  deltaR_.clear();
  charge_.clear();
  ptMu1_.clear();
  ptMu2_.clear();
  etaMu1_.clear();
  etaMu2_.clear();
  phiMu1_.clear();
  phiMu2_.clear();
  ptJPsi_.clear();
  massJPsi_.clear();
  vtxprobJPsi_.clear();
  ///////////////////////////////////////
  difIndex_.clear();
  globalZ_.clear();
  tight_.clear();
  dxyMu3_.clear();
  dxyMu4_.clear();
  dzMu3_.clear();
  dzMu4_.clear();
  SIP3_.clear();
  SIP4_.clear();
  RIsoZ_.clear();
  deltaRZ_.clear();
  chargeZ_.clear();
  ptMu3_.clear();
  ptMu4_.clear();
  etaMu3_.clear();
  etaMu4_.clear();
  phiMu3_.clear();
  phiMu4_.clear();
  massdiMu_.clear();
  massZ_.clear();
  vtxprobZ_.clear();


  ndiMu_ = 0;
  ndiMu_ = 0;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  if (!muonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Muons in event";
    return;
  }

  int nGoodMuons = 0;
  int tmpMu1 = 0;
  int tmpMu2 = 0;
  int tmpMu3 = 0;
  int tmpMu4 = 0;

  edm::View<pat::Muon>::const_iterator muindex_1;
  edm::View<pat::Muon>::const_iterator muindex_2;

  edm::Handle<edm::TriggerResults> trgResultsHandle;
  e.getByToken(trgResultsLabel_, trgResultsHandle);

  const edm::TriggerNames &trgNames = e.triggerNames(*trgResultsHandle);
  int bitEleMuX = -1;
  for (size_t i = 0; i < trgNames.size(); ++i) {
    const string &name = trgNames.triggerName(i);
    if ((name.find("HLT_IsoMu24_v")                  != string::npos))          bitEleMuX = 0;
    else if ((name.find("HLT_Mu30_TkMu11_v")         != string::npos))          bitEleMuX = 1;
    else if ((name.find("HLT_TripleMu_12_10_5_v")    != string::npos))          bitEleMuX = 2;
  }

  if ((bitEleMuX == 0 || bitEleMuX == 1 || bitEleMuX == 2)){

    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {

      // Build transientTrack
      const reco::TransientTrack &tt1 = transientTrackBuilder->build(iMu->bestTrack());

      for (edm::View<pat::Muon>::const_iterator jMu = iMu; jMu != muonHandle->end(); ++jMu) {

	if (iMu == jMu) {
	  tmpMu2++;
	  continue;
	}
	if (iMu->pt() < 3.0 || jMu->pt() < 3.0) continue;
	if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
	if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;

	diMuIndex1_.push_back(tmpMu1);
	diMuIndex2_.push_back(tmpMu2);

	const reco::TransientTrack &tt2 = transientTrackBuilder->build(jMu->bestTrack());
	vector<reco::TransientTrack> t_tks = {tt1,tt2};

	KalmanVertexFitter fitter;
	TransientVertex tmpVertex = fitter.vertex(t_tks);

	if(tmpVertex.isValid()){
	  diMuVtxIsValid_.push_back(1);
	  diMuVtx_.push_back(tmpVertex.position().x());
	  diMuVty_.push_back(tmpVertex.position().y());
	  diMuVtz_.push_back(tmpVertex.position().z());
	  diMuChi2_.push_back(tmpVertex.totalChiSquared());
	  diMuNDF_.push_back(tmpVertex.degreesOfFreedom());
	  diMuVtxProb_.push_back(ChiSquaredProbability(tmpVertex.totalChiSquared(), tmpVertex.degreesOfFreedom()));

	  ////Verificando se a minha solucao vai dar certo
	  if(iMu->isGlobalMuon() && jMu->isGlobalMuon()){
	    global_.push_back(1);
	  }
	  else{
	    global_.push_back(0);
	  }
	  if(iMu->isSoftMuon(vtx) && jMu->isSoftMuon(vtx)){
	    soft_.push_back(1);
	  }
	  else{
	    soft_.push_back(0);
	  }
	  dxyMu1_.push_back(iMu->muonBestTrack()->dxy(pv));
	  dxyMu2_.push_back(jMu->muonBestTrack()->dxy(pv));
	  dzMu1_.push_back(iMu->muonBestTrack()->dz(pv));
	  dzMu2_.push_back(jMu->muonBestTrack()->dz(pv));
	  SIP1_.push_back((fabs(iMu->dB(pat::Muon::PV3D))/iMu->edB(pat::Muon::PV3D)));
	  SIP2_.push_back((fabs(jMu->dB(pat::Muon::PV3D))/jMu->edB(pat::Muon::PV3D)));
	  if(((iMu->pfIsolationR04().sumChargedHadronPt + max( 0., iMu->pfIsolationR04().sumNeutralHadronEt + iMu->pfIsolationR04().sumPhotonEt - 0.5 * iMu->pfIsolationR04().sumPUPt ) ) /  (iMu->pt() )< 0.35) && ((jMu->pfIsolationR04().sumChargedHadronPt + max( 0., jMu->pfIsolationR04().sumNeutralHadronEt + jMu->pfIsolationR04().sumPhotonEt - 0.5 * jMu->pfIsolationR04().sumPUPt ) ) /  (jMu->pt() )< 0.35) ){
	    RIso_.push_back(1);
	  }
	  else{
	    RIso_.push_back(0);
	  }
	  deltaR_.push_back( deltaR(iMu->eta(), iMu->phi(), jMu->eta(), jMu->phi()) );
	  if(iMu->charge() != jMu->charge()){
	    charge_.push_back(1);
	  }
	  else{
	    charge_.push_back(0);
	  }

	  if ( iMu->pt() > jMu->pt()){
	    Mu1.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), 0.1057);
	    Mu2.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), 0.1057);
	  }
	  else{
	    Mu2.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), 0.1057);
	    Mu1.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), 0.1057);
	  }
	  ptMu1_.push_back(Mu1.Pt());
	  ptMu2_.push_back(Mu2.Pt());
	  etaMu1_.push_back(Mu1.Eta());
	  etaMu2_.push_back(Mu2.Eta());
	  phiMu1_.push_back(Mu1.Phi());
	  phiMu2_.push_back(Mu2.Phi());
	  jpsiCand = Mu1 + Mu2;

	  ptJPsi_.push_back(jpsiCand.Pt());
	  massJPsi_.push_back(jpsiCand.M());
	  vtxprobJPsi_.push_back(ChiSquaredProbability(tmpVertex.totalChiSquared(), tmpVertex.degreesOfFreedom()));

	}
	else{
	  diMuVtxIsValid_.push_back(0);
	  diMuVtx_.push_back(-999.);
	  diMuVty_.push_back(-999.);
	  diMuVtz_.push_back(-999.);
	  diMuChi2_.push_back(-999.);
	  diMuNDF_.push_back(-999.);
	  diMuVtxProb_.push_back(-999.);
	  diMu_CosAlpha_.push_back(-999.);
	  diMu_Lxy_.push_back(-999.);
	  diMu_Rxy_.push_back(-999.);
	  diMu_eLxy_.push_back(-999.);
	  diMu_SLxy_.push_back(-999.);
	  diMu_ctau_.push_back(-999.);
	  diMu_ctauErr_.push_back(-999.);
	  ///////////////////////////////////
	  dxyMu1_.push_back(-999.);
	  dxyMu2_.push_back(-999.);
	  dzMu1_.push_back(-999.);
	  dzMu2_.push_back(-999.);
	  SIP1_.push_back(-999.);
	  SIP2_.push_back(-999.);
	  ptMu1_.push_back(-999.);
	  ptMu2_.push_back(-999.);
	  etaMu1_.push_back(-999.);
	  etaMu2_.push_back(-999.);
	  phiMu1_.push_back(Mu1.Phi());
	  phiMu2_.push_back(Mu2.Phi());
	  ptJPsi_.push_back(-999.);
	  massJPsi_.push_back(-999.);
	  vtxprobJPsi_.push_back(-999.);
	  deltaR_.push_back(-999.);
	}
	ndiMu_++;
	tmpMu2++;
      }
      tmpMu1++;
      tmpMu2 = tmpMu1;
    }

    ////////////Imformacoes para a formacao do Z

    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {
      const reco::TransientTrack &tta = transientTrackBuilder->build(muindex_1->bestTrack());
      for (edm::View<pat::Muon>::const_iterator jMu = iMu; jMu != muonHandle->end(); ++jMu) {
	const reco::TransientTrack &ttb = transientTrackBuilder->build(muindex_2->bestTrack());
	if (iMu == jMu) {
	  tmpMu4++;
	  continue;
	}
	/*if(muindex_1 != iMu || muindex_2 != jMu || muindex_2 != iMu || muindex_1 != jMu){
	  difIndex_.push_back(1);
	  }
	  else{
	  difIndex_.push_back(0);
	  }*/
	if(muindex_1 != iMu || muindex_2 != jMu || muindex_2 != iMu || muindex_1 != jMu){
	  const reco::TransientTrack &tt3 = transientTrackBuilder->build(iMu->bestTrack());
	  const reco::TransientTrack &tt4 = transientTrackBuilder->build(jMu->bestTrack());
	  vector<reco::TransientTrack> t_tks_Z = {tta,ttb,tt3,tt4};
	  KalmanVertexFitter fitter_Z;
	  TransientVertex tmpVertex_Z = fitter_Z.vertex(t_tks_Z);
	  if(tmpVertex_Z.isValid()){
	    if(iMu->isGlobalMuon() && jMu->isGlobalMuon()){
	      globalZ_.push_back(1);
	    }
	    else{
	      globalZ_.push_back(0);
	    }
	    if(iMu->isTightMuon(vtx) && jMu->isTightMuon(vtx)){
	      tight_.push_back(1);
	    }
	    else{
	      tight_.push_back(0);
	    }
	    dxyMu3_.push_back(iMu->muonBestTrack()->dxy(pv));
	    dxyMu4_.push_back(jMu->muonBestTrack()->dxy(pv));
	    dzMu3_.push_back(iMu->muonBestTrack()->dz(pv));
	    dzMu4_.push_back(jMu->muonBestTrack()->dz(pv));
	    SIP3_.push_back((fabs(iMu->dB(pat::Muon::PV3D))/iMu->edB(pat::Muon::PV3D)));
	    SIP4_.push_back((fabs(jMu->dB(pat::Muon::PV3D))/jMu->edB(pat::Muon::PV3D)));
	    if(((iMu->pfIsolationR04().sumChargedHadronPt + max( 0., iMu->pfIsolationR04().sumNeutralHadronEt + iMu->pfIsolationR04().sumPhotonEt - 0.5 * iMu->pfIsolationR04().sumPUPt ) ) /  (iMu->pt() )< 0.35) && ((jMu->pfIsolationR04().sumChargedHadronPt + max( 0., jMu->pfIsolationR04().sumNeutralHadronEt + jMu->pfIsolationR04().sumPhotonEt - 0.5 * jMu->pfIsolationR04().sumPUPt ) ) /  (jMu->pt() )< 0.35) ){
	      RIsoZ_.push_back(1);
	    }
	    else{
	      RIsoZ_.push_back(0);
	    }
	    deltaRZ_.push_back( deltaR(iMu->eta(), iMu->phi(), jMu->eta(), jMu->phi()) );
	    if(iMu->charge() != jMu->charge()){
	      chargeZ_.push_back(1);
	    }
	    else{
	      chargeZ_.push_back(0);
	    }
	    if ( iMu->pt() > jMu->pt()){
	      Mu3.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), 0.1057);
	      Mu4.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), 0.1057);
	    }
	    else{
	      Mu4.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), 0.1057);
	      Mu3.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), 0.1057);
	    }
	    ptMu3_.push_back(Mu3.Pt());
	    ptMu4_.push_back(Mu4.Pt());
	    etaMu3_.push_back(Mu3.Eta());
	    etaMu4_.push_back(Mu4.Eta());
	    phiMu3_.push_back(Mu3.Phi());
	    phiMu4_.push_back(Mu4.Phi());
	    LeadMu = Mu3 + Mu4;
	    Z = Mu1 + Mu2 + Mu3 + Mu4;

	    massdiMu_.push_back(LeadMu.M());
	    massZ_.push_back(Z.M());
	    vtxprobZ_.push_back(ChiSquaredProbability(tmpVertex_Z.totalChiSquared(), tmpVertex_Z.degreesOfFreedom()));
	  }//isValid
	  else{
	    dxyMu3_.push_back(-999.);
	    dxyMu4_.push_back(-999.);
	    dzMu3_.push_back(-999.);
	    dzMu4_.push_back(-999.);
	    SIP3_.push_back(-999.);
	    SIP4_.push_back(-999.);
	    deltaRZ_.push_back(-999.);
	    ptMu3_.push_back(-999.);
	    ptMu4_.push_back(-999.);
	    etaMu3_.push_back(-999.);
	    etaMu4_.push_back(-999.);
	    phiMu3_.push_back(-999.);
	    phiMu4_.push_back(-999.);
	    massdiMu_.push_back(-999.);
	    massZ_.push_back(-999.);
	    vtxprobZ_.push_back(-999.);
	  }
	}//if mu index
	else{
	  globalZ_.push_back(0);
	  tight_.push_back(0);
	  RIsoZ_.push_back(0);
	  chargeZ_.push_back(0);
	  dxyMu3_.push_back(-900.);
	  dxyMu4_.push_back(-900.);
	  dzMu3_.push_back(-900.);
	  dzMu4_.push_back(-900.);
	  SIP3_.push_back(-900.);
	  SIP4_.push_back(-900.);
	  deltaRZ_.push_back(-900.);
	  ptMu3_.push_back(-900.);
	  ptMu4_.push_back(-900.);
	  etaMu3_.push_back(-900.);
	  etaMu4_.push_back(-900.);
	  phiMu3_.push_back(-900.);
	  phiMu4_.push_back(-900.);
	  massdiMu_.push_back(-900.);
	  massZ_.push_back(-900.);
	  vtxprobZ_.push_back(-900.);
	}
	ndiLead_++;
	tmpMu4++;
      }
      tmpMu3++;
      tmpMu4 = tmpMu3;
    }
  }
}
