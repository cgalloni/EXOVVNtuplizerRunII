#include "../interface/JpsiMuNtuplizer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
#include <cmath>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include <vector>
#include <TH1.h>

//===================================================================================================================
JpsiMuNtuplizer::JpsiMuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
				  edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				  edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
                                  edm::EDGetTokenT<edm::TriggerResults> triggertoken,
                                  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
				  NtupleBranches* nBranches,
                                  std::map< std::string, bool >&           runFlags    )

  : CandidateNtuplizer ( nBranches )
  , muonToken_	        ( muonToken )
  , verticeToken_          ( verticeToken )
  , packedpfcandidatesToken_(packedpfcandidatesToken) 
  , HLTtriggersToken_   ( triggertoken )
  , triggerObjects_     ( triggerobject )
  , isCutFlow_ ( runFlags["doCutFlow"]  )
  , isTriggerDecisions_ ( runFlags["doTriggerDecisions"] )
  , isTriggerObjects_ ( runFlags["doTriggerObjects"] )

   
{
}

//===================================================================================================================
JpsiMuNtuplizer::~JpsiMuNtuplizer( void )
{

}
std::tuple<Float_t, Float_t> JpsiMuNtuplizer::vertexProb( const std::vector<reco::TransientTrack>& tracks) {
  
  Float_t vprob = -1;
  Float_t vz = -1;
  //  Float_t xy = -1;

  KalmanVertexFitter kalman_fitter;
  TransientVertex vertex;
  vertex = kalman_fitter.vertex(tracks);

  if(vertex.isValid()){

    vprob =  TMath::Prob(vertex.totalChiSquared(), vertex.degreesOfFreedom());

    //    if(vprob==0) vprob = -1;
    vz = vertex.position().z();
  }  

  
  return std::forward_as_tuple(vprob, vz);
  //  return vprob;
}


//==================================================================================================================
/*
struct track_entry_t {
  track_entry_t(int ix, int pid, bool mfit) : trackIx(ix), particleID(pid), massFit(mfit) {}
  int trackIx;
  int particleID;
  bool massFit;
};

void JpsiMuNtuplizer::getAllTracks(std::vector<track_entry_t> *out_vector, int onlyThisVertex) {
  HFDecayTreeTrackIterator trackIt;
  //HFDecayTreeIterator treeIt;
  HFDecayTreeTrackIterator treeIt;

  for (trackIt = fTrackIndices.begin(); trackIt!=fTrackIndices.end(); ++trackIt)
    out_vector->push_back(*trackIt);
}
  
  for (treeIt = fSubVertices.begin(); treeIt!=fSubVertices.end(); ++treeIt) {
    if (!treeIt->fVertexing || !onlyThisVertex)
      treeIt->getAllTracks(out_vector,onlyThisVertex);
  }
  
  }
std::vector<track_entry_t> JpsiMuNtuplizer::getAllTracks(int onlyThisVertex) {
  vector<track_entry_t> tracks;
  getAllTracks(&tracks,onlyThisVertex);
  return tracks;
  }*/
//===================================================================================================================
void JpsiMuNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  //std::cout<<"!!!--->Jpsi Ntuplizer<---!!!"<<std::endl;
  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(muonToken_	, muons_    );
if (isTriggerObjects_){
  event.getByToken(triggerObjects_, triggerObjects);
   }
  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);  
  Float_t vprobJpsi=-1; 
  Float_t vzJpsi=-999; 
  Float_t vprobB=-1;
  Float_t vzB=-999; 
  std::vector<reco::TransientTrack> transient_tracks_muon_pair;
  std::vector<reco::TransientTrack> transient_tracks_trimuon;
 
  std::vector<pat::PackedCandidate> pfcollection; 
  std::vector<Int_t> pfidcollection; 
  pfcollection.clear();
  pfidcollection.clear();
  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    if(TMath::Abs(pf.pdgId())==211 && TMath::Abs(pf.eta()) < 2.3){
      pfcollection.push_back(pf);
      //std::cout<<pf.pt()<<std::endl;
      pfidcollection.push_back(ii);
    }
  }
  Int_t numOfch = (int)pfcollection.size();

  std::vector<pat::Muon> muoncollection;
  std::vector<unsigned int> muoncollection_id;
  muoncollection.clear();
  muoncollection_id.clear();

  int wherepass;
  if ( isCutFlow_ ) {
     wherepass = 0;
     }
  unsigned int nmus = 0;
  unsigned int isJpsiMu_=0;

event.getByToken(HLTtriggersToken_, HLTtriggers_);

// std::cout << "Initialized. Checking Event Trigger: " << std::endl;
// Did the event trigger?
//

  bool isTriggered = false;
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
  std::string finalTriggerName = "";

if (isTriggerDecisions_) {
  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos){
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	finalTriggerName = trigNames.triggerName(i);
      }
    }
  }
  // std::cout << "finalTriggerName = "  << finalTriggerName << std::endl;
  }

// Look for muons that could have come from a J/psi. First we check pT and eta,
// then we check if it's soft 
  reco::VertexCollection::const_iterator firstGoodVertex = vertices_->begin();

  for(unsigned int imuon = 0; imuon < muons_->size(); ++imuon){
    const pat::Muon & mu =(*muons_)[imuon];
    if(mu.pt()<4 ||fabs(mu.eta())>2.4) {
      continue;
         }
    if ( isCutFlow_ ) {
        wherepass = 1; 
       }
    //if( mu.isMuon()==0 || mu.isGlobalMuon()==0 || mu.isTrackerMuon()==0 ) continue;
    //Muon softness
    Int_t isSoft=mu.isSoftMuon(*firstGoodVertex);
    if(isSoft<0.5) {
      continue;
      }
    if ( isCutFlow_ ) {
       wherepass = 2; 
       }
     muoncollection.push_back(mu);
     muoncollection_id.push_back(imuon); 
     nmus++;
     }
    
//  now we see if the muons make a reasonable Jpsi

    for (unsigned int jmu1=0; jmu1<nmus; jmu1++) {
        for (unsigned int jmu2=jmu1+1; jmu2<nmus; jmu2++) {
//we don't care if there aren't enough muons to make a Jpsi, really this just speeds things up a bit
            if(!( nmus >= 2)) {continue;}
            const pat::Muon & mu1=muoncollection[jmu1];
            const pat::Muon & mu2=muoncollection[jmu2];

            if (mu1.charge() + mu2.charge() != 0) continue;

            TLorentzVector tlv_mu1, tlv_mu2;
            tlv_mu1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), mu1.mass());
            tlv_mu2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), mu2.mass());

            TLorentzVector tlv_jpsi = (tlv_mu1+tlv_mu2);
            Float_t jpsi_mass = tlv_jpsi.M();
            Float_t jpsi_pt = tlv_jpsi.Pt();
         
            if(jpsi_mass < 2.6 || jpsi_mass > 3.6) continue;

            if ( isCutFlow_ ) {
               wherepass = 3; 
               }
            // At this point it's passed our local selection cut. Now we compare to trigger.
            // Assuming we're doing trigger matching.
            // Aka, we use TriggerDecision and then compare the muons to TriggerObjects

            if (isTriggerDecisions_) {
               if( !isTriggered ) continue; 
               if ( isCutFlow_ ) {
                  wherepass = 4; 
                  }
               }
            // Now we match to trigger objects
            if (isTriggerObjects_) {

            bool trigMatch1 = false;
            bool trigMatch2 = false;
 
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    if (trigMatch1 && trigMatch2) continue;
                    // ^^^ we don't need to keep looking if we've already found matches!
                    obj.unpackPathNames(trigNames);
                    obj.unpackFilterLabels(event, *HLTtriggers_);
 
                    std::vector<std::string> pathNamesAll  = obj.pathNames(false);
 
                    bool isPathExist = false;
 
                    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
                        if(pathNamesAll[h]==finalTriggerName) isPathExist = true;
                        }
                    
                    if(!isPathExist) continue;
 
                    bool isFilterExist = false;
                
                    for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
                     
                        if(obj.filterLabels()[hh].find("hltJpsiTkVertexFilter") != std::string::npos){
                          isFilterExist = true;
                          }
                        }
                    
                    if(!isFilterExist) continue;
                    
                    Float_t trigger_dR1 = reco::deltaR(obj.eta(), obj.phi(),
                     				mu1.eta(), mu1.phi());
                    Float_t trigger_dR2 = reco::deltaR(obj.eta(), obj.phi(),
                     				mu2.eta(), mu2.phi());
                    
                    if(!trigMatch1 && trigger_dR1 < 0.1) {trigMatch1 = true;}
                    else if (!trigMatch2 && trigger_dR2 < 0.1) {trigMatch2 = true;}
                    }
               if (!trigMatch1 || !trigMatch2 ) continue;
               if ( isCutFlow_ ) {
                  wherepass = 5; 
                  }
               }
/// Ok now we are done with trigger stuff

             //Vertex requirement
             transient_tracks_muon_pair.clear();
             const reco::TrackRef    track1 = mu1.muonBestTrack(); 
             const reco::TrackRef    track2 = mu2.muonBestTrack();
             transient_tracks_muon_pair.push_back((*builder).build(track1.get()));
             transient_tracks_muon_pair.push_back((*builder).build(track2.get()));
             KalmanVertexFitter kalman_fitter;
             TransientVertex Jpsi_vertex;
             Jpsi_vertex = kalman_fitter.vertex(transient_tracks_muon_pair);
             std::tie(vprobJpsi, vzJpsi) = vertexProb(transient_tracks_muon_pair);
 
 
             //Flying part
             double Jpsi_flightSig3D = -999;
             double Jpsi_flightLength3D =  -999; 
             double Jpsi_flightLengthErr3D = -999; 
             
             double Jpsi_flightSig2D =  -999; 
             double Jpsi_flightLength2D = -999; 
             double Jpsi_flightLengthErr2D = -999; 
              
             GlobalVector Jpsi_direction(tlv_jpsi.Px(), tlv_jpsi.Py(), tlv_jpsi.Pz()); //To compute sign of IP
             Jpsi_flightSig3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex, Jpsi_vertex, Jpsi_direction, true).significance();
             Jpsi_flightSig3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,Jpsi_vertex,Jpsi_direction,true).significance();
             Jpsi_flightLength3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,Jpsi_vertex,Jpsi_direction,true).value();
             Jpsi_flightLengthErr3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,Jpsi_vertex,Jpsi_direction,true).error();
             
             Jpsi_flightSig2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,Jpsi_vertex,Jpsi_direction,true).significance();
             Jpsi_flightLength2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,Jpsi_vertex,Jpsi_direction,true).value();
             Jpsi_flightLengthErr2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,Jpsi_vertex,Jpsi_direction,true).error();

       //std::cout<<"-------->VPROB:"<<vprobJpsi<<std::endl;
             unsigned int nmus3=0;
             for(unsigned int kmuon = 0; kmuon < muons_->size(); ++kmuon){
                const pat::Muon mu3 = (*muons_)[kmuon];
           	nmus3++;
           	if(mu3.pt()<2) {
                     continue;
                     }
           	if(nmus3==muoncollection_id[jmu1] || nmus3==muoncollection_id[jmu2]) continue;
                  // we don't want it to be the same muon as mu 1 or 2!
               if ( isCutFlow_ ) {
                  wherepass = 6;
                  }
           	//std::cout<<"1:"<<nmus<<"\t, 2:"<<nmus2<<"\t ,3:"<<nmus3<<std::endl;
           	//std::cout<<"1:"<<mu.pt()<<"\t, 2:"<<mu2.pt()<<"\t ,3:"<<mu3.pt()<<std::endl;
           	TLorentzVector tlv_muon3;
           	TransientVertex B_vertex;
           	tlv_muon3.SetPtEtaPhiM(mu3.pt(), mu3.eta(), mu3.phi(), mu3.mass());
           	TLorentzVector tlv_B = tlv_jpsi + tlv_muon3;
           	double b_mass =tlv_B.M();
           	//std::cout<<"trimuon scenario, bmass: "<< b_mass<<std::endl;
           	//vtx prob part
           	transient_tracks_trimuon.clear();
           	const reco::TrackRef    track3 = mu3.muonBestTrack(); 
           	transient_tracks_trimuon.push_back((*builder).build(track1.get()));
           	transient_tracks_trimuon.push_back((*builder).build(track2.get()));
           	transient_tracks_trimuon.push_back((*builder).build(track3.get()));
           	KalmanVertexFitter B_kalman_fitter;
           	B_vertex = B_kalman_fitter.vertex(transient_tracks_trimuon);
           	//std::cout<<"before getting to the probability thing: "<<std::endl;
           
           	std::tie(vprobB, vzB) = vertexProb(transient_tracks_trimuon);
           	//std::cout<<"after getting to the probability thing: "<<std::endl;
           
           	GlobalVector direction(tlv_B.Px(), tlv_B.Py(), tlv_B.Pz()); //To compute sign of IP
           
           	double flightSig3D = -999;
           	double flightLength3D =  -999; 
           	double flightLengthErr3D = -999; 
           	
           	double flightSig2D =  -999; 
           	double flightLength2D = -999; 
           	double flightLengthErr2D = -999; 
           	flightSig3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex, B_vertex, direction, true).significance();
           	flightSig3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,B_vertex,direction,true).significance();
           	flightLength3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,B_vertex,direction,true).value();
           	flightLengthErr3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,B_vertex,direction,true).error();
           	    
           	flightSig2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,B_vertex,direction,true).significance();
           	flightLength2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,B_vertex,direction,true).value();
           	flightLengthErr2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,B_vertex,direction,true).error();
           	TLorentzVector tlv_mu3;
           	tlv_muon3.SetPtEtaPhiM(mu3.pt(), mu3.eta(), mu3.phi(), mu3.mass());
           
           	Float_t min_dR_pf = 999.;
           	Float_t iso_pt03 =0.;
           	Float_t iso_pt04 =0.;
           	Float_t iso_pt05 =0.;
           	Float_t iso_pt06 =0.;
           	Float_t iso_pt07 =0.;
           
           	//isolation part
           	for(int kkk = 0; kkk < numOfch; kkk ++){
           	  
           	  pat::PackedCandidate _pf = pfcollection[kkk];
           	  if(TMath::Abs(_pf.pdgId())!=211) continue;
           	  TLorentzVector tlv_iso;
           	  tlv_iso.SetPtEtaPhiM(_pf.pt(), _pf.eta(), _pf.phi(), _pf.mass());
           	  Float_t dR_iso = tlv_muon3.DeltaR(tlv_iso);
           	  if(min_dR_pf>dR_iso) min_dR_pf=dR_iso;
           //	  if(dR_iso>0.7) continue; 
           	  if(dR_iso<0.7)  iso_pt07 += tlv_iso.Pt();
           	  if(dR_iso<0.6)  iso_pt06 += tlv_iso.Pt();
           	  if(dR_iso<0.5)  iso_pt05 += tlv_iso.Pt();
           	  if(dR_iso<0.4)  iso_pt04 += tlv_iso.Pt();
           	  if(dR_iso<0.3)  iso_pt03 += tlv_iso.Pt();
           	  //std::cout<<dR_iso<<"<-dr, pt03->"<<iso_pt07<<std::endl;
           	}
          	//isolation part
          	nBranches_->Jpsi_mu3_isopt03.push_back(iso_pt03);
          	nBranches_->Jpsi_mu3_isopt04.push_back(iso_pt04);
          	nBranches_->Jpsi_mu3_isopt05.push_back(iso_pt05);
          	nBranches_->Jpsi_mu3_isopt06.push_back(iso_pt06);
          	nBranches_->Jpsi_mu3_isopt07.push_back(iso_pt07);
          	nBranches_->Jpsi_dr_mu3pf.push_back(min_dR_pf);  
          	//flight part
          	//jpsi part
          	nBranches_->Jpsi_flightSig3D.push_back(Jpsi_flightSig3D);
          	nBranches_->Jpsi_flightLength3D.push_back(Jpsi_flightLength3D);
          	nBranches_->Jpsi_flightLengthErr3D.push_back(Jpsi_flightLengthErr3D);
          	nBranches_->Jpsi_flightSig2D.push_back(Jpsi_flightSig2D);
          	nBranches_->Jpsi_flightLength2D.push_back(Jpsi_flightLength2D);
          	nBranches_->Jpsi_flightLengthErr2D.push_back(Jpsi_flightLengthErr2D);
          
          	//trimuon part
          	nBranches_->Jpsi_trimu_flightSig3D.push_back(flightSig3D);
          	nBranches_->Jpsi_trimu_flightLength3D.push_back(flightLength3D);
          	nBranches_->Jpsi_trimu_flightLengthErr3D.push_back(flightLengthErr3D);
          	nBranches_->Jpsi_trimu_flightSig2D.push_back(flightSig2D);
          	nBranches_->Jpsi_trimu_flightLength2D.push_back(flightLength2D);
          	nBranches_->Jpsi_trimu_flightLengthErr2D.push_back(flightLengthErr2D);
          
          	isJpsiMu_=1;
          	//save muon branches
          	nBranches_->Jpsi_mu1_isLoose.push_back(mu1.isLooseMuon());
          	nBranches_->Jpsi_mu1_isTight.push_back(mu1.isTightMuon(*firstGoodVertex));
          	nBranches_->Jpsi_mu1_isPF.push_back(mu1.isPFMuon());
          	nBranches_->Jpsi_mu1_isGlobal.push_back(mu1.isGlobalMuon());
          	nBranches_->Jpsi_mu1_isTracker.push_back(mu1.isTrackerMuon());
          	nBranches_->Jpsi_mu1_isSoft.push_back(mu1.isSoftMuon(*firstGoodVertex));
          
          	nBranches_->Jpsi_mu2_isLoose.push_back(mu2.isLooseMuon());
          	nBranches_->Jpsi_mu2_isTight.push_back(mu2.isTightMuon(*firstGoodVertex));
          	nBranches_->Jpsi_mu2_isPF.push_back(mu2.isPFMuon());
          	nBranches_->Jpsi_mu2_isGlobal.push_back(mu2.isGlobalMuon());
          	nBranches_->Jpsi_mu2_isTracker.push_back(mu2.isTrackerMuon());
          	nBranches_->Jpsi_mu2_isSoft.push_back(mu2.isSoftMuon(*firstGoodVertex));
          
          
          	nBranches_->Jpsi_mu1_pt.push_back(mu1.pt());
          	nBranches_->Jpsi_mu1_eta.push_back(mu1.eta());
          	nBranches_->Jpsi_mu1_phi.push_back(mu1.phi());
          	nBranches_->Jpsi_mu1_ch.push_back(mu1.charge());
          
          	nBranches_->Jpsi_mu2_pt.push_back(mu2.pt());
          	nBranches_->Jpsi_mu2_eta.push_back(mu2.eta());
          	nBranches_->Jpsi_mu2_phi.push_back(mu2.phi());
          	nBranches_->Jpsi_mu2_ch.push_back(mu2.charge());
          
          	nBranches_->Jpsi_mu3_pt.push_back(mu3.pt());
          	nBranches_->Jpsi_mu3_eta.push_back(mu3.eta());
          	nBranches_->Jpsi_mu3_phi.push_back(mu3.phi());
          	nBranches_->Jpsi_mu3_ch.push_back(mu3.charge());
          	nBranches_->Jpsi_mu3_x.push_back(mu3.vx());
          	nBranches_->Jpsi_mu3_y.push_back(mu3.vy());
          	nBranches_->Jpsi_mu3_z.push_back(mu3.vz());
          
          
          
          	nBranches_->Jpsi_mu3_isLoose.push_back(mu3.isLooseMuon());
          	nBranches_->Jpsi_mu3_isTight.push_back(mu3.isTightMuon(*firstGoodVertex));
          	nBranches_->Jpsi_mu3_isPF.push_back(mu3.isPFMuon());
          	nBranches_->Jpsi_mu3_isGlobal.push_back(mu3.isGlobalMuon());
          	nBranches_->Jpsi_mu3_isTracker.push_back(mu3.isTrackerMuon());
          	nBranches_->Jpsi_mu3_isSoft.push_back(mu3.isSoftMuon(*firstGoodVertex));
          	//std::cout<<firstGoodVertex->position().z()<<"<-PV, Jpsi->"<< Jpsi_vertex.position().z()<< std::endl;
          	//jpsi variables
          	nBranches_->Jpsi_PV_x.push_back(firstGoodVertex->position().x());
          	nBranches_->Jpsi_PV_y.push_back(firstGoodVertex->position().y());
          	nBranches_->Jpsi_PV_z.push_back(firstGoodVertex->position().z());
          	if(Jpsi_vertex.isValid()){
          	  nBranches_->Jpsi_dx.push_back(Jpsi_vertex.position().x()-firstGoodVertex->position().x());
          	  nBranches_->Jpsi_dy.push_back(Jpsi_vertex.position().y()-firstGoodVertex->position().y());
          	  nBranches_->Jpsi_dz.push_back(Jpsi_vertex.position().z()-firstGoodVertex->position().z());	  
          	}else{
          	  nBranches_->Jpsi_dx.push_back(-999);
          	  nBranches_->Jpsi_dy.push_back(-999);
          	  nBranches_->Jpsi_dz.push_back(-999);
          	}
          	nBranches_->Jpsi_pt.push_back(tlv_jpsi.Pt());
          	nBranches_->Jpsi_eta.push_back(tlv_jpsi.Eta());
          	nBranches_->Jpsi_phi.push_back(tlv_jpsi.Phi());
          	nBranches_->Jpsi_mass.push_back(jpsi_mass);
          	nBranches_->Jpsi_vtxprob.push_back(vprobJpsi);
          	nBranches_->Jpsi_vtxz.push_back(vzJpsi);
          	//std::cout<<"this is before the trimuon part"<< std::endl;
          
                  //Trimuon (B) variables
          	nBranches_->Jpsi_trimu_pt.push_back(tlv_B.Pt());
          	nBranches_->Jpsi_trimu_eta.push_back(tlv_B.Eta());
          	nBranches_->Jpsi_trimu_phi.push_back(tlv_B.Phi());
          	nBranches_->Jpsi_trimu_mass.push_back(b_mass);
          
          	//std::cout<<"this is before the trimuon vtx"<< std::endl;
          
          	nBranches_->Jpsi_trimu_vtxprob.push_back(vprobB);
          	//std::cout<<"this is before the trimuon vtx z direction"<< std::endl;
          	//std::cout<<vprobB<<"<-yutas, normal->"<<vprobB2<< std::endl;
          	//std::cout<<vzB<<"<-yutas, normal->"<<B_vertex.position().z()<< std::endl;
          
          	//std::cout<<sizeof(B_vertex.position())<<std::endl;
          	//std::cout<<B_vertex.position()<<std::endl;
          	//std::cout<<B_vertex.position().z()<<std::endl;
          	nBranches_->Jpsi_trimu_vtxz.push_back(vzB);
          	//std::cout<<"this is before the trimuon position"<< std::endl;
          	if(B_vertex.isValid()){
          	  nBranches_->Jpsi_trimu_dx.push_back(B_vertex.position().x()-firstGoodVertex->position().x());
          	  nBranches_->Jpsi_trimu_dy.push_back(B_vertex.position().y()-firstGoodVertex->position().y());
          	  nBranches_->Jpsi_trimu_dz.push_back(B_vertex.position().z()-firstGoodVertex->position().z());
          	}else{
          	   nBranches_->Jpsi_trimu_dx.push_back(-999);
          	   nBranches_->Jpsi_trimu_dy.push_back(-999);
          	   nBranches_->Jpsi_trimu_dz.push_back(-999);
          	}
          	//std::cout<<"this should be the end of the Jpsi part"<< std::endl;
      }
            }
        }
  //
  nBranches_->IsJpsiMu.push_back(isJpsiMu_);
  // within the analyze()
  

  if ( isCutFlow_ ) {
     for (int i = 0; i<=wherepass; i++) {
         nBranches_->cutflow_perevt->Fill(i);
         }
     }
//    nBranches_->cutflow_perevt->Fill(wherepass);


  if(isJpsiMu_==0){
    //std::cout<<"this should be working and i am returning"<<std::endl;
    return;
    }
}
