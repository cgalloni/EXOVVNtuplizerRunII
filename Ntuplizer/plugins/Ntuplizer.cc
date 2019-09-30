#include "../interface/Ntuplizer.h"
#include "../interface/CandidateNtuplizer.h"
#include "../interface/GenJetsNtuplizer.h"
#include "../interface/METsNtuplizer.h"
#include "../interface/PileUpNtuplizer.h"
#include "../interface/GenEventNtuplizer.h"
#include "../interface/GenParticlesNtuplizer.h"
#include "../interface/TriggersNtuplizer.h"
#include "../interface/VerticesNtuplizer.h"
#include "../interface/JpsiMuNtuplizer.h"
#include "../interface/JpsiEleNtuplizer.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

// #include "DataFormats/METReco/interface/PFMET.h"


///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
        beamToken_                  (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
	vtxToken_             	    (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	rhoToken_             	    (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
	packedpfcandidatesToken_    (consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedpfcandidates"))),
	puinfoToken_          	    (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfo"))),
	geneventToken_        	    (consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),     
	lheEventProductToken_       (consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("externallheProducer"))),     
	genparticleToken_     	    (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
	
	

	muonToken_	      	    (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	//mvaValuesMapToken_          (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
	//mvaCategoriesMapToken_      (consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"))),
	ebRecHitsToken_             (consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRecHits"))),

	tauToken_	      	    (consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
				     //tauBoostedTauToken_	    (consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tausBoostedTau"))),

	metToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
	metpuppiToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_puppi"))),
	metmvaToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_mva"))),
	metSigToken_	      	    (consumes<double>(edm::InputTag("METSignificance","METSignificance"))),
	metCovToken_	      	    (consumes<math::Error<2>::type>(edm::InputTag("METSignificance","METCovariance"))),

	jetForMetCorrToken_   	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsForMetCorr"))),

	triggerToken_	      	    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"))),
	triggerObjects_	      	    (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"))),
	triggerPrescales_     	    (consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerprescales"))),
        noiseFilterToken_     	    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"))),
        HBHENoiseFilterLooseResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterLoose"))),
        HBHENoiseFilterTightResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterTight"))),
	HBHENoiseIsoFilterResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseIsoFilter"))),
	 ecalBadCalibFilterUpdateToken_(consumes< bool >(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_ecalBadCalibReducedMINIAODFilter")))

{


  /*=======================================================================================*/
  edm::Service<TFileService> fs;
  TTree* tree = fs->make<TTree>( "tree", "tree" );

  std::map< std::string, bool > runFlags;
  runFlags["runOnMC"] = iConfig.getParameter<bool>("runOnMC");
  runFlags["doGenParticles"] = iConfig.getParameter<bool>("doGenParticles");
  runFlags["doGenEvent"] = iConfig.getParameter<bool>("doGenEvent");
  runFlags["doPileUp"] = iConfig.getParameter<bool>("doPileUp");
  runFlags["doVertices"] = iConfig.getParameter<bool>("doVertices");
  runFlags["doTriggerDecisions"] = iConfig.getParameter<bool>("doTriggerDecisions");
  runFlags["doTriggerObjects"] = iConfig.getParameter<bool>("doTriggerObjects");
  runFlags["doHltFilters"] = iConfig.getParameter<bool>("doHltFilters");
  runFlags["doMissingEt"] = iConfig.getParameter<bool>("doMissingEt");
  runFlags["doMETSVFIT"] = iConfig.getParameter<bool>("doMETSVFIT");
  runFlags["doMVAMET"] = iConfig.getParameter<bool>("doMVAMET");
  runFlags["doJpsiMu"] = iConfig.getParameter<bool>("doJpsiMu");
  runFlags["doJpsiEle"] = iConfig.getParameter<bool>("doJpsiEle");
  runFlags["doGenHist"] = iConfig.getParameter<bool>("doGenHist");
  runFlags["doCutFlow"] = iConfig.getParameter<bool>("doCutFlow");


  
  electronToken_	      	    =consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
    // eleVetoIdMapToken_    	    =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"));
    // eleLooseIdMapToken_   	    =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"));
    // eleMediumIdMapToken_  	    =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"));
    // eleTightIdMapToken_   	    =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"));
    // eleHLTIdMapToken_  	    =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHLTIdMap"));
    // eleHEEPIdMapToken_    	    =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"));
    // eleMVAMediumIdMapToken_     =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAMediumIdMap"));
    // eleMVATightIdMapToken_      =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATightIdMap"));
 

  std::string jecpath = iConfig.getParameter<std::string>("jecpath");
  jecpath = "EXOVVNtuplizerRunII/Ntuplizer/data/" + jecpath;
  //jecpath = std::string("data/") + jecpath;
  std::cout << "jecpath  "<< jecpath  <<std::endl;
  nBranches_ = new NtupleBranches( runFlags, tree );

std::cout << "Naming Histos" << std::endl;
if ( runFlags["doCutFlow"] ){
   nBranches_->cutflow_perevt->GetXaxis()->SetBinLabel(1,"Pre-Cut");
   nBranches_->cutflow_perevt->GetXaxis()->SetBinLabel(2,"Passing p_{T} & #eta cut for #mu_{1} & #mu_{2}");
   nBranches_->cutflow_perevt->GetXaxis()->SetBinLabel(3,"Passing Soft cut for #mu_{1} & #mu_{2}");
   nBranches_->cutflow_perevt->GetXaxis()->SetBinLabel(4,"#mu_{1} & #mu_{2} Make a believable J/#psi");
   nBranches_->cutflow_perevt->GetXaxis()->SetBinLabel(5,"Event Triggered");
   nBranches_->cutflow_perevt->GetXaxis()->SetBinLabel(6,"Event Passed Trigger matching");
   nBranches_->cutflow_perevt->GetXaxis()->SetBinLabel(7,"Passing p_{T} cut for #mu_{3}");
   }
//
 if ( runFlags["doGenHist"] ){
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(1, "#mu");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(2, "#pi^{0}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(3, "#pi^{#pm}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(4, "#rho^{0}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(5, "#rho^{#pm}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(6, "#eta");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(7, "#eta^{`}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(8, "#omega");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(9, "#phi");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(10, "K^{0}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(11, "K^{#pm}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(12, "K^{*0}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(13, "K^{*#pm}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(14, "D^{#pm}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(15, "D^{0}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(16, "#eta_{c}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(17, "#eta_{b}");
    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetBinLabel(18, "#Upsilon");
//    nBranches_->genParticle_Bdau_X_pdgId->SetTitle("Gen Level X PDGID in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_X_pdgId->GetXaxis()->SetTitle("PDGID");
//    nBranches_->genParticle_Bdau_X_pt->SetTitle("Gen Level X p_{T} in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_X_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
//    nBranches_->genParticle_Bdau_X_eta->SetTitle("Gen Level X #eta in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_X_eta->GetXaxis()->SetTitle("#eta");
//    nBranches_->genParticle_Bdau_X_phi->SetTitle("Gen Level X #phi in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_X_phi->GetXaxis()->SetTitle("#phi");
//    nBranches_->genParticle_Bdau_mu1_pt->SetTitle("Gen Level #mu_{J/#psi,1} p_{T} in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_mu1_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
//    nBranches_->genParticle_Bdau_mu1_eta->SetTitle("Gen Level #mu_{J/#psi,1} #eta in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_mu1_eta->GetXaxis()->SetTitle("#eta");
//    nBranches_->genParticle_Bdau_mu1_phi->SetTitle("Gen Level #mu_{J/#psi,1} #phi in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_mu1_phi->GetXaxis()->SetTitle("#phi");
//    nBranches_->genParticle_Bdau_mu2_pt->SetTitle("Gen Level #mu_{J/#psi,2} p_{T} in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_mu2_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
//    nBranches_->genParticle_Bdau_mu2_eta->SetTitle("Gen Level #mu_{J/#psi,2} #eta in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_mu2_eta->GetXaxis()->SetTitle("#eta");
//    nBranches_->genParticle_Bdau_mu2_phi->SetTitle("Gen Level #mu_{J/#psi,2} #phi in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_mu2_phi->GetXaxis()->SetTitle("#phi");
//    nBranches_->genParticle_Bdau_Jpsi_pt->SetTitle("Gen Level J/#psi p_{T} in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_Jpsi_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
//    nBranches_->genParticle_Bdau_Jpsi_eta->SetTitle("Gen Level J/#psi #eta in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_Jpsi_eta->GetXaxis()->SetTitle("#eta");
//    nBranches_->genParticle_Bdau_Jpsi_phi->SetTitle("Gen Level J/#psi #phi in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_Jpsi_phi->GetXaxis()->SetTitle("#phi");
//    nBranches_->genParticle_Bdau_Jpsi_mass->SetTitle("Gen Level J/#psi mass in B->J/#psi+X");
//    nBranches_->genParticle_Bdau_Jpsi_mass->GetXaxis()->SetTitle("Mass (GeV)");
//    nBranches_->genParticle_Bvis_pt->SetTitle("Gen Level Visible B p_{T} in B->J/#psi+X");
//    nBranches_->genParticle_Bvis_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
//    nBranches_->genParticle_Bvis_eta->SetTitle("Gen Level Visible B #eta in B->J/#psi+X");
//    nBranches_->genParticle_Bvis_eta->GetXaxis()->SetTitle("#eta");
//    nBranches_->genParticle_Bvis_phi->SetTitle("Gen Level Visible B #phi in B->J/#psi+X");
//    nBranches_->genParticle_Bvis_phi->GetXaxis()->SetTitle("#phi");
//    nBranches_->genParticle_Bvis_mass->SetTitle("Gen Level Visible B mass in B->J/#psi+X");
//    nBranches_->genParticle_Bvis_mass->GetXaxis()->SetTitle("Mass (GeV)");
    }
  /*=======================================================================================*/
// std::cout << "Histos Named" << std::endl;
 
  /*=======================================================================================*/
  if (runFlags["doMissingEt"]) {
    std::vector<std::string> corrFormulas;
    corrFormulas.push_back(iConfig.getParameter<std::string>("corrMetPx"));
    corrFormulas.push_back(iConfig.getParameter<std::string>("corrMetPy"));

    std::vector<std::string> jecAK4Labels;
    std::vector<std::string> tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK4forMetCorr");
    std::string tmpString = "";
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
       tmpString = jecpath + tmpVec[v];
       std::cout << " tmpString "<< tmpString <<std::endl;
       jecAK4Labels.push_back(edm::FileInPath(tmpString).fullPath());
    }
    
    nTuplizers_["MET"] = new METsNtuplizer( metToken_          , 
					    metpuppiToken_     , 
					    metmvaToken_     , 
                                            jetForMetCorrToken_, 
					    muonToken_         ,
					    rhoToken_	       ,
					    vtxToken_	       ,
					    metSigToken_       ,
					    metCovToken_       ,
					    jecAK4Labels       ,
                                            corrFormulas       ,
					    nBranches_         ,
					    runFlags  );
  }
    
  
  /*=======================================================================================*/  
  std::vector<edm::EDGetTokenT<reco::VertexCollection>> vtxTokens;
  vtxTokens.push_back( vtxToken_  );  


  /*=======================================================================================*/  

						      
  if (runFlags["doVertices"]) {
    std::cout<<"\n\n --->DOING MET<---\n\n"<<std::endl;
    nTuplizers_["vertices"] = new VerticesNtuplizer( vtxTokens   , 
						     beamToken_,
                                                     nBranches_  ,
						     runFlags    );
  }
  if (runFlags["doJpsiMu"]) {
    std::cout<<"\n\n --->GETTING INSIDE HERE<---\n\n"<<std::endl;
    nTuplizers_["JpsiMu"] = new JpsiMuNtuplizer( muonToken_   , 
						 vtxToken_   , 
						 packedpfcandidatesToken_,
                                                 triggerToken_,
                                                 triggerObjects_,
						 nBranches_  ,
                                                 runFlags    );
  }
  if (runFlags["doJpsiEle"]) {
    std::cout<<"\n\n --->GETTING INSIDE THE ELECTRON PART<---\n\n"<<std::endl;
    nTuplizers_["JpsiEle"] = new JpsiEleNtuplizer( electronToken_   , 
					     vtxToken_   , 
					     nBranches_ );
  }

  if (runFlags["doTriggerDecisions"] || runFlags["doTriggerObjects"] || runFlags["doTriggerDecisions"]) {
    nTuplizers_["triggers"] = new TriggersNtuplizer( triggerToken_, 
                                                     triggerObjects_, 
						     triggerPrescales_,
                                                     noiseFilterToken_,
						     HBHENoiseFilterLooseResultToken_,
						     HBHENoiseFilterTightResultToken_,
						     HBHENoiseIsoFilterResultToken_,
						     ecalBadCalibFilterUpdateToken_,
						     nBranches_,
                                                     iConfig,
                                                     runFlags );
  }

 
  /*=======================================================================================*/    
  if ( runFlags["runOnMC"] ){

     
    if (runFlags["doGenParticles"]) {
      std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> genpTokens;
      genpTokens.push_back( genparticleToken_ );

      nTuplizers_["genParticles"] = new GenParticlesNtuplizer( genpTokens, nBranches_, runFlags );
    }

    if (runFlags["doPileUp"]) {
      std::vector<edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > puTokens;
      puTokens.push_back( puinfoToken_ );
      nTuplizers_["PU"] = new PileUpNtuplizer( puTokens, nBranches_, runFlags );
    }

    if (runFlags["doGenEvent"]) {
      std::vector<edm::EDGetTokenT< GenEventInfoProduct > > geneTokens;
      geneTokens.push_back( geneventToken_ );
      std::vector<edm::EDGetTokenT<  LHEEventProduct > > lheTokens;
      lheTokens.push_back( lheEventProductToken_);
      nTuplizers_["genEvent"] = new GenEventNtuplizer( geneTokens, nBranches_ , lheTokens, runFlags);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::~Ntuplizer()
{
	  
  for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it ){
    //std::cout << "deconstructor: Branches: " << it->first << std::endl;
    delete it->second;
  }
   nTuplizers_.clear();
   
   delete nBranches_;
   
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  nBranches_->reset();

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if( vertices->empty() ) return; // skip the event if no PV found
  
  nBranches_->EVENT_event     = iEvent.id().event();
  nBranches_->EVENT_run       = iEvent.id().run();
  nBranches_->EVENT_lumiBlock = iEvent.id().luminosityBlock();  
  //std::cout<<"before the branches loop"<<std::endl; 
  for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it ){
    //std::cout << "Fill Branchines: " << it->first << std::endl;
    (it->second)->fillBranches( iEvent, iSetup );
  }
  nBranches_->fillTree();
  
  nBranches_->reset();    
  
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginJob(){
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endJob() {
std::cout << "Saving Histos" << std::endl;
if ( nBranches_->genParticle_Bdau_X_pt->GetEntries() > 0 ) {
   nBranches_->genParticle_Bdau_X_pdgId->Draw();
   nBranches_->genParticle_Bdau_X_pdgId->Write();
   nBranches_->genParticle_Bdau_X_pt->Draw();
   nBranches_->genParticle_Bdau_X_pt->Write();
   nBranches_->genParticle_Bdau_X_eta->Draw();
   nBranches_->genParticle_Bdau_X_eta->Write();
   nBranches_->genParticle_Bdau_X_phi->Draw();
   nBranches_->genParticle_Bdau_X_phi->Write();
   nBranches_->genParticle_Bdau_mu1_pt->Draw();
   nBranches_->genParticle_Bdau_mu1_pt->Write();
   nBranches_->genParticle_Bdau_mu1_eta->Draw();
   nBranches_->genParticle_Bdau_mu1_eta->Write();
   nBranches_->genParticle_Bdau_mu1_phi->Draw();
   nBranches_->genParticle_Bdau_mu1_phi->Write();
   nBranches_->genParticle_Bdau_mu2_pt->Draw();
   nBranches_->genParticle_Bdau_mu2_pt->Write();
   nBranches_->genParticle_Bdau_mu2_eta->Draw();
   nBranches_->genParticle_Bdau_mu2_eta->Write();
   nBranches_->genParticle_Bdau_mu2_phi->Draw();
   nBranches_->genParticle_Bdau_mu2_phi->Write();
   nBranches_->genParticle_Bdau_Jpsi_mass->Draw();
   nBranches_->genParticle_Bdau_Jpsi_mass->Write();
   nBranches_->genParticle_Bdau_Jpsi_pt->Draw();
   nBranches_->genParticle_Bdau_Jpsi_pt->Write();
   nBranches_->genParticle_Bdau_Jpsi_eta->Draw();
   nBranches_->genParticle_Bdau_Jpsi_eta->Write();
   nBranches_->genParticle_Bdau_Jpsi_phi->Draw();
   nBranches_->genParticle_Bdau_Jpsi_phi->Write();
   nBranches_->genParticle_Bvis_mass->Draw();
   nBranches_->genParticle_Bvis_mass->Write();
   nBranches_->genParticle_Bvis_pt->Draw();
   nBranches_->genParticle_Bvis_pt->Write();
   nBranches_->genParticle_Bvis_eta->Draw();
   nBranches_->genParticle_Bvis_eta->Write();
   nBranches_->genParticle_Bvis_phi->Draw();
   nBranches_->genParticle_Bvis_phi->Write();
   }
if ( nBranches_->cutflow_perevt->GetEntries() > 0 ) {
   nBranches_->cutflow_perevt->Draw();
   nBranches_->cutflow_perevt->Write();
   }
std::cout << "Histos Saved" << std::endl;

//if ( failhist->GetEntries()>0 ){
//   failhist->Draw();
//   failhist->Write();
//   }
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&){
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
