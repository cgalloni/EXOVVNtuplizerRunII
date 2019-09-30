#include "../interface/GenParticlesNtuplizer.h"
 
//===================================================================================================================        
GenParticlesNtuplizer::GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags ) 

   : CandidateNtuplizer( nBranches )
   , genParticlesToken_( tokens[0] )
   , isJpsiMu_( runFlags["doJpsiMu"])
   , isJpsiEle_( runFlags["doJpsiEle"]  )
   , isGenHist_( runFlags["doGenHist"]  )
{

}

//===================================================================================================================        
GenParticlesNtuplizer::~GenParticlesNtuplizer( void )
{
}

//===================================================================================================================        
void GenParticlesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  //chunk to remove those events with no jspi if that analysis is chosen
  if (!isGenHist_) {
  std::vector<int> doJpsi_;
  if(isJpsiEle_) {
    doJpsi_ = nBranches_->IsJpsiEle;
    //std::cout<<"im getting inside the electron part"<<std::endl;
  }else if(isJpsiMu_){
    doJpsi_ = nBranches_->IsJpsiMu;
    //std::cout<<"nbranch thing\t"<<size(isJpsi_)<<"; "<< isJpsi_[0]<<std::endl;
  }
  if(size(doJpsi_)>0) if(doJpsi_[0]==0) return;
  }
   

    event.getByToken(genParticlesToken_ , genParticles_); 

   /* here we want to save  gen particles info*/
   
    std::vector<int> vDau ;
    std::vector<int> vMoth;
    int nMoth = 0;
    int nDau  = 0;  
    //nBranches_->genParticle_N = genParticles_->size(); // the genParticles are filtered below
    for( unsigned p=0; p<genParticles_->size(); ++p ){
      //if( (*genParticles_)[p].status() != 3 ) continue;
      vDau.clear(); vMoth.clear();
      nDau = 0; nMoth = 0;
      
      bool isPrompt( (*genParticles_)[p].statusFlags().isPrompt() );
      bool isDirectPromptTauDecayProduct( (*genParticles_)[p].statusFlags().isDirectPromptTauDecayProduct() );
      bool fromHardProcessFinalState( (*genParticles_)[p].fromHardProcessFinalState() );
      bool isDirectHardProcessTauDecayProductFinalState( (*genParticles_)[p].isDirectHardProcessTauDecayProductFinalState() );
      bool isLepton( abs((*genParticles_)[p].pdgId())>=11 && abs((*genParticles_)[p].pdgId())<=18 );
      bool isQuark( abs((*genParticles_)[p].pdgId())<=6 && abs((*genParticles_)[p].status())<=29 );
      bool isPhoton( abs((*genParticles_)[p].pdgId())==22 && (*genParticles_)[p].pt()>10. );
      bool isGluon( abs((*genParticles_)[p].pdgId())==22 && (*genParticles_)[p].pt()>10. );
      bool isWZH( abs((*genParticles_)[p].pdgId())>=23 && abs((*genParticles_)[p].pdgId())<=25 );
      bool isHeavyMeson( abs((*genParticles_)[p].pdgId())>0 && abs((*genParticles_)[p].pdgId())<=1000 );
      bool isHeavyBaryon( abs((*genParticles_)[p].pdgId())>=1000);
      bool isBSM( (abs((*genParticles_)[p].pdgId())>=30 && abs((*genParticles_)[p].pdgId())<=50) || abs((*genParticles_)[p].pdgId())>=1000000 );
      
      if(!isLepton && !isQuark && !isPhoton && !isGluon && !isWZH && !isHeavyMeson && !isHeavyBaryon && !isBSM && !isDirectPromptTauDecayProduct && !fromHardProcessFinalState && !isDirectHardProcessTauDecayProductFinalState) continue;
      
//      nBranches_->genParticle_px    .push_back((*genParticles_)[p].px()     );
//      nBranches_->genParticle_py    .push_back((*genParticles_)[p].py()     );
//      nBranches_->genParticle_pz    .push_back((*genParticles_)[p].pz()     );
//      nBranches_->genParticle_e     .push_back((*genParticles_)[p].energy() );
      nBranches_->genParticle_pt    .push_back((*genParticles_)[p].pt()     );
      nBranches_->genParticle_eta   .push_back((*genParticles_)[p].eta()    );
      nBranches_->genParticle_phi   .push_back((*genParticles_)[p].phi()    );
      nBranches_->genParticle_mass  .push_back((*genParticles_)[p].mass()   );
      nBranches_->genParticle_status.push_back((*genParticles_)[p].status() );
      nBranches_->genParticle_pdgId .push_back((*genParticles_)[p].pdgId()  );


      // needed for the gen matching
      nBranches_->genParticle_isPrompt.push_back( isPrompt );
      nBranches_->genParticle_isDirectPromptTauDecayProduct.push_back( isDirectPromptTauDecayProduct );

      // needed for the MVA recoil correction
      nBranches_->genParticle_fromHardProcessFinalState.push_back( fromHardProcessFinalState );
      nBranches_->genParticle_isDirectHardProcessTauDecayProductFinalState.push_back( isDirectHardProcessTauDecayProductFinalState );


      for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
        vDau.push_back( (*genParticles_)[p].daughter(d)->pdgId() );
	      nDau++;
      }

if ( isGenHist_ ) {
// Looking At B mesons who decay to Jpsi+X. Catalogue what else they decay to in addition to the Jpsi. Get the particle's pdgid,
// the pT, eta, and phi of it and the two muons from the jpsi, the jpsi's (aka dimuon) pt, eta, phi, and mass, and the B's
// Visible pt, eta, phi, and mass
      if ( (  abs((*genParticles_)[p].pdgId()) == 511
           || abs((*genParticles_)[p].pdgId()) == 521 
           || abs((*genParticles_)[p].pdgId()) == 513 
           || abs((*genParticles_)[p].pdgId()) == 523 
           || abs((*genParticles_)[p].pdgId()) == 515 
           || abs((*genParticles_)[p].pdgId()) == 525 
           || abs((*genParticles_)[p].pdgId()) == 531 
           || abs((*genParticles_)[p].pdgId()) == 533 
           || abs((*genParticles_)[p].pdgId()) == 535 
           || abs((*genParticles_)[p].pdgId()) == 541 
           || abs((*genParticles_)[p].pdgId()) == 543 
           || abs((*genParticles_)[p].pdgId()) == 545 )
         && (*genParticles_)[p].status() == 2 ) {
         TLorentzVector mu1, mu2, mu3;
         std::vector<int> mu3pdgid;
         for (unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ) {
//             std::cout << p << " B Daughter: " << (*genParticles_)[p].daughter(d)->pdgId() 
//                       << " status: " << (*genParticles_)[p].daughter(d)->status() << std::endl;
             if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 12 
               || abs((*genParticles_)[p].daughter(d)->pdgId()) == 14 
               || abs((*genParticles_)[p].daughter(d)->pdgId()) == 16 ) {continue;} 
             if ( (*genParticles_)[p].daughter(d)->pdgId() == 443 ) { 
// Loop over jpsi daughters & get the two mus. if there aren't two mus then skip it
                for ( unsigned int jd=0; jd<(*genParticles_)[p].daughter(d)->numberOfDaughters(); ++jd ) {
//                    std::cout<<"J/psi Daughter: " << (*genParticles_)[p].daughter(d)->daughter(jd)->pdgId() << std::endl;
                    if ( (*genParticles_)[p].daughter(d)->daughter(jd)->pdgId() == 13 ) {
                       mu1.SetPtEtaPhiM((*genParticles_)[p].daughter(d)->daughter(jd)->pt(),
                                        (*genParticles_)[p].daughter(d)->daughter(jd)->eta(),
                                        (*genParticles_)[p].daughter(d)->daughter(jd)->phi(), 
                                        (*genParticles_)[p].daughter(d)->daughter(jd)->mass());
                       } 
                    else if ( (*genParticles_)[p].daughter(d)->daughter(jd)->pdgId() == -13 ) {
                       mu2.SetPtEtaPhiM((*genParticles_)[p].daughter(d)->daughter(jd)->pt(), 
                                        (*genParticles_)[p].daughter(d)->daughter(jd)->eta(),
                                        (*genParticles_)[p].daughter(d)->daughter(jd)->phi(), 
                                        (*genParticles_)[p].daughter(d)->daughter(jd)->mass());
                       }
                    }
                 } 
             }
          if (mu1.Pt() > 0 && mu2.Pt() > 0) {
//             std::cout<<"Muon pT 1: " << mu1.Pt() << " Muon pT 2: " << mu2.Pt() << std::endl;
             for (unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ) {
//                  std::cout << (*genParticles_)[p].daughter(d)->pdgId() <<
//                  " Particle Status: " << (*genParticles_)[p].daughter(d)->status() << std::endl;
//                  std::cout<<"Particle: "<< (*genParticles_)[p].daughter(d)->pdgId()<<std::endl;
//	            for (unsigned int m=0; m<(*genParticles_)[p].daughter(d)->numberOfMothers(); ++m) {
//	                std::cout<<"Particle Mother: "<< (*genParticles_)[p].daughter(d)->mother(m)->pdgId()<<std::endl;
//	                }
//                  for (unsigned int xd=0; xd<(*genParticles_)[p].numberOfDaughters(); ++xd) {
//                      std::cout << "X daughter: " << (*genParticles_)[p].daughter(d)->daughter(xd)->pdgId() 
//                                << " Status: " << (*genParticles_)[p].daughter(d)->daughter(xd)->status() << std::endl;
//                      }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 12 
                    || abs((*genParticles_)[p].daughter(d)->pdgId()) == 14 
                    || abs((*genParticles_)[p].daughter(d)->pdgId()) == 16 ) {continue;} 
                  if ( (*genParticles_)[p].daughter(d)->pdgId() == 443 ) {continue;}
//                  nBranches_->genParticle_Bdau_X_pdgId.push_back(
//                                    (*genParticles_)[p].daughter(d)->pdgId());
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 13) { 
                  //mu+
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(0); 
                  }
                  if ( (*genParticles_)[p].daughter(d)->pdgId() == 111) { 
                  //pi0
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(1); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 211) { 
                  //pi+
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(2); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 113) { 
                  //rho0
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(3); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 213) { 
                  //rho+
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(4); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 221) { 
                  //eta
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(5); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 331) { 
                  //eta'
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(6); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 223) { 
                  //omega
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(7); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 333) { 
                  //phi
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(8); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 311) { 
                  //K0
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(9); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 321) { 
                  //K+
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(10); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 313) { 
                  //K*0
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(11); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 323) { 
                  //K*+
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(12); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 411) { 
                  //D+
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(13); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 421) { 
                  //D0
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(14); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 441) { 
                  //eta_c
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(15); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 551) { 
                  //eta_b
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(16); 
                  }
                  if ( abs((*genParticles_)[p].daughter(d)->pdgId()) == 553) { 
                  //Gamma
                  nBranches_->genParticle_Bdau_X_pdgId->Fill(17); 
                  }
                  TLorentzVector temp;
                  temp.SetPtEtaPhiM((*genParticles_)[p].daughter(d)->pt(),
                                    (*genParticles_)[p].daughter(d)->eta(),
                                    (*genParticles_)[p].daughter(d)->phi(),
                                    (*genParticles_)[p].daughter(d)->mass());
                  mu3 += temp;
                  }
             if (mu3.Pt() > 0) {
                nBranches_->genParticle_Bdau_X_pt->Fill(mu3.Pt()); 
                nBranches_->genParticle_Bdau_X_eta->Fill(mu3.Eta()); 
                nBranches_->genParticle_Bdau_X_phi->Fill(mu3.Phi()); 
                nBranches_->genParticle_Bdau_mu1_pt->Fill(mu1.Pt()); 
                nBranches_->genParticle_Bdau_mu1_eta->Fill(mu1.Eta()); 
                nBranches_->genParticle_Bdau_mu1_phi->Fill(mu1.Phi()); 
                nBranches_->genParticle_Bdau_mu2_pt->Fill(mu2.Pt()); 
                nBranches_->genParticle_Bdau_mu2_eta->Fill(mu2.Eta()); 
                nBranches_->genParticle_Bdau_mu2_phi->Fill(mu2.Phi()); 
                nBranches_->genParticle_Bdau_Jpsi_mass->Fill((mu1+mu2).M()); 
                nBranches_->genParticle_Bdau_Jpsi_pt->Fill((mu1+mu2).Pt()); 
                nBranches_->genParticle_Bdau_Jpsi_eta->Fill((mu1+mu2).Eta()); 
                nBranches_->genParticle_Bdau_Jpsi_phi->Fill((mu1+mu2).Phi()); 
                nBranches_->genParticle_Bvis_mass->Fill((mu1+mu2+mu3).M()); 
                nBranches_->genParticle_Bvis_pt->Fill((mu1+mu2+mu3).Pt()); 
                nBranches_->genParticle_Bvis_eta->Fill((mu1+mu2+mu3).Eta()); 
                nBranches_->genParticle_Bvis_phi->Fill((mu1+mu2+mu3).Phi()); 
//                nBranches_->genParticle_Bdau_X_pt.push_back(mu3.Pt());
//                nBranches_->genParticle_Bdau_X_eta.push_back(mu3.Eta());
//                nBranches_->genParticle_Bdau_X_phi.push_back(mu3.Phi());
//                nBranches_->genParticle_Bdau_mu1_pt.push_back(mu1.Pt());
//                nBranches_->genParticle_Bdau_mu1_eta.push_back(mu1.Eta());
//                nBranches_->genParticle_Bdau_mu1_phi.push_back(mu1.Phi());
//                nBranches_->genParticle_Bdau_mu2_pt.push_back(mu2.Pt());
//                nBranches_->genParticle_Bdau_mu2_eta.push_back(mu2.Eta());
//                nBranches_->genParticle_Bdau_mu2_phi.push_back(mu2.Phi());
//                nBranches_->genParticle_Bdau_Jpsi_pt.push_back((mu1+mu2).Pt());
//                nBranches_->genParticle_Bdau_Jpsi_eta.push_back((mu1+mu2).Eta());
//                nBranches_->genParticle_Bdau_Jpsi_phi.push_back((mu1+mu2).Phi());
//                nBranches_->genParticle_Bdau_Jpsi_mass.push_back((mu1+mu2).M());
//                nBranches_->genParticle_Bvis_pt.push_back((mu1+mu2+mu3).Pt());
//                nBranches_->genParticle_Bvis_eta.push_back((mu1+mu2+mu3).Eta());
//                nBranches_->genParticle_Bvis_phi.push_back((mu1+mu2+mu3).Phi());
//                nBranches_->genParticle_Bvis_mass.push_back((mu1+mu2+mu3).M());
                }
             }
         }
}

      if(abs((*genParticles_)[p].pdgId())==15 && (*genParticles_)[p].statusFlags().isPrompt() && (*genParticles_)[p].status()==2){

        if(nDau>1){

          bool flag_radioactive_gamma = false;
          bool flag_radioactive_tau = false;
          
          for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
            Int_t taupdgId = abs((*genParticles_)[p].daughter(d)->pdgId());
            
            if(taupdgId==22) flag_radioactive_gamma = true;   
            if(taupdgId==15) flag_radioactive_tau = true; 
              
          }

          if(!(flag_radioactive_gamma && flag_radioactive_tau)){

            TLorentzVector tau;
            tau.SetPtEtaPhiM(0,0,0,0);
            Int_t decaymode = -1;
          
            for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
              Float_t taupt = (*genParticles_)[p].daughter(d)->pt();
              Float_t taueta = (*genParticles_)[p].daughter(d)->eta();
              Float_t tauphi = (*genParticles_)[p].daughter(d)->phi();
              Float_t taumass = (*genParticles_)[p].daughter(d)->mass();
              Int_t taupdgId = abs((*genParticles_)[p].daughter(d)->pdgId());
              
              TLorentzVector taudau;
              taudau.SetPtEtaPhiM(taupt, taueta, tauphi, taumass);
              if(!(taupdgId >= 11 && taupdgId<=16)){
	              tau += taudau;
	              decaymode = 4;
              }
              if(taupdgId==11){ // electron decay
	              decaymode = 2;
	              tau += taudau;
	            }
              if(taupdgId==13){ // muon decay
		            decaymode = 3;
	              tau += taudau;
	            }
          
            }

            nBranches_->genParticle_tauvispt  .push_back( (float)tau.Pt()  );
            nBranches_->genParticle_tauviseta  .push_back( (float)tau.Eta()  );
            nBranches_->genParticle_tauvisphi  .push_back( (float)tau.Phi()  );
            nBranches_->genParticle_tauvismass  .push_back( (float)tau.M()  );
            nBranches_->genParticle_taudecay  .push_back( decaymode  );
          }
          else{
            nBranches_->genParticle_tauvispt  .push_back( -99.  );
            nBranches_->genParticle_tauviseta  .push_back( -99.  );
            nBranches_->genParticle_tauvisphi  .push_back( -99.  );
            nBranches_->genParticle_tauvismass  .push_back( -99.  );
            nBranches_->genParticle_taudecay  .push_back( 0  ); // self decay (tau -> tau)
          }
          
        }else{
          nBranches_->genParticle_tauvispt  .push_back( -99.  );
          nBranches_->genParticle_tauviseta  .push_back( -99.  );
          nBranches_->genParticle_tauvisphi  .push_back( -99.  );
          nBranches_->genParticle_tauvismass  .push_back( -99.  );
          
          nBranches_->genParticle_taudecay  .push_back( -1  ); // self decay (tau -> tau)
          
        }
      }
      

      for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
        vMoth.push_back( (*genParticles_)[p].mother(m)->pdgId() );
    	  nMoth++;
      }

      nBranches_->genParticle_nDau  .push_back( nDau  );
      nBranches_->genParticle_nMoth .push_back( nMoth );      
      nBranches_->genParticle_mother.push_back( vMoth );
      nBranches_->genParticle_dau   .push_back( vDau  );      

    }

    nBranches_->genParticle_N = nBranches_->genParticle_pt.size(); // save number of save genParticles
    
}

