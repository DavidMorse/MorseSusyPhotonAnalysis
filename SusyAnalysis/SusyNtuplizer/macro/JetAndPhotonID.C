//----------JETS------------
      std::vector<susy::PFJet*>  pfJets;
    
      if(printLevel>0) cout << "Find pfJets in the event." << endl;
  
      std::map<TString,susy::PFJetCollection>::iterator pfJets_it = event->pfJets.find("ak5");

      if(pfJets_it != event->pfJets.end()){

	susy::PFJetCollection& jetColl = pfJets_it->second;

	for(std::vector<susy::PFJet>::iterator it = jetColl.begin(); it != jetColl.end(); it++) {
	  // set up corrections for PFJets
	  std::map<TString,Float_t>::iterator s_it_L1FastL2L3 = it->jecScaleFactors.find("L1FastL2L3");
	  std::map<TString,Float_t>::iterator s_it_L2L3       = it->jecScaleFactors.find("L2L3");
	  if (s_it_L1FastL2L3 == it->jecScaleFactors.end() || s_it_L2L3 == it->jecScaleFactors.end()) {
	    if (s_it_L1FastL2L3 == it->jecScaleFactors.end())cout << "L1FastL2L3 JEC is not available for this jet!!!" << endl;
	    if (s_it_L2L3 == it->jecScaleFactors.end()) cout << "L2L3 JEC is not available for this jet!!!" << endl;
	    continue;
	  }
	  float scale = s_it_L1FastL2L3->second;
	  float scaleMatch = s_it_L1FastL2L3->second/s_it_L2L3->second;

	  TLorentzVector corrP4 = scale * it->momentum;
	  TLorentzVector corrP4Match = scaleMatch * it->momentum;

	  float JetThreshold=30.;
	  if(   corrP4.Pt()>=JetThreshold
		&& std::abs(corrP4.Eta()) <= 2.5
		&& it->nConstituents>1
		//use uncorrected jet E for fractions
		&& it->neutralHadronEnergy/it->momentum.E()<0.99
		&& it->neutralEmEnergy/it->momentum.E()<0.99 
		&& it->passPuJetIdLoose(susy::kFull)
		) {
	    if(std::abs(corrP4.Eta()) >= 2.4){
	      //cout<<"Jet Eta:            "<<it->momentum.Eta()<<endl;
	      //set jet momentum to L1FastL2L3 corrected one
	      it->momentum=corrP4;
	      //cout<<"Jet corrP4 Eta:     "<<it->momentum.Eta()<<endl;
	      pfJets.first.push_back(&*it);
	    }
	    else if( it->chargedMultiplicity>0 
		     && it->chargedHadronEnergy/it->momentum.E()>0
		     && it->chargedEmEnergy/it->momentum.E()<0.99){
	      it->momentum=corrP4;
	      pfJets.first.push_back(&*it);
	    }
	  }
	}// for jet
      }// if not end

      std::sort(pfJets.begin(),pfJets.end(),EtGreater<susy::PFJet>);
  
      //------end JETS--------------
   if(printLevel>0) cout << "Find loose and tight photons in the event." << endl;
      //----------
      //find photons, sort by energy
      double leadEtCut=40.0,
	trailEtCut    =25.0,
	maxEta        =1.4442,
	maxR9         =1.0,
	maxHoverE     =0.05,
	maxSihih      =0.012,
	maxCombIso    =6.0,
	maxFakeCombIso=20.0,
	maxMETdphi    =2.8,
	//minJetMETdphi     =0.5;
	minJetMETdphi     =-999.;
      //----------
      std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
      if(phoMap != event->photons.end()) {
            
	std::vector<susy::Photon*>   pho_Cands;

	//loop over photon collection 
	for(std::vector<susy::Photon>::iterator it_Pho = phoMap->second.begin(); it_Pho != phoMap->second.end(); it_Pho++) {
	  if(printLevel>0) cout<<"looping over photon collection"<<endl;

	  bool phoCand=false;
	  //----------------set up cuts-------------------
	  float chargedHadronIso=it_Pho->chargedHadronIso;
	  float neutralHadronIso=it_Pho->neutralHadronIso;
	  float photonIso=it_Pho->photonIso;
	  float PhoEt=it_Pho->momentum.Et();
	  // fiducial cuts. Look for only barrel now
	  bool etaCut = (std::abs(it_Pho->caloPosition.Eta()) < maxEta);

	  // Spike cleaning
	  bool isSpike = (it_Pho->r9>maxR9 || it_Pho->sigmaIetaIeta<=0.001 || it_Pho->sigmaIphiIphi<=0.001);

	  // Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	  bool EtCut = (PhoEt > trailEtCut);
	  bool EtCutLead = (PhoEt > leadEtCut);

	  // H/E
	  bool heCut = (it_Pho->hadTowOverEm < maxHoverE);
	  // sigma_ietaieta
	  bool sIetaCut = (it_Pho->sigmaIetaIeta < maxSihih);

	  //PF Isolation
	  bool Iso2012IdCutLoose = false,Iso2012IdCutMed = false,Iso2012IdCutTight = false;
	  if(fabs(it_Pho->caloPosition.Eta())<1.0){
	    Iso2012IdCutLoose=((chargedHadronIso<2.6+0.012*Rho25) && (neutralHadronIso<3.5+PhoEt*0.04+0.03*Rho25) && (photonIso<1.3+PhoEt*0.005+.148*Rho25));
	    Iso2012IdCutMed=((chargedHadronIso<1.5+0.012*Rho25) && (neutralHadronIso<1.0+PhoEt*0.04+0.03*Rho25) && (photonIso<0.7+PhoEt*0.005+.148*Rho25));
	    Iso2012IdCutTight=((chargedHadronIso<0.7+0.012*Rho25) && (neutralHadronIso<0.4+PhoEt*0.04+0.03*Rho25) && (photonIso<0.5+PhoEt*0.005+.148*Rho25));
	  }
	  else{
	    Iso2012IdCutLoose=((chargedHadronIso<2.6+0.010*Rho25) && (neutralHadronIso<3.5+PhoEt*0.04+0.057*Rho25) && (photonIso<1.3+PhoEt*0.005+.130*Rho25));
	    Iso2012IdCutMed=((chargedHadronIso<1.5+0.010*Rho25) && (neutralHadronIso<1.0+PhoEt*0.04+0.057*Rho25) && (photonIso<0.7+PhoEt*0.005+.130*Rho25));
	    Iso2012IdCutTight=((chargedHadronIso<0.7+0.010*Rho25) && (neutralHadronIso<0.4+PhoEt*0.04+0.057*Rho25) && (photonIso<0.5+PhoEt*0.005+.130*Rho25));
	  }

	  bool PhoCutPF =  ( etaCut && EtCut && !isSpike && heCut && pixelCut && Iso2012IdCutLoose && sIetaCut);

	  //Now fill pho_Cands
	  if( PhoCutPF){
	      pho_Cands.push_back(&*it_Pho);
	  }
     
	  if(printLevel>0) cout<<"End of Photon Loop"<<endl;
	}//for(it_Pho)

	//sort pho_Cands and fake_Cands by Pt
	std::sort(pho_Cands.begin(), pho_Cands.end(), EtGreater<susy::Photon>);
	if(printLevel>0)cout<<"phoCands size= "<<pho_Cands.size()<<endl;
