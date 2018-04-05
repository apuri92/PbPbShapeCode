#define TriggerDef_cxx
#include "pPbFragmentation/pPbFragmentation.h"

using namespace std;

void BaseClass::SetTrigger_chains(){
	
	if (_dataset == 3){
	// pp 2015, 5.02TeV
		if(_isMB)
		{
			_nTriggers=1;
			trigger_chains.push_back("HLT_mb_sptrk");
			
			trigger_thresholds.push_back(0);
			
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(0); (jet_pt_trig[0]).push_back(10000.);
		}
		else if (!_doForward)
		{
			_nTriggers=7;
			_trigger_collection="a4tcemsubjesFS";
			trigger_chains.push_back("HLT_j20");
			trigger_chains.push_back("HLT_j30_L1TE5");
			trigger_chains.push_back("HLT_j40_L1TE10");
            trigger_chains.push_back("HLT_j50_L1J12");
			trigger_chains.push_back("HLT_j60_L1J15");
			trigger_chains.push_back("HLT_j75_L1J20");
			trigger_chains.push_back("HLT_j85");
           
		
			trigger_thresholds.push_back(20);
			trigger_thresholds.push_back(30);
			trigger_thresholds.push_back(40);
			trigger_thresholds.push_back(50);
			trigger_thresholds.push_back(60);
			trigger_thresholds.push_back(75);
			trigger_thresholds.push_back(85);
		
            
			jet_pt_trig.resize(_nTriggers);
            

			(jet_pt_trig[0]).push_back(26.0); (jet_pt_trig[0]).push_back(35.0);
			(jet_pt_trig[1]).push_back(35.0); (jet_pt_trig[1]).push_back(44.5);
			(jet_pt_trig[2]).push_back(44.5); (jet_pt_trig[2]).push_back(59.0);
			(jet_pt_trig[3]).push_back(59.0); (jet_pt_trig[3]).push_back(70.0);
			(jet_pt_trig[4]).push_back(70.0); (jet_pt_trig[4]).push_back(79.0);
			(jet_pt_trig[5]).push_back(79.0); (jet_pt_trig[5]).push_back(89.0);
			(jet_pt_trig[6]).push_back(89.0); (jet_pt_trig[6]).push_back(10000.);
			
			if (_addmuons){
			    _nTriggers=13;
				trigger_chains.push_back("HLT_mu4");
				trigger_chains.push_back("HLT_mu4_mu4noL1");
				trigger_chains.push_back("HLT_mu4noL1");
				trigger_chains.push_back("HLT_mu6");
				trigger_chains.push_back("HLT_mu8");
				trigger_chains.push_back("HLT_mu14");
				
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
						
				jet_pt_trig.resize(_nTriggers);
				(jet_pt_trig[7]).push_back(0); (jet_pt_trig[7]).push_back(10000.);
				(jet_pt_trig[8]).push_back(0); (jet_pt_trig[8]).push_back(10000.);
				(jet_pt_trig[9]).push_back(0); (jet_pt_trig[9]).push_back(10000.);
				(jet_pt_trig[10]).push_back(0); (jet_pt_trig[10]).push_back(10000.);
				(jet_pt_trig[11]).push_back(0); (jet_pt_trig[11]).push_back(10000.);
				(jet_pt_trig[12]).push_back(0); (jet_pt_trig[12]).push_back(10000.);
				
			}	
		}
		else
		{
			_nTriggers=6;
			_trigger_collection="a4tcemsubjesFS";
			trigger_chains.push_back("HLT_j10_320eta490");
			trigger_chains.push_back("HLT_j15_320eta490");
			trigger_chains.push_back("HLT_j25_320eta490");
			trigger_chains.push_back("HLT_j35_320eta490");
			trigger_chains.push_back("HLT_j45_320eta490");
			trigger_chains.push_back("HLT_j55_320eta490");
		
			trigger_thresholds.push_back(10);
			trigger_thresholds.push_back(15);
			trigger_thresholds.push_back(25);
			trigger_thresholds.push_back(35);
			trigger_thresholds.push_back(45);
			trigger_thresholds.push_back(55);
		
		
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(20); (jet_pt_trig[0]).push_back(25);
			(jet_pt_trig[1]).push_back(25); (jet_pt_trig[1]).push_back(35);
			(jet_pt_trig[2]).push_back(35); (jet_pt_trig[2]).push_back(45);
			(jet_pt_trig[3]).push_back(45); (jet_pt_trig[3]).push_back(55);
			(jet_pt_trig[4]).push_back(55); (jet_pt_trig[4]).push_back(65);
			(jet_pt_trig[5]).push_back(65); (jet_pt_trig[5]).push_back(10000.);	
		}
	}
	
	//PbPb 2015, 5.02TeV
	if (_dataset == 4){
		//Read prescales sets
		//Order of Y bins starting from bin #3: "HLT_j40_ion_L1TE20", "HLT_j50_ion_L1TE20", "HLT_j60_ion_L1TE50", "HLT_j75_ion_L1TE50", "HLT_j100_ion_L1TE50","HLT_j20", "HLT_j30_L1TE5", "HLT_j40_L1TE10",  "HLT_j50_L1J12","HLT_j60_L1J15","HLT_j75_L1J20","HLT_j85","HLT_mb_sptrk_ion_L1ZDC_A_C_VTE50","HLT_noalg_mb_L1TE50","HLT_mb_sptrk"]
		cout << "Setting triggers....";
		//First trigger for PbPb FF
		_first_trigger = 1; // <=> j40
		
		
		if(_isMB)
		{
			if(_isOverlay){
				_nTriggers=4;
				_trigger_collection="a4ionemsubjesFS";
				trigger_chains.push_back("HLT_mb_sptrk_ion_L1ZDC_A_C_VTE50_OVERLAY");
				trigger_chains.push_back("HLT_noalg_L1TE50_OVERLAY");
				trigger_chains.push_back("HLT_noalg_mb_L1TE1500.0ETA49_OVERLAY");
				trigger_chains.push_back("HLT_noalg_mb_L1TE6500.0ETA49_OVERLAY");
					
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
			
				trigger_PS.push_back(1.);
				trigger_PS.push_back(1.);
				trigger_PS.push_back(1.);
				trigger_PS.push_back(1.);
			
				jet_pt_trig.resize(_nTriggers);
				(jet_pt_trig[0]).push_back(0); (jet_pt_trig[0]).push_back(10000.);
				(jet_pt_trig[1]).push_back(0); (jet_pt_trig[1]).push_back(10000.);
				(jet_pt_trig[2]).push_back(0); (jet_pt_trig[2]).push_back(10000.);
				(jet_pt_trig[3]).push_back(0); (jet_pt_trig[3]).push_back(10000.);
			}
			else{	
				_nTriggers=2;
				_trigger_collection="a4ionemsubjesFS";
				trigger_chains.push_back("HLT_noalg_mb_L1TE50");
				trigger_chains.push_back("HLT_mb_sptrk_ion_L1ZDC_A_C_VTE50");
					
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
			
				trigger_PS.push_back(21.857);
				trigger_PS.push_back(21.51);
			
				jet_pt_trig.resize(_nTriggers);
				(jet_pt_trig[0]).push_back(0); (jet_pt_trig[0]).push_back(10000.);
				(jet_pt_trig[1]).push_back(0); (jet_pt_trig[1]).push_back(10000.);
			}	
		}
		else
		{
			_nTriggers=5;
			_trigger_collection="a4ionemsubjesFS";
			trigger_chains.push_back("HLT_j30_ion_L1TE20");//0
			trigger_chains.push_back("HLT_j40_ion_L1TE20");//1
			trigger_chains.push_back("HLT_j50_ion_L1TE20");//2
			trigger_chains.push_back("HLT_j60_ion_L1TE50");//3
			trigger_chains.push_back("HLT_j75_ion_L1TE50");//4
			//trigger_chains.push_back("HLT_j100_ion_L1TE50");//5
					
			trigger_thresholds.push_back(30);
			trigger_thresholds.push_back(40);
			trigger_thresholds.push_back(50);
			trigger_thresholds.push_back(60);
			trigger_thresholds.push_back(75);
			//trigger_thresholds.push_back(100);
			
			
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(45);   (jet_pt_trig[0]).push_back(58); //TODO: need to be reevaluated
			(jet_pt_trig[1]).push_back(58); (jet_pt_trig[1]).push_back(68); 
			(jet_pt_trig[2]).push_back(68); (jet_pt_trig[2]).push_back(82); 
			(jet_pt_trig[3]).push_back(82); (jet_pt_trig[3]).push_back(91); 
			(jet_pt_trig[4]).push_back(91); (jet_pt_trig[4]).push_back(10000.);
			//(jet_pt_trig[5]).push_back(116.); (jet_pt_trig[5]).push_back(10000.);
			
			if (_addmuons){
				_nTriggers=10;
				trigger_chains.push_back("HLT_mu4");
				trigger_chains.push_back("HLT_mu4_mu4noL1");
				trigger_chains.push_back("HLT_mu4noL1");
				trigger_chains.push_back("HLT_mu6");
				trigger_chains.push_back("HLT_mu8");
				
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
				trigger_thresholds.push_back(0);
						
				jet_pt_trig.resize(_nTriggers);
				(jet_pt_trig[5]).push_back(0); (jet_pt_trig[5]).push_back(10000.);
				(jet_pt_trig[6]).push_back(0); (jet_pt_trig[6]).push_back(10000.);
				(jet_pt_trig[7]).push_back(0); (jet_pt_trig[7]).push_back(10000.);
				(jet_pt_trig[8]).push_back(0); (jet_pt_trig[8]).push_back(10000.);
				(jet_pt_trig[9]).push_back(0); (jet_pt_trig[9]).push_back(10000.);
				
			}
		}
	}
	
	//pPb @ 8TeV
	if (_dataset == 5){
		
		cout << "Setting triggers....";
		_first_trigger = 1;
		
		
		if(_isMB)
		{
			_nTriggers=1;
			_trigger_collection="a4ionemsubjesFS";
			trigger_chains.push_back("HLT_mb_sptrk_L1MBTS_1");
					
			trigger_thresholds.push_back(0);
			
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(0); (jet_pt_trig[0]).push_back(10000.);
		}
		else
		{
			_nTriggers=10;
			_trigger_collection="a4ionemsubjesFS";
			// 313063-313603
			trigger_chains.push_back("HLT_j40_ion_L1J5");//0
			trigger_chains.push_back("HLT_j50_ion_L1J10");//1
			trigger_chains.push_back("HLT_j60_ion_L1J20");//2
			trigger_chains.push_back("HLT_j75_ion_L1J20");//3
			trigger_chains.push_back("HLT_j100_ion_L1J20");//4
			//313629-314170
			trigger_chains.push_back("HLT_j40_L1J5");//1
			trigger_chains.push_back("HLT_j50_L1J10");//2
			trigger_chains.push_back("HLT_j60");//3
			trigger_chains.push_back("HLT_j75_L1J20");//4
			trigger_chains.push_back("HLT_j100_L1J20");//5
			
			
			trigger_thresholds.push_back(41);
			trigger_thresholds.push_back(51);
			trigger_thresholds.push_back(61);
			trigger_thresholds.push_back(76);
			trigger_thresholds.push_back(101);
			
			trigger_thresholds.push_back(40);
			trigger_thresholds.push_back(50);
			trigger_thresholds.push_back(60);
			trigger_thresholds.push_back(75);
			trigger_thresholds.push_back(100);
			
			
			jet_pt_trig.resize(_nTriggers);
			//TODO: these ranges are basically crap, need to be changed for any real analysis
			(jet_pt_trig[0]).push_back(45);  (jet_pt_trig[0]).push_back(55); 
			(jet_pt_trig[1]).push_back(55);  (jet_pt_trig[1]).push_back(65); 
			(jet_pt_trig[2]).push_back(65);  (jet_pt_trig[2]).push_back(75); 
			(jet_pt_trig[3]).push_back(75);  (jet_pt_trig[3]).push_back(110); 
			(jet_pt_trig[4]).push_back(110); (jet_pt_trig[4]).push_back(10000.);
			
			(jet_pt_trig[5]).push_back(45);  (jet_pt_trig[5]).push_back(55); 
			(jet_pt_trig[6]).push_back(55);  (jet_pt_trig[6]).push_back(65); 
			(jet_pt_trig[7]).push_back(65);  (jet_pt_trig[7]).push_back(75); 
			(jet_pt_trig[8]).push_back(75);  (jet_pt_trig[8]).push_back(110); 
			(jet_pt_trig[9]).push_back(110); (jet_pt_trig[9]).push_back(10000.);
		}
	}
	
}

void BaseClass::SetTrigger_hist(TH2D* h){
	for(int i=1; i<=h->GetNbinsX();i++){
		//h->GetXaxis()->SetBinLabel(i,Form("j%.0f",trigger_thresholds.at(i-1)));
		h->GetXaxis()->SetBinLabel(i,Form("%s",trigger_chains.at(i-1).c_str()) );
	}
}
