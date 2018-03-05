#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <string>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
TriggerObjects TriggerTools::matchHLT(const double eta, const double phi, 
				      const std::vector<TriggerRecord>& triggerRecords,
				      const trigger::TriggerEvent& triggerEvent)
{
  const double dRMax = 0.2;

  TriggerObjects matchBits;
  for(unsigned int irec=0; irec<triggerRecords.size(); irec++) {     
    for(unsigned int iobj=0; iobj<triggerRecords[irec].objectMap.size(); iobj++) {
      const std::string   filterName = triggerRecords[irec].objectMap[iobj].first;
      const unsigned int  filterBit  = triggerRecords[irec].objectMap[iobj].second;
      
      edm::InputTag filterTag(filterName,"","HLT");
      // filterIndex must be less than the size of trgEvent or you get a CMSException: _M_range_check
      if(triggerEvent.filterIndex(filterTag) < triggerEvent.sizeFilters()) {
        const trigger::TriggerObjectCollection& toc(triggerEvent.getObjects());      
        const trigger::Keys& keys(triggerEvent.filterKeys(triggerEvent.filterIndex(filterTag)));
        
        for(unsigned int hlto=0; hlto<keys.size(); hlto++) {
          trigger::size_type hltf = keys[hlto];
          const trigger::TriggerObject& tobj(toc[hltf]);
          if(reco::deltaR(eta,phi,tobj.eta(),tobj.phi()) < dRMax) {
            matchBits[filterBit] = 1;
          }
        }
      }
    }
  }
  
  return matchBits;
}

TriggerObjects TriggerTools::matchHLT(const double eta, const double phi,  
				      const edm::Event &iEvent,
                                      const std::vector<TriggerRecord> &triggerRecords,
				      const edm::TriggerNames &triggerNames,
				      const edm::TriggerResults &triggerResults,
                                      pat::TriggerObjectStandAloneCollection &triggerObjects)
{

  const double dRMax = 0.2;
  TriggerObjects matchBits;  

  for (unsigned int irec=0; irec<triggerRecords.size(); irec++) {
    for (unsigned int iobj=0; iobj<triggerRecords[irec].objectMap.size(); iobj++) {
      std::vector<std::string> dumb; dumb.push_back(triggerRecords[irec].objectMap[iobj].first);
      const unsigned int filterBit = triggerRecords[irec].objectMap[iobj].second;

      for (pat::TriggerObjectStandAlone tobj : triggerObjects) {
	tobj.unpackPathNames(triggerNames);
	if (tobj.pathNames().size()==0) continue;
	tobj.unpackFilterLabels(iEvent, triggerResults);
	//std::cout << tobj.filterLabels().size() << std::endl;
	if (tobj.hasFilterLabel(dumb.at(0))) {
	  if (reco::deltaR(eta, phi, tobj.eta(), tobj.phi())<dRMax) {
	    matchBits[filterBit] = 1;
	  }
	}
	
      }
    }
  }

  
//  if (triggerObjects.at(0).pathNames().size()==0) return matchBits;
//  triggerObjects.at(0).unpackPathNames(triggerNames);
//  std::cout << "----" << std::endl;
//  std::cout << triggerObjects.at(0).pathNames().size() << std::endl;
//
//  for (pat::TriggerObjectStandAlone tobj : triggerObjects) {
//    std::cout << "QQ" << std::endl;
//    tobj.unpackFilterLabels(dumb);
//    std::cout << "QQQ" << std::endl;
//    std::cout << tobj.filterLabels().size() << std::endl;
//  }

  //for(unsigned int irec=0; irec<triggerRecords.size(); irec++) {     
  //  for(unsigned int iobj=0; iobj<triggerRecords[irec].objectMap.size(); iobj++) {
  //    const std::string   filterName = triggerRecords[irec].objectMap[iobj].first;
  //    const unsigned int  filterBit  = triggerRecords[irec].objectMap[iobj].second;
  //    std::vector<std::string> dumb;
  //    dumb.push_back(filterName);
  //    dumb.push_back("");
  //    std::cout << "hai" << std::endl;
  //    for(pat::TriggerObjectStandAlone tobj : triggerObjects) {
  //	tobj.unpackFilterLabels(dumb);
  //	std::cout << tobj.filterLabels().size() << std::endl;
  //	std::cout << tobj.filterLabels()[0] << std::endl;
  //      if(tobj.hasFilterLabel(filterName)) {
  //        if(reco::deltaR(eta,phi,tobj.eta(),tobj.phi()) < dRMax) {
  //          matchBits[filterBit] = 1; 
  //        }
  //      }
  //    }
  //  }
  //}

  return matchBits;
}

/*TriggerObjects TriggerTools::matchHLT(const double eta, const double phi,  
				      const trigger::Results &triggerResults,
                                      const pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  const double dRMax = 0.2;

  //TriggerObjects matchBits;

  for (pat::TriggerObjectStandAlone tobj : triggerObjects) {
    if(reco::deltaR(eta,phi,tobj.eta(),tobj.phi()) < dRMax) {    
      tobj.unpackFilterLabels(
    }
}

  for(unsigned int irec=0; irec<triggerRecords.size(); irec++) {     
    //for(unsigned int iobj=0; iobj<triggerRecords[irec].objectMap.size(); iobj++) {
    //const std::string   filterName = triggerRecords[irec].objectMap[iobj].first;
    //const unsigned int  filterBit  = triggerRecords[irec].objectMap[iobj].second;
      //std::cout << "hai" << std::endl;
      for(pat::TriggerObjectStandAlone tobj : triggerObjects) {
        if(tobj.hasFilterLabel(filterName)) {

            matchBits[filterBit] = 1; 
          }
        }
      }
    }
  }

  return matchBits;
}
*/
