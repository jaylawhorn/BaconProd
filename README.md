# BaconProd
from ksung25, csa14 version


##To include custom MVA MET


  git-cms-merge-topic sabrandt:MVAForBaconv2
  
##To include PUPPI
  From instructions on PUPPI Twiki (https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI#Recipe_for_74X_developing)
  
  
  git cms-addpkg CommonTools/PileupAlgos
  
  
  git cms-merge-topic nhanvtran:puppi-etadep-741-v1
  currently we are keeping puppi version in as a placeholder
  if checking out puppi gives too many cmssw packages, you can probably just copy these from me:  
  /afs/cern.ch/user/s/sabrandt/work/public/MVAMet/baconation/CMSSW_7_4_3_patch1/src/RecoJets
  /afs/cern.ch/user/s/sabrandt/work/public/MVAMet/baconation/CMSSW_7_4_3_patch1/src/CommonTools/
  
  
Use config file crab/makingBaconPuppiMVAMets_MC.py to include MVA and PUPPI calculations
