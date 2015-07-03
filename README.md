# BaconProd
from ksung25, csa14 version


##To include custom MVA MET


  git-cms-merge-topic sabrandt:MVAForBaconv2
  
##To include PUPPI
  From instructions on PUPPI Twiki (https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI#Recipe_for_74X_developing)
  
  
  git cms-addpkg CommonTools/PileupAlgos
  
  
  git cms-merge-topic nhanvtran:puppi-etadep-741-v1
  
  
  
Use config file crab/makingBaconPuppiMVAMets_MC.py to include MVA and PUPPI calculations
