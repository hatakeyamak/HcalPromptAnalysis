rm Output_Histo.root
root -b -q 'fillHisto_HcalDigis.C ( "TFileServiceOutputTree.root"," Output_Histo.root" )'
