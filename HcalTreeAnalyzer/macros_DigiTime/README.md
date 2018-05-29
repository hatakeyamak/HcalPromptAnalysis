
# Example

## Making histograms

```
root.exe -b -q 'ana_DigiTiming.C+("list_trees_pi50_MCfull_CMSSW_10_2_0_pre3.txt",   "hcal_timestudy_pi50_histograms.root")'
root.exe -b -q 'ana_DigiTiming.C+("list_trees_minbias_MCfull_CMSSW_10_2_0_pre3.txt","hcal_timestudy_minbias_histograms.root")'
root.exe -b -q 'ana_DigiTiming.C+("list_trees_ttbar_MCfull_CMSSW_10_2_0_pre3.txt",  "hcal_timestudy_ttbar_histograms.root")'
root.exe -b -q 'ana_DigiTiming.C+("list_trees_2017Data_HLTPhysics_CMSSW_10_2_0_pre3.txt",  "hcal_timestudy_2017data_HLTPhysics_histograms.root")'
```

## Plotting histograms
```
root.exe
root> .x plot_digi_SOI.C
```
