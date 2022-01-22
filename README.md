# Mu2e-II Snowmass 2021 sensitivity tools

Analysis tools for the Mu2e-II Snowmass 2021 sensitivity group

## Building

```
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
cd /mu2e/app/users/${USER}/
mkdir mu2eii
cd mu2eii
git clone --single-branch --branch Mu2eII_SM21 git@github.com:Mu2e/Offline.git
git clone git@github.com:Mu2e/TrkAna.git
git clone git@github.com:michaelmackenzie/Mu2eII_Snowmass.git mu2eii
# Some Stntuple classes currently have issues due to Offline advancing, use a version patched for now
git clone /mu2e/app/users/mmackenz/muse/Stntuple/ Stntuple
./mu2eii/scripts/build_config_muse
muse setup
muse build -j12
```

# Histogramming TrkAna trees

Histogramming is implemented using an object that takes in a TrkAna tree and writes out sets of histograms
of the event/track parameters with different event/track selections.

```
root.exe -q -b mu2eii/scripts/process_trkana.C
```

## Optimizing the signal window

Currently this runs using the SU2020 histograms with the code written by Pasha and imported from Mu2e/su2020, with the code
imported in the `stat` package.

This evaluates the mean discovery R for scans in the track p vs t space. The histogram directory is taken from
the .rootrc file, where these should be migrated from /mu2e/data/projects/su2020/hist/ to /mu2e/data/projects/mu2eii_snowmass/hist/.
This path can be updated to a new histogram path if needed.
```
root.exe mu2eii/mumem_sensitivity.C
root> mumem ana(13); //See stat/constants.cc to find the Mode options
root> ana.scan_pmin(104.1, 104.9, 680, 1700, 6);
root> ana.scan_tmin(104.0, 104.9, 750, 1700, 10);
```

## Evaluating the median expected limit and discovery

This uses the `calc` package, similarly imported from Mu2e/su2020.

This has been updated to run in a compiled mode, and a new version has been implemented to create the Poisson model from a datacard
that lists the sources of signals and backgrounds, the mean expected conttributions of each, and the systematic uncertainties on each,
where each uncertainty is either 100% or 0% correlated with the other sources.

Original script based model implementation:
```
root.exe -q -b mu2eii/Mu2eII_model.C
```

Datacard based model implementation (using the datacard `mu2eii/datacard.txt` with estimates using SU2020 histograms and Mu2e-II parameters):
```
root.exe -q -b mu2eii/model_builder.C
```

