#! /bin/sh

echo 'Processing SingleMuPt1'
cmsRun CAfitterMu1_cfg.py
mv histotrip.root file_ok_MuPt1.root

echo 'Processing SingleMuPt10'
cmsRun CAfitterMu10_cfg.py
mv histotrip.root file_ok_MuPt10.root

echo 'Processing SingleMuPt100'
cmsRun CAfitterMu100_cfg.py
mv histotrip.root file_ok_MuPt100.root

echo 'Processing SingleMuPt1000'
cmsRun CAfitterMu1000_cfg.py
mv histotrip.root file_ok_MuPt1000.root

echo 'Processing TTbar'
cmsRun CAfitterTTbar_cfg.py
mv histotrip.root file_ok_ttbar.root
