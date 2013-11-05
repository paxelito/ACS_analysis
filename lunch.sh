python ./acsStatesAnalysis.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1 2 1e-18;
python ./acsAttractorAnalysisInTime.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1 2 101;
python ./acsDynStatInTime.py -p ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1 -c 2;
python ./acsAttractorAnalysis.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1 2;
python ./acsSpeciesActivities.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1;