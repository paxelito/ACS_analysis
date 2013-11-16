#nice python ./acsStatesAnalysis.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1/1e-15 2 1e-15;
nice python ./acsAttractorAnalysis.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1/1e-15 2;
nice python ./acsSpeciesActivities.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1/1e-15;
nice python ./acsAttractorAnalysisInTime.py ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1/1e-15 2 1001;
nice python ./acsDynStatInTime.py -p ~/Documents/SIMS_PROTOCELL/protocellMembrane/max6_CH1/1e-15 -c 2;