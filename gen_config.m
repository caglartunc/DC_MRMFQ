config.params.am1 = 0.01; % MCG LB TP increase rate (10 sec avg)
config.params.rm1 = 0.01; % MCG LB TP decrease rate (5 sec avg) --> avg TP: (20+50)/3 = 23.333 Mbps
config.params.as1 = 1/100; % SCG HB TP increase rate 
config.params.rs1 = 1/100; % SCG HB TP decrease rate --> outage prob: 0.5

config.params.mu_m1 = 0.01; % avg rate: 10 Mbps
config.params.mu_m2 = 0.05; % avg rate: 50 Mbps
config.params.mu_s1 = 0; % avg rate: 0 Mbps
config.params.mu_s2 = 0.8; % avg rate: 500 Mbps

config.params.packet_drift = 0.1; % (s) drift of the buffer during packet arrival
config.params.ba = 0.1; % (1/Mb) 1/avg packet size arriving at PDCP: 1 Mb (normalized with config.params.packet_drift)
config.params.rate = 0.1; % (Mb/ms) data rate requirement in Mbps
config.params.ta = config.params.ba*config.params.rate/config.params.packet_drift; % (1/s) avg packet arrival rate at PDCP buffer, so the rate is: 8Mb*100/8 = 100 Mb

config.params.T0 = 0;
config.params.B = 100; % PDCP buffer size in Mb
config.params.Td = 0; % deactivation threshold for PDCP buffer
config.paracms.Ta = 7.5; % activation threshold for PDCP buffer
config.params.T = unique(sort([0 config.params.Td config.params.Ta 5 10:5:config.params.B])); % splitting the buffer into multiple regimes for numeric stability

config.params.delta = 0.01; % delta for plotting PDFs/CDFs
config.output_dir = ['DATA/scg_',num2str(config.params.mu_s2*1e3),'mbps/scg_',num2str(1/config.params.as1),'_',num2str(1/config.params.rs1),'/',num2str(config.params.rate*1000),'Mb/'];
config.file_name = 'data_'+strjoin(string([config.params.Td config.params.Ta config.params.B]),'_')+'.mat';

config.numRegimes = length(config.params.T)-1;
config.numStates = 12;

config.save_data = 0;
config.plot_data = 0;