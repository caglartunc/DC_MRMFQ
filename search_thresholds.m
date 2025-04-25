% Script to generate the search space for activation and deactivation thresholds

clear all
close all
clc

warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

B = 100;
granularity = 0.5;

rate_vec = [200];
scg_trans_vec = [50,100];
scg_data_rate_vec = [400,800];


td_vec = 0:granularity:10;
tic
caseIter = 1;
numCases = length(rate_vec) * length(scg_trans_vec) * length(scg_data_rate_vec);
h = waitbar(0,'Please wait while simulation is in progress...');
run gen_config.m
for rateIdx = 1:length(rate_vec)
    for transIdx = 1:length(scg_trans_vec)
        for scgIdx = 1:length(scg_data_rate_vec)
                waitbar(caseIter/numCases,h, ['Progress: ' num2str(caseIter/numCases*100) '%... Estimated time remaining: ' num2str(toc*numCases/caseIter/60) ' minutes.'])
            for tdIdx = 1:length(td_vec)
                Td = td_vec(tdIdx);

                ta_vec = Td:granularity:td_vec(end);

                for taIdx = 1:length(ta_vec)
                    Ta = ta_vec(taIdx);
                    config.params.rate = rate_vec(rateIdx)/1000;
                    config.params.ta = config.params.ba*config.params.rate/config.params.packet_drift; % (1/s) avg packet arrival rate at PDCP buffer, so the rate is: 8Mb*100/8 = 100 Mb
                    config.params.as1 = 1/scg_trans_vec(transIdx);
                    config.params.rs1 = 1/scg_trans_vec(transIdx);
                    config.params.mu_s2 = scg_data_rate_vec(scgIdx)/1000;

                    config.params.Td = Td; % deactivation threshold for PDCP buffer
                    config.params.Ta = Ta; % activation threshold for PDCP buffer
                    config.params.T = unique(sort([0 config.params.Td config.params.Ta 5 10:5:config.params.B]));
                    config.numRegimes = length(config.params.T)-1;
                    config.output_dir = ['DATA/scg_',num2str(config.params.mu_s2*1e3),'mbps/scg_',num2str(1/config.params.as1),'_',num2str(1/config.params.rs1),'/',num2str(config.params.rate*1000),'Mb/'];
                    config.file_name = 'data_'+strjoin(string([config.params.Td config.params.Ta config.params.B]),'_')+'.mat';
                    if exist(config.output_dir+config.file_name, 'file') == 0 % if file does not exist, run main. Else skip.
                        run main_test.m
                    else
                        disp("File " + config.output_dir+config.file_name + " exists.")
                    end
                end
            end
            caseIter = caseIter + 1;
            toc
        end
    end
end