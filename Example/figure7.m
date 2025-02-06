close all;
clear;
addpath('..\routine');
rng('default')

sigma = 1;    theta = 5;    N = 5000;    s = 1;
P = sort([50  80  130]); % period of PGP signals
SNR = 0;    delta = 10^(-SNR/20);
lob = [0.3 3];   upb = [2 10];
period = (20:160)';
est_periods = zeros(s, length(P), 2); 

filename = sprintf('3components.xlsx');
folder_path = '..\data'; 
full_path = fullfile(folder_path, filename);

original = readtable(full_path);
signals = original.signals;
trend = original.trend;
seeasonal_elements = [original.period50, original.period80, original.period130];

for i = 1:s
    SRPGPmodels = [];    noisySignal = signals;   
    for j = 1 : length(P)
        SRPGPmodel = fit_SRPGP(period, noisySignal, @regpoly2, @period_sin_gauss_cov, lob, upb);
        noisySignal = noisySignal - SRPGPmodel.Z;
        SRPGPmodels = [SRPGPmodels; SRPGPmodel];
    end
    [est_periods(i,:,2), index] = match_periods([SRPGPmodels(:).period], P);
    SRPGPmodels = [];    noisySignal = signals;  
    for j = 1 : length(P)
        SRPGPmodel = fit_SRPGP(est_periods(i,j,2), noisySignal, @regpoly2, @period_sin_gauss_cov, lob, upb);
        noisySignal = noisySignal - SRPGPmodel.Z;
        SRPGPmodels = [SRPGPmodels; SRPGPmodel];
    end
    estimated_noise = noisySignal - SRPGPmodels(end).trend;
    estimated_signals = [SRPGPmodels(end).trend horzcat(SRPGPmodels.Z)];
    original_signals = [trend seeasonal_elements];
end

subplot(3,2,1);     
plot(signals, 'b-', 'DisplayName', 'Signal Data');  
title('Signal'); 
xlim([0 500]);  
subplot(3,2,2);     
plot(trend, 'b-', 'DisplayName', 'Original Trend');    
hold on;    
plot(SRPGPmodels(end).trend, 'r--', 'DisplayName', 'Estimated Trend');  
title('Trend');

subplot(3,2,3);     
plot(original_signals(:, 2), 'b-', 'DisplayName', 'Original Period A');    
hold on;    
plot(estimated_signals(:, 2), 'r--', 'DisplayName', 'Estimated Period A');  
title('Periodic Component A');
xlim([0 200]);  

subplot(3,2,4);     
plot(original_signals(:, 3), 'b-', 'DisplayName', 'Original Period B');    
hold on;    
plot(estimated_signals(:, 3), 'r--', 'DisplayName', 'Estimated Period B');  
title('Periodic Component B');
xlim([0 300]);   

subplot(3,2,5);     
plot(original_signals(:, 4), 'b-', 'DisplayName', 'Original Period C');    
hold on;    
plot(estimated_signals(:, 4), 'r--', 'DisplayName', 'Estimated Period C');  
title('Periodic Component C');
xlim([0 500]);   

subplot(3,2,6);     
plot(original.noise, 'b-', 'DisplayName', 'Original Data');    
hold on;    
plot(estimated_noise, 'r--', 'DisplayName', 'Estimated Signals using SSTD');  
title('Noise');
xlim([0 500]);   

lgd = legend('Location', 'southoutside', 'Orientation', 'horizontal');
lgd.Position = [0.25, 0.02, 0.5, 0.01];


function [period, order] = match_periods(est_period, real_period)
    bak_period = est_period;
    period = nan(size(real_period));
    for i = 1 : length(real_period)
        index = est_period==real_period(i);
        if sum(index)>0
            period(i) = est_period(index);
            est_period = est_period(~index);
        end
    end
    index = isnan(period);
    if sum(index)>0
        period(index) = sort(est_period);
    end
    [~, order] = ismember(period, bak_period);
end