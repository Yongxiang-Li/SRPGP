close all;   clear;   clc;
addpath('..\routine');
rng('default')

% Read signals from Excel
N = 6000;
raw_data = csvread('..\data\VIC.csv',1,0);
ele_signal = raw_data(481:N,3);
temp_data = csvread('..\data\Melbourne.csv',1,0);
temp = temp_data(481:N,3);

% Set parameters of decomposition
number = 350;
period = 300:number;
lob = [0.1 1];
upb = [3 10];

% d: degree of Bspline (1 for linear and 2 for quadratic) 
% p: degree freedom of basis matrix
bs.d = 3;    bs.p = 30;
lb = 0; ub = 1;
bs.knots = [lb*ones(1,bs.d) equalspace(lb,ub,bs.p-bs.d+1) ub*ones(1,bs.d)]';  % contrl point Px
bs.func = @(x)BsplineMatrix(bs.d, bs.p, bs.knots, x); % Bspline basis function
regr = @(x)temp_regr(temp,bs.func,x);

tic
GPFit = fit_SRPGP(period, ele_signal, regr, @period_sin_gauss_cov, lob, upb);
toc

signal1 = ele_signal - GPFit.Z;
period1 = 1:100;
GPFit1 = fit_SRPGP(period1, signal1, regr, @period_sin_gauss_cov, lob, upb);
remainder = signal1 - GPFit1.Z - GPFit1.H*GPFit1.alpha;

startDate = datetime(2000,1,11,0,0,0); 
dateArray = startDate + minutes(0:30:(5518*30)); 


figure;
% Original Data
ax1 = subplot(7,1,1); 
plot(ele_signal); 
title('(a) Original Data');
ylim([3000 7000]);axis tight; 
xticks(1:1440:length(ele_signal)); 
xticklabels(datestr(dateArray(1:1440:end), 'mmmm'));

% Trend
ax2 = subplot(7,1,2); 
plot(GPFit1.H(:,3:end)*GPFit1.alpha(3:end)); 
title('(b) Trend');
xlim([0 5519]); 
ylim([4200 5500]); 
xticks(1:1440:length(GPFit1.H(:,3:end)*GPFit1.alpha(3:end))); 
xticklabels(datestr(dateArray(1:1440:end), 'mmmm'));

% Weekly Feature
ax3 = subplot(7,1,3); 
plot(GPFit.Z); 
title('(c) Weekly Pattern'); 
xlim([0 2688]); 
ylim([-1000 1000]); 
xticks(1:720:length(GPFit.Z)); 
xticklabels(datestr(dateArray(1:720:end), 'mmm. dd'));

% Daily Feature
ax4 = subplot(7,1,4); 
plot(GPFit1.Z); 
title('(d) Daily Pattern'); 
xlim([0 768]); 
ylim([-440 300]); 
xticks(1:192:length(GPFit1.Z)); 
xticklabels(datestr(dateArray(1:192:end), 'mmm. dd'));

% Temperatures
ax5 = subplot(7,1,5);
plot(GPFit1.H(:,1).*GPFit1.alpha(1));
title('(e) Temperature');
xlim([0 5519]);
xticks(1:1440:length(GPFit1.H(:,1).*GPFit1.alpha(1)));
xticklabels(datestr(dateArray(1:1440:end), 'mmmm'));

% Squared Temperatures
ax6 = subplot(7,1,6);
plot(GPFit1.H(:,2).*GPFit1.alpha(2));
title('(f) Squared Temperature');
xlim([0 5519]);
xticks(1:1440:length(GPFit1.H(:,2).*GPFit1.alpha(2)));
xticklabels(datestr(dateArray(1:1440:end), 'mmmm'));

% Remainder
ax7 = subplot(7,1,7);
plot(remainder);
title('(g) Remainder');
ylim([-1000 1000]);
xlim([0 5519]);
xticks(1:1440:length(remainder));
xticklabels(datestr(dateArray(1:1440:end), 'mmmm'));

figure;
subplot(2,1,1);  plot(GPFit.Z);    title('Weekly Pattern'); xlim([0 1008]); ylim([-1000 1000]);xticks(1:144:length(GPFit.Z)); 
weekdays = {'Sun.', 'Mon.', 'Tue.', 'Wed.', 'Thu.', 'Fri.', 'Sat.'};
selectedDates = dateArray(1:144:end);
selectedWeekdays = weekdays(weekday(selectedDates));
xticklabels(selectedWeekdays);
subplot(2,1,2);  plot(GPFit1.Z);    title('Daily Pattern'); xlim([0 144]); ylim([-440 300]);xticks(1:12:length(GPFit1.Z)); xticklabels(datestr(dateArray(1:12:end), 'hh'))

figure;
weekdays = {'Tue. ', 'Wed. ', 'Thu. ', 'Fri. ', 'Sat. ','Sun. ', 'Mon. '};
tickNames = cellstr(datestr(dateArray(1:12:end),'HH:MM'));
t = tiledlayout(1,1);
ax1 = axes(t);
plot(GPFit.Z+GPFit1.Z)
xlim([0 336]); 
ylim([-1300 900]);
title('Daily + Weekly Pattern');

xticks(1:12:length(GPFit1.Z)); xticklabels(tickNames)
set(gca,'xaxislocation','top','yaxislocation','left');
set(gca,'Linewidth',2,'Fontsize',10,'box','on');
set(gca, 'XTickLabelRotation',40); 

ax2 = axes(t);
plot(GPFit.Z+GPFit1.Z)
xlim([0 336]); 
ylim([-1300 900]);

xticks(1:48:length(GPFit1.Z)); xticklabels(weekdays)
set(gca,'xaxislocation','bottom','yaxislocation','right');
set(gca,'yticklabel',[]);
% set(gca,'Linewidth',2,'Fontsize',10,'box','on');
