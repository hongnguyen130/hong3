clc; clear; close all;


%% ===== Load du lieu =====
% Pareto tu Python (mean + std)
data30 = load('pareto_mean_P30.mat');   % co RE, SEE_mean, SEE_std
data40 = load('pareto_mean_P40.mat');
data50 = load('pareto_mean_P50.mat');

% AO tu MATLAB
ao = load('AO_points_vs_Pmax.mat');  
% co:
% RE_AO_mean (1x3)
% SEE_AO_mean (1x3)
% Pmax_set

%% ===== Too figure =====
figure('Position',[200 200 700 500]); hold on; grid on;

%% ===== Ve Pareto NSGA-II =====
% --- Pmax = 30 ---
plot(data30.RE, data30.SEE_mean, 'b-', 'LineWidth', 2);
fill([data30.RE; flipud(data30.RE)], ...
     [data30.SEE_mean - data30.SEE_std; flipud(data30.SEE_mean + data30.SEE_std)], ...
     'b', 'FaceAlpha', 0.15, 'EdgeColor','none');

% --- Pmax = 40 ---
plot(data40.RE, data40.SEE_mean, 'r-', 'LineWidth', 2);
fill([data40.RE; flipud(data40.RE)], ...
     [data40.SEE_mean - data40.SEE_std; flipud(data40.SEE_mean + data40.SEE_std)], ...
     'r', 'FaceAlpha', 0.15, 'EdgeColor','none');

% --- Pmax = 50 ---
plot(data50.RE, data50.SEE_mean, 'g-', 'LineWidth', 2);
fill([data50.RE; flipud(data50.RE)], ...
     [data50.SEE_mean - data50.SEE_std; flipud(data50.SEE_mean + data50.SEE_std)], ...
     'g', 'FaceAlpha', 0.15, 'EdgeColor','none');

%% ===== Ve điem AO =====
% Marker lon đe noi bat

plot(ao.RE_AO_mean(1), ao.SEE_AO_mean(1), 'bo', 'MarkerSize', 10, 'LineWidth',2);
plot(ao.RE_AO_mean(2), ao.SEE_AO_mean(2), 'ro', 'MarkerSize', 10, 'LineWidth',2);
plot(ao.RE_AO_mean(3), ao.SEE_AO_mean(3), 'go', 'MarkerSize', 10, 'LineWidth',2);

%% ===== Nhan & legend =====
xlabel('Eavesdropper Rate R_E','FontSize',12);
ylabel('Secrecy Energy Efficiency (SEE)','FontSize',12);

title('NSGA-II Pareto vs AO-Dinkelbach','FontSize',13);

legend(...
    'Pmax=30 (NSGA-II)','', ...
    'Pmax=40 (NSGA-II)','', ...
    'Pmax=50 (NSGA-II)','', ...
    'AO Pmax=30','AO Pmax=40','AO Pmax=50', ...
    'Location','best');

set(gca,'FontSize',11);

%% ===== Luu hinh =====
%saveas(gcf,'NSGAII_vs_AO.png');