clear all; close all; clc;


% Specify the file name
filename = 'Figs. 7(a)-(b).xlsx';

%% Plot for "1D" sheet

data_1D = xlsread(filename, 'eta_vs_EB(m=1)');

% Extract x-axis, y-axis, and matrix values
y = data_1D(2:end, 1)/1e3;        % X-axis data from the first column (ignoring the first cell)
x = data_1D(1, 2:end);        % Y-axis data from the first row (ignoring the first cell)
values_1D = data_1D(2:end, 2:end); % Values from the remaining cells

% Create figure for Monoclastic
figure1 = figure('Name', '\eta vs m for 1D Bending');
colormap('jet');

% Create axes
axes1 = axes('Parent', figure1, 'FontSize', 18, 'FontName', 'Times New Roman');
view(axes1, [0.5 90]);
box(axes1, 'on');
hold(axes1, 'on');

% Plot the surface
surf(y, x, transpose(values_1D), 'Parent', axes1);
shading interp;

% Title and labels
title('RED', 'interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman', 'Rotation', 0);
xlabel('$\bar{E}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('$\eta$', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Set ticks and limits
xticks([0.2e-3 0.4e-3 0.6e-3 0.8e-3 1e-3]);
xlim([0.1e-3 1e-3]);
yticks([3 4 5 6 7 8 9 10]);
ylim([3 10]);

maxValue = max(values_1D(:));
minValue = min(values_1D(:));

clim([0 1]);
colorbar;
%cb = colorbar;
%cb.Ticks = linspace(minValue, maxValue, 5); 
%cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);







% Specify the file name
filename = 'Figs. 7(a)-(b).xlsx';

%% Plot for "1D" sheet

data_1D = xlsread(filename, 'eta_vs_EB(m=1.25)');

% Extract x-axis, y-axis, and matrix values
y = data_1D(2:end, 1)/1e3;        % X-axis data from the first column (ignoring the first cell)
x = data_1D(1, 2:end);        % Y-axis data from the first row (ignoring the first cell)
values_1D = data_1D(2:end, 2:end); % Values from the remaining cells

% Create figure for Monoclastic
figure1 = figure('Name', '\eta vs m for 1D Bending');
colormap('jet');

% Create axes
axes1 = axes('Parent', figure1, 'FontSize', 18, 'FontName', 'Times New Roman');
view(axes1, [0.5 90]);
box(axes1, 'on');
hold(axes1, 'on');

% Plot the surface
surf(y, x, transpose(values_1D), 'Parent', axes1);
shading interp;

% Title and labels
title('RED', 'interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman', 'Rotation', 0);
xlabel('$\bar{E}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('$\eta$', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Set ticks and limits
xticks([0.2e-3 0.4e-3 0.6e-3 0.8e-3 1e-3]);
xlim([0.1e-3 1e-3]);
yticks([3 4 5 6 7 8 9 10]);
ylim([3 10]);

maxValue = max(values_1D(:));
minValue = min(values_1D(:));

clim([0 1]);
colorbar;
%cb = colorbar;
%cb.Ticks = linspace(minValue, maxValue, 5); 
%cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);









% Specify the file name
filename = 'Fig. 7(c).xlsx';

%% Plot for "1D" sheet

data_1D = xlsread(filename, 'm = 1');

% Extract x-axis, y-axis, and matrix values
x = data_1D(2:end, 1);        % X-axis data from the first column (ignoring the first cell)
y = data_1D(1, 2:end);        % Y-axis data from the first row (ignoring the first cell)
values_1D = data_1D(2:end, 2:end)'; % Values from the remaining cells

% Create figure for Monoclastic
figure1 = figure('Name', '\eta vs \omega for 1D Bending');
colormap('jet');

% Create axes
axes1 = axes('Parent', figure1, 'FontSize', 18, 'FontName', 'Times New Roman');
view(axes1, [0.5 90]);
box(axes1, 'on');
hold(axes1, 'on');

% Plot the surface
surf(x, y, values_1D, 'Parent', axes1);
shading interp;

% Title and labels
title('RED', 'interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman', 'Rotation', 0);
xlabel('$\Omega / \Omega_n$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('$\eta$', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Set ticks and limits
%yticks([3 4 5 6 7 8 9 10]);
%ylim([3 10]);
%xticks([0.5 1 1.5 2]);
%ylim([0.5 1.5]);

maxValue = max(values_1D(:));
minValue = min(values_1D(:));

%clim([0 0.7]);
colorbar;
cb = colorbar;
cb.Ticks = linspace(minValue, maxValue, 5); 
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);





% Specify the file name
filename = 'Fig. 7(d).xlsx';

%% Plot for "1D" sheet

data_1D = xlsread(filename, 'm=1');

% Extract x-axis, y-axis, and matrix values
x = data_1D(2:end, 1);        % X-axis data from the first column (ignoring the first cell)
y = data_1D(1, 2:end);        % Y-axis data from the first row (ignoring the first cell)
values_1D = data_1D(2:end, 2:end)'; % Values from the remaining cells

% Create figure for Monoclastic
figure1 = figure('Name', '\eta vs \theta_0 for 1D Bending');
colormap('jet');

% Create axes
axes1 = axes('Parent', figure1, 'FontSize', 18, 'FontName', 'Times New Roman');
view(axes1, [0.5 90]);
box(axes1, 'on');
hold(axes1, 'on');

% Plot the surface
surf(x, y, values_1D, 'Parent', axes1);
shading interp;

% Title and labels
title('RED', 'interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman', 'Rotation', 0);
xlabel('$\Omega / \Omega_n$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('$\theta_0 (^o)$', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Set ticks and limits
%yticks([3 4 5 6 7 8 9 10]);
ylim([0 58.5]);
%xticks([0.5 1 1.5 2]);
%xlim([0.1 1]);

maxValue = max(values_1D(:));
minValue = min(values_1D(:));

%clim([0 0.7]);
colorbar;
cb = colorbar;
cb.Ticks = linspace(minValue, maxValue, 5); 
cb.TickLabels = arrayfun(@(x) sprintf('%.2f', x), cb.Ticks, 'UniformOutput', false);
