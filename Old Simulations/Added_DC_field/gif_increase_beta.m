%% Parameters
t = linspace(0, 1/1000, 1000000); % Time vector
f_drive = 1000; % Driving frequency in Hz
w_drive = 2 * pi * f_drive; % Angular frequency
Hd = 25e-3; % 25 mT amplitude
beta = 1000; % DC field in mT
cm = 1; % Concentration * magnetic moment
H_AC = Hd * cos(w_drive .* t); % Pure AC field without DC offset

% Define the range of beta values
H_DCs = linspace(10^(-3), 40 * 10^(-3), 10);
num_DC = length(H_DCs);

% Create a folder to store images
if ~exist('temp_images', 'dir')
    mkdir('temp_images');
end

% Loop over beta values and create plots
for i = 1:num_DC
    H_DC = H_DCs(i);
    H = H_AC + H_DC;
    H_mag = sqrt(H_DC^2 + H_AC.^2); % Magnitude of the field
    beta_H_mag = beta .* H_mag; % Used in Langevin function
    
    % Magnetization M(t)
    M_mag = cm * coth(beta_H_mag) - 1 ./ beta_H_mag; 
    
    % Analytical solution for dMz/dt
    
    Hz_unit = H_AC ./ H_mag; % z-component of the unit vector of the field
    dH_magdt = (-Hd.^2 .* w_drive .* cos(w_drive .*t) .* sin(w_drive .*t))./H_mag; %the time differential of the magnetic field's magnitude
    dHzdt = -Hd .* w_drive .* sin(w_drive .*t); %the time differential of the z component of the field (d(H_AC)/dt)

    dHz_unit_dt = (dHzdt .* H_mag - dH_magdt .* H_AC) ./ (H_mag .^2); %unit vector of H in z direction

    Mz = cm .* (coth(beta_H_mag) - 1 ./ beta_H_mag) .* Hz_unit;
    dM_mag_dt = cm .* (-csch(beta .* sqrt(H_DC^2 + H_AC.^2)).^2 + 1 ./ (beta.^2 .* (H_DC^2 + H_AC.^2))) ...
        .* (-beta .* Hd.^2 .* w_drive .* sin(w_drive .* t) .* cos(w_drive .* t) ./ sqrt(H_DC^2 + H_AC.^2));
    dMzdt_calc = dM_mag_dt .* Hz_unit + dHz_unit_dt .* M_mag;
    
    % Plot dMz/dt
    fig = figure('Visible', 'off');
    plot(t, dMzdt_calc, 'DisplayName', sprintf('H_DC = %.2f', H_DC));
    title(sprintf('dMz/dt for H_DC = %.2f', H_DC));
    xlabel('Time');
    ylabel('dMz/dt');
    legend;
    xlim([0, max(t)]); % Adjust x-axis limits as needed
    ylim([-max(dMzdt_calc), max(dMzdt_calc)]); % Adjust y-axis limits as needed
    
    % Save plot as image
    filename = sprintf('temp_images/dMzdt_Beta_%.2f.png', H_DC);
    saveas(fig, filename);
    close(fig);
end

%% Create GIF
gif_filename = 'dMzdt_evolution.gif';

for i = 1:num_DC
    H_DC = H_DCs(i);
    filename = sprintf('temp_images/dMzdt_Beta_%.2f.png', H_DC);
    
    % Read image
    img = imread(filename);
    
    % Convert image to grayscale if necessary
    if size(img, 3) == 3 % RGB image
        img = rgb2gray(img); % Convert RGB to grayscale
    end
    
    % Write to GIF
    if i == 1
        imwrite(img, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
    else
        imwrite(img, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

% Clean up
rmdir('temp_images', 's');
