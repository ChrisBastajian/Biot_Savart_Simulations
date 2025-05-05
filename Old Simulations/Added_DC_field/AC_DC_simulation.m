%%Initializing paramaters:
t = linspace(0, 1/(1000), 1000000);
f = 1./t;
beta = 1000; %beta = u0m*T % 
f_drive = 1000; % 1 kHz
w_drive = 2 * pi * f_drive;
Hd = 25 * (10^(-3)); %25 mT amplitude
cm = 1; %concentration * magnetic moment

H_AC = Hd * cos(w_drive .*t); %pure AC field without DC offset in x 
H_DC = 0; % in mT

%%This is the magnetic field applied to the sample:
H = H_AC + H_DC;
H_mag = ((H_DC.^2) + (H_AC.^2)).^(0.5); %the magnitude of the field (what is actually used)
beta_H_mag = beta .* H_mag; %this is what is used in langevin function

%plot(t, H_mag);

%%This is the magnetization M(t):
M_mag = cm * coth(beta_H_mag) - 1./(beta_H_mag); 

%%This is the numerical result of dM/dt
dMdt = gradient(M_mag)./gradient(t); 

%%This is the calculated dMdt by hand:
dM_mag_dt = cm * (-csch(beta .* sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2)).^2 + 1 ./ (beta.^2 .* (H_DC.^2 + (Hd .* cos(w_drive .* t)).^2))) .* (-beta .* Hd.^2 .* w_drive .* sin(w_drive .* t) .* cos(w_drive .* t) ./ sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2));% the time differential of the total magnitude of the magnetization

%%What we really want is dMzdt:
Hz_unit = H_AC ./H_mag; %the z component of the unit vector of the entire field with respect to time
Mz = cm * (coth(beta_H_mag) - 1./(beta_H_mag)) .*Hz_unit;

dMzdt_num = gradient(Mz)./ gradient(t);

%analytically solving, Mz = M_mag * Hz_unit in the z direction
dH_magdt = (-Hd.^2 .* w_drive .* cos(w_drive .*t) .* sin(w_drive .*t))./H_mag; %the time differential of the magnetic field's magnitude
dHzdt = -Hd .* w_drive .* sin(w_drive .*t); %the time differential of the z component of the field (d(H_AC)/dt)
dHz_unit_dt = (dHzdt .* H_mag - dH_magdt .* H_AC) ./ (H_mag .^2); %unit vector of H in z direction

dMzdt_calc = dM_mag_dt .* Hz_unit + dHz_unit_dt .* M_mag;

%%Plots:
figure;

% Plot 1: Magnetic Field
subplot(3,2,1);
plot(t, H, 'DisplayName', 'H field');
title('Magnetic Field Applied');
xlabel('Time');
ylabel('H (T)');
legend;

% Plot 2: Magnetization M(t)
subplot(3,2,2);
plot(t, M_mag, 'DisplayName', 'Magnetization M(t)');
title('Magnetization');
xlabel('Time');
ylabel('M');
legend;

% Plot 3: Numerical and Analytical dM/dt
subplot(3,2,3);
plot(t, dMdt, 'DisplayName', 'Numerical dM/dt');
hold on;
plot(t, dM_mag_dt, '--', 'DisplayName', 'Analytical dM/dt');
title('dM/dt Comparison');
xlabel('Time');
ylabel('dM/dt');
legend;
hold off;

% Plot 4: Numerical and Analytical dMz/dt
subplot(3,2,4);
plot(t, dMzdt_num, 'DisplayName', 'Numerical dMz/dt');
hold on;
plot(t, dMzdt_calc, '--', 'DisplayName', 'Analytical dMz/dt');
title('dMz/dt Comparison');
xlabel('Time');
ylabel('dMz/dt');
legend;
hold off;

% Plot 5: Magnetization Magnitude and z Component
subplot(3,2,5);
plot(t, M_mag, 'DisplayName', 'Magnetization Magnitude');
hold on;
plot(t, Mz, '--', 'DisplayName', 'Magnetization in z-direction');
title('Magnetization Comparison');
xlabel('Time');
ylabel('Magnetization');
legend;
hold off;

annotation('textbox', [0.7, 0, 0, 0.3], 'String', ...
    {['\beta = ', num2str(beta)], ...
     ['f_{drive} = ', num2str(f_drive)], ...
     ['Hd = ', num2str(Hd)], ...
     ['H_{DC} = ', num2str(H_DC)], ...
     ['c*m = ', num2str(cm)]}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');


%%Rewriting the full expression with no intermediate variable substitutions
%%then simplifying

dMzdt_calc_full = cm * ((-csch(beta .* sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2)).^2 + 1 ./ (beta.^2 .* (H_DC.^2 + (Hd .* cos(w_drive .* t)).^2))) .* ...
    (-beta .* Hd.^2 .* w_drive .* sin(w_drive .* t) .* cos(w_drive .* t) ./ sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2)) .* ...
    (Hd .* cos(w_drive .* t)) ./ sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2)) + ...
    ((-Hd .* w_drive .* sin(w_drive .* t) .* sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2) - ...
    ((-Hd.^2 .* w_drive .* cos(w_drive .* t) .* sin(w_drive .* t)) ./ sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2)) .* Hd .* cos(w_drive .* t)) ./ ...
    (H_DC.^2 + (Hd .* cos(w_drive .* t)).^2)) .* ...
    (cm * (coth(beta .* sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2)) - 1 ./ (beta .* sqrt(H_DC.^2 + (Hd .* cos(w_drive .* t)).^2))));

plot(t, dMzdt_calc_full, LineWidth= 2);
 %a  = dMzdt_calc_full ./dMzdt_calc %to check if they're equal

dMzdt_calc_simplified = cm * (-csch(beta_H_mag).^2 + 1 ./ (beta_H_mag.^2)) .* ...
    (-beta .* Hd.^2 .* w_drive .* sin(w_drive .* t) .* cos(w_drive .* t) ./ H_mag) .* Hz_unit + ...
    ((-Hd .* w_drive .* sin(w_drive .* t) .* H_mag - ((-Hd.^2 .* w_drive .* cos(w_drive .* t) .* sin(w_drive .* t)) ./ H_mag) .* Hd .* cos(w_drive .* t)) ./ (H_mag .^2)) .* ...
    (cm * (coth(beta_H_mag) - 1 ./ beta_H_mag));

b = dMzdt_calc_simplified ./dMzdt_calc_full;  %to check if they're equal

