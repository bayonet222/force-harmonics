% Program for initial sizing of a PM Machine
% Casper Klop, 2021

% Design parameters
P = 10e6;                       % Nominal Power [W]
S = 10;                         % Rated speed [RPM]
N_ph = 3;                       % Number of phases
N_spp = 0.4;                    % Number of slots/pole/phase ?? paper gives 0.4, book only gives Nspp>1
N_m = 320;                      % Number of poles

N_sp = N_spp*N_m;               % Number of slots per phase
N_s = N_sp * N_ph;              % Number of slots
N_sm = N_spp*N_ph;              % Number of slots per pole
Temp = 140;                     % Operating temperature copper

% Geometric constraints
L = 1.7;                        % Axial length [m]
R_ro = 4.85;                    % Outside rotor radius [m]
R_so = 5.0;                     % Outside stator radius [m]
g = 0.010;                      % Air gap [m]

l_m = 0.03;                     % Magnet length [m]
w_s = 0.041;                     % Slot width [m]

% Material values
rho_20 = 1.7241e-8;         	% Resistivity copper [Omega.m]
beta_copper = 4.3e-3;           % Temp coefficient [/C]
rho = rho_20 * (1+beta_copper*(Temp-20));

mu_R = 1.05;                    % Magnetic recoil permeability
mu_0 = 4*pi*10^-7;              % Vaccuum permeability [H/m]
B_r = 1.2;                      % PM remanence [T]

rho_s = 7600;                   % Density of stator [kg/m3]
Gamma_s = 0.6;                  % Core loss density (Bmax, fe) [W/kg] ??

% Voltage output
V_DC = 100e3;                   % DC Output voltage of generator [VDC]
V_LL_rms = V_DC/1.35;           % Equivalent LL voltage, assuming diode bridge [Vrms]
V_ph_rms = V_LL_rms/sqrt(3);    % Equivalent phase voltage [Vrms]
V_coil = V_ph_rms/(N_sp/2);     % Voltage per coil [Vrms]
E_max = V_ph_rms;               % Phase back emf [Vrms]

% Additional assumptions
J_max = 3;                      % Max slot current density [A/mm2]
B_max = 1.5;                    % Max flux density [T]

k_st = 0.8;                     % Lamination stacking factor ?? 
k_cp = 0.76;                    % Conductor packing factor ??
alpha_m = 0.7;                  % Magnet fraction

pf = 0.85;                      % Power factor [-]

% Derived values
R_si = R_ro + g + l_m;          % Inside stator radius [m]

f_m = S/60;                     % Mechanical frequency [Hz]
omega_m = 2*pi*f_m;             % Mechanical rot. speed [rad/s]
omega_e = N_m/2 * omega_m;      % Electrical speed [rad/s]
f_e = omega_e/2/pi;             % Electrical frequency [Hz]

T = P/omega_m;                  % Equivalent torque

% Derived pitches
theta_p = 2*pi/N_m;             % Pole pitch [rad]
tau_p = R_si * theta_p;         % Pole pitch [m]

%alpha_cp = floor(N_spp)/N_spp; % Coil-pole fraction

theta_s = 2*pi/N_s;             % Slot pitch [rad]
tau_s = R_si * theta_s;         % Slot pitch (air gap) [m]

theta_c = 2*theta_s;            % Coil pitch, single-layered winding
%tau_c = alpha_cp * tau_p;      % Coil pitch
tau_c = theta_c * R_si;         % Coil pitch
alpha_cp = tau_c/tau_p;         % Coil-pitch fraction

% Magnetic calculations
C_phi = 2*alpha_m/(1+alpha_m);  % Flux concentration factor
P_c = l_m/g/C_phi;              % Permeance coefficient
g_c = g + l_m/mu_R;             % Effective air gap
k_c = 1/(1-1/(tau_s/w_s*(5*g_c/w_s+1)));    % Carter coefficient
k_ml = 1 + 4*l_m/pi/mu_R/alpha_m/tau_p* log(1 + pi*g/(tau_p*(1-alpha_m)));
% Magnet leakage factor

A_g = tau_p*L*(1+alpha_m)/2;                % Air gap area [m2]
B_g = C_phi*B_r/(1+mu_R*k_c*k_ml/P_c);      % Air gap flux density [T]. Average value, corresponding to COMSOL
%B_g = 0.85;
phi_g = B_g * A_g;                          % Air gap flux flowing to tooth [Wb]. 
phi_g_rms = phi_g/sqrt(2);                  % Flux rms

% We could define this as peak phi_g, divide sqrt(2) to get the rms
% and then get Vrms by multiplying with omega_e

% Geometric calculations
w_bi = phi_g/B_max/k_st/L;      % Back iron width [m]       CHANGE FROM HANSELMAN AS FLUX DOES NOT SPLIT EQUALLY
w_th = w_bi/N_sm;               % Tooth width [m]           SAME

R_sb = R_so - w_bi;             % Stator back iron radius [m]
w_sb = R_sb*theta_s - w_th;     % Slot bottom width [m]
w_si = R_si*theta_s - w_th;     % Slot inner width [m]

d_s = R_sb - R_si;              % Total slot depth [m]
A_s = d_s*(w_sb+w_si)/2;        % Conductor area [m2]
V_st = (pi*(R_so^2-R_si^2) - N_s*A_s)*L*k_st;   % Stator volume [m3]

% Electric calculations
k_d = 1;                        % Distribution factor (Concentrated) ??
gamma_s = pi * N_m/N_s;         % Slot pitch angle
epsilon = pi - gamma_s;         % Coil-span angle
k_e = cos(0.5 * epsilon);       % Chording factor 
k_w = k_d * k_e;                % Winding factor
k_s = 1;                        % SKEWED? (Neglected for now)

%n_s = floor(E_max/k_w/k_s/(phi_g/tau_p)/R_ro/N_sp/omega_m);     % Turns per slot
%e_max = N_m*k_w*k_s*(phi_g/tau_p)*R_ro*N_spp*n_s*omega_m;     % Actual peak emf

n_s = floor(E_max/k_w/k_s/phi_g_rms/omega_e/(N_sp/2));  % Turns per slot
e_max = N_sp/2*k_w*k_s*phi_g_rms*n_s*omega_e;         % Actual peak emf

i = P/3/e_max/pf;               % Coil current equals phase current [Arms] 
I_s = n_s * i;                  % Total slot current [Arms]
I_ph = i;                       % Phase current [A] (optimistic estimate)
J_c = I_s/k_cp/A_s/1e6;         % Current density [A/mm2] (rms)

A_wire = A_s * k_cp/n_s;        % Wire area

if J_c < J_max
    disp(['Current density OK: ', num2str(J_c),' A/mm2'])
else
    disp(['Current density too high: ', num2str(J_c),' A/mm2'])
end 

Z_Y = 2 * V_coil / i;           % Equivalent load in Y-connection

R_s = rho*n_s^2*L/k_cp/A_s;             % Slot resistance
R_e = rho*n_s^2*pi*tau_c/2/k_cp/A_s;    % End turn resistance
R_coil = 2 * (R_s + R_e);               % Coil resistance
R_seg = 4 * R_coil;                     % Segment (4 coils) resistance
R_ph = N_sp*(R_s+R_e);                  % Phase resistance
% AC effects?

g_e = k_c * g;                          % Effective air gap length

L_g = n_s^2*mu_R*mu_0*L*tau_p*k_d/2/(l_m+mu_R*g_e); % Air gap L/slot
L_s = n_s^2*mu_0*L*(d_s^2/3/A_s);                   % Slot L
L_e = n_s^2*mu_0*tau_c/16*log(tau_c^2*pi/2/A_s);    % End turn leakage
L_seg = 4 * (L_g + L_s + L_e);                      % Segment (4coils) ind
L_ph = N_sp * (L_g + L_s + L_e);                    % Phase inductance

% Performance calculations
P_r = N_ph * I_ph^2 * R_ph;     % Ohmic motor loss [W]
P_cl = rho_s * V_st * Gamma_s;  % Core loss [W]
P_s = 0.02 * P;                 % Stray loss assumption

eff = P/(P+P_r+P_cl+P_s);       % Efficiency

% For comparison with Solveig
K_s = A_s*J_c*1e6*N_s*k_w*k_cp/pi/(2*R_si);
T_r = sqrt(2)*pi/4 * B_g*K_s*4*R_si^2*L;

F_d = P / 2 / omega_m / pi / R_si^2 / L;            % Force density
if F_d > 60000
    disp(['Force density too high: ', num2str(F_d), ' N/m2'])
elseif F_d < 30000
    disp(['Force density too low: ', num2str(F_d), ' N/m2'])
else
    disp(['Force density OK: ', num2str(F_d), ' N/m2'])
end

if K_s * J_c > 200000
    disp(['KJ too high: ', num2str(K_s * J_c/100), ' A^2/mm^2/cm'])
else
    disp(['KJ OK: ', num2str(K_s * J_c/100), ' A^2/mm^2/cm'])
end

% Checking teeth flux density due to armature
B_th = 1.5 * mu_0 * I_s * sqrt(2) / (g + l_m);
disp(['Armature reaction in teeth: ', num2str(B_th), ' T (aim for ~.3 T)'])