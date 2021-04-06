% 4-node ephaptic model parameters
    
p.L = 100e-4; % cell length (cm, which is 10^4 microns
p.r = 11e-4; % cell radius (cm, which is 10^4 microns)
p.cm = 1.0; % membrane capacitance per unit area (uF/cm^2)
p.rho_myo = 0.150; % myoplasmic resistivity (Kohm-cm)
p.rho_ext = 0.150; % extracellular resistivity (Kohm-cm)
p.w = 8e-7; % 45e-7; % cleft width (cm, which is 10^7 nanometers)
p.V_rest = -85.0; % resting potential (mV)
p.V_fi = +15.0; % fully depolarized potential (mV)

% Note: p.Rm is split into two resistors, each half the value
% of the default value of p.Rg:
p.Rm = 0.5*p.rho_myo*p.L/(pi*p.r^2); % Kohms
% p.Re = 0.01*p.Rm; % NOT USED. Make the system effectively monodomain
p.Cm = p.cm*2*pi*p.r*p.L; % uF
p.Rg = 395.0; % Kohms

% Default ephaptic parameters:
p.Rcl = p.rho_ext*p.w/(pi*p.r^2); % Kohms (not present in the abbrev. model)
p.Rr = p.rho_ext/(8*pi*p.w); %Kohms
p.Ccl = p.cm*pi*p.r^2; % uF    


