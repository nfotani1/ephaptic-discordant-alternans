% Optimal strong ephaptic case:
%   1. All the sodium channels on the ends
%   2. Low cleft capacitance
%   3. Very low gap junction coupling
%   4. Optimal cleft width (12 nm)

% The following puts all the sodium channels on the cell surfaces
% facing the clefts:
total_cell_surf_area = 2*pi*p.r^2 + 2*pi*p.r*p.L;
p.g_fi_cleft = 0.5*total_cell_surf_area/(pi*p.r^2);
p.g_fi_tm = 0;

p.Ccl = 0.1*p.Ccl; % Change the cleft capacitor capacitances

p.Rg = 1000*395.0; % Very low gap junction coupling
% p.Rg = 395.0; % Normal gap junction coupling

% Adjust the cleft width (which affects the resistor Rr)
p.w = 12*1.e-7; % (12 nm) Strong ephaptic regime
% p.w = 40*1.e-7; % (40 nm) Weak ephaptic regime
p.Rr = p.rho_ext/(8*pi*p.w); %Kohms