% Default ephaptic discretized cylindrical model parameters
    
% Geometric factors:
p.L = 100e-4; % cell length (cm, which is 10^4 microns
p.r = 11e-4; % cell radius (cm, which is 10^4 microns)
p.w = 8e-7; % 45e-7; % cleft width (cm, which is 10^7 nanometers)

% Intrinsic membrane capacitances:
p.cm = 1.0; % -lateral- membrane capacitance per unit area (uF/cm^2)
p.cm_cleft = p.cm; % cleft membrane capacitance / unit area (uF/cm^2)
p.Ccl_factor = 1.0; % change this factor to change cleft capacitances without
% other components in the cleft membrane (default: 1.0).

% Intrinsic resistances and conductances:
p.rho_myo = 0.150; % myoplasmic resistivity (Kohm-cm)
p.rho_cleft = 0.150; % cleft-region resistivity (kOhm-cm), from Kucera et al. 
% Note: not clear if Kucera's derivation of R_radial is valid.
p.re = 395/p.L; % resistance per unit length (Kohm/cm) in the region 
% outside the lateral membrane.  Note: I'm defining this to be equal to 
% the default gap junction resistance (395 kOhms) per cell length.
% Note #2: Now using p.re to characterize resistance outside the lateral 
% membrane instead of p.rho_cleft.
p.sigma_g = (1/395.0)/(pi*p.r^2);% Gap junction coupling conductance 
% per unit area (kOhms^(-1)*cm^(-2)). Default chosen so that the total
% gap junction resistance integrated over the cleft-facing portion of the
% cell = 395 kOhms.

% Membrane potential range:
p.V_rest = -85.0; % resting potential (mV)
p.V_fi = +15.0; % fully depolarized potential (mV)

