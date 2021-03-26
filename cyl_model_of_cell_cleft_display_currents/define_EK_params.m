% Echebarria-Karma 3-variable model parameters

% To (nearly) duplicate Run #1 from 2018, uncomment Run #1
% lines in two places below. (More ephaptic case.)
% To exactly duplicate Run #2 from 2018, uncomment Run #2
% lines in two places below. (Less ephaptic case.)
fprintf('define_EK_params: Spatial distribution of I_Na channels is adjustable.\n');
p.tau_so = 15; 
p.tau_fi = 0.8;
% p.tau_si = 4.1; % more ephaptic: 2018 Run #1
p.tau_si = 4.0; % less ephaptic: 2018 Run #2

% Not sure what this comment means anymore:
% To obtain standing discordant alternans nodes, make these
% both 4.0 and change tau_h2 = 2.0 in advance q_EK:
% tau_si_main = 4.0; % typical values between 4 and 4.2
% tau_si_different = 4.0;

p.tau_h1 = 4.8;  p.tau_f1 = 100; p.tau_f2 = 30;
p.tau_h2 = 6; % typical values between 2 and 20

% Comment one of the following out:

% From 2018 Run #1:
% This puts all the sodium channels on the cell surfaces
% facing the clefts:
total_cell_surf_area = 2*pi*p.r^2 + 2*pi*p.r*p.L;
p.g_fi_cleft = 0.5*total_cell_surf_area/(pi*p.r^2);
p.g_fi_tm = 0;

% %From 2018 Run #2:
% This spreads the sodium channels uniformly over
% the surfaces of the cells:
% p.g_fi_cleft = 1;
% p.g_fi_tm = 1;

