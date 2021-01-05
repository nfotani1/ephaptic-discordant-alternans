% Construct an simplified circuit describing how voltages
% delivered by firing the PIECEWISE-LINEAR pre-cleft 
% ion channel component
% creates the voltage needed to fire the post-cleft ion channel component.
%
% Derived from: oo_monod_2_abbrev_eph_model_equiv_circuit_1_v2.m

% Allow Matlab to see the definition of the circuit object in another
% directory:
addpath('circuit-simulation');

clear;

% ----Model parameters:----

% Default model parameters:
define_bidomain_regular_ephaptic_params;
define_EK_params;

% Make the following changes from the default:

% The following puts all the sodium channels on the cell surfaces
% facing the clefts:
total_cell_surf_area = 2*pi*p.r^2 + 2*pi*p.r*p.L;
p.g_fi_cleft = 0.5*total_cell_surf_area/(pi*p.r^2);
p.g_fi_tm = 0;

p.tau_si = 4.1; % 4.1 promotes discordant alteranns

p.Ccl = 0.1*p.Ccl; % Change the cleft capacitor capacitances

% Adjust the cleft width (which affects the resistor Rr)
p.w = 12*1.e-7;
p.Rr = p.rho_ext/(8*pi*p.w); %Kohms

% Adjust the coupling coefficient (which affects the resistor Rg):
coupling_fract_of_normal = 0.001; %
p.Rg = 395/coupling_fract_of_normal; % KOhms

%---End of model parameters---

% Numerical parameters:
Delta_x = p.L; %Spacing between the nodes
Delta_t = 0.0001; % Timestep size
Nt = 1000; % Number of timestepe
alpha = 1; % 0.5 = Trapezoidal method; 1 = Backward Euler
% Add these to the list parameters passed to nonlinear-element functions:
p.Dt = Delta_t; p.Dx = Delta_x;

% Define the circuit object.

% Version using V3 from the 2-abbreviated monodomain model;
% Use the exact, nonlinear model for I_cleft:
N_nodes = 4; % or 5, for the voltage_source + internal resistor version.
ground = N_nodes; % The last node is the ground node (since Matlab doesn't allow a node 0).
ckt = circuit(N_nodes,alpha,p.Dt);
ckt.add_resistor(2,ground,p.Rr);
ckt.add_capacitor(3,ground,p.Cm,-85.0);
ckt.add_capacitor(3,2,p.Ccl,-85.0);  % comment or uncomment
ckt.add_voltage_source(1,ground,@E0);
ckt.add_capacitor(1,2,p.Ccl,-85.0); % comment or uncomment

% This one or the next one should always be commented:

ckt.add_I_nonlin(1,2,@I_cleft_EK_piecewise_linear_2,@advance_EK,...
    @DI_cleft_EK_piecewise_linear_2,[1 1]);

% ckt.add_I_nonlin(1,2,@I_cleft_EK,@advance_EK,...
%     @DI_cleft_EK,[1 1]);
        
% Define the arrays to store the data for the small number of nodes
% and nonlinear currents we intend to study:
V_equiv_ckt = zeros(N_nodes,Nt);  % nodal voltages
t_ckt = zeros(1,Nt); % time

% Initialize the simulation:

ckt.prepare_matrices;
ckt.calc_initial_conditions(p,0,[-85;0;0;0;0;0]); % w/ the 2 Ccl caps.
% ckt.calc_initial_conditions(p,0,[-85;0;0;0;0]); %w/o the 2 Ccl caps.

% Main timestep loop:
for it = 1:Nt
    
    if (mod(it,1000)==0)
        fprintf('Timestep no. %i\n',it);
    end
    
    ckt.advance_circuit(p);
    % Record the nodal voltages each timestep:
    V_equiv_ckt(:,it) = ckt.V;
    t_ckt(it) = ckt.t;
    
end % End of timestep loop

plot(t_ckt,V_equiv_ckt(1:3,:)'); legend('V1','V2','V3');

rmpath('circuit-simulation');



