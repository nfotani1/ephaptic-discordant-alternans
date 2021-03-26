% Cylindrical, spatially-detailed model of the role of ephaptic 
% coupling in action potential propagation.
% Each cell, the cleft region between cells, and the extracellular
% space is represented by several grid points.
% This version of the code creates colormap movies of the potential
% at each grid point, and also plots the currents flowing in
% resistors.

tic;
% clear;

addpath('../circuit-simulation');

% Define default model parameters:
define_default_ephaptic_discretized_cyl_params;
define_EK_params;

% Default numerical parameters:
N_cells = 30; % 10; % Number of cells
Nt = 10000; %100000; % Number of timesteps to run
Delta_t = 1.e-4; % Timestep size
it_plot = 500; % Plotting interval
alpha = 0.5; % 0.5 = Trapezoidal method; 1 = Backward Euler

% --Make the following changes from the default--

p.tau_si = 4.1; % 4.1 promotes discordant alteranns

% strong_ephaptic_case_param_changes;
% weak_ephaptic_low_velocity_case_param_changes;
% weak_ephaptic_high_velocity_case_param_changes;
change_default_params;
movie_name = sprintf('%s_%s_%s_%s',s_INa,s_Ccl,s_Rg,s_w);

% Each of the cells, numbered 0 to N_cells-1, contains
% nr interior nodes in the r direction, and
% nz interior nodes in the z direction.
% There are also additional nodes just inside the membrane
% in both the r and z directions.
nr = 4; % number of interior nodes inside the cell in the radial direction
nz = 5; % number of interior nodes inside the cell in the axial direction
Delta_z = p.L/nz; %Spacing between the nodes in the z-direction.
Delta_r = p.r/(nr-0.5); % Spacing of the nodes in the r-direction.

% Add these to the list parameters passed to nonlinear-element functions:
p.Dt = Delta_t; p.Dz = Delta_z; p.Dr = Delta_r;
p.Dx = p.Dz; % Dx in the I_membrane_EK and DI_membrane_EK routines is Dz here.

% --End of changes to default parameters--

% Check numerical stability criterion for the diffusion term:
% Fix for cylindrical code:
% fprintf('(Delta t)/((Rm+Re)*Cm) = %f\n',Delta_t/((p.Rm+p.Re)*p.Cm));

% Define the circuit object.

% Total number of nodes associated with nodes 0 to N_cells-2:
% nr*nz "interior" nodes,
% nz nodes just inside the lateral boundary,
% nz nodes just ouside the lateral boundary, in the extracellular space
% nr nodes just inside the left boundary,
% nr nodes just outside the left boundary in the cleft region,
% 1 node at the intersection of the left cleft region and extracellular space,
% and nr nodes just inside the right boundary.
N_nodes_per_cell = nr*nz + 2*nz + 2*nr + 1 + nr;
% Note, the nr nodes just outside the right boundary are considered
% part of the next cell, except for the rightmost cell (i.e., cell
% no. N_cells-1), when these nodes are considered part of this cell.
% Thus, Cell no. N_cells-1 has N_nodes_per_cell + nr + 1 nodes.

% So, the total number of cells is,
N_nodes = N_nodes_per_cell*N_cells + nr + 1;

ckt = circuit(N_nodes,alpha,p.Dt);
% All the cells in the 0th cell are numbered first, then
% all the cells in the 1st cell are numbered, etc.
%
% Within each cell, the node numbering order is;
% The interior nodes are first, with nr being the "fast" index,
% The nz nodes just inside the lateral boundary are next,
% The nz nodes just outside the lateral boundary are next,
% The nr nodes just inside the left boundary are next,
% The nr nodes just outside the left boundary are next,
% The 1 node at the intersection of the nodes just outside
% the left boundary and the nodes just outside the lateral
% boundary are next,
% The nr nodes just inside the right boundary are next,
% and, if this the last cell,
% The nr nodes just outside the right boundary are last,

% Define an arrays of function handles for the cleft membrane
% nonlinear currents.  (Arrays are needed, because the currents
% depend on the value of r.)
I_cleft_EK_for_ir = cell(nr,1);
DI_cleft_EK_for_ir = cell(nr,1);
for ir = 1:nr
    if (ir==1)
        area = pi*(0.5*p.Dr)^2;
    else
        area = pi*((ir-0.5)*p.Dr)^2 - pi*((ir-1.5)*p.Dr)^2;
    end
    % area
    I_cleft_EK_for_ir{ir} = @(V,q,p) J_cleft_EK(V,q,p)*area;
    DI_cleft_EK_for_ir{ir} = @(V,q,p) DJ_cleft_EK(V,q,p)*area;
end

% Add all the components to the circuit, one set of nodes
% at a time:

diagnostics = false;

if (diagnostics)
    figure(1); clf;
    plot_cell_boundaries_v2 (nr,nz,p.Dr,p.Dz,p.w,N_cells);
end

% Set up a table containing information needed to calculate
% branch currents later:
current_table = cell(0,4);
nct = 0; % initialize number of currents to plot

for cell_no = 0:(N_cells-1)
    
    if (cell_no==0)
        cap_voltage = -10.0;
    else
        cap_voltage = -85.0;
    end
    
    nodes_added_so_far = 0; % Keep track of the number of nodes wired so far in this cell.
    
    node = @(i) (cell_no*N_nodes_per_cell + i);
    
    % Wire up the interior nodes inside the cell:
    % Intracellular interior radial resistors:
    for iz = 1:nz
        for ir = 1:(nr-1)
            area = 2*pi*(ir-0.5)*p.Dr*p.Dz; length = p.Dr;
            ckt.add_resistor(node((iz-1)*nr+ir),node((iz-1)*nr+ir+1),p.rho_myo*length/area);
            if (diagnostics)
                comment = sprintf('\\rho_{myo} = %f, area = 2\\pi%2.1f\\Deltar\\Deltaz, length = \\Deltar',...
                    p.rho_myo,(ir-0.5));
                plot_branches (node((iz-1)*nr+ir),node((iz-1)*nr+ir+1),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
            end
            nct = nct+1;
            current_table(nct,:) = {'resistor',node((iz-1)*nr+ir),node((iz-1)*nr+ir+1),p.rho_myo*length/area};
        end
    end
    % Intracellular interior axial resistors:
    for iz = 1:(nz-1)
        for ir = 1:nr
            if (ir==1)
                area = pi*(0.5*p.Dr)^2; s_area = sprintf('\\pi(0.5\\Deltar)^2');
            else
                area = pi*((ir-0.5)*p.Dr)^2 - pi*((ir-1.5)*p.Dr)^2;
                s_area = sprintf('\\pi(%2.1f\\Deltar)^2-\\pi(%2.1f\\Deltar)^2',ir-0.5,ir-1.5);
            end
            length = p.Dz;
            ckt.add_resistor(node((iz-1)*nr+ir),node(iz*nr+ir),p.rho_myo*length/area);
            if (diagnostics)
                comment = sprintf('\\rho_{myo} = %f, area = %s, length = \\Deltaz',...
                    p.rho_myo,s_area);
                plot_branches (node((iz-1)*nr+ir),node(iz*nr+ir),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
            end
            nct = nct+1;
            current_table(nct,:) = {'resistor',node((iz-1)*nr+ir),node(iz*nr+ir),p.rho_myo*length/area};
        end
    end
    nodes_added_so_far = nodes_added_so_far + nr*nz;
    
    % Resistor-capacitor-ion channel combinations connected to the
    % nodes just inside the lateral membrane.
    % Note:
    % Nodes just inside the lateral membrane are numbered starting at lat_mem_in_base+1.
    % Nodes just outside the lateral membrane are numbered starting at lat_mem_out_base+1.
    lat_mem_in_base = nodes_added_so_far;
    lat_mem_out_base = nodes_added_so_far + nz;
    for iz = 1:nz
        area = 2*pi*(nr-0.75)*p.Dr*p.Dz; length = 0.5*p.Dr;
        ckt.add_resistor(node(iz*nr),node(lat_mem_in_base+iz),p.rho_myo*length/area);
        if (diagnostics)
            comment = sprintf('\\rho_{myo} = %f, area = 2\\pi%2.2f\\Deltar\\Deltaz, length = 0.5\\Deltar',...
                p.rho_myo,nr-0.75);
            plot_branches (node(iz*nr),node(lat_mem_in_base+iz),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
            nct = nct+1;
            current_table(nct,:) = {'resistor',node(iz*nr),node(lat_mem_in_base+iz),p.rho_myo*length/area};
        area = 2*pi*(nr-0.5)*p.Dr*p.Dz;
        ckt.add_capacitor(node(lat_mem_in_base+iz),node(lat_mem_out_base+iz),p.cm*area,cap_voltage);
        ckt.add_I_nonlin(node(lat_mem_in_base+iz),node(lat_mem_out_base+iz),...
            @I_membrane_EK,@advance_EK,@DI_membrane_EK,[1 1]);
        if (diagnostics)
            comment = sprintf('cm = %f, I_{mem}, area = 2\\pi%2.1f\\Deltar\\Deltaz',p.cm,nr-0.5);
            plot_branches (node(lat_mem_in_base+iz),node(lat_mem_out_base+iz),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
    end
    
    % Connect extracellular resistors adjacent lateral membrane:
    % Define the extracellular resistance to be the default
    for iz = 1:(nz-1)
        length = p.Dz;
        ckt.add_resistor(node(lat_mem_out_base+iz),node(lat_mem_out_base+iz+1),p.re*length);
        if (diagnostics)
            comment = sprintf('re = %f, length = \\Deltaz',p.re);
            plot_branches (node(lat_mem_out_base+iz),node(lat_mem_out_base+iz+1),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
        nct = nct+1;
        current_table(nct,:) = {'resistor',node(lat_mem_out_base+iz),node(lat_mem_out_base+iz+1),p.re*length};
    end
    
    nodes_added_so_far = nodes_added_so_far + 2*nz;
    
    % Resistor-capacitor-ion channel combinations connected to the
    % nodes just inside the cleft-facing membrane on the left side of each
    % cell.
    % Note:
    % Nodes just inside the cleft membrane, are numbered
    %   left_cleft_mem_in_base+1 to left_cleft_mem_in_base+nr.
    % Nodes just outside the cleft membrane, are numbered
    %   left_cleft_mem_out_base+1 to left_cleft_mem_in_base+nr.
    left_cleft_mem_in_base = nodes_added_so_far;
    left_cleft_mem_out_base = nodes_added_so_far + nr;
    for ir = 1:nr
        if (ir==1)
            area = pi*(0.5*p.Dr)^2; s_area = sprintf('\\pi(0.5\\Deltar)^2');
        else
            area = pi*((ir-0.5)*p.Dr)^2 - pi*((ir-1.5)*p.Dr)^2;
            s_area = sprintf('\\pi(%2.1f\\Deltar)^2 - \\pi(%2.1f\\Deltar)^2',ir-0.5,ir-1.5);
        end
        length = 0.5*p.Dz; s_length = sprintf('0.5\\Deltaz');
        ckt.add_resistor(node(ir),node(left_cleft_mem_in_base+ir),p.rho_myo*length/area);
        if (diagnostics)
            comment = sprintf('\\rho_{myo} = %f, area = %s, length = %s',p.rho_myo,s_area,s_length);
            plot_branches (node(ir),node(left_cleft_mem_in_base+ir),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
        nct = nct+1;
        current_table(nct,:) = {'resistor',node(ir),node(left_cleft_mem_in_base+ir),p.rho_myo*length/area};
        %         if (cell_no==0)
        %             init_cap_voltage = 0.0;
        %         else
        %             init_cap_voltage = -85.0;
        %         end
        ckt.add_capacitor(node(left_cleft_mem_in_base+ir),node(left_cleft_mem_out_base+ir),...
            p.Ccl_factor*p.cm_cleft*area,cap_voltage);
        ckt.add_I_nonlin(node(left_cleft_mem_in_base+ir),node(left_cleft_mem_out_base+ir),...
            I_cleft_EK_for_ir{ir},@advance_EK,DI_cleft_EK_for_ir{ir},[1 1]);
        if (diagnostics)
            comment = sprintf('cm_{cleft} = %f, I_{cleft}(%i), area = %s',p.Ccl_factor*p.cm_cleft,ir,s_area);
            plot_branches (node(left_cleft_mem_in_base+ir),node(left_cleft_mem_out_base+ir),...
                nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
    end
    
    % Connect resistors in the cleft region to the left of the cell:
    for ir = 1:(nr-1)
        area = 2*pi*(ir-0.5)*p.Dr*p.w; s_area = sprintf('2\\pi%2.1f\\Deltar w',ir-0.5);
        length = p.Dr; s_length = sprintf('\\Deltar');
        ckt.add_resistor(node(left_cleft_mem_out_base+ir),node(left_cleft_mem_out_base+ir+1),...
            p.rho_cleft*length/area);
        if (diagnostics)
            comment = sprintf('\\rho_{cleft} = %f, area = %s, length = %s',p.rho_cleft,s_area,s_length);
            plot_branches (node(left_cleft_mem_out_base+ir),node(left_cleft_mem_out_base+ir+1),...
                nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
        nct = nct+1;
        current_table(nct,:) = {'resistor',node(left_cleft_mem_out_base+ir),node(left_cleft_mem_out_base+ir+1),...
            p.rho_cleft*length/area};
        
    end
    
    nodes_added_so_far = nodes_added_so_far + 2*nr;
    
    % Connect up the node that lies at the intersection of the
    % left cleft region and the extracellular space:
    left_intersect_node = nodes_added_so_far + 1;
    area = 2*pi*(nr-0.75)*p.Dr*p.w; s_area = sprintf('2\\pi%2.2f\\Deltar w',nr-0.75);
    length = 0.5*p.Dr; s_length = sprintf('0.5\\Deltar');
    ckt.add_resistor(node(left_cleft_mem_out_base+nr),node(left_intersect_node),...
        p.rho_cleft*length/area);
    if (diagnostics)
        comment = sprintf('\\rho_{cleft} = %f, area = %s, length = %s',p.rho_cleft,s_area,s_length);
        plot_branches (node(left_cleft_mem_out_base+nr),node(left_intersect_node),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
    end
    nct = nct+1;
    current_table(nct,:) = {'resistor',node(left_cleft_mem_out_base+nr),node(left_intersect_node),...
        p.rho_cleft*length/area};
    length = 0.5*p.Dz + 0.5*p.w; s_length = sprintf('0.5\\Deltaz + 0.5w');
    ckt.add_resistor(node(left_intersect_node),node(lat_mem_out_base+1),p.re*length);
    if (diagnostics)
        comment = sprintf('re = %f, length = %s',p.re,s_length);
        plot_branches (node(left_intersect_node),node(lat_mem_out_base+1),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
    end
    nct = nct+1;
    current_table(nct,:) = {'resistor',node(left_intersect_node),node(lat_mem_out_base+1),p.re*length};
    
    nodes_added_so_far = nodes_added_so_far + 1;
    
    % Resistor-capacitor-ion channel combinations connected to the
    % nodes just inside the cleft-facing membrane on the right side of each
    % cell.
    % Note:
    % Nodes just inside the cleft membrane, are numbered
    %   right_cleft_mem_in_base+1 to right_cleft_mem_in_base+nr.
    % Nodes just outside the cleft membrane, are numbered
    %   right_cleft_mem_out_base+1 to right_cleft_mem_in_base+nr.
    right_cleft_mem_in_base = nodes_added_so_far;
    % Nodes in the cleft region to the right of a given cell belong to the
    % next cell; except if this is the last cell, then consider these nodes
    % to be associated with this cell:
    if (cell_no==N_cells-1)
        right_cleft_mem_out_base = nodes_added_so_far + nr;
    else
        right_cleft_mem_out_base = N_nodes_per_cell + left_cleft_mem_out_base;
    end
    for ir = 1:nr
        if (ir==1)
            area = pi*(0.5*p.Dr)^2; s_area = sprintf('\\pi(0.5\\Deltar)^2');
        else
            area = pi*((ir-0.5)*p.Dr)^2 - pi*((ir-1.5)*p.Dr)^2;
            s_area = sprintf('\\pi(%2.1f\\Deltar)^2-\\pi(%2.1f\\Deltar)^2',ir-0.5,ir-1.5);
        end
        length = 0.5*p.Dz; s_length = sprintf('0.5\\Deltaz');
        ckt.add_resistor(node(nr*(nz-1)+ir),node(right_cleft_mem_in_base+ir),p.rho_myo*length/area);
        if (diagnostics)
            comment = sprintf('\\rho_{myo} = %f, area = %s, length = %s',p.rho_myo,s_area,s_length);
            plot_branches (node(nr*(nz-1)+ir),node(right_cleft_mem_in_base+ir),nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
        nct = nct+1;
        current_table(nct,:) = {'resistor',node(nr*(nz-1)+ir),node(right_cleft_mem_in_base+ir),p.rho_myo*length/area};
        ckt.add_capacitor(node(right_cleft_mem_in_base+ir),node(right_cleft_mem_out_base+ir),...
            p.Ccl_factor*p.cm_cleft*area,cap_voltage);
        ckt.add_I_nonlin(node(right_cleft_mem_in_base+ir),node(right_cleft_mem_out_base+ir),...
            I_cleft_EK_for_ir{ir},@advance_EK,DI_cleft_EK_for_ir{ir},[1 1]);
        if (diagnostics)
            comment = sprintf('cm_{cleft} = %f, I_{cleft}(%i), area = %s',p.Ccl_factor*p.cm_cleft,ir,s_area);
            plot_branches (node(right_cleft_mem_in_base+ir),node(right_cleft_mem_out_base+ir),...
                nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
    end
    
    % If this is the last cell, connect resistors in the cleft region
    % to the right of the cell.  (Otherwise, it will be done as part of
    % of connecting up the next cell.)
    if (cell_no==N_cells-1)
        for ir = 1:(nr-1)
            area = 2*pi*(ir-0.5)*p.Dr*p.w;
            s_area = sprintf('2\\pi%2.1f\\Deltar w',ir-0.5);
            length = p.Dr; s_length = sprintf('\\Deltar');
            ckt.add_resistor(node(right_cleft_mem_out_base+ir),node(right_cleft_mem_out_base+ir+1),...
                p.rho_cleft*length/area);
            if (diagnostics)
                comment = sprintf('\\rho_{cleft} = %f, area = %s, length = %s',p.rho_cleft,s_area,s_length);
                plot_branches (node(right_cleft_mem_out_base+ir),node(right_cleft_mem_out_base+ir+1),...
                    nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
            end
            nct = nct+1;
            current_table(nct,:) = {'resistor',node(right_cleft_mem_out_base+ir),node(right_cleft_mem_out_base+ir+1),...
                p.rho_cleft*length/area};
        end
    end
    
    % Connect up the gap junction resistors to the right of each cell,
    % unless this is the last cell:
    if (cell_no~=N_cells-1)
        for ir = 1:nr
            if (ir==1)
                area = pi*(0.5*p.Dr)^2; s_area = sprintf('\\pi(0.5\\Deltar)^2');
            else
                area = pi*((ir-0.5)*p.Dr)^2 - pi*((ir-1.5)*p.Dr)^2;
                s_area = sprintf('\\pi(%2.1f\\Deltar)^2 - \\pi(%2.1f\\Deltar)^2',ir-0.5,ir-1.5);
            end
            ckt.add_resistor(node(right_cleft_mem_in_base+ir),...
                node(N_nodes_per_cell+left_cleft_mem_in_base+ir),1/(p.sigma_g*area));
            if (diagnostics)
                comment = sprintf('\\sigma_g = %f, area = %s',p.sigma_g,s_area);
                plot_branches (node(right_cleft_mem_in_base+ir),...
                    node(N_nodes_per_cell+left_cleft_mem_in_base+ir),...
                    nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
            end
            nct = nct+1;
            current_table(nct,:) = {'resistor',node(right_cleft_mem_in_base+ir),...
                node(N_nodes_per_cell+left_cleft_mem_in_base+ir),1/(p.sigma_g*area)};
        end
    end
    
    if (cell_no==N_cells-1)
        nodes_added_so_far = nodes_added_so_far + 2*nr;
    else
        nodes_added_so_far = nodes_added_so_far + nr;
    end
    
    % Connect up the node that lies at the intersection of the
    % right cleft region and the extracellular space.
    % Define this node if this is the last cell. Otherwise
    % use its future definition in the next cell:
    if (cell_no==N_cells-1)
        right_intersect_node = nodes_added_so_far + 1;
    else
        right_intersect_node = N_nodes_per_cell + left_intersect_node;
    end
    if (cell_no==N_cells-1)
        area = 2*pi*(nr-0.75)*p.Dr*p.w;
        s_area = sprintf('2\\pi%2.2f\\Deltar w',nr-0.75);
        length = 0.5*p.Dr; s_length = sprintf('0.5\\Deltar');
        ckt.add_resistor(node(right_cleft_mem_out_base+nr),node(right_intersect_node),...
            p.rho_cleft*length/area);
        if (diagnostics)
            comment = sprintf('\\rho_{cleft} = %f, area = %s, length = %s',p.rho_cleft,s_area,s_length);
            plot_branches (node(right_cleft_mem_out_base+nr),node(right_intersect_node),...
                nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
        end
        nct = nct+1;
        current_table(nct,:) = {'resistor',node(right_cleft_mem_out_base+nr),node(right_intersect_node),...
            p.rho_cleft*length/area};
    end
    length = 0.5*p.Dz + 0.5*p.w; s_length = sprintf('0.5\\Deltaz + 0.5w');
    ckt.add_resistor(node(right_intersect_node),node(lat_mem_out_base+nz),p.re*length);
    if (diagnostics)
        comment = sprintf('re = %f, length = %s',p.re,s_length);
        plot_branches (node(right_intersect_node),node(lat_mem_out_base+nz),...
            nr,nz,p.Dr,p.Dz,p.w,N_cells,comment);
    end
    nct = nct + 1;
    current_table(nct,:) = {'resistor',node(right_intersect_node),node(lat_mem_out_base+nz),p.re*length};
    if (cell_no==N_cells-1)
        nodes_added_so_far = nodes_added_so_far + 1;
    end
    
    N_nodes_in_this_cell = nodes_added_so_far;
    
    % Check:
    if (cell_no~=N_cells-1)
        if (N_nodes_per_cell~=N_nodes_in_this_cell)
            disp('N_nodes_per_cell does not match N_nodes_in_this_cell');
        end
    end
    
end

hold on;
plot_cell_boundaries_v2 (nr,nz,p.Dr,p.Dz,p.w,N_cells);
hold off;

ckt.prepare_matrices;

% internal_node_guess = [-85,-85,-85,-85,-85,-85,0,0];
% V_augmented_guess = [repmat(internal_node_guess,1,N_cells-1),-85,...
%     zeros(1,3*N_cells-2)]';
% ckt.calc_initial_conditions(p,0.0,V_augmented_guess);
ckt.calc_initial_conditions(p,0.0);

vidObj = VideoWriter(movie_name,'MPEG-4');
open(vidObj);

% Main timestep loop:
for it = 1:Nt
    
    if (mod(it-1,20)==0)
        figure(2);
        hold off;
        colorplot_of_array (ckt.V,nr,nz,p.Dr,p.Dz,p.w,N_cells);
        hold on;
        
        % plot_cell_boundaries_v2 (nr,nz,p.Dr,p.Dz,p.w,N_cells); hold on;
        r_from = zeros(1,nct); r_to = zeros(1,nct);
        z_from = zeros(1,nct); z_to = zeros(1,nct);
        current = zeros(1,nct);
        for ict = 1:nct
            if (strcmp(current_table{ict,1},'resistor'))
                node_from = current_table{ict,2};
                node_to = current_table{ict,3};
                resistance = current_table{ict,4};
                [r_from(ict),z_from(ict)] = node_coordinates (node_from,nr,nz,p.Dr,p.Dz,p.w,N_cells);
                [r_to(ict),z_to(ict)] = node_coordinates (node_to,nr,nz,p.Dr,p.Dz,p.w,N_cells);
                if (node_from == N_nodes)
                    V_from = 0;
                else
                    V_from = ckt.V(node_from);
                end
                if (node_to == N_nodes)
                    V_to = 0;
                else
                    V_to = ckt.V(node_to);
                end
                current(ict) = (V_from-V_to)/resistance;
                if (current(ict)<0)
                    temp = r_from(ict); r_from(ict) = r_to(ict); r_to(ict) = temp;
                    temp = z_from(ict); z_from(ict) = z_to(ict); z_to(ict) = temp;
                    current(ict) = - current(ict);
                end
            end
        end
        max_current = max(current);
        norm_current = current/max_current;
        quiver(z_from,-r_from,...
            norm_current.*(z_to-z_from),norm_current.*(r_from-r_to),0,'k',...
            'LineWidth',1);
        
        axis([0,0.06,-1.5e-3,0]);
        xlabel('z (cm)'); ylabel('-r (cm)');
        title(sprintf('Time = %6.3f ms, I_{max} = %6.3e \\muA, V_{min} = %6.3f mV, V_{max} = %6.3f mV\n',...
            (it-1)*p.Dt,max_current,min(ckt.V),max(ckt.V)));
        
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
        
        drawnow;
    end
    
    %     if (mod(it,100)==0)
    %         figure(2); hold off;
    %         plot_cell_boundaries_v2 (nr,nz,p.Dr,p.Dz,p.w,N_cells);
    %         hold on;
    %         plot_dyn_variable_array (ckt.V,nr,nz,p.Dr,p.Dz,p.w,N_cells);
    %         hold off;
    %         title(sprintf('Timestep = %i,V_{min} = %e, V_{max} = %e\n',it-1,min(ckt.V),max(ckt.V)));
    %         % pause;
    %     end
    
    if (mod(it,1000)==0)
        fprintf('Timestep no. %i\n',it);
    end
    ckt.advance_circuit(p);
    
    
    % Add a stimulus every so often:
    % if (mod(it,it_stim)==1); ckt.V(1) = ckt.V(1) + 100.0; end
    
end

close(vidObj);

rmpath('../circuit-simulation');

toc;


