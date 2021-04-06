% Object-oriented version of the bidomain, regular ephaptic model
% in 1D:

tic;
clear;

% Define model parameters:
define_bidomain_regular_ephaptic_params;
define_EK_params;
p.tau_si = 4.1;

% Numerical parameters:
N_cells = 320; % Number of nodes <- Run with 3 or 20
Nt = 30000; % Number of timesteps to run
it_stim = 13800; %1310; % Pacing interval (in timesteps)
Delta_x = p.L; %Spacing between the nodes
Delta_t = 0.01; % Timestep size
it_plot = 500; % Plotting interval
alpha = 1; % 0.5 = Trapezoidal method; 1 = Backward Euler

% Remove redundant nodes: end up with just two per cell.
% Nodal voltages to plot:
% Plot nodes 1 and 2 for all cells:
cell_no_hist = [0:(N_cells-1),0:(N_cells-2)];
local_node_hist = [1*ones(1,N_cells),2*ones(1,N_cells-1)];
node_hist = cell_no_hist*2 + local_node_hist; % Nodes to plot vs. time.

% Nonlinear currents to plot:
INL_cell_no = 10*ones(1,3);
INL_local_component_no = 1:3;
INL_component_no = INL_cell_no*3 + INL_local_component_no;

% Define history arrays:
V_hist = zeros(length(node_hist),Nt+1);
INL_hist = zeros(length(INL_component_no),Nt+1);
t_hist = zeros(1,Nt+1);

% Add these to the list parameters passed to nonlinear-element functions:
p.Dt = Delta_t; p.Dx = Delta_x;

coupling_fract_of_normal = logspace(-3,0,10);
w = (2:2:26)*1.e-7;

figure(5); clf;

velocity = zeros(length(coupling_fract_of_normal),length(w));

for icfn = 1:length(coupling_fract_of_normal)
    for iw = 1:length(w)
        

% Check numerical stability criterion for the diffusion term:
% fprintf('(Delta t)/((Rm+Re)*Cm) = %f\n',Delta_t/((p.Rm+p.Re)*p.Cm));

% Define the circuit object.
% It has N_cells cells.
% The first (N_cells-1) cells has all 8 internal nodes.
% The last cell just has the first 2 internal nodes:
N_nodes = 2*(N_cells-1)+2;
ckt = circuit(N_nodes,alpha,p.Dt);
% Numbering of nodes:
% N_cells cells numbered from 0 to N_cells-1.
% Associated with cell cell_no are the following "internal" nodes:
% Node 2*cell_no+1: Intracellular node.
% Node 2*cell_no+2: Cell node in cleft.
% Node 2*(N_cells-1)+2 -> Node 0: Ground.
% Construct the circuit.
%
% Notes:
%
% (1) By "internal node," I mean the numbers 0 through 7, as
% defined above, in the definition of the node numbers. Internal node
% 8 and 9 are nodes 0 and 1 of the next cell.
%
% (2) Matlab's "anonymous function" construct is used to define
% the function "nodes," below, which calculates the
% node associated with internal node i associated with cell cell_no.
% The function also checks if the node is the ground node (cell_no=0, i=0)
% and sets it to the last node number in the system, N_nodes.
% (The latter is the current convention of the "circuit" class).
%
        p.w = w(iw);
        p.Rg = 395/coupling_fract_of_normal(icfn); % KOhms
        p.Rr = p.rho_ext/(8*pi*p.w); %Kohms
        
        % Add conventional transmembrane components to all cells:
        for cell_no = 0:(N_cells-1)
            node = @(i) (2*cell_no+i)*(i~=0) + N_nodes*(i==0);
            ckt.add_capacitor(node(1),node(0),p.Cm,-85.0);
            ckt.add_I_nonlin(node(1),node(0),...
                @I_membrane_EK,@advance_EK,@DI_membrane_EK,[1 1]);
            if (cell_no~=N_cells-1)
                ckt.add_resistor(node(1),node(3),p.Rg); % <- Try removing (gap junction resistor) lowest priority
                ckt.add_capacitor(node(1),node(2),0.1*p.Ccl,-85.0); % <-Removing these two capacitors is also interesting
                ckt.add_I_nonlin(node(1),node(2),...
                    @I_cleft_EK,@advance_EK,@DI_cleft_EK,[1 1]);
                ckt.add_resistor(node(2),node(0),p.Rr);
                ckt.add_capacitor(node(3),node(2),0.1*p.Ccl,-85.0); % <-This is the other capacitor to remove
                ckt.add_I_nonlin(node(3),node(2),...
                    @I_cleft_EK,@advance_EK,@DI_cleft_EK,[1 1]);
            end
        end
        
        ckt.prepare_matrices;
        
        % internal_node_guess = [-85,-85,-85,-85,-85,-85,0,0];
        % V_augmented_guess = [repmat(internal_node_guess,1,N_cells-1),-85,...
        %     zeros(1,3*N_cells-2)]';
        % ckt.calc_initial_conditions(p,0.0,V_augmented_guess);
        ckt.calc_initial_conditions(p,0.0);
        V_hist(:,1) = ckt.V(node_hist);
        INL_hist(:,1) = ckt.I_NL(INL_component_no);
        t_hist(1) = ckt.t;
        
        
        V1_hist = zeros(ceil(Nt/it_plot),N_cells);
        i_hist = 0;
        
        cell_A = 10;
        cell_B = N_cells-10;
        V_threshold = 0.0;
        
        wave_detected_in_cell_A = false;
        wave_detected_in_cell_B = false;
        wave_crossing_time_cell_A = nan;
        wave_crossing_time_cell_B = nan;
        
        % Main timestep loop:
        it = 0;
        while ((~wave_detected_in_cell_A)|| (~wave_detected_in_cell_B))
            
            it = it + 1;
            if (mod(it,1000)==0)
                fprintf('Timestep no. %i\n',it);
            end
            ckt.advance_circuit(p);
            % Add a stimulus every so often:
            if (mod(it,it_stim)==1); ckt.V(1) = ckt.V(1) + 200.0; end
            
%             if mod(it,100)==0
%                 fprintf('%f %f\n',ckt.V(2*cell_A+1),ckt.V(2*cell_B+1));
%             end
            
            if (ckt.V(2*cell_A+1)>V_threshold && ~wave_detected_in_cell_A)
                wave_crossing_time_cell_A = ckt.t;
                wave_detected_in_cell_A = true;
            end
            
            if (ckt.V(2*cell_B+1)>V_threshold && ~wave_detected_in_cell_B)
                wave_crossing_time_cell_B = ckt.t;
                wave_detected_in_cell_B = true;
            end
            
            
            V_hist(:,it+1) = ckt.V(node_hist);
            INL_hist(:,it+1) = ckt.I_NL(INL_component_no);
            t_hist(it+1) = ckt.t;
            
            % Plot every so often:
            if (mod(it,it_plot)==0)
                V1 = ckt.V((2*0+1):2:(2*(N_cells-1)+1));
                i_hist = i_hist+1;
                V1_hist(i_hist,:) = V1';
                %         figure(1); clf;
                %         subplot(2,1,1);
                %         plot(V1,'r','LineWidth',2);
                %         ylabel('V_1');
                %         title(sprintf('Timestep = %i',it));
                %         set(gca,'FontSize',16);
                %         subplot(2,1,2);
                %         plot(ckt.V(2:2:end),'b','LineWidth',2);
                %         xlabel('Cell no.');
                %         ylabel('V_2');
                %         set(gca,'FontSize',16);
                %         drawnow;
            end
            
        end % end of timestep loop
        
        % Calculate velocity of this wave (cm/s):
        velocity(icfn,iw) = 1000*p.L*(cell_B-cell_A)/ ...
            (wave_crossing_time_cell_B-wave_crossing_time_cell_A);
        
        figure(5);
        if (iw~=1)
            semilogy(1.e7*w((iw-1):iw),velocity(icfn,(iw-1):iw),'b','LineWidth',2);
        end
        semilogy(1.e7*w(iw),velocity(icfn,iw),'b*'); hold on;  
        drawnow;
        
        toc;
        
        figure(2);
        imagesc((1:(N_cells-1))*p.Dx,(1:i_hist)*it_plot*p.Dt,V1_hist); axis xy;
        title('Membrane potential V (mV)');
        xlabel('x (cm)');
        ylabel('Time (ms)');
        colorbar;
        set(gca,'FontSize',18);
        drawnow;
        
%         % Plot nodal voltages:
%         figure(3);
%         clf;
%         
%         subplot(2,1,1);
%         plot(t_hist,V_hist(1:N_cells,:),'LineWidth',2);
%         ylabel('Local node 1 voltages (mV)');
%         set(gca,'FontSize',16);
%         
%         subplot(2,1,2);
%         plot(t_hist,V_hist((N_cells+1):end,:),'LineWidth',2);
%         ylabel('Local node 2 voltages (mV)');
%         xlabel('Time (ms)');
%         set(gca,'FontSize',16);
%         
%         drawnow;
        
        %
        % for j = 1:length(node_hist)
        %     % subplot(5,2,j)
        %     plot(t_hist,V_hist(j,:),...
        %         'LineWidth',2,...
        %         'DisplayName',sprintf('Cell %i, node %i',cell_no_hist(j),local_node_hist(j)));
        %     legend;
        %     ylabel('Nodal voltages (mV)');
        %     set(gca,'FontSize',16);
        %
        % end
        
%         % Plot nonlinear ion-channel currents:
%         figure(4);
%         clf
%         for j = 1:length(INL_component_no)
%             plot(t_hist,INL_hist(j,:),...
%                 'LineWidth',2,...
%                 'DisplayName',sprintf('Cell %i, comp %i',INL_cell_no(j),INL_local_component_no(j)));
%             hold on;
%         end
%         legend;
%         xlabel('Time (ms)'); ylabel('Nonlin. currents (\muA)');
%         set(gca,'FontSize',16);
%         hold off;
%         drawnow;
        
        % pause;
        
    end
end

figure(5); clf;
for icfn = 1:length(coupling_fract_of_normal)
semilogy(1.e7*w,velocity(icfn,:),'-*','LineWidth',2,...
    'DisplayName',sprintf('%5.3f',coupling_fract_of_normal(icfn)));
hold on;
end
legend;
xlabel('width w (nm)');
ylabel('Conduction velocity (cm/s)');
set(gca,'FontSize',16);
hold off;
