% Object-oriented version of the 4-node ephaptic model in 1D.

tic;
clear;

addpath('../circuit-simulation');

% Define default model parameters:
define_default_4_node_ephaptic_params;
define_EK_params;

% Default numerical parameters:
N_cells = 320; % Number of nodes <- Run with 3 or 20
Nt = 40*50000; % Number of timesteps to run
it_stim = 13100; %1310; % Pacing interval (in timesteps)
Delta_x = p.L; %Spacing between the nodes
Delta_t = 0.01; % Timestep size
it_plot = 500; % Plotting interval
alpha = 1; % 0.5 = Trapezoidal method; 1 = Backward Euler

% --Make the following changes from the default--

p.tau_si = 4.1; % 4.1 promotes discordant alteranns

% For now, just make one of these true:
velocity_measurement_mode = false;
plot_extended_diagnostics = true;

% Comment one of these:
% % All sodium channels on ends:
% total_cell_surf_area = 2*pi*p.r^2 + 2*pi*p.r*p.L;
% p.g_fi_cleft = 0.5*total_cell_surf_area/(pi*p.r^2);
% p.g_fi_tm = 0;
% Or sodium channels are uniformly distributed:
p.g_fi_cleft = 1;
p.g_fi_tm = 1;

% Reduce cleft capacitance by a factor:
fraction_normal_Ccl = 1.0;
p.Ccl = fraction_normal_Ccl*p.Ccl;

% Use the following range of gap junction couplings:
% gap_junction_coupling_fraction_range = logspace(-3,0,10);
gap_junction_coupling_fraction_range = 0.1;

% Cleft width range
% w = (2:2:30)*1.e-7;
% w = 26*1.e-7;
w = 12.e-7;

% --End of changes to default parameters--

velocity = zeros(length(gap_junction_coupling_fraction_range),length(w));

for igc = 1:length(gap_junction_coupling_fraction_range)
    p.Rg = 395.0/gap_junction_coupling_fraction_range(igc);
    for iw = 1:length(w)
        p.w = w(iw);
        p.Rr = p.rho_ext/(8*pi*p.w); %Kohms
        
        velocity_calculated = false;
        
        % For velocity measurement diagnostics:
        % Measure velocity between cells A and B:
        if (velocity_measurement_mode)
            cell_A = 10;
            cell_B = N_cells - 10;
            wave_detected_in_cell_A = false;
            wave_detected_in_cell_B = false;
            Nt = 1000000; % basically infinity. 
        end
        
        % Nodal voltages t0 plot:
        % Plot nodes 1 and 3 for all cells:
%         cell_no_hist = [0:(N_cells-1),0:(N_cells-2)];
%         local_node_hist = [1*ones(1,N_cells),3*ones(1,N_cells-1)];
%         node_hist = cell_no_hist*4 + local_node_hist; % Nodes to plot vs. time.
        local_node_hist = [2,3,4,5,6];
        cell_no_hist = round(N_cells/2)*ones(1,length(local_node_hist));
        node_hist = cell_no_hist*4 + local_node_hist; % Nodes to plot vs. time.
        
        % Nonlinear currents to plot:
        INL_local_component_no = 1:4;
        INL_cell_no = round(N_cells/2)*ones(1,length(INL_local_component_no));
        INL_component_no = INL_cell_no*3 + INL_local_component_no;
        
        % Define history arrays:
        V_hist = zeros(length(node_hist),Nt+1);
        INL_hist = zeros(length(INL_component_no),Nt+1);
        t_hist = zeros(1,Nt+1);
        
        % Add these to the list parameters passed to nonlinear-element functions:
        p.Dt = Delta_t; p.Dx = Delta_x;
        
        % Set up arrays to store information for each
        % wavefront arrival event:
        N_waves = round(Nt/it_stim)+1; % approx no. of waves.
        wf.time = nan(N_waves,N_cells);
        wf.f = nan(N_waves,N_cells);
        wf.h = nan(N_waves,N_cells);
        wf.DI = nan(N_waves,N_cells);
        wf.i_wave = zeros(1,N_cells);
        V_wavefront_threshold = -65.0;
        
        %Set up arrays to store information for each wavefront arrival event:
        wb.time = nan(N_waves,N_cells);
        wb.f = nan(N_waves,N_cells);
        wb.h = nan(N_waves,N_cells);
        wb.APD = nan(N_waves,N_cells);
        wb.i_wave = zeros(1,N_cells);
        V_waveback_threshold = -65.0;
        
        % Define the circuit object.
        N_nodes = 4*N_cells+2;
        ckt = circuit(N_nodes,alpha,p.Dt);
        % Numbering of nodes:
        % N_cells cells numbered from 0 to N_cells-1.
        % Associated with cell cell_no are the following "internal" nodes:
        % Node 4*cell_no+1: Node in cleft to the left of the cell.
        % Node 4*cell_no+2: Node just inside the left end of the cell.
        % Node 4*cell_no+3: Node in the center of the cell.
        % Node 4*cell_no+4: Node just inside the right end of the cell.
        % Two extra nodes:
        % Node 4*(N_cells-1)+1: Node in the cleft to the right of the last cell.
        % Node 4*(N_cells-1)+2 -> Node 0: Ground.
        %
        % Construct the circuit.
        %
        % Notes:
        %
        % (1) By "internal node," I mean the numbers 1 through 4, as
        % defined above, in the definition of the node numbers. 
        %
        % (2) Matlab's "anonymous function" construct is used to define
        % the function "nodes," below, which calculates the
        % node associated with internal node i associated with cell cell_no.
        % The function also checks if the node is the ground node (cell_no=0, i=0)
        % and sets it to the last node number in the system, N_nodes.
        % (The latter is the current convention of the "circuit" class).
        %
        % Add conventional transmembrane components to all cells:
        for cell_no = 0:(N_cells-1)
            node = @(i) (4*cell_no+i)*(i~=0) + N_nodes*(i==0);
            ckt.add_resistor(node(1),node(0),p.Rr);
            ckt.add_capacitor(node(2),node(1),p.Ccl,-85.0);
            ckt.add_I_nonlin(node(2),node(1),...
                @I_cleft_EK,@advance_EK,@DI_cleft_EK,[1 1]);
            ckt.add_resistor(node(3),node(2),p.Rm);
            ckt.add_capacitor(node(3),node(0),p.Cm,-85.0);
            ckt.add_I_nonlin(node(3),node(0),...
                @I_membrane_EK,@advance_EK,@DI_membrane_EK,[1 1]);
            ckt.add_resistor(node(3),node(4),p.Rm);
            ckt.add_capacitor(node(4),node(5),p.Ccl,-85.0); 
            ckt.add_I_nonlin(node(4),node(5),...
                @I_cleft_EK,@advance_EK,@DI_cleft_EK,[1 1]);            
            if (cell_no~=N_cells-1)
                ckt.add_resistor(node(4),node(6),p.Rg);
            end
            if (cell_no==N_cells-1)
                ckt.add_resistor(node(5),node(0),p.Rr);
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
        
        % Main timestep loop:
        for it = 1:Nt
            
            % Keep track of previous V for use by the wavefront
            % arrival recorder:
            ckt_V_old = ckt.V;
            
            if (mod(it,1000)==0)
                fprintf('Timestep no. %i\n',it);
            end
            ckt.advance_circuit(p);
            % Add a stimulus every so often:
            if (mod(it,it_stim)==1); ckt.V(3) = ckt.V(3) + 100.0; end
            
            % Wavefront arrival recorder.
            % If a wavefront has arrived in any of the cells,
            % designated in Wavefront_ix_record, record data:
            %
            % First, find the indices of those cells in which
            % a wavefront arrival event has occurred:
            wavefront_cells = find((ckt_V_old(3:4:end)<V_wavefront_threshold) ...
                & (ckt.V(3:4:end)>=V_wavefront_threshold));
            % For each of those cells, record data:
            for i_cell = wavefront_cells' % Note: i_cell = cell_no+1.
                wf.i_wave(i_cell) = wf.i_wave(i_cell) + 1;
                wave_no = wf.i_wave(i_cell);
                wf.time(wave_no,i_cell) = ckt.t;
                wf.h(wave_no,i_cell) = ckt.q_array(2+3*(i_cell-1),1);
                wf.f(wave_no,i_cell) = ckt.q_array(2+3*(i_cell-1),2);
                if wave_no > 1
                    wf.DI(wave_no,i_cell) = wf.time(wave_no,i_cell)- wb.time(wave_no-1,i_cell);% probably put this later(next loop)
                end
                if (velocity_measurement_mode)
                    if (i_cell==cell_A && ~wave_detected_in_cell_A)
                        wave_detected_in_cell_A = true;
                        wave_crossing_time_cell_A = ckt.t;
                    end
                    if (i_cell==cell_B)
                        wave_detected_in_cell_B = true;
                        if (wave_detected_in_cell_A==false)
                            disp('Wave reached cell B before cell A');
                        end
                        wave_crossing_time_cell_B = ckt.t;
                        % Velocity in cm/s:
                        velocity(igc,iw) = 1000*(cell_B-cell_A)*Delta_x / ...
                            (wave_crossing_time_cell_B - wave_crossing_time_cell_A);
                        fprintf('Calc. velocity = %4.2f cm/s for coupling fraction %4.2e and w = %i nm.\n',...
                             velocity(igc,iw),gap_junction_coupling_fraction_range(igc),...
                             w(iw)*1.e+7);
                        velocity_calculated = true;
                    end
                end
            end
            if (velocity_calculated); break; end
            
            % Here is my attempt at the waveback calculation:
            waveback_cells = find((ckt_V_old(3:4:end)>=V_waveback_threshold)...
                & (ckt.V(3:4:end)<V_waveback_threshold));
            for i_cell = waveback_cells'
                wb.i_wave(i_cell) = wb.i_wave(i_cell) + 1;
                waveb_no = wb.i_wave(i_cell);
                wb.time(waveb_no,i_cell) = ckt.t;
                wb.h(waveb_no,i_cell) = ckt.q_array(2+3*(i_cell-1),1);
                wb.f(waveb_no,i_cell) = ckt.q_array(2+3*(i_cell-1),2);
                wb.APD(waveb_no,i_cell) = wb.time(waveb_no,i_cell) - wf.time(waveb_no,i_cell);
            end
            V_hist(:,it+1) = ckt.V(node_hist);
            INL_hist(:,it+1) = ckt.I_NL(INL_component_no);
            t_hist(it+1) = ckt.t;
            
            % Plot every so often:
            
            
            if (mod(it,it_plot)==0)
                V1 = ckt.V((4*0+3):4:(4*(N_cells-1)+3));
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
            
        end
        
        toc;
        
        if (plot_extended_diagnostics)
            
            figure(2);
            imagesc((1:(N_cells-1))*p.Dx,(1:i_hist)*it_plot*p.Dt,V1_hist); axis xy;
            title('Membrane potential V (mV)');
            xlabel('x (cm)');
            ylabel('Time (ms)');
            colorbar;
            set(gca,'FontSize',18);
            
            % Plot nodal voltages:
            figure(3);
            clf;
            
            plot(t_hist,V_hist,'LineWidth',2);
            set(gca,'FontSize',16);
            legend('Node 2','Node 3','Node 4','Node 5','Node 6');

%             
%             subplot(2,1,1);
%             plot(t_hist,V_hist(1:N_cells,:),'LineWidth',2);
%             ylabel('Local node 1 voltages (mV)');
%             set(gca,'FontSize',16);
%             
%             subplot(2,1,2);
%             plot(t_hist,V_hist((N_cells+1):end,:),'LineWidth',2);
%             ylabel('Local node 2 voltages (mV)');
%             xlabel('Time (ms)');
%             set(gca,'FontSize',16);
            
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
            
            % Plot nonlinear ion-channel currents:
            figure(4);
            clf
            for j = 1:length(INL_component_no)
                plot(t_hist,INL_hist(j,:),...
                    'LineWidth',2,...
                    'DisplayName',sprintf('Cell %i, comp %i',INL_cell_no(j),INL_local_component_no(j)));
                hold on;
            end
            legend;
            xlabel('Time (ms)'); ylabel('Nonlin. currents (\muA)');
            set(gca,'FontSize',16);
            hold off;
            
            % Plot wavefronts
            figure(2);hold on
            %clf;
            plot((1:320)*p.Dx,wf.time','g','LineWidth',2);
            plot((1:320)*p.Dx,wb.time' ,'r','LineWidth',2);
            %plot(wf.DI','c','LineWidth',2);
            set(gca,'FontSize',16); hold off;
            
            
            figure(6);
            clf;
            [nr,nc] = size(wf.time);
            ix = 30:(N_cells-30);
            for ir = 1:nr
                vel = 1000*(20*p.L)./(wf.time(ir,ix+10)-wf.time(ir,ix-10));
                plot(ix*p.L,vel,'LineWidth',2);
                xlabel('x (cm)'); ylabel('Conduction velocity (cm/s)');
                set(gca,'FontSize',16);
                hold on;
                pause;
            end
            hold off;
            
            figure(7);
            [nr,nc] = size(wf.time);
            ix_range = 30:(N_cells-30);
            h = zeros(size(ix_range));
            for ir = 1:nr
                for i = 1:length(ix_range)
                    ix = ix_range(i);
                    h(i) = mean(wf.f(ir,(ix-10):(ix+10)));
                end
                plot(ix_range*p.L,h,'LineWidth',2);
                xlabel('x (cm)'); ylabel('f');
                set(gca,'FontSize',16);
                hold on;
                pause;
            end
            hold off;
            
            figure(8);
            [nr,nc] = size(wf.DI);
            ix_range = 30:(N_cells-30);
            DI = zeros(size(ix_range));
            for ir = 2:nr
                for i = 1:length(ix_range)
                    ix = ix_range(i);
                    DI(i) = mean(wf.DI(ir,(ix-10):(ix+10)));
                end
                plot(ix_range*p.L,DI,'LineWidth',2);
                xlabel('x (cm)'); ylabel('DI');
                set(gca,'FontSize',16);
                hold on;
                pause;
            end
            hold off;
            
            % Plot wavefront velocities vs. DI from wavefront
            % data from the simulation.
            figure(9);
            [nr,nc] = size(wf.DI);
            ix_range = 30:(N_cells-30);
            DI = zeros(size(ix_range));
            vel = zeros(size(ix_range));
            x =[]; y =[];
            for ir = 2:nr
                for i = 1:length(ix_range)
                    ix = ix_range(i);
                    DI(i) = mean(wf.DI(ir,(ix-10):(ix+10)));
                    vel(i) = 1000*(20*p.L)./(wf.time(ir,ix+10)-wf.time(ir,ix-10));
                end
                x = [x,DI]; y = [y,vel];
                plot(DI,vel,'LineWidth',2);
                xlabel('DI(ms)'); ylabel('Vel (cm/s)');
                set(gca,'FontSize',16);
                hold on;
                pause;
            end
            hold off;
            
            % Fit the velocity vs. DI curve.
            % Keep only those indices for which both x and y are defined:
            ind = find((~isnan(x)) & (~isnan(y)));
            xx = x(ind); yy = y(ind);
            plot(xx,yy,'.');
            xlabel('DI (ms)'); ylabel('Velocity c (ms)');
            modelfun = @(b,x) b(1) - b(2)*exp(-x/b(3));
            % Find an initial guess [b1,b2,b3] for our curve-fitter:
            [y1,i1] = min(yy); x1 = xx(i1);
            [y2,i2] = max(yy); x2 = xx(i2);
            x_mid = 0.5*(x1+x2);
            [x_rel,i3] = min(abs(xx-x_mid)); x3 = x_mid+x_rel; y3 = yy(i3);
            b1 = y2;
            b3 = (x3-x1)/log((y1-b1)/(y3-b1));
            b2 = -(y1-b1)*exp(x1/b3);
            % Now curvefit:
            % beta = nlinfit(xx,yy,modelfun,[47,1,1])
            beta = nlinfit(xx,yy,modelfun,[b1,b2,b3])
            
            c_vs_DI = @(DI) beta(1)-beta(2)*exp(-DI/beta(3));
            dcdDI = @(DI) 1000*(beta(2)/beta(3))*exp(-DI/beta(3)); % in cm/s^2
            
            DI = 10:0.01:90;
            hold on;
            plot(DI,c_vs_DI(DI),'r');
            hold off;
            
            % Plot APD vs. DI from wavefront and waveback
            % data from the simulation.
            figure(10); clf;
            [nr,nc] = size(wf.DI);
            ix_range = 30:(N_cells-30);
            DI = zeros(size(ix_range));
            APD = zeros(size(ix_range));
            x =[]; y =[];
            for ir = 2:nr
                for i = 1:length(ix_range)
                    ix = ix_range(i);
                    DI(i) = mean(wf.DI(ir,(ix-10):(ix+10)));
                    APD(i) = mean(wb.APD(ir,(ix-10):(ix+10)));
                end
                x = [x,DI]; y = [y,APD];
                plot(DI,APD,'LineWidth',2);
                xlabel('DI(ms)'); ylabel('APD (ms)');
                set(gca,'FontSize',16);
                hold on;
                pause;
            end
            hold off;
            
            % Fit the APD vs. DI curve.
            % Keep only those indices for which both x and y are defined:
            ind = find((~isnan(x)) & (~isnan(y)));
            xx = x(ind); yy = y(ind);
            plot(xx,yy,'.');
            xlabel('DI (ms)'); ylabel('APD (ms)');
            modelfun = @(b,x) b(1) - b(2)*exp(-x/b(3));
            % Find an initial guess [b1,b2,b3] for our curve-fitter:
            [y1,i1] = min(yy); x1 = xx(i1);
            [y2,i2] = max(yy); x2 = xx(i2);
            x_mid = 0.5*(x1+x2);
            [x_rel,i3] = min(abs(xx-x_mid)); x3 = x_mid+x_rel; y3 = yy(i3);
            b1 = y2;
            b3 = (x3-x1)/log((y1-b1)/(y3-b1));
            b2 = -(y1-b1)*exp(x1/b3);
            % Now curvefit:
            beta = nlinfit(xx,yy,modelfun,[b1,b2,b3])
            
            APD_vs_DI = @(DI) beta(1)-beta(2)*exp(-DI/beta(3));
            dAPDdDI = @(DI) (beta(2)/beta(3))*exp(-DI/beta(3)); % in (cm/s)/ms
            
            DI = 10:0.01:90;
            hold on;
            plot(DI,APD_vs_DI(DI),'r');
            
            axis equal;
            
            fprintf('Apply Echebaria-Karma theory to this simulation.\n');
            
            % Calculate APD and DI where the slope of the
            % APD vs. DI curve is 1.
            
            fprintf('Here are the parameters needed by this theory:\n\n');
            DI_bar = - beta(3)*log(beta(3)/beta(2));
            fprintf('DI_bar = %f ms\n',DI_bar);
            APD_bar = APD_vs_DI(DI_bar);
            fprintf('APD_bar = %f ms\n',APD_bar);
            plot(DI_bar,APD_bar,'*');
            BCL = DI_bar + APD_bar;
            fprintf('BCL = %f ms\n',BCL);
            plot([BCL,0],[0,BCL],'g');
            
            set(gca,'FontSize',16');
            c = c_vs_DI(DI_bar); % cm/s
            fprintf('c = %f cm/s\n',c);
            c_prime = dcdDI(DI_bar); %cm/s^2
            fprintf('c'' = %f cm/s^2\n',c_prime);
            Lambda = c^2/(2*c_prime); % cm
            fprintf('Lambda = %f cm\n',Lambda);
            
            % Calculate Echebarria-Karma (2002 PRL) parameters:
            % Width of transition region between longer and
            % shorter APD regions in a simulation designed to
            % show this region:
            xi = 0.016; % cm
            fprintf('xi = %f cm\n',xi);
            
            D = 1000*xi^2/APD_bar; % cm^2/s
            fprintf('D = %f cm^2/s\n',D);
            
            w = 2*D / c ; % cm
            fprintf('w = %f cm\n',w);
            fprintf('\n');
            
            % Twice the spacing between discordant alterans
            % nodes, from Echebarria & Karma (2002):
            
            fprintf('Dispersion is weak if 4*Lambda >> xi^4/w^3.\n');
            fprintf('Dispersion is strong if 4*Lambda << xi^4/w^3.\n');
            fprintf('\nHere, we have that:\n\n');
            
            fprintf('4*Lambda = %f\n',4*Lambda);
            fprintf('xi^4/w^3 = %f\n\n',xi^4/w^3);
            
            fprintf('Echebarria-Karma theory predicts that twice the\n');
            fprintf('spacing between discordant alternans nodes should be:\n\n');
            lambda_1 = 2*pi*sqrt(w*Lambda); % cm
            fprintf('lambda = %f cm is dispersion is weak, and\n',lambda_1);
            lambda_2 = (4*pi/sqrt(3))*(2*xi*xi*Lambda)^(1/3); % cm
            fprintf('lambda = %f cm is dispersion is strong.\n\n',lambda_2);
            
            rmpath('../circuit-simulation');
            
        end
        
    end
end

if (velocity_measurement_mode)
    figure(21); clf;
    for igc = 1:length(gap_junction_coupling_fraction_range)
        semilogy(1.e7*w,velocity(igc,:),'-*','LineWidth',2,...
            'DisplayName',sprintf('%5.3f',gap_junction_coupling_fraction_range(igc)));
        hold on;
    end
    legend;
    xlabel('width w (nm)');
    ylabel('Conduction velocity (cm/s)');
    set(gca,'FontSize',16)
    hold off;
end

