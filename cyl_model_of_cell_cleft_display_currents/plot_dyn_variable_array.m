function plot_dyn_variable_array (q,nr,nz,Dr,Dz,w,N_cells)
% Plot the array q as a function of r and z.

N_nodes_per_cell = nr*nz + 2*nz + 2*nr + 1 + nr;

epsilon = 0.1; % define distance to plot to correspond to
% "just inside" the membrane.

r_plot = zeros(1,length(q)); z_plot = zeros(1,length(q));

for node_no = 1:length(q)
    node_within_cell = round(mod(node_no-1,N_nodes_per_cell)+1);
    cell_no = round((node_no-node_within_cell)/N_nodes_per_cell);
    if (cell_no==N_cells) % Add right cleft nodes to the last cell.
        cell_no = cell_no-1;
        node_within_cell = node_within_cell + N_nodes_per_cell;
    end
    Z = cell_no * (w+nz*Dz);
    if (node_within_cell <= nr*nz)
        ir = round(mod(node_within_cell-1,nr)+1);
        iz = round( ((node_within_cell-1)-(ir-1))/nr + 1 );
        r = (ir-1)*Dr;
        z = Z + 0.5*w + (iz-0.5)*Dz;
    elseif (node_within_cell <= nr*nz + 2*nz)
        rel_node = node_within_cell - nr*nz;
        if (rel_node <= nz)
            r = (nr-0.5-epsilon)*Dr;
        else
            r = (nr-0.5+epsilon)*Dr;
        end
        iz = round(mod(rel_node-1,nz)+1);
        z = Z + 0.5*w + (iz-0.5)*Dz;
    elseif (node_within_cell <= nr*nz + 2*nz + 2*nr)
        rel_node = node_within_cell - nr*nz - 2*nz;
        ir = round(mod(rel_node-1,nr)+1);
        r = (ir-1)*Dr;
        if (rel_node<=nr)
            z = Z + 0.5*w + epsilon*Dz;
        else
            z = Z;
        end
    elseif (node_within_cell == nr*nz + 2*nz + 2*nr + 1)
        r = (nr-0.5+epsilon)*Dr;
        z = Z;
    elseif (node_within_cell <= nr*nz + 2*nz + 2*nr + 1 + nr)
        rel_node = node_within_cell - nr*nz - 2*nz - 2*nr - 1;
        r = (rel_node-1)*Dr;
        z = Z + 0.5*w + nz*Dz - epsilon*Dz;
    elseif (node_within_cell <= nr*nz + 2*nz + 2*nr + 1 + nr + nr)
        rel_node = node_within_cell - nr*nz - 2*nz - 2*nr - 1 - nr;
        r = (rel_node-1)*Dr;
        z = Z + w + nz*Dz;
    elseif (node_within_cell == nr*nz + 2*nz + 2*nr + 1 + nr + nr + 1)
        r = (nr-0.5+epsilon)*Dr;
        z = Z + w + nz*Dz;
    else
        fprintf('plot_dyn_variable_array error: node(%i) is not a legal node\n',node_no);
    end
    r_plot(node_no) = r; z_plot(node_no) = z;
end

c = colormap;
max_q = max(q); min_q = min(q);

for node_no = 1:length(q)
    if (max_q-min_q==0)
        i = 1;
    else
        i = round(1+255*(q(node_no)-min_q)/(max_q-min_q));
    end
    plot(z_plot(node_no),-r_plot(node_no),'ko','MarkerFaceColor',c(i,:),'MarkerSize',10);
    hold on;
end
hold off;
    

% hold on;
% plot(z_plot,-r_plot,'r*'); plot(z_plot,-r_plot,'r');
% xlabel('z'); ylabel('r');
% % pause;
% plot(z_plot,-r_plot,'b*'); plot(z_plot,-r_plot,'b-');
% hold off;

    
        


        
        

        
