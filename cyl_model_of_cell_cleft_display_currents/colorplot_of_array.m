function colorplot_of_array (q,nr,nz,Dr,Dz,w,N_cells)
% Create a colorplot of the array q as a function of r and z.

% Map the data, which is nonuniformly distributed in space, to a
% rectangular array. Define the each of the rows and columns to correspond
% to specific values of r and z, respectively.
% Modify the routine plot_dyn_variable_array.m to do this.

N_nodes_per_cell = nr*nz + 2*nz + 2*nr + 1 + nr;

epsilon = 0.1; eps_ext = 1.0;
w = 0.5*Dz; % Exaggerate the cleft region width

n_rows = nr+2;
n_cols = (N_cells-1)*(nz+3)+(nz+4);
q_array = zeros(n_rows+1,n_cols+1); % array to store data, including a pad row and col

% Define the correspondence between rows & columns in q_array
% and their r and z values:
r_points = [Dr*(0:(nr-1)),Dr*(nr-0.5),Dr*(nr-0.5+epsilon)];
r_points = [r_points,r_points(end)+Dr]; % pad cell needed for pcolor;
z_points_in_cell = [0, 0.5*w+epsilon*Dz , 0.5*w+0.5*Dz+Dz*(0:(nz-1)) ,...
    0.5*w+0.5*Dz+Dz*(nz-1)+0.5*Dz-epsilon*Dz];
z_points = [];
for i_cell = 0:(N_cells-1)
    z_points = [z_points,i_cell*(w+nz*Dz)+z_points_in_cell];
end
z_points = [z_points,N_cells*(w+nz*Dz)];
z_points = [z_points,N_cells*(w+nz*Dz)+Dz]; % pad cell needed for pcolor

% Define the widths of the elements containing each node:
r_widths = [0.5*Dr,Dr*ones(1,nr-1),epsilon*Dr,eps_ext*Dr];
z_widths_for_one_cell = [w,epsilon*Dz,Dz*ones(1,nz),epsilon*Dz];
z_widths = [];
for i_cell = 0:(N_cells-1)
    z_widths = [z_widths,z_widths_for_one_cell];
end
z_widths = [z_widths,w];
% Use these widths to define the boundaries of the elements surrounding
% each node:
r_bdys = cumsum(r_widths); r_bdys = [0,r_bdys];
z_bdys = cumsum(z_widths); z_bdys = [0,z_bdys];

% Note that: N_nodes_per_cell = (nr+2)*(nz+3) - 5.  
% So five elements of the array in each cell does not correspond
% to a data point.  For the time being, just allow these to be 0.

% Redefine various length parameters so the calculations of (r+1) and 
% (z+1) indices into q_array:
% epsilon = 1; Dr = 1; Dz = 1; half_Dz = 1; half_w = 1;

for node_no = 1:length(q)
    node_within_cell = round(mod(node_no-1,N_nodes_per_cell)+1);
    cell_no = round((node_no-node_within_cell)/N_nodes_per_cell);
    if (cell_no==N_cells) % Add right cleft nodes to the last cell.
        cell_no = cell_no-1;
        node_within_cell = node_within_cell + N_nodes_per_cell;
    end
    Z = cell_no * (1+1+(nz-1)+1+1);
    if (node_within_cell <= nr*nz)
        ir = round(mod(node_within_cell-1,nr)+1);
        iz = round( ((node_within_cell-1)-(ir-1))/nr + 1 );
        r = (ir-1);
        z = Z + 1+1 +(iz-1);
    elseif (node_within_cell <= nr*nz + 2*nz)
        rel_node = node_within_cell - nr*nz;
        if (rel_node <= nz)
            r = nr;
        else
            r = nr+1;
        end
        iz = round(mod(rel_node-1,nz)+1);
        z = Z + 1+1 + (iz-1);
    elseif (node_within_cell <= nr*nz + 2*nz + 2*nr)
        rel_node = node_within_cell - nr*nz - 2*nz;
        ir = round(mod(rel_node-1,nr)+1);
        r = (ir-1);
        if (rel_node<=nr)
            z = Z + 1;
        else
            z = Z;
        end
    elseif (node_within_cell == nr*nz + 2*nz + 2*nr + 1)
        r = nr+1;
        z = Z;
    elseif (node_within_cell <= nr*nz + 2*nz + 2*nr + 1 + nr)
        rel_node = node_within_cell - nr*nz - 2*nz - 2*nr - 1;
        r = (rel_node-1);
        z = Z + 1+1 + (nz-1) + 1;
    elseif (node_within_cell <= nr*nz + 2*nz + 2*nr + 1 + nr + nr)
        rel_node = node_within_cell - nr*nz - 2*nz - 2*nr - 1 - nr;
        r = (rel_node-1);
        z = Z + 1+1 + (nz-1) + 1+1;
    elseif (node_within_cell == nr*nz + 2*nz + 2*nr + 1 + nr + nr + 1)
        r = nr+1;
        z = Z + 1+1 + (nz-1) + 1+1;
    else
        fprintf('colorplot_of_array error: node(%i) is not a legal node\n',node_no);
    end
    q_array(r+1,z+1) = q(node_no);
end

% Fill in reasonable values for those elements of the array that
% do not correspond to nodes:
for cell_no = 0:(N_cells-1)
    Z = cell_no * (1+1+(nz-1)+1+1);
    
    r = nr; z = Z;
    q_array(r+1,z+1) = q_array(r+2,z+1);
    
    r = nr; z = Z+1;
    q_array(r+1,z+1) = ( (0.5*Dz)*q_array(r,z+1) + (0.5*Dr)*q_array(r+1,z+2) )...
        / ( 0.5*Dz + 0.5*Dr );
    
    r = nr+1; z = Z+1;
    q_array(r+1,z+1) = ( (0.5*Dz)*q_array(r+1,z) + (0.5*w)*q_array(r+1,z+2) ) ...
        / ( 0.5*Dz + 0.5*w );
    
    r = nr; z = Z + 1+1+(nz-1)+1;
    q_array(r+1,z+1) = ( (0.5*Dz)*q_array(r,z+1) + (0.5*Dr)*q_array(r+1,z) )...
        / ( 0.5*Dz + 0.5*Dr );
    
    r = nr+1; z = Z + 1+1+(nz-1)+1;
    q_array(r+1,z+1) = ( (0.5*Dz)*q_array(r+1,z+2) + (0.5*w)*q_array(r+1,z) ) ...
        / ( 0.5*Dz + 0.5*w );
    
    if (cell_no==N_cells-1)
        r = nr; z = Z + 1+1+(nz-1)+1+1;
        q_array(r+1,z+1) = q_array(r+2,z+1);
    end
    
end
        
% pcolor(z_points,-r_points,q_array); caxis([-90,30]); axis xy;
% shading interp;
pcolor(z_bdys,-r_bdys,q_array); caxis([-90,30]); axis xy;
colorbar;
