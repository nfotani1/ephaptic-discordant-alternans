function [r_coord,z_coord] = node_coordinates (node_no,nr,nz,Dr,Dz,w,N_cells)
% Find the (r,z) coordinates for node number node_no in the cell-cleft
% cylindrical system.

N_nodes_per_cell = nr*nz + 2*nz + 2*nr + 1 + nr;

epsilon = 0.1; % define distance to plot to correspond to
% "just inside" the membrane.
eps_ext = 1.0; % distance "just outside" the membrane
w = 0.5*Dz; % Exaggerate the cleft region width

node_within_cell = round(mod(node_no-1,N_nodes_per_cell)+1);
cell_no = round((node_no-node_within_cell)/N_nodes_per_cell);
if (cell_no==N_cells) % Add right cleft nodes to the last cell.
    cell_no = cell_no-1;
    node_within_cell = node_within_cell + N_nodes_per_cell;
end
Z = 0.5*w+cell_no * (2*epsilon*Dz+w+nz*Dz);
if (node_within_cell <= nr*nz)
    ir = round(mod(node_within_cell-1,nr)+1);
    iz = round( ((node_within_cell-1)-(ir-1))/nr + 1 );
    r = (ir-1)*Dr;
    z = Z + 0.5*w + epsilon*Dz + (iz-0.5)*Dz;
elseif (node_within_cell <= nr*nz + 2*nz)
    rel_node = node_within_cell - nr*nz;
    if (rel_node <= nz)
        r = (nr-0.5-0.5*epsilon)*Dr;
    else
        r = (nr-0.5+epsilon+0.5*eps_ext)*Dr;
    end
    iz = round(mod(rel_node-1,nz)+1);
    z = Z + 0.5*w +epsilon*Dz + (iz-0.5)*Dz;
elseif (node_within_cell <= nr*nz + 2*nz + 2*nr)
    rel_node = node_within_cell - nr*nz - 2*nz;
    ir = round(mod(rel_node-1,nr)+1);
    r = (ir-1)*Dr;
    if (rel_node<=nr)
        z = Z + 0.5*w + 0.5*epsilon*Dz;
    else
        z = Z;
    end
elseif (node_within_cell == nr*nz + 2*nz + 2*nr + 1)
    r = (nr-0.5+epsilon+0.5*eps_ext)*Dr;
    z = Z;
elseif (node_within_cell <= nr*nz + 2*nz + 2*nr + 1 + nr)
    rel_node = node_within_cell - nr*nz - 2*nz - 2*nr - 1;
    r = (rel_node-1)*Dr;
    z = Z + 0.5*w + epsilon*Dz + nz*Dz + 0.5*epsilon*Dz;
elseif (node_within_cell <= nr*nz + 2*nz + 2*nr + 1 + nr + nr)
    rel_node = node_within_cell - nr*nz - 2*nz - 2*nr - 1 - nr;
    r = (rel_node-1)*Dr;
    z = Z + (2*epsilon*Dz+w+nz*Dz);
elseif (node_within_cell == nr*nz + 2*nz + 2*nr + 1 + nr + nr + 1)
    r = (nr-0.5+epsilon+0.5*eps_ext)*Dr;
    z = Z + (2*epsilon*Dz+w+nz*Dz);
else
    fprintf('node_coordinates error: %i is not a legal node number.\n',node_no);
end
r_coord = r; z_coord = z;







