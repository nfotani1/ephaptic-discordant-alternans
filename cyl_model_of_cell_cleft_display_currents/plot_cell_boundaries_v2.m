function plot_cell_boundaries_v2 (nr,nz,Dr,Dz,w,N_cells)
% Plot outlines of the cells in r-z geometry.

w = 0.5*Dz; % make cleft region artificially wide, for clarity.
epsilon = 0.1; % Width of internal boundary cells in units of Dz or Dr.
eps_ext = 1.0; % Width of the extracellular space in units of Dr.

axial_cell_length = (nz+2*epsilon)*Dz;
axial_cell_spacing = w + axial_cell_length;

prior_hold = ishold;
plot([0,N_cells*axial_cell_spacing+w],[0,0],'c--'); 
hold on;
for i_cell = 0:(N_cells-1)
    z_left = w + i_cell*axial_cell_spacing;
    z_right = z_left + axial_cell_length;
    r_center = 0;
    r_outer = -(nr-0.5+epsilon)*Dr;
    plot([z_left,z_left],[r_center,r_outer],'c');
    plot([z_right,z_right],[r_center,r_outer],'c');
    plot([z_left,z_right],[r_outer,r_outer],'c');
end
axis([0,N_cells*axial_cell_spacing+w,-(nr-0.5+epsilon+eps_ext)*Dr,0]);
% Return hold status back to what it was
if (prior_hold); hold on; else; hold off; end
