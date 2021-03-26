function plot_branches (node_1,node_2,nr,nz,Dr,Dz,w,N_cells,comment)
% Previously called plot_nodes.  The new name is more accurate.
% For code debugging purposes, plot nodes 1 and 2.

node(1) = node_1; node(2) = node_2;
r_plot = zeros(1,2); z_plot = zeros(1,2);
for i = 1:2
    [r_coord,z_coord] = node_coordinates (node(i),nr,nz,Dr,Dz,w,N_cells);
    r_plot(i) = r_coord; z_plot(i) = z_coord;
end

hold on;
plot(z_plot,-r_plot,'r*'); plot(z_plot,-r_plot,'r');
title(comment);
xlabel('z'); ylabel('r');
drawnow;
plot(z_plot,-r_plot,'b*'); plot(z_plot,-r_plot,'b-');
hold off;
