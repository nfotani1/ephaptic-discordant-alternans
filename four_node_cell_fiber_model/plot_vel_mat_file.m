figure(8); clf;
% load vel_normal_Ccl; 
load vel_16xCcl_finish;
velocity1 = velocity;
load vel_16xCcl;
velocity = velocity + velocity1;

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