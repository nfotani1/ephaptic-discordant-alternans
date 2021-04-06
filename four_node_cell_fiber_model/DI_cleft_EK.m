function DI = DI_cleft_EK(V,q,p)
% Assume V is a column vector of component voltages,
% and q is an array, each whose rows contains
% auxiliary variables for one component:

u = (V-p.V_rest)/(p.V_fi-p.V_rest);
du_dV = 1/(p.V_fi-p.V_rest);

h = q(:,1); f = q(:,2);

u2= u/0.2;
u26 = u2.^6;
m_inf = (u>0).*(1-1./(1+u26));
dm_inf_du = (u>0)./((1+u26).*(1+u26)).*(6*u2.^5)/0.2;

u4 = u/0.4;
u44 = u4.^4;
d_inf = (u>0).*(1-1./(1+u44));
dd_inf_du = (u>0)./(1+u44).^2.*(4*u4.^3)/0.4;

% J_fi = p.g_fi_cleft*h.*m_inf.*(u-1.3)/p.tau_fi;
dJ_fi_du = p.g_fi_cleft*h.*( dm_inf_du.*(u-1.3) + m_inf )/p.tau_fi;

% J_si = f.*d_inf.*(u-1.4)/p.tau_si;
dJ_si_du = f.*( dd_inf_du.*(u-1.4) + d_inf)/p.tau_si;

% J_so = (1-exp(-4*u))/p.tau_so;
dJ_so_du = 4*exp(-4*u)/p.tau_so;

dJ_du = ( dJ_fi_du + dJ_si_du + dJ_so_du ) ;
dJ_dV = dJ_du * du_dV;

dJ_fi_dh = p.g_fi_cleft*m_inf.*(u-1.3)/p.tau_fi;
dJ_si_df = d_inf.*(u-1.4)/p.tau_si;

dJ_dh = dJ_fi_dh;
dJ_df = dJ_si_df;
DJ = [dJ_dV,dJ_dh,dJ_df];
DI = (pi*(p.r)^2)*DJ * (p.V_fi-p.V_rest) * p.cm; 

