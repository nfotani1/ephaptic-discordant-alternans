function J = J_cleft_EK(V,q,p)
% Same as I_cleft_EK, but returns the current -density-
% (current per unit area) instead of the current assuming
% a particular area.  
% Another difference: uses the capacitance per unit
% area given by p.cm_cleft, not by p.cm.
u = (V-p.V_rest)/(p.V_fi-p.V_rest);

h = q(:,1); f = q(:,2);

m_inf = (u>0).*(u/0.2).^6./(1+(u/0.2).^6);
d_inf = (u>0).*(u/0.4).^4./(1+(u/0.4).^4);

J_fi = p.g_fi_cleft*h.*m_inf.*(u-1.3)/p.tau_fi;
J_si = f.*d_inf.*(u-1.4)./p.tau_si;
J_so = (1-exp(-4*u))/p.tau_so;

J = J_fi + J_si + J_so; % (1/ms, as defined in the EPJ article)
J = J * (p.V_fi-p.V_rest); % (mV/ms = uA/uF, since V/s = A/F)
J = J * p.cm_cleft; % current per unit membrane area (uA/cm^2)