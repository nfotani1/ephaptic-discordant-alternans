function q_new = advance_EK(V,q,p)

% Advance auxiliary variables in the Echebarria-Karma 3-variable model from EPJ.

u = (V-p.V_rest)/(p.V_fi-p.V_rest);

h = q(:,1);
f = q(:,2);

tau_h = p.tau_h1 + p.tau_h2*exp(-20*(u-0.1).^2);
tau_f = p.tau_f2 + (p.tau_f1-p.tau_f2)*u.^3;

h_inf = (u>0)./(1+(u/0.1).^6) + (u<=0);
f_inf = (u>0)./(1+(u/0.1).^4) + (u<=0);

h = h + p.Dt*(h_inf - h)./tau_h;
f = f + p.Dt*(f_inf - f)./tau_f;

q_new = [h, f];


