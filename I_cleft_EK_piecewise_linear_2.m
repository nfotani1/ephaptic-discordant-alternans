function I = I_cleft_EK_piecewise_linear_2(V,q,p)

V_thresh = -72.0; E2 = 44.7; R2_int = 20675.0;
% V_thresh = -72.0; E2 = 45; R2_int = 20229.0;

I = (V>V_thresh).* ((V-E2)/R2_int);