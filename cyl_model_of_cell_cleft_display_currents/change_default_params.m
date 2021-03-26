% Make changes to the default parameters in this file.
%
% We currently believe that the optimal strong ephaptic case is:
%   1. All the sodium channels on the ends
%   2. Low cleft capacitance
%   3. Very low gap junction coupling
%   4. Optimal cleft width (12 nm)

% Uncomment one of each of the following 4 pairs:

if (strcmp(s_INa,'INa_ends'))
    % All the sodium channels on the ends:
    total_cell_surf_area = 2*pi*p.r^2 + 2*pi*p.r*p.L;
    p.g_fi_cleft = 0.5*total_cell_surf_area/(pi*p.r^2);
    p.g_fi_tm = 0;
elseif (strcmp(s_INa,'INa_unif'))
    % Or spread the sodium channels uniformly:
    p.g_fi_cleft = 1.0;
    p.g_fi_tm = 1.0;
end

if (strcmp(s_Ccl,'low_Ccl'))
    % Make the cleft capacitance low:
    p.Ccl_factor = 0.1;
elseif (strcmp(s_Ccl,'normal_Ccl'))
    % Or leave it at its default value;
    p.Ccl_factor = 1.0;
end

if (strcmp(s_Rg,'low_gap_juncts'))
    % Use very low gap junction coupling:
    p.sigma_g = 0.001*(1/395.0)/(pi*p.r^2);
elseif (strcmp(s_Rg,'normal_gap_juncts'))
    % Or use the normal value:
    p.sigma_g = (1/395.0)/(pi*p.r^2);
elseif (strcmp(s_Rg,'custom_gap_juncts'))
    % Define gap_alpha in this case, before calling this routine:
    p.sigma_g = gap_alpha*(1/395.0)/(pi*p.r^2);
end

% (Note: may need to re-determine for the discretized 
% cyclindrical simulation)
if (strcmp(s_w,'w_12'))
    % Use the optimal cleft width for strong ephaptic effect:
    p.w = 12*1.e-7; % (12 nm)
elseif (strcmp(s_w,'w_40'))
    % Or use the default width:
    p.w = 40*1.e-7; % (40 nm)
end