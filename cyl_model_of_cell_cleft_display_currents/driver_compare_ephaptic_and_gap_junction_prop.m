% Compare the spatial patterns between:
% (a) ephaptic-dominated wave propagation, and
% (b) gap-junction dominated wave propagation,
% for the same propagation speed.

clear;

s_INa_opt = {'INa_ends','INa_unif'};
s_Ccl_opt = {'low_Ccl','normal_Ccl'};
s_Rg_opt = {'low_gap_juncts','normal_gap_juncts'};
s_w_opt = {'w_12','w_40'};

% % % Ephaptic-dominated propagation:
% % s_INa = 'INa_ends';
% % s_Ccl = 'low_Ccl';
% % s_Rg = 'low_gap_juncts';
% % s_w = 'w_12';
% % cyl_cell_cleft_model_w_currents;

% % Gap-junction-dominated propagation w/ Na channels on ends and low Ccl:
% s_INa = 'INa_ends';
% s_Ccl = 'low_Ccl';
% s_Rg = 'custom_gap_juncts'; gap_alpha = 0.045;
% s_w = 'w_40';
% cyl_cell_cleft_model_w_currents;

% % % Gap-junction-dominated propagation w/ Na channels on ends and normal Ccl:
% % s_INa = 'INa_ends';
% % s_Ccl = 'normal_Ccl';
% % s_Rg = 'custom_gap_juncts'; gap_alpha = 0.045;
% % s_w = 'w_40';
% % cyl_cell_cleft_model_w_currents;

% Gap-junction-dominated propagation w/ Na channels uniformly distributed:
s_INa = 'INa_unif';
s_Ccl = 'low_Ccl';
s_Rg = 'custom_gap_juncts'; gap_alpha = 0.095;
s_w = 'w_40';
cyl_cell_cleft_model_w_currents;