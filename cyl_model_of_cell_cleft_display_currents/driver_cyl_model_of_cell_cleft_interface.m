% Driver that will run cyl_model_of_cell_cleft_interface 
% several times, with different values of the key parameters

clear;
s_INa_opt = {'INa_ends','INa_unif'};
s_Ccl_opt = {'low_Ccl','normal_Ccl'};
s_Rg_opt = {'low_gap_juncts','normal_gap_juncts'};
s_w_opt = {'w_12','w_40'};

for i_INa = 1:2
    for i_Ccl = 1:2
        for i_Rg = 1:2
            for i_w = 1:2
                s_INa = s_INa_opt{i_INa};
                s_Ccl = s_Ccl_opt{i_Ccl};
                s_Rg = s_Rg_opt{i_Rg};
                s_w = s_w_opt{i_w};
                cyl_model_of_cell_cleft_interface;
            end
        end
    end
end
                