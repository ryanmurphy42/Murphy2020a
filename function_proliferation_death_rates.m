function [G_m,G_prolif_m,G_death_m] = ....
    function_proliferation_death_rates(N_cells,m_springs_per_cell, prolif_law,death_law,beta_prolif_m, gamma_death_m,x_current_m,l_crit_death)

%Calculate proliferation and death propensities 
%Prolif law: 201- constant, 202-linear, 203-logistic
%Death law: 201- constant, 202-linear, 203-logistic


G_prolif_m=zeros(N_cells*m_springs_per_cell,1);

if prolif_law == 201 %constant
    G_prolif_m =  beta_prolif_m;
elseif prolif_law == 202 %linear
    spring_sizes = diff(x_current_m);
    N_springs = length(spring_sizes);
    for ii=1:N_springs
        G_prolif_m(ii) = beta_prolif_m(ii)*m_springs_per_cell*spring_sizes(ii);
    end
elseif prolif_law == 203 %logistic
    G_prolif_m =  beta_prolif_m;
end


G_death_m=zeros(N_cells*m_springs_per_cell,1);
if death_law == 201 %constant
    G_death_m =  gamma_death_m;
elseif death_law == 202 %linear
    spring_sizes = diff(x_current_m);
    N_springs = length(spring_sizes);
    for ii=1:N_springs
        if spring_sizes(ii) < l_crit_death
            G_death_m(ii) = gamma_death_m(ii)*(l_crit_death - m_springs_per_cell*spring_sizes(ii));
        else
            G_death_m(ii) = 0;
        end
    end
elseif death_law == 203 %logistic
    spring_sizes = diff(x_current_m);
    N_springs = length(spring_sizes);
    for ii=1:N_springs
        G_death_m(ii) = gamma_death_m(ii)*( 1/(m_springs_per_cell*spring_sizes(ii)));
    end
end


G_m = [G_prolif_m; G_death_m];


end