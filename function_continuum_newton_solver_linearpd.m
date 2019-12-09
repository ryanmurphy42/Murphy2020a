function [q,k,a,beta,gamma,s,node_velocity,bc_switch_made,w,res] = function_continuum_newton_solver_linearpd(q,k,a,beta,gamma,l_crit_death,s,nodesx,...
    dt,dx,eta,max_iters,err_tol,k_hom, a_hom, beta_hom, gamma_hom,one_or_two_pop,bc_switch_made)

%Newton-Raphson method for linear proliferation and death mechanisms.

w = 1; %Initial iterate
res = 1; %Initial residual (set to be above err_tol to start iterating)

%Set the values from the previous timestep
q_old =q; 
k_old=k;
a_old=a;
beta_old=beta;
gamma_old=gamma;
s_old = s;

%Spatial discretisation
x_discretisation=0:dx:((nodesx-1)*dx);

%% Newton-Raphson iteration

while w <= max_iters && res > err_tol
    
    %% Solve for q
    
    q_update = q; % pass solution
    
    [L_q,D_q,U_q] = function_jacobian_q_prolifdeathlinear(q,k,a,beta,gamma,l_crit_death,nodesx,dt,dx,eta);
    rhs_q = -function_discretised_q_prolifdeathlinear(q,q_old,k,a,beta,gamma,l_crit_death,nodesx,dt,dx,eta);
    
    % solve solution using thomas algorthim
    dq = tridia3(nodesx,L_q,D_q,U_q,rhs_q);
    q = q_update + dq;
    
    %% Calculate the node velocity
    
    [node_velocity,node_velocity_sign] = function_node_velocity(q,k,a,nodesx,dt,dx,eta);
    
    
    %% Solve for k
    
    if k_hom == 1
        dk=0;
    else
        k_update = k; % pass solution
        
        [L_k,D_k,U_k] =function_jacobian_k(q,k,a,nodesx,dt,dx,eta,node_velocity_sign);
        rhs_k = -function_discretised_func_k(q,k,k_old,a,nodesx,dt,dx,eta,node_velocity_sign);
        
        % solve solution using thomas algorthim
        dk = tridia3(nodesx,L_k,D_k,U_k,rhs_k);
        k = k_update + dk;
        
    end
      
    %% Solve for a
    
    if a_hom == 1
        da=0;
    else
        a_update = a; % pass solution
        
        [L_a,D_a,U_a] = function_jacobian_a(q,k,a,nodesx,dt,dx,eta,node_velocity_sign);
        rhs_a = -function_discretised_func_a(q,k,a,a_old,nodesx,dt,dx,eta,node_velocity_sign);
        
        % solve solution using thomas algorthim
        da = tridia3(nodesx,L_a,D_a,U_a,rhs_a);
        a = a_update + da;
        
    end
    
    
    %% Solve for beta
    
    if beta_hom == 1
        dbeta=0;
    else
        beta_update = beta; % pass solution

        [L_beta,D_beta,U_beta] = function_jacobian_beta(q,k,a,beta,nodesx,dt,dx,eta,node_velocity_sign);
        rhs_beta = -function_discretised_func_beta(q,k,a,beta,beta_old,nodesx,dt,dx,eta,node_velocity_sign);
        
        % solve solution using thomas algorthim
        dbeta = tridia3(nodesx,L_beta,D_beta,U_beta,rhs_beta);
        beta = beta_update + dbeta;
        
    end
    
    
    %% Solve for gamma
    
    if gamma_hom == 1
        dgamma=0;
    else
        gamma_update = gamma; % pass solution

        [L_gamma,D_gamma,U_gamma] = function_jacobian_gamma(q,k,a,gamma,nodesx,dt,dx,eta,node_velocity_sign);
        rhs_gamma = -function_discretised_func_gamma(q,k,a,gamma,gamma_old,nodesx,dt,dx,eta,node_velocity_sign);
        
        % solve solution using thomas algorthim
        dgamma = tridia3(nodesx,L_gamma,D_gamma,U_gamma,rhs_gamma);
        gamma = gamma_update + dgamma;
        
    end
    
    %% Interface position
    
   
    if one_or_two_pop == 2
        %If there are initially two populations
        if bc_switch_made == 0
            %If there are still 2 populations
            %Update the interface position with the velocity at the closest node
            %find i
            [~,snode] = min(abs(x_discretisation - s));
            [node_velocity,~] = function_node_velocity(q,k,a,nodesx,dt,dx,eta);
            ds = dt*node_velocity(snode);
            s = s_old + ds;
        elseif bc_switch_made == 1
            ds=0;
        end
    else
        ds=0;
    end
    
    
    %% Calculate residual
    
    res = norm([dq;dk;da;dbeta;dgamma;ds],inf);
    
    %% Increase iteration
    
    w = w + 1;
    
end

% If number of iterations is bigger than max iterations print this.
if w > max_iters
    fprintf('\nmax iters reached\n')
    res
end


end