function g = function_discretised_q_prolifdeathlinear(q,q_old,k,a,beta, gamma, l_crit_death,nodesx,dt,dx,eta)

% Discretised cell density equation - linear proliferation and death mechanisms

r = dt/(eta*dx^2);

g = zeros(nodesx,1);

for i = 2:nodesx-1
    if (1/q(i)) < l_crit_death
        g(i) = - q(i) + q_old(i) ...
            -r/2*(k(i-1)*(1/q(i-1)-a(i-1)) - 2*k(i)*(1/q(i)-a(i)) + k(i+1)*(1/q(i+1)-a(i+1)))...
            -r/2*(k(i-1)*(1/q_old(i-1)-a(i-1)) - 2*k(i)*(1/q_old(i)-a(i)) + k(i+1)*(1/q_old(i+1)-a(i+1))) ...
            + (dt/2)*( (beta(i)+gamma(i) - gamma(i)*q(i)*l_crit_death) + (beta(i)+gamma(i) - gamma(i)*q_old(i)*l_crit_death)   ); 
    else
        g(i) = - q(i) + q_old(i) ...
            -r/2*(k(i-1)*(1/q(i-1)-a(i-1)) - 2*k(i)*(1/q(i)-a(i)) + k(i+1)*(1/q(i+1)-a(i+1)))...
            -r/2*(k(i-1)*(1/q_old(i-1)-a(i-1)) - 2*k(i)*(1/q_old(i)-a(i)) + k(i+1)*(1/q_old(i+1)-a(i+1))) ...
            + (dt/2)*((beta(i)) + (beta(i))  ); 
    end
end

g(1) = k(2)*(1/q(2)-a(2)) - k(1)*(1/q(1)-a(1));
g(end) = k(end)*(1/q(end)-a(end)) - k(end-1)*(1/q(end-1) - a(end-1));


end