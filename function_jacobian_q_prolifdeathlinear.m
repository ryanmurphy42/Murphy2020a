function [L,D,U] = function_jacobian_q_prolifdeathlinear(q,k,a,beta,gamma,l_crit_death,nodes,dt,dx,eta)

%Jacobian for cell density with linear proliferation and death mechanism

L = zeros(nodes,1);
D = zeros(nodes,1);
U = zeros(nodes,1);

r = dt/(eta*dx^2);


for i = 2:nodes-1
    if (1/q(i)) < l_crit_death
        L(i) = r/2*k(i-1)*1/q(i-1)^2;
        D(i) = -1 - r*k(i)*1/q(i)^2 - (dt/2)*gamma(i)*l_crit_death ;
        U(i) = r/2*k(i+1)*1/q(i+1)^2;
    else
        L(i) = r/2*k(i-1)*1/q(i-1)^2;
        D(i) = -1 - r*k(i)*1/q(i)^2;
        U(i) = r/2*k(i+1)*1/q(i+1)^2;
    end
end

D(1) = k(1)/q(1)^2;
U(1) = -k(2)/q(2)^2;
D(end) = -k(end)/q(end)^2;
L(end) = k(end-1)/q(end-1)^2;

end