function [L,D,U] = function_jacobian_beta(q,k,a,beta,nodesx,dt,dx,eta,node_velocity_sign)

%Jacobian for beta

L = zeros(nodesx,1);
D = zeros(nodesx,1);
U = zeros(nodesx,1);

r = dt/(eta*dx^2);

for i = 2:nodesx-1
    if node_velocity_sign(i) >=0 %use backwards difference
        L(i) = -(r/2)*(1/q(i))*(( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))  )*(-1);
        D(i) = -1 -(r/2)*(1/q(i))*(( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))  )*(1);
    else %forwards difference
        D(i) = -1 -(r/2)*(1/q(i))*( ( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))  )*(-1);
        U(i) = -(r/2)*(1/q(i))*(( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))  )*(1);
    end
    
end

%% Boundary condition at x=0

i=1;
D(i)=1;



%% Boundary condition at x=L

i=nodesx;
D(i)=1;



end