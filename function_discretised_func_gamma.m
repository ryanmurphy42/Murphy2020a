function g = function_discretised_func_gamma(q,k,a,gamma,gamma_old,nodesx,dt,dx,eta,node_velocity_sign)

% Discretised gamma equation

g = zeros(nodesx,1);

r=dt/(eta*dx^2);

for i = 2:nodesx-1
	if node_velocity_sign(i) >=0 %use backwards difference 
		g(i) = - gamma(i) + gamma_old(i) ...
			-r/2*(1/q(i))*( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))*(gamma(i)-gamma(i-1)) ...
			-r/2*(1/q(i))*( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))*(gamma_old(i)-gamma_old(i-1));
	else %use forwards difference
	    g(i) = - gamma(i) + gamma_old(i) ...
			-r/2*(1/q(i))*( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))*(gamma(i+1)-gamma(i)) ...
			-r/2*(1/q(i))*( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))*(gamma_old(i+1)-gamma_old(i));
	end
end

%% Boundary condition at x=0
%fixed bc
g(1)=gamma(1) - gamma_old(1);



%% Boundary condition at x=L

%fixed bc
g(nodesx)=gamma(nodesx) - gamma_old(nodesx);



end