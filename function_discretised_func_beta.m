function g = function_discretised_func_beta(q,k,a,beta,beta_old,nodesx,dt,dx,eta,node_velocity_sign)

% Discretised beta equation

g = zeros(nodesx,1);

r=dt/(eta*dx^2);

for i = 2:nodesx-1
	if node_velocity_sign(i) >=0 %use backwards difference 
		g(i) = - beta(i) + beta_old(i) ...
			-r/2*(1/q(i))*( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))*(beta(i)-beta(i-1)) ...
			-r/2*(1/q(i))*( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))*(beta_old(i)-beta_old(i-1));
	else %use forwards difference
	    g(i) = - beta(i) + beta_old(i) ...
			-r/2*(1/q(i))*( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))*(beta(i+1)-beta(i)) ...
			-r/2*(1/q(i))*( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))*(beta_old(i+1)-beta_old(i));
	end
end

%% Boundary condition at x=0
%fixed bc
g(1)=beta(1) - beta_old(1);



%% Boundary condition at x=L

%fixed bc
g(nodesx)=beta(nodesx) - beta_old(nodesx);



end