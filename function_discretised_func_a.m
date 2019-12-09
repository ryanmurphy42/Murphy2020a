function g = function_discretised_func_a(q,k,a,a_old,nodes,dt,dx,eta,node_velocity_sign)

% Discretised resting cell length equation

g = zeros(nodes,1);

r=dt/(eta*dx^2);

for i = 2:nodes-1
	if node_velocity_sign(i) >=0 %use backwards difference 
		g(i) = - a(i) + a_old(i) ...
			-r/2*(1/q(i))*( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))*(a(i)-a(i-1)) ...
			-r/2*(1/q(i))*( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))*(a_old(i)-a_old(i-1));
	else %use forwards difference
	    g(i) = - a(i) + a_old(i) ...
			-r/2*(1/q(i))*( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))*(a(i+1)-a(i)) ...
			-r/2*(1/q(i))*( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))*(a_old(i+1)-a_old(i));
	end
end

%% Boundary condition at x=0
%fixed bc
g(1)=a(1) - a_old(1);

%% Boundary condition at x=L

%fixed bc
g(nodes)=a(nodes) - a_old(nodes);


end