function g = function_discretised_func_k(q,k,k_old,a,nodesx,dt,dx,eta,node_velocity_sign)

% Discretised cell stiffness equation

g = zeros(nodesx,1);

r=dt/(eta*dx^2);

for i = 2:nodesx-1
	if node_velocity_sign(i) >=0 %use backwards difference 
		g(i) = - k(i) + k_old(i) ...
			-r/2*(1/q(i))*( (k(i)*(1/q(i)-a(i))) -  (k(i-1)*(1/q(i-1)-a(i-1))))*(k(i)-k(i-1)) ...
			-r/2*(1/q(i))*( (k_old(i)*(1/q(i)-a(i))) -  (k_old(i-1)*(1/q(i-1)-a(i-1))))*(k_old(i)-k_old(i-1));
	else %use forwards difference
	    g(i) = - k(i) + k_old(i) ...
			-r/2*(1/q(i))*( (k(i+1)*(1/q(i+1)-a(i+1))) -  (k(i)*(1/q(i)-a(i))))*(k(i+1)-k(i)) ...
			-r/2*(1/q(i))*( (k_old(i+1)*(1/q(i+1)-a(i+1))) -  (k_old(i)*(1/q(i)-a(i))))*(k_old(i+1)-k_old(i));
	end
end

g(1)=k(1) - k_old(1);
g(nodesx)=k(nodesx) - k_old(nodesx);

end