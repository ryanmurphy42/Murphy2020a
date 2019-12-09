function F = function_p1d1_fsolve_initial_positions_gaussian_v1(N,L,gaussian_n,gaussian_sigma,gaussian_mu,   x )
%Determine the initial positions with Gaussian initial condition

F(1) = x(1) - 0;
for i=2:N
    F(i) = (1/(x(i) - x(i-1)) + 1/(x(i) - x(i+1)))/(x(i-1)/2 - x(i+1)/2) + ...
        (gaussian_n*(x(i)-gaussian_mu)/(sqrt(2*pi)*gaussian_sigma^3))*(exp(-((x(i)-gaussian_mu)^2)/(2*gaussian_sigma^2))) ;
end
F(N+1)= x(N+1) - L;



end

