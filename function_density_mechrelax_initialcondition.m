function [ x_ss ] = function_density_mechrelax_initialcondition( N, L, k, a )
%   Cell positions at mechanical relaxation of N cells, in fixed domain of length L, with cell stiffness vector k and resting cell length vector a

x_ss = zeros(N+1,1);

x_ss(1)=0;

sum_a = sum(a);
sum_1_k = sum(1./k);

for ii = 2:N
    x_ss(ii) =  x_ss(ii-1) + (a(ii-1) + (L - sum_a)/(k(ii-1)*sum_1_k));
end

x_ss(N+1) = L;

end

