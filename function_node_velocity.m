function [node_velocity,node_velocity_sign] = function_node_velocity(q,k,a,nodesx,dt,dx,eta)


%% Node velocity internal nodes
node_velocity = zeros(nodesx,1);
node_velocity_sign = zeros(nodesx,1);
for i=2:(nodesx-1)
    node_velocity(i) = (1/dx)*(1/q(i))*(1/eta)*0.5*(     k(i+1)*((1/(q(i+1)))-a(i+1) ) -  k(i-1)*((1/(q(i-1)))-a(i-1) ) ) ;
    node_velocity_sign(i) = sign(node_velocity(i));
end

%% Node velocity first node
i=1;
node_velocity(i) = (1/dx)*(1/q(i))*(1/eta)*(     k(i+1)*((1/(q(i+1)))-a(i+1) ) -  k(i)*((1/(q(i)))-a(i) ) ) ;
node_velocity_sign(i) = sign(node_velocity(i));

%% Node velocity final node
i=nodesx;
node_velocity(i) = (1/dx)*(1/q(i))*(1/eta)*(     k(i)*((1/(q(i)))-a(i) ) -  k(i-1)*((1/(q(i-1)))-a(i-1) ) ) ;
node_velocity_sign(i) = sign(node_velocity(i));

end