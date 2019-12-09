function [ x ] = tridia3( nodesx, A,B,C, D )
    %tridia Thomas Algorithm for solving tridiagonal matrices problems
    %   Uses row eliminations as per Gaussian Elimination. The form of the
    %   matrix allows for the faster solving.

    % A is the sub diagonal of the tridiagonal matrix
    % B is the diagonal of the tridiagonal matrix
    % C is the superdiagonal of the tridiagonal matrix
    % D is the know RHS
    % x are the unknows which are solved for

    %% Read in the size of the inputs (only allows for n*n)
    rows=nodesx;

    %% Read in the sub-diagonal, diagonal and super-diagonal
    
    a=A;
    b=B;
    c=C;
    f=D;
    
    %% Initialise the variables to default to zero

    alpha = zeros(rows,1);
    beta = zeros(rows,1);
    y = zeros(rows,1);
    x = zeros(rows,1);

    %% Determine the intial values
    alpha(1) = b(1); 
    beta(1) = c(1)/alpha(1);
    y(1) = f(1)/alpha(1);

    %% Loop to determine the remaining values

    for i= 2:rows
       alpha(i) = b(i) - a(i)*beta(i-1);
       beta(i) = c(i)/alpha(i);
       y(i) = (f(i)-a(i)*y(i-1))/alpha(i);
    end

    %% Backwards substitution to determine the unknowns
    x(rows) = y(rows);

    for j = rows-1:-1:1
      x(j) = y(j) - beta(j)*x(j+1);
    end

end

