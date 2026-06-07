function f_eq = eqm_d2q9(Rho,U)

% establish the D2Q9 lattice
e = [ 0  0;
      1  0;
      0  1;
     -1  0;
      0 -1;
      1  1;
     -1  1;
     -1 -1;
      1 -1 ];

% weights for each lattice direction
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]';

% compute e_i dot U for all 9 directions at once
e_dot_u = e*U;

% compute U dot U once
u_squared = U(1)^2 + U(2)^2;

% compute all 9 equilibrium PDFs at once
f_eq = w*Rho .* (1 + 3*e_dot_u + (9/2)*e_dot_u.^2 - (3/2)*u_squared);

end
