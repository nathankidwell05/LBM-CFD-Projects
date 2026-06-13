function [Rho, U] = moment_rho_u_d2q9(f)

% create lattice model
e = [ 0  0;
      1  0;
      0  1;
     -1  0;
      0 -1;
      1  1;
     -1  1;
     -1 -1;
      1 -1 ];

% compute density
Rho = sum(f);

% velocity components in the X and Y direction
Ux = sum(f .* e(:,1)) / Rho;
Uy = sum(f .* e(:,2)) / Rho;

U = [Ux ; Uy];

end