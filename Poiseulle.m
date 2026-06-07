% row, column
clear all;
clc;
%% Parameters
% D2Q9 Lattice
Ksi = [0 1 0 -1 0 1 -1 -1 1;...
    0 0 1 0 -1 1 1 -1 -1 ];
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
c_s = 1/sqrt(3);
% Other LBM Related Parameters
Tau = 0.8;
Rho_in = 2;
N_x = 100;
N_y = 20;
%% Initialization
Rho=ones(1,N_y,N_x)*Rho_in;
U=zeros(2,N_y,N_x);
f=zeros(9,N_y,N_x);
% Alternative  initialization for PDF
for j=1:N_y
    for i=1:N_x
        f(:,j,i)=eqm_d2q9(Rho(:,j,i),U(:,j,i));
    end
end

f_new=f;
f_eq=f;


%% Solving in a loop
T_max = 1000;
for t=1:T_max
    % Streaming + Boundary Conditions
    for j = 1:N_y
        for i = 1:N_x
            if j==1 % this is the top boundary
                if i==1 % Top-Left corner node
                    f_new(1,j,i)=f(1,j,i);
                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(7,j,i)=f(7,j+1,i+1);

                    f_new(2,j,i)=f_new(4,j,i);
                    f_new(5,j,i)=f_new(3,j,i);
                    f_new(6,j,i)=Rho_in/2 - f_new(1,j,i)/2 - f_new(3,j,i)-f_new(4,j,i)-f_new(7,j,i);
                    f_new(8,j,i)=f_new(6,j,i);
                    f_new(9,j,i)=f_new(7,j,i);
                elseif i==N_x % Top-Right corner node
                    f_new(1,j,i)=f(1,j,i);
                    f_new(2,j,i)=f(2,j,i-1);
                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(6,j,i)=f(6,j+1,i-1);


                    f_new(4,j,i)=f_new(2,j,i);
                    f_new(5,j,i)=f_new(3,j,i);
                    f_new(7,j,i)=f_new(7,j,i-1);
                    f_new(8,j,i)=f_new(6,j,i);
                    f_new(9,j,i)=f_new(7,j,i);



                else % All other nodes on top boundary
                    f_new(1,j,i)=f(1,j,i);
                    f_new(2,j,i)=f(2,j,i-1);
                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(6,j,i)=f(6,j+1,i-1);
                    f_new(7,j,i)=f(7,j+1,i+1);


                    f_new(5,j,i)=f_new(3,j,i);
                    f_new(8,j,i)=f_new(6,j,i)+(f_new(2,j,i)-f_new(4,j,i))/2;
                    f_new(9,j,i)=f_new(7,j,i)-(f_new(2,j,i)-f_new(4,j,i))/2;
                end


            elseif j==N_y % This is the bottom boundary
                if i==1 % Bottom-Left corner node
                    f_new(1,j,i)=f(1,j,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(5,j,i)=f(5,j-1,i);
                    f_new(8,j,i)=f(8,j-1,i+1);


                    f_new(2,j,i)=f_new(4,j,i);
                    f_new(3,j,i)=f_new(5,j,i);
                    f_new(6,j,i)=f_new(8,j,i);
                    f_new(7,j,i)=Rho_in/2-f_new(1,j,i)/2-f_new(4,j,i)-f_new(5,j,i)-f_new(8,j,i);
                    f_new(9,j,i)=f_new(7,j,i);

                elseif i==N_x % Bottom-Right corner node
                    f_new(1,j,i)=f(1,j,i);
                    f_new(2,j,i)=f(2,j,i-1);
                    f_new(5,j,i)=f(5,j-1,i);
                    f_new(9,j,i)=f(9,j-1,i-1);

                    f_new(3,j,i)=f_new(5,j,i);
                    f_new(4,j,i)=f_new(2,j,i);
                    f_new(6,j,i)=f_new(6,j,i-1);
                    f_new(7,j,i)=f_new(9,j,i);
                    f_new(8,j,i)=f_new(6,j,i);


                else % all other nodes on bottom boundary
                    f_new(1,j,i)=f(1,j,i);
                    f_new(2,j,i)=f(2,j,i-1);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(5,j,i)=f(5,j-1,i);
                    f_new(8,j,i)=f(8,j-1,i+1);
                    f_new(9,j,i)=f(9,j-1,i-1);

                    f_new(3,j,i)=f_new(5,j,i);
                    f_new(6,j,i)=f_new(8,j,i)+(f_new(4,j,i)-f_new(2,j,i))/2;
                    f_new(7,j,i)=f_new(9,j,i)-(f_new(4,j,i)-f_new(2,j,i))/2;
                end

            elseif i==1 % this is the left boundary
                f_new(1,j,i)=f(1,j,i);

                f_new(3,j,i)=f(3,j+1,i);
                f_new(4,j,i)=f(4,j,i+1);
                f_new(5,j,i)=f(5,j-1,i);

                f_new(7,j,i)=f(7,j+1,i+1);
                f_new(8,j,i)=f(8,j-1,i+1);

                U_in=1-(f_new(1,j,i)+f_new(3,j,i)+f_new(5,j,i)...
                    +2*(f_new(4,j,i)+f_new(7,j,i)+f_new(8,j,i)))/Rho_in;
                f_new(2,j,i)=f_new(4,j,i)+Rho_in*U_in*(2/3);
                f_new(6,j,i)=f_new(8,j,i)+(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6;
                f_new(9,j,i)=f_new(7,j,i)-(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6;
            elseif i==N_x
                f_new(1,j,i)=f(1,j,i);
                f_new(2,j,i)=f(2,j,i-1);
                f_new(3,j,i)=f(3,j+1,i);
                f_new(5,j,i)=f(5,j-1,i);
                f_new(6,j,i)=f(6,j+1,i-1);
                f_new(9,j,i)=f(9,j-1,i-1);

                f_new(4,j,i)=f_new(4,j,i-1);
                f_new(7,j,i)=f_new(7,j,i-1);
                f_new(8,j,i)=f_new(8,j,i-1);
            else % this is all interior nodes
                f_new(1,j,i)=f(1,j,i);
                f_new(2,j,i)=f(2,j,i-1);
                f_new(3,j,i)=f(3,j+1,i);
                f_new(4,j,i)=f(4,j,i+1);
                f_new(5,j,i)=f(5,j-1,i);
                f_new(6,j,i)=f(6,j+1,i-1);
                f_new(7,j,i)=f(7,j+1,i+1);
                f_new(8,j,i)=f(8,j-1,i+1);
                f_new(9,j,i)=f(9,j-1,i-1);

            end
        end
    end
    % Collision
    % Moment Calculation
    for j=1:N_y
        for i=1:N_x
            [Rho(:,j,i),U(:,j,i)] = moment_rho_u_d2q9(f_new(:,j,i));
        end
    end
    % f_eq calculation
    for j=1:N_y
        for i=1:N_x
            f_eq(:,j,i) = eqm_d2q9(Rho(:,j,i),U(:,j,i));
        end
    end
    % Collision + Update PDF
    f=f_new-((f_new-f_eq)/Tau);
end
%% Post-Processing
figure
quiver(flipud(squeeze(U(1,:,:))),flipud(squeeze(U(2,:,:))),10)
axis equal tight

figure
contourf(flipud(squeeze(Rho)),30)
axis equal tight