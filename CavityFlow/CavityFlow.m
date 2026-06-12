% row, column
clear all;
clc; 
%% Parameters
Ksi = [0 1 0 -1 0 1 -1 -1 1;...
       0 0 1 0 -1 1 1 -1 -1 ]; % D2Q9 Lattice
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]; % Weight
c_s = 1/sqrt(3); % Speeed of Sound
Re = 100; % Reynolds Number
N_x = 100;
N_y = N_x;
L = N_x-1;
Tau = 0.8; % Relaxation Time
nu = (Tau-0.5)/3; % Kinematic Viscoscity
U_lid = Re*(nu/L); % Speed of the moving lid
M =[1  1  1  1  1  1  1  1  1;
   -4 -1 -1 -1 -1  2  2  2  2;
    4 -2 -2 -2 -2  1  1  1  1;
    0  1  0 -1  0  1 -1 -1  1;
    0 -2  0  2  0  1 -1 -1  1;
    0  0  1  0 -1  1  1 -1 -1;
    0  0 -2  0  2  1  1 -1 -1;
    0  1 -1  1 -1  0  0  0  0;
    0  0  0  0  0  1 -1  1 -1];
m_eq=zeros(9,N_y,N_x);
s_e=1.19; %fudge factor
s_epsilon=1.4; %fudge factor
s_q=1.2; %fudge factor
s_v=1/Tau; %controls viscoscity
case_id = 'base';   % change this for each run
switch case_id

    case 'base'
        S = diag([0,s_e,s_epsilon,0,s_q,0,s_q,s_v,s_v]);

    case 'e'
        S = diag([0.5,s_e,s_epsilon,0,s_q,0,s_q,s_v,s_v]);

    case 'f'
        S = diag([0,s_e,s_epsilon,0.5,s_q,0,s_q,s_v,s_v]);

    case 'g'
        S = diag([0,s_e,s_epsilon,0,s_q,0.5,s_q,s_v,s_v]);

    case 'h'
        S = diag([0,s_e,s_epsilon,0,s_q,0,s_q,s_v,0.5*s_v]);

    case 'i'
        S = diag([0,0.8,s_epsilon,0,s_q,0,s_q,s_v,s_v]);

    case 'j'
        S = diag([0,s_e,0.8,0,s_q,0,s_q,s_v,s_v]);

    case 'k'
        S = diag([0,s_e,s_epsilon,0,0.8,0,0.8,s_v,s_v]);

    case 'l'
        S = diag([0,s_e,s_epsilon,0,1.2,0,1.8,s_v,s_v]);
    case 'z'
        S = diag([0,1.19,1.4,0,1.2,0,1.8,s_v,s_v]);
end

disp('Using S diagonal:')
disp(diag(S)')
M_inverse=inv(M);
F_C=1; % 0---SRT; 1---MRT
%% Initialization
Rho=ones(1,N_y,N_x);
U=zeros(2,N_y,N_x);
f=zeros(9,N_y,N_x);


% Alternative  initialization for PDF
for j=1:N_y
    for i=1:N_x
        f(:,j,i)=eqm_d2q9(Rho(:,j,i),U(:,j,i));
    end
end


%% Solving in a loop
T_max = 3000;
tic
for t=1:T_max
     f_new = zeros(9,N_y,N_x);
% Streaming + Boundary Conditions
for j = 1:N_y
    for i = 1:N_x
        if j==1 % this is the top boundary
            if i==1 % Top-Left corner nodes
                    % known
                    f_new(1,j,i)=f(1,j,i);
                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(7,j,i)=f(7,j+1,i+1);
                    % bounce back
                    f_new(2,j,i)=f_new(4,j,i);
                    f_new(5,j,i)=f_new(3,j,i);
                    f_new(9,j,i)=f_new(7,j,i);
                    % density extrapolation
                    Rho_TL = (Rho(:,j+1,i)+Rho(:,j,i+1))/2;
                    % resolve f6 from density summation
                    f_new(6,j,i)=(Rho_TL-f_new(1,j,i)-f_new(3,j,i)...
                        -f_new(4,j,i)-f_new(7,j,i)-f_new(2,j,i)...
                        -f_new(5,j,i)-f_new(9,j,i))/2;
                    f_new(8,j,i)=f_new(6,j,i);
                    

            elseif i==N_x % Top-Right corner node
                        % known
                        f_new(1,j,i)=f(1,j,i);
                        f_new(2,j,i)=f(2,j,i-1);
                        f_new(3,j,i)=f(3,j+1,i);
                        f_new(6,j,i)=f(6,j+1,i-1);
                        % bounce back
                        f_new(4,j,i)=f_new(2,j,i);
                        f_new(5,j,i)=f_new(3,j,i);
                        f_new(8,j,i)=f_new(6,j,i);
                        % density extrapolation
                        Rho_TR = (Rho(:,j+1,i)+Rho(:,j,i-1))/2;
                        % resolve f9 from density summation
                        f_new(9,j,i)=(Rho_TR-f_new(1,j,i)-f_new(2,j,i)...
                            -f_new(3,j,i)-f_new(6,j,i)-f_new(4,j,i)-...
                            f_new(5,j,i)-f_new(8,j,i))/2;
                        f_new(7,j,i)=f_new(9,j,i);



            else % all other nodes on top boundary 
                    f_new(1,j,i)=f(1,j,i);
                    f_new(2,j,i)=f(2,j,i-1);
                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(6,j,i)=f(6,j+1,i-1);
                    f_new(7,j,i)=f(7,j+1,i+1);
                    
                    % must resolve density first
                    Rho_top = f_new(1,j,i)+f_new(2,j,i)+f_new(4,j,i)+...
                              2*(f_new(3,j,i)+f_new(6,j,i)+f_new(7,j,i));

                    f_new(5,j,i)=f_new(3,j,i);
                    f_new(8,j,i)=f_new(6,j,i)-(Rho_top*U_lid+f_new(4,j,i)-f_new(2,j,i))/2;
                    f_new(9,j,i)=(Rho_top*U_lid+f_new(4,j,i)-f_new(2,j,i)+2*f_new(7,j,i))/2;
            end
            

        elseif j==N_y % this is the bottom boundary
            if i==1 % Bottom-Left corner node
                    % known
                    f_new(1,j,i)=f(1,j,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(5,j,i)=f(5,j-1,i);
                    f_new(8,j,i)=f(8,j-1,i+1);
                    % bounce back
                    f_new(2,j,i)=f_new(4,j,i);
                    f_new(3,j,i)=f_new(5,j,i);
                    f_new(6,j,i)=f_new(8,j,i);
                    % density extrapolation
                    Rho_BL = (Rho(:,j-1,i)+Rho(:,j,i+1))/2;
                    % resolve f7 from density summation
                    f_new(7,j,i)=(Rho_BL-f_new(1,j,i)-f_new(2,j,i)...
                        -f_new(4,j,i)-f_new(5,j,i)-f_new(6,j,i)-...
                        f_new(8,j,i)-f_new(3,j,i))/2;
                    f_new(9,j,i)=f_new(7,j,i);

            elseif i==N_x % Bottom-Right corner node
                    % known
                    f_new(1,j,i)=f(1,j,i);
                    f_new(2,j,i)=f(2,j,i-1);
                    f_new(5,j,i)=f(5,j-1,i);
                    f_new(9,j,i)=f(9,j-1,i-1);
                    % bounce back
                    f_new(3,j,i)=f_new(5,j,i);
                    f_new(4,j,i)=f_new(2,j,i);
                    f_new(7,j,i)=f_new(9,j,i);
                    % density extrapolation
                    Rho_BR = (Rho(:,j-1,i)+Rho(:,j,i-1))/2;
                    % resolve f6 from density summation
                    f_new(6,j,i)=(Rho_BR-f_new(1,j,i)-f_new(2,j,i)...
                        -f_new(4,j,i)-f_new(5,j,i)-f_new(7,j,i)-...
                        f_new(9,j,i)-f_new(3,j,i))/2;
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
                    
                    f_new(2,j,i)=f_new(4,j,i);
                    f_new(6,j,i)=f_new(8,j,i)+(f_new(5,j,i)-f_new(3,j,i))/2;
                    f_new(9,j,i)=f_new(7,j,i)-(f_new(5,j,i)-f_new(3,j,i))/2;
        elseif i==N_x % this is the right boundary
               f_new(1,j,i)=f(1,j,i);
               f_new(2,j,i)=f(2,j,i-1);
               f_new(3,j,i)=f(3,j+1,i);
               f_new(5,j,i)=f(5,j-1,i);
               f_new(6,j,i)=f(6,j+1,i-1);
               f_new(9,j,i)=f(9,j-1,i-1);
               
               f_new(4,j,i)=f_new(2,j,i);
               f_new(7,j,i)=f_new(9,j,i)+(f_new(5,j,i)-f_new(3,j,i))/2;
               f_new(8,j,i)=f_new(6,j,i)-(f_new(5,j,i)-f_new(3,j,i))/2;
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
    % f_eq calculation or m_eq depending on SRT or MRT
    if F_C==0
        for j=1:N_y
            for i=1:N_x
                f_eq(:,j,i) = eqm_d2q9(Rho(:,j,i),U(:,j,i));
            end
        end
    elseif F_C==1 % MRT
        m_eq(1,:,:)=Rho;
        m_eq(2,:,:)=Rho.*(-2+3*sum(U.*U,1));
        m_eq(3,:,:)=Rho.*(1-3*sum(U.*U,1));
        m_eq(4,:,:)=Rho.*U(1,:,:);
        m_eq(5,:,:)=-Rho.*U(1,:,:);
        m_eq(6,:,:)=Rho.*U(2,:,:);
        m_eq(7,:,:)=-Rho.*U(2,:,:);
        m_eq(8,:,:)=Rho.*(U(1,:,:).^2-U(2,:,:).^2);
        m_eq(9,:,:)=Rho.*U(1,:,:).*U(2,:,:);
    else
        error('Wrong Collison Model')
    end
    % Collision + Update PDF
    if F_C==0
        f = f_new - ((f_new-f_eq)/Tau);
    elseif F_C==1
        m = pagemtimes(M,f_new);
        f = pagemtimes(M_inverse, m - pagemtimes(S,(m-m_eq)));
    else
        error('Wrong Collision Model')
    end

    if any(isnan(f(:))) || any(isinf(f(:)))
        fprintf('Failed at Tau = %.4f, t = %d\n', Tau, t)
        break
    end
end

% final update for plotting/error
for j = 1:N_y
    for i = 1:N_x
        [Rho(:,j,i),U(:,j,i)] = moment_rho_u_d2q9(f(:,j,i));
    end
end

runtime = toc;
fprintf('Runtime: %.4f seconds\n', runtime)
runtime=toc;
fprintf('Runtime: %.4f seconds\n', runtime)
%% Post-Processing
Ux = squeeze(U(1,:,:));
Uy = squeeze(U(2,:,:));
plot_streamlines(U,N_x,N_y)
figure
quiver(flipud(Ux),flipud(Uy),10)
axis equal tight
title('Velocity Field')
xlabel('x')
ylabel('y')

figure
contourf(flipud(squeeze(Rho)),30)
axis equal tight
colorbar
title('Density Field')
xlabel('x')
ylabel('y')

%% Ghia Comparison at Re = 100
% normalize coordinates form 0-1
x_star = (0:N_x-1)/(N_x-1);
y_star = (N_y-1:-1:0)/(N_y-1);
% middle for vertical centerline
i_mid = round(N_x/2);
% middle for horizontal centerline
j_mid = round(N_y/2);
% normalize data using lid velo
u_sim = Ux(:,i_mid)/U_lid;
v_sim = Uy(j_mid,:)/U_lid;
% sort y-data for plotting
[y_sim_plot, idx] = sort(y_star);
u_sim_plot = u_sim(idx);

% Ghia data
y_ghia = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 ...
          0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];

u_ghia = [0 -0.03717 -0.04192 -0.04775 -0.06434 -0.1015 -0.15662 ...
          -0.2109 -0.20581 -0.13641 0.00332 0.23151 0.68717 ...
          0.73722 0.78871 0.84123 1];

x_ghia = [0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 ...
          0.5 0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1];

v_ghia = [0 0.09233 0.10091 0.1089 0.12317 0.16077 0.17507 ...
          0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 ...
          -0.08864 -0.07391 -0.05906 0];

figure
plot(u_sim_plot, y_sim_plot, 'LineWidth', 1.5)
hold on
plot(u_ghia, y_ghia, 'ro')
hold off
xlabel('u^*')
ylabel('y^*')
title('Vertical Centerline: y^* vs u^*')
legend('Simulation','Ghia Data')
grid on

figure
plot(x_star, v_sim, 'LineWidth', 1.5)
hold on
plot(x_ghia, v_ghia, 'ro')
hold off
xlabel('x^*')
ylabel('v^*')
title('Horizontal Centerline: v^* vs x^*')
legend('Simulation','Ghia Data')
grid on

%% Error Calculation

% get sim vals at the same locations as given data
u_sim_at_ghia = interp1(y_sim_plot, u_sim_plot, y_ghia);
v_sim_at_ghia = interp1(x_star, v_sim, x_ghia);

% avg abs error
u_error = mean(abs(u_sim_at_ghia - u_ghia));
v_error = mean(abs(v_sim_at_ghia - v_ghia));

fprintf('Average error for u*: %.6f\n', u_error)
fprintf('Average error for v*: %.6f\n', v_error)