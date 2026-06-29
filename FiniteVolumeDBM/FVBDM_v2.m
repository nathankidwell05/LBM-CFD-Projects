clear
clc
close all

%% Finite-volume discrete Boltzmann model: lid-driven cavity
% D2Q9, cell-centered finite volumes, first-order upwind face fluxes,
% BGK collision, and non-equilibrium extrapolation at solid walls.

% Physical and numerical parameters
Re = 100;
L = 1;
U_lid = 0.1;
rho_ref = 1;
N_x = 81;
N_y = N_x;
max_iterations = 100000;
report_interval = 1000;
convergence_tolerance = 1e-8;

c_s_sq = 1/3;
nu = U_lid*L/Re;
tau = nu/c_s_sq;
dx = L/N_x;
dy = L/N_y;
dt = min(0.8*tau,0.4*min(dx,dy));

Ksi = [0 1 0 -1 0 1 -1 -1 1; ...
       0 0 1 0 -1 1 1 -1 -1];
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]';

fprintf('Re = %g, grid = %d x %d, U_lid = %g\n',Re,N_x,N_y,U_lid)
fprintf('nu = %.4g, tau = %.4g, dt = %.4g, CFL = %.4g\n', ...
        nu,tau,dt,dt/min(dx,dy))

if U_lid/sqrt(c_s_sq) >= 0.3
    error('Lid Mach number must stay below 0.3 for this weakly compressible model.')
end
if dt/tau >= 1
    error('Choose dt/tau < 1 for a stable explicit BGK update.')
end

%% Initialization
rho = rho_ref*ones(1,N_y,N_x);
U = zeros(2,N_y,N_x);
f = equilibrium_field(rho,U,Ksi,w);
residual_history = nan(ceil(max_iterations/report_interval),1);
history_index = 0;
U_previous = U;

tic
for iteration = 1:max_iterations
    % Macroscopic fields and equilibrium populations.
    [rho,U] = moments_field(f,Ksi);
    f_eq = equilibrium_field(rho,U,Ksi,w);

    % Non-equilibrium extrapolation gives the wall population. The ghost
    % population follows from linear face interpolation: f_wall=(f+f_g)/2.
    U_top = zeros(2,1,N_x);
    U_top(1,:,:) = U_lid;
    f_wall_top = equilibrium_field(rho(:,1,:),U_top,Ksi,w) + ...
                 (f(:,1,:)-f_eq(:,1,:));
    f_ghost_top = 2*f_wall_top-f(:,1,:);

    U_bottom = zeros(2,1,N_x);
    f_wall_bottom = equilibrium_field(rho(:,end,:),U_bottom,Ksi,w) + ...
                    (f(:,end,:)-f_eq(:,end,:));
    f_ghost_bottom = 2*f_wall_bottom-f(:,end,:);

    U_left = zeros(2,N_y,1);
    f_wall_left = equilibrium_field(rho(:,:,1),U_left,Ksi,w) + ...
                  (f(:,:,1)-f_eq(:,:,1));
    f_ghost_left = 2*f_wall_left-f(:,:,1);

    U_right = zeros(2,N_y,1);
    f_wall_right = equilibrium_field(rho(:,:,end),U_right,Ksi,w) + ...
                   (f(:,:,end)-f_eq(:,:,end));
    f_ghost_right = 2*f_wall_right-f(:,:,end);

    % First-order upwind finite-volume flux divergence.
    divergence = zeros(size(f));
    for k = 1:9
        f_cell = squeeze(f(k,:,:));
        east = [f_cell(:,2:end),reshape(f_ghost_right(k,:,:),N_y,1)];
        west = [reshape(f_ghost_left(k,:,:),N_y,1),f_cell(:,1:end-1)];
        north = [reshape(f_ghost_top(k,:,:),1,N_x);f_cell(1:end-1,:)];
        south = [f_cell(2:end,:);reshape(f_ghost_bottom(k,:,:),1,N_x)];

        cx = Ksi(1,k);
        cy = Ksi(2,k);
        if cx >= 0
            flux_east = cx*f_cell;
        else
            flux_east = cx*east;
        end
        if cx <= 0
            flux_west = -cx*f_cell;
        else
            flux_west = -cx*west;
        end
        if cy >= 0
            flux_north = cy*f_cell;
        else
            flux_north = cy*north;
        end
        if cy <= 0
            flux_south = -cy*f_cell;
        else
            flux_south = -cy*south;
        end

        divergence(k,:,:) = (flux_east+flux_west)/dx + ...
                            (flux_north+flux_south)/dy;
    end

    f = f-dt*divergence-(dt/tau)*(f-f_eq);

    if any(~isfinite(f(:)))
        error('Solution became non-finite at iteration %d.',iteration)
    end

    if mod(iteration,report_interval)==0
        [~,U_check] = moments_field(f,Ksi);
        residual = max(abs(U_check(:)-U_previous(:)))/U_lid;
        history_index = history_index+1;
        residual_history(history_index) = residual;
        fprintf('iteration %6d: residual = %.3e\n',iteration,residual)
        U_previous = U_check;
        if residual < convergence_tolerance
            break
        end
    end
end
runtime = toc;

[rho,U] = moments_field(f,Ksi);
Ux = squeeze(U(1,:,:));
Uy = squeeze(U(2,:,:));
speed = hypot(Ux,Uy);
residual_history = residual_history(1:history_index);

fprintf('Completed %d iterations in %.2f seconds.\n',iteration,runtime)
fprintf('Speed range: [%.4g, %.4g]\n',min(speed(:)),max(speed(:)))
fprintf('Density range: [%.6g, %.6g]\n',min(rho(:)),max(rho(:)))

%% Figures
results_directory = fileparts(mfilename('fullpath'));
if ~isfolder(results_directory)
    mkdir(results_directory)
end

x = ((1:N_x)-0.5)*dx;
y = ((1:N_y)-0.5)*dy;
[X,Y] = meshgrid(x,y);
Y_plot = flipud(Y);

figure('Color','w')
contourf(X,Y_plot,speed/U_lid,30,'LineColor','none')
axis equal tight
colorbar
colormap(turbo)
xlabel('x/L')
ylabel('y/L')
title(sprintf('Normalized velocity magnitude, Re = %g',Re))
exportgraphics(gcf,fullfile(results_directory,'velocity_magnitude.png'),'Resolution',200)

figure('Color','w')
quiver(flipud(Ux),flipud(Uy),10)
axis equal tight
xlabel('x')
ylabel('y')
title(sprintf('Velocity Field, Re = %g',Re))
exportgraphics(gcf,fullfile(results_directory,'velocity_vectors.png'),'Resolution',200)

figure('Color','w')
contourf(X,Y_plot,squeeze(rho),30,'LineColor','none')
axis equal tight
colorbar
colormap(turbo)
xlabel('x/L')
ylabel('y/L')
title(sprintf('Density field, Re = %g',Re))
exportgraphics(gcf,fullfile(results_directory,'density_field.png'),'Resolution',200)

% Ghia et al. benchmark data at Re=100.
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

i_mid = round(N_x/2);
j_mid = round(N_y/2);
u_center = flipud(Ux(:,i_mid))/U_lid;
v_center = Uy(j_mid,:)/U_lid;

y_profile = [0,y,1];
u_profile = [0;u_center;1];
x_profile = [0,x,1];
v_profile = [0,v_center,0];
u_at_ghia = interp1(y_profile,u_profile,y_ghia,'linear');
v_at_ghia = interp1(x_profile,v_profile,x_ghia,'linear');
u_mean_absolute_error = mean(abs(u_at_ghia-u_ghia));
v_mean_absolute_error = mean(abs(v_at_ghia-v_ghia));
fprintf('Ghia centerline mean absolute error: u = %.4e, v = %.4e\n', ...
        u_mean_absolute_error,v_mean_absolute_error)

figure('Color','w')
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
nexttile
plot(u_center,y,'LineWidth',1.6)
hold on
plot(u_ghia,y_ghia,'o','MarkerSize',5)
grid on
xlabel('u/U_{lid}')
ylabel('y/L')
title('Vertical centerline')
legend('FVDBM','Ghia et al.','Location','best')
nexttile
plot(x,v_center,'LineWidth',1.6)
hold on
plot(x_ghia,v_ghia,'o','MarkerSize',5)
grid on
xlabel('x/L')
ylabel('v/U_{lid}')
title('Horizontal centerline')
legend('FVDBM','Ghia et al.','Location','best')
exportgraphics(gcf,fullfile(results_directory,'centerline_comparison.png'),'Resolution',200)

figure('Color','w')
semilogy((1:history_index)*report_interval,residual_history,'LineWidth',1.6)
grid on
xlabel('Iteration')
ylabel('Velocity residual')
title('Convergence history')
exportgraphics(gcf,fullfile(results_directory,'convergence_history.png'),'Resolution',200)

save(fullfile(results_directory,'FVBDM_results.mat'), ...
     'Re','L','U_lid','N_x','N_y','x','y','rho','U','residual_history','iteration', ...
     'u_mean_absolute_error','v_mean_absolute_error')

function f_eq = equilibrium_field(rho,U,Ksi,w)
    eu = reshape(Ksi(1,:),9,1,1).*U(1,:,:) + ...
         reshape(Ksi(2,:),9,1,1).*U(2,:,:);
    u_sq = U(1,:,:).^2+U(2,:,:).^2;
    f_eq = reshape(w,9,1,1).*rho.* ...
           (1+3*eu+4.5*eu.^2-1.5*u_sq);
end

function [rho,U] = moments_field(f,Ksi)
    rho = sum(f,1);
    momentum_x = sum(reshape(Ksi(1,:),9,1,1).*f,1);
    momentum_y = sum(reshape(Ksi(2,:),9,1,1).*f,1);
    U = cat(1,momentum_x./rho,momentum_y./rho);
end
