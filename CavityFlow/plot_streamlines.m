function plot_streamlines(U,N_x,N_y)

% Extract velocity components
Ux = squeeze(U(1,:,:));   % x velocity
Uy = squeeze(U(2,:,:));   % y velocity

% Create grid
[X,Y] = meshgrid(1:N_x,1:N_y);

% Velocity magnitude
vel_mag = sqrt(Ux.^2 + Uy.^2);

% Plot
figure
contourf(X,Y,vel_mag,30,'LineColor','none')
hold on
streamslice(X,Y,Ux,Uy)

axis equal tight
set(gca,'YDir','normal')
xlabel('x')
ylabel('y')
title('Velocity Streamlines')
colorbar

end