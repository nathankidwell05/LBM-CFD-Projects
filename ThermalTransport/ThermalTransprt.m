% row, column
clear;
clc; 
%% Parameters
% D2Q9 Lattice
Ksi = [0 1 0 -1 0 1 -1 -1 1;...
       0 0 1 0 -1 1 1 -1 -1 ];
w = [0 1/6 1/6 1/6 1/6 1/12 1/12 1/12 1/12];
c_s = 1/sqrt(3);
% Other LBM Related Parameters 
Tau = 0.8; 
N_x = 101;
N_y = N_x;
%% Initialization
T_One=1/2;
T_Two=1/3;
T_Three=1;
T_inf=4/5;
k=1/8;
h=1;
R=8.314;
Rho=ones(1,N_y,N_x);
T=ones(1,N_y,N_x);
f=zeros(9,N_y,N_x);
% Circle Cylinder
R_cyl=round(0.2*N_x);
center_y = N_y/2;
center_x = N_x/2;
%% Domain ID
% Domain=1 ---- Fluid
% Domain=0 ---- Solid
Domain_ID=zeros(N_y,N_x);
for j=1:N_y
    for i=1:N_x
        if test_circle(i-1,j-1,R_cyl,center_x,center_y)
            Domain_ID(j,i)=0;
        else
            Domain_ID(j,i)=1;
        end
    end
end
contourf(Domain_ID,30)
axis equal tight
%% Zone ID
% Zone ID=0 ---Dead Zone
% Zone ID=1 ---Boundary Nodes
% Zone ID=2 ---Fluid Nodes
% Zone ID=3 ---All other nodes in fluid domain, including the nodes on the
% outer boundaries
Zone_ID=zeros(N_y,N_x);
for j=1:N_y
    for i=1:N_x
        if Domain_ID(j,i) == 0 % solid domain
            if Domain_ID(j,i+1)==1 || Domain_ID(j,i-1)==1 || Domain_ID(j+1,i)==1 || Domain_ID(j-1,i)==1|| Domain_ID(j+1,i+1)==1 || Domain_ID(j+1,i-1)==1 || Domain_ID(j-1,i+1)==1 || Domain_ID(j-1,i-1)==1
                Zone_ID(j,i)=1;
            else
                Zone_ID(j,i)=0;
            end
        else % fluid domain
            if j==1 || j==N_y || i==1 || i==N_x
                Zone_ID(j,i)=3;
            else
                if Domain_ID(j,i+1)==0 || Domain_ID(j,i-1)==0 || Domain_ID(j+1,i)==0 || Domain_ID(j-1,i)==0 || Domain_ID(j+1,i+1)==0 || Domain_ID(j+1,i-1)==0 || Domain_ID(j-1,i+1)==0 || Domain_ID(j-1,i-1)==0
                    
                    Zone_ID(j,i)=2; % fluid boundary near solid
                    
                else
                    Zone_ID(j,i)=3; % normal fluid node

                end
            end
        end
    end
end
contourf(Zone_ID,30)
axis equal tight

%% Initialization
% Alternative  initialization for PDF
for j=1:N_y
    for i= 1:N_x
        f(:,j,i) = w'* Rho(1,j,i) * R * T(1,j,i);
    end
end
f_new=f;
f_eq=f;
%% Solving in a loop
T_max = 20000;
for t=1:T_max
    f_new=f;
    T_bottom = (k*squeeze(T(1,N_y-1,:)) + h*T_inf) / (k + h);
    % Streaming + Boundary Conditions
    for j = 1:N_y
        for i = 1:N_x
            if Zone_ID(j,i)==0
                ;
            elseif Zone_ID(j,i)==1
                    ;
            elseif Zone_ID(j,i)==2 % Fluid Nodes, where the circular boundary scheme is implemented
                    %% Direction 1
                    f_new(1,j,i)=f(1,j,i);
                    %% Direction 2
                    if Zone_ID(j,i+1)==1
                        x1=i; % boundary node x
                        y1=j-1; % boundary node y
                        x2=i-1; % fluid node x
                        y2=j-1; % fluid node y
                        C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        f_eq_alpha_f = w(2)*Rho(1,j,i)*R*T(1,j,i);
                        f_neq_f=f(2,j,i) - f_eq_alpha_f;
                        if delta >= 0.25
                            f_neq = f_neq_f;
                        else
                            f_eq_alpha_ff = w(2)*Rho(1,j,i-1)*R*T(1,j,i-1);
                            f_neq_ff = f(2,j,i-1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                            f_eq_abar = w(4)*Rho(1,j,i)*R*T_Three;
                            f_new(4,j,i)=f_eq_abar + f_neq;
                    else
                        f_new(4,j,i)=f(4,j,i+1);
                    end
            
                     %% Direction 3
                     if Zone_ID(j-1,i)==1
                         x1=i-1;
                         y1=j-2;
                         x2=i-1;
                         y2=j-1;
                         C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                         delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                         f_eq_alpha_f = w(3)*Rho(1,j,i)*R*T(1,j,i);
                         f_neq_f=f(3,j,i) - f_eq_alpha_f;
                         if delta >= 0.25
                             f_neq = f_neq_f;
                         else
                             f_eq_alpha_ff = w(3)*Rho(1,j+1,i)*R*T(1,j+1,i);
                             f_neq_ff = f(3,j+1,i) - f_eq_alpha_ff;
                             f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                         end
                         f_eq_abar = w(5)*Rho(1,j,i)*R*T_Three;
                         f_new(5,j,i)=f_eq_abar + f_neq;
                     else
                         f_new(5,j,i)=f(5,j-1,i);
                     end

                    %% Direction 4
                    if Zone_ID(j,i-1)==1
                        x1=i-2;
                        y1=j-1;
                        x2=i-1;
                        y2=j-1;
                        C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        f_eq_alpha_f = w(4)*Rho(1,j,i)*R*T(1,j,i);
                        f_neq_f=f(4,j,i) - f_eq_alpha_f;
                        if delta >= 0.25
                            f_neq = f_neq_f;
                        else
                            f_eq_alpha_ff = w(4)*Rho(1,j,i+1)*R*T(1,j,i+1);
                            f_neq_ff = f(4,j,i+1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(2)*Rho(1,j,i)*R*T_Three;
                        f_new(2,j,i)=f_eq_abar + f_neq;
                    else
                        f_new(2,j,i)=f(2,j,i-1);
                    end

                    %% Direction 5
                    if Zone_ID(j+1,i)==1
                        x1=i-1;
                        y1=j;
                        x2=i-1;
                        y2=j-1;
                        C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        f_eq_alpha_f = w(5)*Rho(1,j,i)*R*T(1,j,i);
                        f_neq_f=f(5,j,i) - f_eq_alpha_f;
                        if delta >= 0.25
                            f_neq = f_neq_f;
                        else
                            f_eq_alpha_ff = w(5)*Rho(1,j-1,i)*R*T(1,j-1,i);
                            f_neq_ff = f(5,j-1,i) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(3)*Rho(1,j,i)*R*T_Three;
                        f_new(3,j,i)=f_eq_abar + f_neq;
                    else
                        f_new(3,j,i)=f(3,j+1,i);
                    end
                    %% Direction 6
                    if Zone_ID(j-1,i+1)==1
                        x1=i;
                        y1=j-2;
                        x2=i-1;
                        y2=j-1;
                        C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        f_eq_alpha_f = w(6)*Rho(1,j,i)*R*T(1,j,i);
                        f_neq_f=f(6,j,i) - f_eq_alpha_f;
                        if delta >= 0.25
                            f_neq = f_neq_f;
                        else
                            f_eq_alpha_ff = w(6)*Rho(1,j+1,i-1)*R*T(1,j+1,i-1);
                            f_neq_ff = f(6,j+1,i-1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(8)*Rho(1,j,i)*R*T_Three;
                        f_new(8,j,i)=f_eq_abar + f_neq;
                    else
                        f_new(8,j,i)=f(8,j-1,i+1);
                    end
                     %% Direction 7
                    if Zone_ID(j-1,i-1)==1
                        x1=i-2;
                        y1=j-2;
                        x2=i-1;
                        y2=j-1;
                        C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        f_eq_alpha_f = w(7)*Rho(1,j,i)*R*T(1,j,i);
                        f_neq_f=f(7,j,i) - f_eq_alpha_f;
                        if delta >= 0.25
                            f_neq = f_neq_f;
                        else
                            f_eq_alpha_ff = w(7)*Rho(1,j+1,i+1)*R*T(1,j+1,i+1);
                            f_neq_ff = f(7,j+1,i+1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(9)*Rho(1,j,i)*R*T_Three;
                        f_new(9,j,i)=f_eq_abar + f_neq;
                    else
                        f_new(9,j,i)=f(9,j-1,i-1);
                    end
                      %% Direction 8
                    if Zone_ID(j+1,i-1)==1
                        x1=i-2;
                        y1=j;
                        x2=i-1;
                        y2=j-1;
                        C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        f_eq_alpha_f = w(8)*Rho(1,j,i)*R*T(1,j,i);
                        f_neq_f=f(8,j,i) - f_eq_alpha_f;
                        if delta >= 0.25
                            f_neq = f_neq_f;
                        else
                            f_eq_alpha_ff = w(8)*Rho(1,j-1,i+1)*R*T(1,j-1,i+1);
                            f_neq_ff = f(8,j-1,i+1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(6)*Rho(1,j,i)*R*T_Three;
                        f_new(6,j,i)=f_eq_abar + f_neq;
                    else
                        f_new(6,j,i)=f(6,j+1,i-1);
                    end
                       %% Direction 9
                    if Zone_ID(j+1,i+1)==1
                        x1=i;
                        y1=j;
                        x2=i-1;
                        y2=j-1;
                        C_w=find_the_wall_point(x1,y1,x2,y2,R_cyl,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        f_eq_alpha_f = w(9)*Rho(1,j,i)*R*T(1,j,i);
                        f_neq_f=f(9,j,i) - f_eq_alpha_f;
                        if delta >= 0.25
                            f_neq = f_neq_f;
                        else
                            f_eq_alpha_ff = w(9)*Rho(1,j-1,i-1)*R*T(1,j-1,i-1);
                            f_neq_ff = f(9,j-1,i-1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(7)*Rho(1,j,i)*R*T_Three;
                        f_new(7,j,i)=f_eq_abar + f_neq;
                    else
                        f_new(7,j,i)=f(7,j+1,i+1);
                    end
            else % Zone ID = 3
                if j==1 % this is the top boundary
                    if i==1 % Top-Left corner node
                        f_new(1,j,i)=f(1,j,i);
                        f_new(3,j,i)=f(3,j+1,i);
                        f_new(4,j,i)=f(4,j,i+1);
                        f_new(7,j,i)=f(7,j+1,i+1);
                        T_W=(T_One+T_Two)/2;
                        f_new(2,j,i)=f_new(4,j,i);
                        f_new(5,j,i)=f_new(3,j,i);
                        f_new(6,j,i)=(Rho(1,j,i)*R*T_W-f_new(1,j,i)-2*f_new(3,j,i)-2*f_new(4,j,i)-2*f_new(7,j,i))/2;
                        f_new(8,j,i)=f_new(6,j,i);
                        f_new(9,j,i)=f_new(7,j,i);
                    elseif i==N_x % Top-Right corner node
                        f_new(1,j,i)=f(1,j,i);
                        f_new(2,j,i)=f(2,j,i-1);
                        f_new(3,j,i)=f(3,j+1,i);
                        f_new(6,j,i)=f(6,j+1,i-1);
                        T_W=T_One;
                        f_new(4,j,i)=f_new(2,j,i);
                        f_new(5,j,i)=f_new(3,j,i);
                        f_new(8,j,i)=f_new(6,j,i);
                        f_new(7,j,i)=(Rho(1,j,i)*R*T_W-f_new(1,j,i)-2*f_new(2,j,i)-2*f_new(3,j,i)-2*f_new(6,j,i))/2;
                        f_new(9,j,i)=f_new(7,j,i);
                    else % this is the top boundary
                        f_new(1,j,i)=f(1,j,i);
                        f_new(2,j,i)=f(2,j,i-1);
                        f_new(3,j,i)=f(3,j+1,i);
                        f_new(4,j,i)=f(4,j,i+1);
                        f_new(6,j,i)=f(6,j+1,i-1);
                        f_new(7,j,i)=f(7,j+1,i+1);

                        T_W=T_One;
                        f_new(5,j,i)=f_new(3,j,i);
                        f_new(8,j,i)=(Rho(1,j,i)*R*T_W-f_new(1,j,i)-f_new(2,j,i)-2*f_new(3,j,i)-f_new(4,j,i)-f_new(6,j,i)-f_new(7,j,i))/2;
                        f_new(9,j,i)=f_new(8,j,i);
                    end

                elseif j==N_y % this is the bottom boundary
                    if i==1 % Bottom-Left corner node
                        f_new(1,j,i)=f(1,j,i);
                        f_new(4,j,i)=f(4,j,i+1);
                        f_new(5,j,i)=f(5,j-1,i);
                        f_new(8,j,i)=f(8,j-1,i+1);
                        T_W=T_Two;
                        f_new(2,j,i)=f_new(4,j,i);
                        f_new(3,j,i)=f_new(5,j,i);
                        f_new(6,j,i)=f_new(8,j,i);
                        f_new(7,j,i)=(Rho(1,j,i)*R*T_W-f_new(1,j,i)-2*f_new(4,j,i)-2*f_new(5,j,i)-2*f_new(8,j,i))/2;
                        f_new(9,j,i)=f_new(7,j,i);

                    elseif i==N_x % Bottom-Right corner node
                        f_new(1,j,i)=f(1,j,i);
                        f_new(2,j,i)=f(2,j,i-1);
                        f_new(5,j,i)=f(5,j-1,i);
                        f_new(9,j,i)=f(9,j-1,i-1);
                        T_W=T_bottom(i);
                        f_new(3,j,i)=f_new(5,j,i);
                        f_new(4,j,i)=f_new(2,j,i);
                        f_new(7,j,i)=f_new(9,j,i);
                        f_new(6,j,i)=(Rho(1,j,i)*R*T_W-f_new(1,j,i)-2*f_new(2,j,i)-2*f_new(5,j,i)-2*f_new(9,j,i))/2;
                        f_new(8,j,i)=f_new(6,j,i);


                    else % all other nodes on bottom boundary
                        f_new(1,j,i)=f(1,j,i);
                        f_new(2,j,i)=f(2,j,i-1);
                        f_new(4,j,i)=f(4,j,i+1);
                        f_new(5,j,i)=f(5,j-1,i);
                        f_new(8,j,i)=f(8,j-1,i+1);
                        f_new(9,j,i)=f(9,j-1,i-1);
                        T_W=T_bottom(i);
                        f_new(3,j,i)=-f_new(5,j,i)+2*w(3)*Rho(1,j,i)*R*T_W;
                        f_new(6,j,i)=-f_new(8,j,i)+2*w(6)*Rho(1,j,i)*R*T_W;
                        f_new(7,j,i)=-f_new(9,j,i)+2*w(7)*Rho(1,j,i)*R*T_W;
                    end

                elseif i==1 % this is the left boundary
                    f_new(1,j,i)=f(1,j,i);
                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(5,j,i)=f(5,j-1,i);
                    f_new(7,j,i)=f(7,j+1,i+1);
                    f_new(8,j,i)=f(8,j-1,i+1);

                    T_W=T_Two;
                    f_new(2,j,i)=f_new(4,j,i);
                    f_new(6,j,i)=(Rho(1,j,i)*R*T_W-f_new(1,j,i)-2*f_new(4,j,i)-f_new(5,j,i)-f_new(7,j,i)-f_new(8,j,i)-f_new(3,j,i))/2;
                    f_new(9,j,i)=f_new(6,j,i);
             
                elseif i==N_x % this is the right boundary
                    f_new(1,j,i)=f(1,j,i);
                    f_new(2,j,i)=f(2,j,i-1);
                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(5,j,i)=f(5,j-1,i);
                    f_new(6,j,i)=f(6,j+1,i-1);
                    f_new(9,j,i)=f(9,j-1,i-1);

                    % insulated wall: dT/dx = 0, so T_wall = nearby interior temperature
                    T_W = T(1,j,i-1);

                    f_new(4,j,i)=f_new(2,j,i);

                    f_new(7,j,i)= ...
                        (Rho(1,j,i)*R*T_W ...
                        - f_new(1,j,i) ...
                        - 2*f_new(2,j,i) ...
                        - f_new(3,j,i) ...
                        - f_new(5,j,i) ...
                        - f_new(6,j,i) ...
                        - f_new(9,j,i))/2;

                    f_new(8,j,i)=f_new(7,j,i);
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
    end


            
    
    % Collision - MRT

    T = sum(f_new,1)./R./Rho;

    for j = 1:N_y
        for i = 1:N_x
            f_eq(:,j,i) = w' * Rho(1,j,i) * R * T(1,j,i);
        end
    end

    M = [1  1  1  1  1  1  1  1  1;
        -4 -1 -1 -1 -1  2  2  2  2;
        4 -2 -2 -2 -2  1  1  1  1;
        0  1  0 -1  0  1 -1 -1  1;
        0 -2  0  2  0  1 -1 -1  1;
        0  0  1  0 -1  1  1 -1 -1;
        0  0 -2  0  2  1  1 -1 -1;
        0  1 -1  1 -1  0  0  0  0;
        0  0  0  0  0  1 -1  1 -1];

    Minv = inv(M);

    sT = 1/Tau;

    S = diag([0, sT, sT, 1.0, 1.0, 1.0, 1.2, 1.2, 1.0]);

    for j = 1:N_y
        for i = 1:N_x
            if Zone_ID(j,i) ~= 0 && Zone_ID(j,i) ~= 1

                m = M * f_new(:,j,i);
                m_eq = M * f_eq(:,j,i);

                m_post = m - S * (m - m_eq);

                f(:,j,i) = Minv * m_post;

            end
        end
    end
end
%% Post-Processing

% Temperature contour
figure
contourf(flipud(squeeze(T(1,:,:))),30)
colorbar
axis equal tight
title('Temperature Distribution')

%% Load benchmark data
load('Project3_Benchmark Data.mat')

%% Horizontal mid-plane
j_mid = round(N_y/2);

T_hori = squeeze(T(1,j_mid,:))';
x_star = (0:N_x-1)/(N_x-1);

% Put hole region at T_Three
for i = 1:N_x
    if Domain_ID(j_mid,i)==0 
        T_hori(i) = T_Three;
    end
end

T_star_hori = (T_hori - T_Two)/(T_Three - T_Two);

%% Vertical mid-plane
i_mid = round(N_x/2);

T_vert = squeeze(T(1,:,i_mid));
y_star_vert = 1 - (0:N_y-1)/(N_y-1);

% Put hole region at T_Three
for j = 1:N_y
    if Domain_ID(j,i_mid)==0 || Zone_ID(j,i_mid)==1
        T_vert(j) = T_Three;
    end
end

T_star_vert = (T_vert - T_Two)/(T_Three - T_Two);

% Sort y data so plot goes from y*=0 to y*=1
[y_star_vert, idx] = sort(y_star_vert);
T_star_vert = T_star_vert(idx);

%% Interpolate simulation onto benchmark points
T_sim_hori_interp = interp1(x_star, T_star_hori, x_benchmark, 'linear', 'extrap');
T_sim_vert_interp = interp1(y_star_vert, T_star_vert, y_benchmark, 'linear', 'extrap');

%% Plot horizontal mid-plane
figure
plot(x_star, T_star_hori, 'b-', 'LineWidth', 1.5)
hold on
plot(x_benchmark, T_benchmark_hori, 'ro', 'MarkerSize', 5)
xlabel('x^*')
ylabel('T^*')
legend('Simulation', 'Benchmark')
title('Temperature on Horizontal Mid-Plane')
grid on

%% Plot vertical mid-plane
figure
plot(y_star_vert, T_star_vert, 'b-', 'LineWidth', 1.5)
hold on
plot(y_benchmark, T_benchmark_vert, 'ro', 'MarkerSize', 5)
xlabel('y^*')
ylabel('T^*')
legend('Simulation', 'Benchmark')
title('Temperature on Vertical Mid-Plane')
grid on

%% L2 error
L2_hori = sqrt(mean((T_sim_hori_interp(:) - T_benchmark_hori(:)).^2));
L2_vert = sqrt(mean((T_sim_vert_interp(:) - T_benchmark_vert(:)).^2));

fprintf('\n----- Results -----\n');
fprintf('Tau = %.3f\n', Tau);
fprintf('N_x = %d\n', N_x);
fprintf('L2 error horizontal = %.6f\n', L2_hori);
fprintf('L2 error vertical   = %.6f\n', L2_vert);