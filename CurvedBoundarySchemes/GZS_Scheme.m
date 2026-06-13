clear;
clc; 
%% Parameters
% D2Q9 Lattice
Ksi = [0 1 0 -1 0 1 -1 -1 1;...
       0 0 1 0 -1 1 1 -1 -1 ];
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
c_s = 1/sqrt(3);
% Other LBM Related Parameters 
Tau = 0.60; 
Rho_in = 2;
R=7;
D=2*R;
N_x = 35*D;
N_y = 9*D;
% Circle Cylinder
center_y = N_y/2;
center_x = 5*D;
%% Re related parameters
Re = 40;
nu = c_s^2*(Tau - 0.5);
U_in = Re*nu/(2*R);
Ma = U_in/c_s;
%% Domain ID
% Domain=1 ---- Fluid
% Domain=0 ---- Solid
Domain_ID=zeros(N_y,N_x);
for j=1:N_y
    for i=1:N_x
        if test_circle(i-1,j-1,R,center_x,center_y)
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
T_max = 1500;
for t=1:T_max
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
                        C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        U_f=squeeze(U(:,j,i));
                        U_ff=squeeze(U(:,j,i-1));
                        U_b1=((delta-1)*U_f)/(delta);
                        f_eq_alpha_f=w(2)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,2)'*U_f)/c_s^2 + ((Ksi(:,2)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                        f_neq_f=f(2,j,i) - f_eq_alpha_f;
                        if delta >= 0.75
                            U_b = U_b1;
                            f_neq = f_neq_f;
                            
                        else
                            U_b2=((delta-1)*U_ff)/(1+delta);
                            U_b = delta*U_b1+(1-delta)*U_b2;
                            f_eq_alpha_ff = w(2)*squeeze(Rho(:,j,i-1))*(1 + (Ksi(:,2)'*U_ff)/c_s^2 + ((Ksi(:,2)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                            f_neq_ff = f(2,j,i-1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                            f_eq_abar = w(4)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,4)'*U_b)/c_s^2 + ((Ksi(:,4)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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
                         C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                         delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                         U_f=squeeze(U(:,j,i));
                         U_ff=squeeze(U(:,j+1,i));
                         U_b1=((delta-1)*U_f)/(delta);
                         f_eq_alpha_f=w(3)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,3)'*U_f)/c_s^2 + ((Ksi(:,3)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                         f_neq_f=f(3,j,i) - f_eq_alpha_f;
                         if delta >= 0.75
                             U_b = U_b1;
                             f_neq = f_neq_f;

                         else
                             U_b2=((delta-1)*U_ff)/(1+delta);
                             U_b = delta*U_b1+(1-delta)*U_b2;
                             f_eq_alpha_ff = w(3)*squeeze(Rho(:,j+1,i))*(1 + (Ksi(:,3)'*U_ff)/c_s^2 + ((Ksi(:,3)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                             f_neq_ff = f(3,j+1,i) - f_eq_alpha_ff;
                             f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                         end
                         f_eq_abar = w(5)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,5)'*U_b)/c_s^2 + ((Ksi(:,5)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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
                        C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        U_f=squeeze(U(:,j,i));
                        U_ff=squeeze(U(:,j,i+1));
                        U_b1=((delta-1)*U_f)/(delta);
                        f_eq_alpha_f=w(4)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,4)'*U_f)/c_s^2 + ((Ksi(:,4)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                        f_neq_f=f(4,j,i) - f_eq_alpha_f;
                        if delta >= 0.75
                            U_b = U_b1;
                            f_neq = f_neq_f;

                        else
                            U_b2=((delta-1)*U_ff)/(1+delta);
                            U_b = delta*U_b1+(1-delta)*U_b2;
                            f_eq_alpha_ff = w(4)*squeeze(Rho(:,j,i+1))*(1 + (Ksi(:,4)'*U_ff)/c_s^2 + ((Ksi(:,4)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                            f_neq_ff = f(4,j,i+1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(2)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,2)'*U_b)/c_s^2 + ((Ksi(:,2)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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
                        C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        U_f=squeeze(U(:,j,i));
                        U_ff=squeeze(U(:,j-1,i));
                        U_b1=((delta-1)*U_f)/(delta);
                        f_eq_alpha_f=w(5)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,5)'*U_f)/c_s^2 + ((Ksi(:,5)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                        f_neq_f=f(5,j,i) - f_eq_alpha_f;
                        if delta >= 0.75
                            U_b = U_b1;
                            f_neq = f_neq_f;

                        else
                            U_b2=((delta-1)*U_ff)/(1+delta);
                            U_b = delta*U_b1+(1-delta)*U_b2;
                            f_eq_alpha_ff = w(5)*squeeze(Rho(:,j-1,i))*(1 + (Ksi(:,5)'*U_ff)/c_s^2 + ((Ksi(:,5)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                            f_neq_ff = f(5,j-1,i) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(3)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,3)'*U_b)/c_s^2 + ((Ksi(:,3)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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
                        C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        U_f=squeeze(U(:,j,i));
                        U_ff=squeeze(U(:,j+1,i-1));
                        U_b1=((delta-1)*U_f)/(delta);
                        f_eq_alpha_f=w(6)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,6)'*U_f)/c_s^2 + ((Ksi(:,6)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                        f_neq_f=f(6,j,i) - f_eq_alpha_f;
                        if delta >= 0.75
                            U_b = U_b1;
                            f_neq = f_neq_f;

                        else
                            U_b2=((delta-1)*U_ff)/(1+delta);
                            U_b = delta*U_b1+(1-delta)*U_b2;
                            f_eq_alpha_ff = w(6)*squeeze(Rho(:,j+1,i-1))*(1 + (Ksi(:,6)'*U_ff)/c_s^2 + ((Ksi(:,6)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                            f_neq_ff = f(6,j+1,i-1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(8)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,8)'*U_b)/c_s^2 + ((Ksi(:,8)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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
                        C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        U_f=squeeze(U(:,j,i));
                        U_ff=squeeze(U(:,j+1,i+1));
                        U_b1=((delta-1)*U_f)/(delta);
                        f_eq_alpha_f=w(7)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,7)'*U_f)/c_s^2 + ((Ksi(:,7)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                        f_neq_f=f(7,j,i) - f_eq_alpha_f;
                        if delta >= 0.75
                            U_b = U_b1;
                            f_neq = f_neq_f;

                        else
                            U_b2=((delta-1)*U_ff)/(1+delta);
                            U_b = delta*U_b1+(1-delta)*U_b2;
                            f_eq_alpha_ff = w(7)*squeeze(Rho(:,j+1,i+1))*(1 + (Ksi(:,7)'*U_ff)/c_s^2 + ((Ksi(:,7)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                            f_neq_ff = f(7,j+1,i+1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(9)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,9)'*U_b)/c_s^2 + ((Ksi(:,9)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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
                        C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        U_f=squeeze(U(:,j,i));
                        U_ff=squeeze(U(:,j-1,i+1));
                        U_b1=((delta-1)*U_f)/(delta);
                        f_eq_alpha_f=w(8)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,8)'*U_f)/c_s^2 + ((Ksi(:,8)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                        f_neq_f=f(8,j,i) - f_eq_alpha_f;
                        if delta >= 0.75
                            U_b = U_b1;
                            f_neq = f_neq_f;

                        else
                            U_b2=((delta-1)*U_ff)/(1+delta);
                            U_b = delta*U_b1+(1-delta)*U_b2;
                            f_eq_alpha_ff = w(8)*squeeze(Rho(:,j-1,i+1))*(1 + (Ksi(:,8)'*U_ff)/c_s^2 + ((Ksi(:,8)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                            f_neq_ff = f(8,j-1,i+1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(6)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,6)'*U_b)/c_s^2 + ((Ksi(:,6)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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
                        C_w=find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y);
                        delta=sqrt((C_w(1)-x2)^2+(C_w(2)-y2)^2)/sqrt((x1-x2)^2+(y1-y2)^2);
                        U_f=squeeze(U(:,j,i));
                        U_ff=squeeze(U(:,j-1,i-1));
                        U_b1=((delta-1)*U_f)/(delta);
                        f_eq_alpha_f=w(9)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,9)'*U_f)/c_s^2 + ((Ksi(:,9)'*U_f)^2)/(2*c_s^4) - (U_f'*U_f)/(2*c_s^2));
                        f_neq_f=f(9,j,i) - f_eq_alpha_f;
                        if delta >= 0.75
                            U_b = U_b1;
                            f_neq = f_neq_f;

                        else
                            U_b2=((delta-1)*U_ff)/(1+delta);
                            U_b = delta*U_b1+(1-delta)*U_b2;
                            f_eq_alpha_ff = w(9)*squeeze(Rho(:,j-1,i-1))*(1 + (Ksi(:,9)'*U_ff)/c_s^2 + ((Ksi(:,9)'*U_ff)^2)/(2*c_s^4) - (U_ff'*U_ff)/(2*c_s^2));
                            f_neq_ff = f(9,j-1,i-1) - f_eq_alpha_ff;
                            f_neq = delta*f_neq_f+(1-delta)*f_neq_ff;
                        end
                        f_eq_abar = w(7)*squeeze(Rho(:,j,i))*(1 + (Ksi(:,7)'*U_b)/c_s^2 + ((Ksi(:,7)'*U_b)^2)/(2*c_s^4) - (U_b'*U_b)/(2*c_s^2));
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



                    else % all other nodes on top boundary: periodic

                        f_new(1,j,i)=f(1,j,i);

                        f_new(2,j,i)=f(2,j,i-1);
                        f_new(3,j,i)=f(3,j+1,i);
                        f_new(4,j,i)=f(4,j,i+1);

                        % Periodic wrap: j-1 becomes N_y
                        f_new(5,j,i)=f(5,N_y-1,i);

                        f_new(6,j,i)=f(6,j+1,i-1);
                        f_new(7,j,i)=f(7,j+1,i+1);

                        % Periodic wrap: j-1 becomes N_y
                        f_new(8,j,i)=f(8,N_y-1,i+1);
                        f_new(9,j,i)=f(9,N_y-1,i-1);
                    end

                elseif j==N_y % this is the bottom boundary
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


                    else % all other nodes on bottom boundary: periodic

                        f_new(1,j,i)=f(1,j,i);

                        f_new(2,j,i)=f(2,j,i-1);

                        % Periodic wrap: j+1 becomes 1
                        f_new(3,j,i)=f(3,2,i);

                        f_new(4,j,i)=f(4,j,i+1);
                        f_new(5,j,i)=f(5,j-1,i);

                        % Periodic wrap: j+1 becomes 1
                        f_new(6,j,i)=f(6,2,i-1);
                        f_new(7,j,i)=f(7,2,i+1);

                        f_new(8,j,i)=f(8,j-1,i+1);
                        f_new(9,j,i)=f(9,j-1,i-1);
                    end

                elseif i==1 % this is the left boundary
                    f_new(1,j,i)=f(1,j,i);

                    f_new(3,j,i)=f(3,j+1,i);
                    f_new(4,j,i)=f(4,j,i+1);
                    f_new(5,j,i)=f(5,j-1,i);

                    f_new(7,j,i)=f(7,j+1,i+1);
                    f_new(8,j,i)=f(8,j-1,i+1);

                    U_in = Re*nu/(2*R);

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
    end


            
    
% Collision
    % Moment Calculation
    % Improvement 1: There is no need to calculate the moment near a dead
    % zone or boundary node - no meaningful density of velocity to be
    % computed.
    for j = 1:N_y
        for i = 1:N_x
            if Zone_ID(j,i) ~= 0 && Zone_ID(j,i) ~= 1
                [Rho(:,j,i),U(:,j,i)] = moment_rho_u_d2q9(f_new(:,j,i));
            end
        end
    end
    % Equilibrium Calculation
    % Improvement 2: Same reasoning as moment calculation
    for j = 1:N_y
        for i = 1:N_x
            if Zone_ID(j,i) ~= 0 && Zone_ID(j,i) ~= 1
                f_eq(:,j,i) = eqm_d2q9(Rho(:,j,i), U(:,j,i));
                f_new(:,j,i) = f_new(:,j,i) - (1/Tau) * (f_new(:,j,i) - f_eq(:,j,i));
            end
        end
    end

f = f_new;  % update f for next timestep
end
%% Post-Processing
figure
quiver(flipud(squeeze(U(1,:,:))),flipud(squeeze(U(2,:,:))),10)
axis equal tight

figure
contourf(flipud(squeeze(Rho)),30)
axis equal tight
%% More most processing 
%% L/R and theta_s
%% Wake length L/R
j_center = round(center_y) + 1;
u_center = squeeze(U(1,j_center,:));

x_B = center_x + R;
i_back = round(x_B) + 1;

x_C = NaN;

for i = i_back:N_x-1
    if u_center(i) < 0 && u_center(i+1) >= 0

        x1 = i - 1;
        x2 = i;

        u1 = u_center(i);
        u2 = u_center(i+1);

        x_C = x1 - u1*(x2 - x1)/(u2 - u1);
        break
    end
end

if isnan(x_C)
    L = NaN;
    L_over_R = NaN;
else
    L = x_C - x_B;
    L_over_R = L/R;
end

dL = 0.5;
dL_over_R = dL/R;
percent_error_L = abs(dL_over_R/L_over_R)*100;


%% Separation angle theta_s

theta_list = [];
utheta_list = [];

for j = 2:N_y-1
    for i = 2:N_x-1

        if Zone_ID(j,i) == 2

            x = i - 1;
            y = j - 1;

            x_rel = x - center_x;
            y_rel = center_y - y;

            theta = atan2d(y_rel,x_rel);

            if theta > 0 && theta < 180

                ux = U(1,j,i);
                uy = -U(2,j,i);

                utheta = -ux*sind(theta) + uy*cosd(theta);

                theta_list = [theta_list; theta];
                utheta_list = [utheta_list; utheta];

            end
        end
    end
end

[theta_sort, idx] = sort(theta_list);
utheta_sort = utheta_list(idx);

theta_s = NaN;

for k = 1:length(theta_sort)-1
    if utheta_sort(k)*utheta_sort(k+1) < 0

        theta1 = theta_sort(k);
        theta2 = theta_sort(k+1);

        u1 = utheta_sort(k);
        u2 = utheta_sort(k+1);

        theta_s = theta1 - u1*(theta2 - theta1)/(u2 - u1);
        break
    end
end

dtheta = atan2d(0.5,R);
percent_error_theta = abs(dtheta/theta_s)*100;


%% Print results

fprintf('\n----- Measurement Results -----\n');
fprintf('Re = %.3f\n', Re);
fprintf('Tau = %.3f\n', Tau);
fprintf('nu = %.3f\n', nu);
fprintf('U_in = %.3f\n', U_in);
fprintf('Ma = %.3f\n', Ma);
fprintf('N_x = %d, N_y = %d\n', N_x, N_y);
fprintf('R = %.3f\n', R);

fprintf('\nWake length:\n');
fprintf('x_B = %.3f\n', x_B);
fprintf('x_C = %.3f\n', x_C);
fprintf('L = %.3f\n', L);
fprintf('L/R = %.3f\n', L_over_R);
fprintf('Estimated L/R measurement uncertainty = %.3f\n', dL_over_R);
fprintf('Estimated L/R percent measurement error = %.3f %%\n', percent_error_L);

fprintf('\nSeparation angle:\n');
fprintf('theta_s = %.3f degrees\n', theta_s);
fprintf('Estimated theta uncertainty = %.3f degrees\n', dtheta);
fprintf('Estimated theta percent measurement error = %.3f %%\n', percent_error_theta);