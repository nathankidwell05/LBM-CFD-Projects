clear all 
clc

%% Parameter Setup
% Domain/Geometry
L=100;
H=L;
N_nd_x=101;
N_nd_y=N_nd_x;
N_cv_x=N_nd_x-1;
N_cv_y=N_nd_y-1;
dx=L/N_cv_x;
dy=H/N_cv_y;
A=dx*dy;
n1=[1;0];
n2=[0;-1];
n3=[-1;0];
n4=[0;1];
% DBM Related
Ksi = [0 1 0 -1 0 1 -1 -1 1;...
       0 0 1 0 -1 1 1 -1 -1 ]; % Lattice Velocities
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
c_s = 1/sqrt(3); % Speed of Sound
Tau=0.4; % Relaxation Time
Re=100;
U_lid=Re*Tau*c_s^2/L;
dt=0.1;
%% Initialization
Rho_ref=2;
Rho_nd=ones(1,N_nd_y,N_nd_x)*Rho_ref;
Rho_cv=ones(1,N_cv_y,N_cv_x)*Rho_ref;
U_nd=zeros(2,N_nd_y,N_nd_x);
U_cv=zeros(2,N_cv_y,N_cv_x);

f_nd=ones(9,N_nd_y,N_nd_x);
for j=1:N_nd_y
    for i=1:N_nd_x
       f_nd(:,j,i) = eqm_d2q9(squeeze(Rho_nd(1,j,i)),squeeze(U_nd(:,j,i)));
    end
end

f_cv=ones(9,N_cv_y,N_cv_x);
for j=1:N_cv_y
    for i=1:N_cv_x
       f_cv(:,j,i) = eqm_d2q9(squeeze(Rho_cv(1,j,i)),squeeze(U_cv(:,j,i)));
    end
end

f_eq_nd=f_nd;
f_eq_cv=f_cv;
f_neq_nd=f_nd;
f_new_cv=f_cv;

F=zeros(9,N_cv_y,N_cv_x);
%% Solving
Timer=50000;
report_interval=1000;
residual_history=nan(1,ceil(Timer/report_interval));
history_index=0;
U_previous=U_cv;
tic
for t=1:Timer
% Nodal PDF on the Boundaries
for j=1:N_nd_y
    for i=1:N_nd_x
        if j==1 % Top Boundary
            if i==1 % Top Left Corner
            %f_eq
            Rho_nd(1,j,i)=Rho_cv(1,1,1);
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[U_lid;0]);
            %f_neq
            f_cv_corner=f_cv(:,1,1);
           
            f_eq_cv_corner=f_eq_cv(:,1,1);
           
            f_neq_node=(f_cv_corner-f_eq_cv_corner);

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;
            elseif i==N_nd_x% Top Right Corner
            %f_eq
            Rho_nd(1,j,i)=Rho_cv(1,1,end);
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[U_lid;0]);
            %f_neq
            f_cv_corner=f_cv(:,1,end);
           
            f_eq_cv_corner=f_eq_cv(:,1,end);
           
            f_neq_node=(f_cv_corner-f_eq_cv_corner);

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;
            else % All other nodes on the top boundary
            %f_eq
            Rho_nd(1,j,i)=(Rho_cv(1,1,i-1)+Rho_cv(1,1,i))/2;
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[U_lid;0]);
            %f_neq
            f_cv_1=f_cv(:,1,i-1);
            f_cv_2=f_cv(:,1,i);

            f_eq_cv_1=f_eq_cv(:,1,i-1);
            f_eq_cv_2=f_eq_cv(:,1,i);
            f_neq_node=((f_cv_1-f_eq_cv_1)+(f_cv_2-f_eq_cv_2))/2;

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;
            end
        elseif j==N_nd_y % Bottom Boundary
            if i==1 % Bottom Left Corner
            %f_eq
            Rho_nd(1,j,i)=Rho_cv(1,end,1);
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[0;0]);
            %f_neq
            f_cv_corner=f_cv(:,end,1);
           
            f_eq_cv_corner=f_eq_cv(:,end,1);
           
            f_neq_node=(f_cv_corner-f_eq_cv_corner);

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;
            elseif i==N_nd_x% Bottom Right Corner
            %f_eq
            Rho_nd(1,j,i)=Rho_cv(1,end,end);
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[0;0]);
            %f_neq
            f_cv_corner=f_cv(:,end,end);
           
            f_eq_cv_corner=f_eq_cv(:,end,end);
           
            f_neq_node=(f_cv_corner-f_eq_cv_corner);

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;
            else % All other nodes on the bottom boundary
            %f_eq
            Rho_nd(1,j,i)=(Rho_cv(1,end,i-1)+Rho_cv(1,end,i))/2;
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[0;0]);
            %f_neq
            f_cv_1=f_cv(:,end,i-1);
            f_cv_2=f_cv(:,end,i);

            f_eq_cv_1=f_eq_cv(:,end,i-1);
            f_eq_cv_2=f_eq_cv(:,end,i);
            f_neq_node=((f_cv_1-f_eq_cv_1)+(f_cv_2-f_eq_cv_2))/2;

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;
            end
        elseif i==1 % Left Boundary
            %f_eq
            Rho_nd(1,j,i)=(Rho_cv(1,j,1)+Rho_cv(1,j-1,1))/2;
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[0;0]);
            %f_neq
            f_cv_1=f_cv(:,j-1,1);
            f_cv_2=f_cv(:,j,1);

            f_eq_cv_1=f_eq_cv(:,j-1,1);
            f_eq_cv_2=f_eq_cv(:,j,1);
            f_neq_node=((f_cv_1-f_eq_cv_1)+(f_cv_2-f_eq_cv_2))/2;

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;
        elseif i==N_nd_x % Right Boundary
            %f_eq
            Rho_nd(1,j,i)=(Rho_cv(1,j,end)+Rho_cv(1,j-1,end))/2;
            f_eq_nd(:,j,i)=eqm_d2q9(Rho_nd(1,j,i),[0;0]);
            %f_neq
            f_cv_1=f_cv(:,j-1,end);
            f_cv_2=f_cv(:,j,end);

            f_eq_cv_1=f_eq_cv(:,j-1,end);
            f_eq_cv_2=f_eq_cv(:,j,end);
            f_neq_node=((f_cv_1-f_eq_cv_1)+(f_cv_2-f_eq_cv_2))/2;

            % f=f_eq+f_neq
            f_nd(:,j,i)=f_eq_nd(:,j,i)+f_neq_node;

        else % All Interior Nodes
            % Do Nothing
        end
    end
end
% Flux
for j=1:N_cv_y
    for i=1:N_cv_x
        if j==1 % Top Boundary
            if i==1 % Top Left Corner
                %% Flux on face 1 - EAST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n1>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i+1);
                    else
                        f_u(k,1)=f_cv(k,j,i+1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc1=f_u;
                F1=fc1.*(Ksi'*n1)*dy;
            %% Flux on face 2 - SOUTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n2>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j+1,i);
                    else
                        f_u(k,1)=f_cv(k,j+1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc2=f_u;
                F2=fc2.*(Ksi'*n2)*dx;
            %% Flux on face 3 - WEST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,1,1);
                f_nd2=f_nd(:,2,1);
                f_ghost=f_nd1+f_nd2-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n3>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc3=f_u;
                F3=fc3.*(Ksi'*n3)*dy;
            %% Flux on face 4 - NORTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,j,i);
                f_nd2=f_nd(:,j,i+1);
                f_ghost=f_nd1+f_nd2-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n4>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc4=f_u;
                F4=fc4.*(Ksi'*n4)*dx;
                % Total Flux
                F(:,j,i)=(F1+F2+F3+F4)/A;
            elseif i==N_cv_x % Top Right Corner
                 %% Flux on face 1 - EAST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,1,end);
                f_nd2=f_nd(:,2,end);
                f_ghost=f_nd1+f_nd2-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n1>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc1=f_u;
                F1=fc1.*(Ksi'*n1)*dy;
            %% Flux on face 2 - SOUTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n2>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j+1,i);
                    else
                        f_u(k,1)=f_cv(k,j+1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc2=f_u;
                F2=fc2.*(Ksi'*n2)*dx;
            %% Flux on face 3 - WEST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n3>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i-1);
                    else
                        f_u(k,1)=f_cv(k,j,i-1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc3=f_u;
                F3=fc3.*(Ksi'*n3)*dy;
            %% Flux on face 4 - NORTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,j,i);
                f_nd2=f_nd(:,j,i+1);
                f_ghost=f_nd1+f_nd2-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n4>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc4=f_u;
                F4=fc4.*(Ksi'*n4)*dx;
                % Total Flux
                F(:,j,i)=(F1+F2+F3+F4)/A;
            else % All other CVs on the Top
                %% Flux on face 1 - EAST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n1>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i+1);
                    else
                        f_u(k,1)=f_cv(k,j,i+1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc1=f_u;
                F1=fc1.*(Ksi'*n1)*dy;
            %% Flux on face 2 - SOUTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n2>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j+1,i);
                    else
                        f_u(k,1)=f_cv(k,j+1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc2=f_u;
                F2=fc2.*(Ksi'*n2)*dx;
            %% Flux on face 3 - WEST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n3>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i-1);
                    else
                        f_u(k,1)=f_cv(k,j,i-1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc3=f_u;
                F3=fc3.*(Ksi'*n3)*dy;
            %% Flux on face 4 - NORTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,j,i);
                f_nd2=f_nd(:,j,i+1);
                f_ghost=f_nd1+f_nd2-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n4>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc4=f_u;
                F4=fc4.*(Ksi'*n4)*dx;
                % Total Flux
                F(:,j,i)=(F1+F2+F3+F4)/A;
            end
        elseif j==N_cv_y % Bottom Boundary
            if i==1 % Bottom Left Corner
                 %% Flux on face 1 - EAST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n1>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i+1);
                    else
                        f_u(k,1)=f_cv(k,j,i+1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc1=f_u;
                F1=fc1.*(Ksi'*n1)*dy;
                %% Flux on face 2 - SOUTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,end,i);
                f_nd2=f_nd(:,end,i+1);
                f_ghost=(f_nd1+f_nd2)-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n2>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc2=f_u;
                F2=fc2.*(Ksi'*n2)*dx;
                %% Flux on face 3 - WEST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,end,1);
                f_nd2=f_nd(:,end-1,1);
                f_ghost=(f_nd1+f_nd2)-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n3>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc3=f_u;
                F3=fc3.*(Ksi'*n3)*dy;
                %% Flux on face 4 - NORTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n4>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j-1,i);
                    else
                        f_u(k,1)=f_cv(k,j-1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc4=f_u;
                F4=fc4.*(Ksi'*n4)*dx;
                % Total Flux
                F(:,j,i)=(F1+F2+F3+F4)/A;
            elseif i==N_cv_x % Bottom Right Corner
                %% Flux on face 1 - EAST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,end-1,end);
                f_nd2=f_nd(:,end,end);
                f_ghost=(f_nd1+f_nd2)-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n1>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc1=f_u;
                F1=fc1.*(Ksi'*n1)*dy;
                %% Flux on face 2 - SOUTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,end,i);
                f_nd2=f_nd(:,end,i+1);
                f_ghost=(f_nd1+f_nd2)-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n2>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc2=f_u;
                F2=fc2.*(Ksi'*n2)*dx;
                %% Flux on face 3 - WEST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n3>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i-1);
                    else
                        f_u(k,1)=f_cv(k,j,i-1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc3=f_u;
                F3=fc3.*(Ksi'*n3)*dy;
                %% Flux on face 4 - NORTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n4>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j-1,i);
                    else
                        f_u(k,1)=f_cv(k,j-1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                    % Flux Scheme
                    fc4=f_u;
                    F4=fc4.*(Ksi'*n4)*dx;
                    % Total Flux
                    F(:,j,i)=(F1+F2+F3+F4)/A;
                end
            else % All other CVs on the Bottom
                %% Flux on face 1 - EAST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n1>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i+1);
                    else
                        f_u(k,1)=f_cv(k,j,i+1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc1=f_u;
                F1=fc1.*(Ksi'*n1)*dy;
                %% Flux on face 2 - SOUTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                f_nd1=f_nd(:,end,i);
                f_nd2=f_nd(:,end,i+1);
                f_ghost=(f_nd1+f_nd2)-f_cv(:,j,i);
                for k=1:9
                    if Ksi(:,k)'*n2>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_ghost(k,1,1);
                    else
                        f_u(k,1)=f_ghost(k,1,1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc2=f_u;
                F2=fc2.*(Ksi'*n2)*dx;
                %% Flux on face 3 - WEST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n3>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i-1);
                    else
                        f_u(k,1)=f_cv(k,j,i-1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc3=f_u;
                F3=fc3.*(Ksi'*n3)*dy;
                %% Flux on face 4 - NORTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n4>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j-1,i);
                    else
                        f_u(k,1)=f_cv(k,j-1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc4=f_u;
                F4=fc4.*(Ksi'*n4)*dx;
                % Total Flux
                F(:,j,i)=(F1+F2+F3+F4)/A;
            end

        elseif i==1 % Left Boundary
             %% Flux on face 1 - East
             f_u=zeros(9,1);
             f_d=zeros(9,1);
            for k=1:9
                if Ksi(:,k)'*n1>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_cv(k,j,i+1);
                else
                    f_u(k,1)=f_cv(k,j,i+1);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc1=f_u;
            F1=fc1.*(Ksi'*n1)*dy;
            %% Flux on face 2 - SOUTH
            f_u=zeros(9,1);
            f_d=zeros(9,1);
            for k=1:9
                if Ksi(:,k)'*n2>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_cv(k,j+1,i);
                else
                    f_u(k,1)=f_cv(k,j+1,i);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc2=f_u;
            F2=fc2.*(Ksi'*n2)*dx;
            %% Flux on face 3 - WEST
            f_u=zeros(9,1);
            f_d=zeros(9,1);
            f_nd1=f_nd(:,j,1);
            f_nd2=f_nd(:,j+1,1);
            f_ghost=(f_nd1+f_nd2)-f_cv(:,j,i);
            for k=1:9
                if Ksi(:,k)'*n3>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_ghost(k,1,1);
                else
                    f_u(k,1)=f_ghost(k,1,1);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc3=f_u;
            F3=fc3.*(Ksi'*n3)*dy;
            %% Flux on face 4 - NORTH
            f_u=zeros(9,1);
            f_d=zeros(9,1);
            for k=1:9
                if Ksi(:,k)'*n4>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_cv(k,j-1,i);
                else
                    f_u(k,1)=f_cv(k,j-1,i);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc4=f_u;
            F4=fc4.*(Ksi'*n4)*dx;
            % Total Flux
            F(:,j,i)=(F1+F2+F3+F4)/A;
        elseif i==N_cv_x % Right Boundary
            %% Flux on face 1 - East
            f_u=zeros(9,1);
            f_d=zeros(9,1);
            f_nd1=f_nd(:,j,end);
            f_nd2=f_nd(:,j+1,end);
            f_ghost=(f_nd1+f_nd2)-f_cv(:,j,i);
            for k=1:9
                if Ksi(:,k)'*n1>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_ghost(k,1,1);
                else
                    f_u(k,1)=f_ghost(k,1,1);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc1=f_u;
            F1=fc1.*(Ksi'*n1)*dy;
            %% Flux on face 2 - SOUTH
            f_u=zeros(9,1);
            f_d=zeros(9,1);
            for k=1:9
                if Ksi(:,k)'*n2>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_cv(k,j+1,i);
                else
                    f_u(k,1)=f_cv(k,j+1,i);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc2=f_u;
            F2=fc2.*(Ksi'*n2)*dx;
            %% Flux on face 3 - WEST
            f_u=zeros(9,1);
            f_d=zeros(9,1);
            for k=1:9
                if Ksi(:,k)'*n3>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_cv(k,j,i-1);
                else
                    f_u(k,1)=f_cv(k,j,i-1);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc3=f_u;
            F3=fc3.*(Ksi'*n3)*dy;
            %% Flux on face 4 - NORTH
            f_u=zeros(9,1);
            f_d=zeros(9,1);
            for k=1:9
                if Ksi(:,k)'*n4>=0
                    f_u(k,1)=f_cv(k,j,i);
                    f_d(k,1)=f_cv(k,j-1,i);
                else
                    f_u(k,1)=f_cv(k,j-1,i);
                    f_d(k,1)=f_cv(k,j,i);
                end
            end
            % Flux Scheme
            fc4=f_u;
            F4=fc4.*(Ksi'*n4)*dx;
            % Total Flux
            F(:,j,i)=(F1+F2+F3+F4)/A;
        else % All other CVs internally
            %% Flux on face 1 - EAST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n1>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i+1);
                    else
                        f_u(k,1)=f_cv(k,j,i+1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc1=f_u;
                F1=fc1.*(Ksi'*n1)*dy;
            %% Flux on face 2 - SOUTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n2>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j+1,i);
                    else
                        f_u(k,1)=f_cv(k,j+1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc2=f_u;
                F2=fc2.*(Ksi'*n2)*dx;
            %% Flux on face 3 - WEST
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n3>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j,i-1);
                    else
                        f_u(k,1)=f_cv(k,j,i-1);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc3=f_u;
                F3=fc3.*(Ksi'*n3)*dy;
            %% Flux on face 4 - NORTH
                f_u=zeros(9,1);
                f_d=zeros(9,1);
                for k=1:9
                    if Ksi(:,k)'*n4>=0
                        f_u(k,1)=f_cv(k,j,i);
                        f_d(k,1)=f_cv(k,j-1,i);
                    else
                        f_u(k,1)=f_cv(k,j-1,i);
                        f_d(k,1)=f_cv(k,j,i);
                    end
                end
                % Flux Scheme
                fc4=f_u;
                F4=fc4.*(Ksi'*n4)*dx;
                % Total Flux
                F(:,j,i)=(F1+F2+F3+F4)/A;
        end
    end
end
% Moment Calculations
for j=1:N_cv_y
    for i=1:N_cv_x
        [Rho_cv(1,j,i), U_cv(:,j,i)] = moment_rho_u_d2q9(squeeze(f_cv(:,j,i)));
    end
end
% Equilibrium Calculations
for j=1:N_cv_y
    for i=1:N_cv_x
        f_eq_cv(:,j,i) = eqm_d2q9(squeeze(Rho_cv(1,j,i)),squeeze(U_cv(:,j,i)));
    end
end
% Collision
f_cv=((Tau-dt)/Tau)*f_cv+(dt/Tau)*f_eq_cv-dt*F;

if any(~isfinite(f_cv(:)))
    error('Solution became non-finite at iteration %d.',t)
end

if mod(t,report_interval)==0
    history_index=history_index+1;
    residual_history(history_index)=max(abs(U_cv(:)-U_previous(:)))/U_lid;
    U_previous=U_cv;
    fprintf('iteration %d of %d, residual = %.3e\n',t,Timer,residual_history(history_index))
end
end
runtime=toc;

%% Post Processing
for j=1:N_cv_y
    for i=1:N_cv_x
        [Rho_cv(1,j,i), U_cv(:,j,i)] = moment_rho_u_d2q9(squeeze(f_cv(:,j,i)));
    end
end

figure
quiver(flipud(squeeze(U_cv(1,:,:))),flipud(squeeze(U_cv(2,:,:))),10)
axis equal tight
title('Velocity Field')
xlabel('x')
ylabel('y')
exportgraphics(gcf,'velocity_vectors.png','Resolution',200)

figure
contourf(flipud(squeeze(Rho_cv)),30)
axis equal tight
colorbar
title('Density Field')
xlabel('x')
ylabel('y')
exportgraphics(gcf,'density_field.png','Resolution',200)

Ux=squeeze(U_cv(1,:,:));
Uy=squeeze(U_cv(2,:,:));
VelocityMagnitude=sqrt(Ux.^2+Uy.^2);

figure
contourf(flipud(VelocityMagnitude/U_lid),30,'LineColor','none')
axis equal tight
colorbar
title('Normalized Velocity Magnitude')
xlabel('x')
ylabel('y')
exportgraphics(gcf,'velocity_magnitude.png','Resolution',200)

% Ghia Comparison at Re = 100
x_star=((1:N_cv_x)-0.5)/N_cv_x;
y_star=((1:N_cv_y)-0.5)/N_cv_y;
i_mid=round(N_cv_x/2);
j_mid=round(N_cv_y/2);
u_center=flipud(Ux(:,i_mid))/U_lid;
v_center=Uy(j_mid,:)/U_lid;

y_ghia=[0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 ...
        0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1];
u_ghia=[0 -0.03717 -0.04192 -0.04775 -0.06434 -0.1015 -0.15662 ...
        -0.2109 -0.20581 -0.13641 0.00332 0.23151 0.68717 ...
        0.73722 0.78871 0.84123 1];
x_ghia=[0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 ...
        0.5 0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1];
v_ghia=[0 0.09233 0.10091 0.1089 0.12317 0.16077 0.17507 ...
        0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 ...
        -0.08864 -0.07391 -0.05906 0];

u_mean_absolute_error=mean(abs(interp1([0 y_star 1],[0;u_center;1],y_ghia)-u_ghia));
v_mean_absolute_error=mean(abs(interp1([0 x_star 1],[0 v_center 0],x_ghia)-v_ghia));

figure
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
nexttile
plot(u_center,y_star,'LineWidth',1.5)
hold on
plot(u_ghia,y_ghia,'ro')
hold off
grid on
xlabel('u/U_{lid}')
ylabel('y/L')
title('Vertical Centerline')
legend('FVDBM','Ghia Data','Location','best')
nexttile
plot(x_star,v_center,'LineWidth',1.5)
hold on
plot(x_ghia,v_ghia,'ro')
hold off
grid on
xlabel('x/L')
ylabel('v/U_{lid}')
title('Horizontal Centerline')
legend('FVDBM','Ghia Data','Location','best')
exportgraphics(gcf,'centerline_comparison.png','Resolution',200)

figure
semilogy((1:history_index)*report_interval,residual_history(1:history_index),'LineWidth',1.5)
grid on
xlabel('Iteration')
ylabel('Velocity Residual')
title('Convergence History')
exportgraphics(gcf,'convergence_history.png','Resolution',200)

fprintf('Runtime: %.3f seconds\n',runtime)
fprintf('Density range: [%g, %g]\n',min(Rho_cv(:)),max(Rho_cv(:)))
fprintf('Speed range: [%g, %g]\n',min(VelocityMagnitude(:)),max(VelocityMagnitude(:)))
fprintf('Ghia MAE: u = %.6g, v = %.6g\n',u_mean_absolute_error,v_mean_absolute_error)
save('FVBDM_results.mat','Re','L','U_lid','N_cv_x','N_cv_y','Rho_cv','U_cv', ...
     'residual_history','runtime','u_mean_absolute_error','v_mean_absolute_error')
