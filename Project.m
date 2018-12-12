%% Final Project Part 1
%% Cody Charland & Ian Cooke
% We are using the Simple Skycrane.
clear all; close all;
saving = 0; % Will generate new data, save = 0 will use saved data
load('skycrane_finalproj_KFdata.mat');
% Nominal System Params
rho = 0.020; % [kg/m^3]
g = 3.711; % [m/s^2]
beta = pi/4; % [rad]
C_D = 0.2; % [none]
m_f = 390; % [kg]
w_f = 1; % [m]
h_f = 0.5; % [m]
d_f = 1; % [m]
m_b = 1510; % [kg]
w_b = 3.2; % [m]
h_b = 2.5; % [m]
d_b = 2.9; % [m]
h_cm = 0.9421; % [m]
w_cm = w_b/2;
A_side = (h_b*d_b) + (h_f*d_f); % [m^2]
A_bot = (w_b*d_b) + (w_f*d_f); % [m^2]
T_1 = 0.5*g*(m_b + m_f)/cos(beta);
T_2 = T_1;
I_n = 1/12*(m_b*(w_b^2 + h_b^2) + m_f*(w_f^2 + h_f^2));
theta = 0; % [rad]
alpha = 0; % [rad]
xi_dot = 0; % [m/s]
z_dot = 0;
theta_dot= 0;
dt = 0.1;
Gamma = [zeros(1,3);1,0,0;zeros(1,3);0,1,0;zeros(1,3);0,0,1];


% CT Jacobians
A_tilde = [0,1,0,0,0,0;...
    0, -1/2*rho*C_D/(m_b+m_f)*(A_side*(cos(theta - alpha)*(2*xi_dot^2 + z_dot^2)/sqrt(xi_dot^2 + z_dot^2) - sin(theta - alpha)*z_dot*xi_dot/sqrt(xi_dot^2 + z_dot^2)) + A_bot*(sin(theta - alpha)*(2*xi_dot^2 + z_dot^2)/sqrt(xi_dot^2 + z_dot^2) + cos(theta - alpha)*z_dot*xi_dot/sqrt(xi_dot^2 + z_dot^2))), 0, -1/2*rho*C_D/(m_b+m_f)*(A_side*(cos(theta - alpha)*z_dot*xi_dot/sqrt(xi_dot^2 + z_dot^2) + sin(theta - alpha)*xi_dot^2/sqrt(xi_dot^2 + z_dot^2)) + A_bot*(sin(theta - alpha)*z_dot*xi_dot)/sqrt(xi_dot^2 + z_dot^2) - cos(theta - alpha)*xi_dot^2/sqrt(xi_dot^2 + z_dot^2)), 1/(m_b+m_f)*(T_1*(cos(beta)*cos(theta) - sin(beta)*sin(theta)) + T_2*(cos(beta)*cos(theta) + sin(beta)*sin(theta))) - 1/2*rho*C_D/(m_b+m_f)*((-A_side*sin(theta - alpha) + A_bot*cos(theta - alpha))*xi_dot*sqrt(xi_dot^2 + z_dot^2)), 0;...
    0,0,0,1,0,0;...
    0, -1/2*rho*C_D/(m_b+m_f)*(A_side*(cos(theta - alpha)*z_dot*xi_dot/sqrt(xi_dot^2 + z_dot^2) - sin(theta - alpha)*z_dot^2/sqrt(xi_dot^2 + z_dot^2)) + A_bot*(sin(theta-alpha)*z_dot*xi_dot/sqrt(xi_dot^2 + z_dot^2) + cos(theta - alpha)*z_dot^2/sqrt(xi_dot^2 + z_dot^2))), 0, -1/2*rho*C_D/(m_b+m_f)*(A_side*(cos(theta-alpha)*(2*z_dot^2 + xi_dot^2)/sqrt(xi_dot^2 + z_dot^2) + sin(theta - alpha)*xi_dot*z_dot/sqrt(xi_dot^2 + z_dot^2)) + A_bot*(sin(theta - alpha)*(2*z_dot^2 + xi_dot^2)/sqrt(xi_dot^2 + z_dot^2) - cos(theta - alpha)*xi_dot*z_dot/sqrt(xi_dot^2 + z_dot^2))), 1/(m_b+m_f)*(T_1*(-cos(beta)*sin(theta) - sin(beta)*cos(theta)) + T_2*(-cos(beta)*sin(theta) + sin(beta)*cos(theta))) - 1/2*rho*C_D/(m_b+m_f)*((-A_side*sin(theta - alpha) + A_bot*cos(theta - alpha))*z_dot/sqrt(xi_dot^2 + z_dot^2)), 0;...
    zeros(1,5),1;...
    zeros(1,6)];
[nan_i, nan_j] = find(isnan(A_tilde));
for i = 1:length(nan_i)
    A_tilde(nan_i(i), nan_j(i)) = 0;
end

B_tilde = [0 0;...
    
(cos(beta)*sin(theta)+sin(beta)*cos(theta))/(m_b+m_f) (cos(beta)*sin(theta)-sin(beta)*cos(theta))/(m_b+m_f);...

0 0;...

(cos(beta)*cos(theta)-sin(beta)*sin(theta))/(m_b+m_f) (cos(beta)*cos(theta)+sin(beta)*sin(theta))/(m_b+m_f);...

0 0; ...

1/I_n*(cos(beta)*w_b/2-sin(beta)*h_cm) 1/I_n*(-w_b/2*cos(beta)+sin(beta)*h_cm);];




a1 = -1/2*C_D*rho/(m_b+m_f)*((-A_side*sin(theta-alpha)*z_dot/(xi_dot^2+z_dot^2)+A_bot*cos(theta-alpha)*z_dot/(xi_dot^2+z_dot^2))*xi_dot^2*sqrt(xi_dot^2+z_dot^2)+(A_side*cos(theta-alpha)+A_bot*sin(theta-alpha))*sqrt(xi_dot^2+z_dot^2)+(A_side*cos(theta-alpha)+A_bot*sin(theta-alpha))*xi_dot^2/sqrt(xi_dot^2+z_dot^2));

a2 = -1/2*C_D*rho/(m_b+m_f)*((-A_side*sin(theta-alpha)*xi_dot/(xi_dot^2+z_dot^2)-A_bot*cos(theta-alpha)*xi_dot/(xi_dot^2+z_dot^2))*xi_dot*sqrt(xi_dot^2+z_dot^2)+z_dot*xi_dot/sqrt(xi_dot^2+z_dot^2)*(A_side*cos(theta-alpha)+A_bot*sin(theta-alpha)));

a3 = (T_1*(cos(beta)*cos(theta)-sin(beta)*sin(theta))+T_2*(cos(beta)*cos(theta)+sin(beta)*sin(theta)))/(m_b+m_f)-1/2*rho*C_D/(m_b+m_f)*(-A_side*sin(theta-alpha)+A_bot*cos(theta-alpha))*xi_dot*sqrt(xi_dot^2+z_dot^2);



C_tilde = [1 0 0 0 0 0;...
    
0 0 1 0 0 0;...

0 0 0 0 0 1;...

0 a1 0 a2 a3 0];

[nan_i, nan_j] = find(isnan(C_tilde));
for i = 1:length(nan_i)
    C_tilde(nan_i(i), nan_j(i)) = 0;
end



D_tilde = [ 0 0; 0 0; 0 0; (cos(beta)*sin(theta)+sin(beta)*cos(theta))/(m_b+m_f) (cos(beta)*sin(theta)+sin(beta)*cos(theta))/(m_b+m_f)];


%% b) Compute DT Jacobians using Euler's method and discuss Controllability, Stability, and Observability

F_tilde = eye(6) + dt*A_tilde;
G_tilde = dt*B_tilde;
H_tilde = C_tilde;
M_tilde = dt*D_tilde;
Omega_tilde = dt*Gamma;

system = struct;
system.F_tilde = F_tilde;
system.G_tilde = G_tilde;
system.H_tilde = H_tilde;
system.M_tilde = M_tilde;
system.Omega_tilde = Omega_tilde;
system.Klin = Klin;
system.Qtrue = Qtrue;
system.Rtrue = Rtrue;

%% Build Ground Truth Data
x_NL = [0,0,20,0,0,0]';
P = diag([1 1 1 1 .05 .05]);
x_Nom = x_NL;
%NL_data = x_NL;
u_Nom = [0.5*g*(m_b + m_f)/cos(beta); 0.5*g*(m_b + m_f)/cos(beta);];
if saving
    numTrials = 100;
    for k = 1:numTrials
        [NL_data(:,:,k),y_(:,:,k),u_(:,:,k)] = MakeTruth(0,50,500,system,x_NL,u_Nom);
        %     for j = 1:6
        %         subplot(3,2,j)
        %         plot(0:0.1:50,NL_data(j,:,k))
        %         hold on
        %     end
    end
    % Save off Data
    
    filename = 'truthData.mat';
    save(filename, 'NL_data', 'y_', 'u_', 'numTrials');
    
else
    load('truthData.mat');
end
figure
p= mean(var(NL_data,0,3),2)
%% Build Nominal Trajectory
x_nom_traj = repmat(x_Nom,1,501);
u_nom_traj = repmat([T_1;T_2],1,501);
y_nom_traj = repmat([x_nom_traj(1); x_nom_traj(3); x_nom_traj(6); 0],1,501);
%% Linear KF
P = diag(p);
P(1,1) = 2;
P(2,2) = 1/20;
P(3,3) = 1/12;
P(4,4) = 1.5;
P(5,5) = .005;
delta_x_hat = [0;0;0;0;0;0;];
delta_u = [0;0];
Q = Qtrue;
Q_KF = diag([25 .0008 .0002 1 .1 .1]);


R = Rtrue
R(3,3) = .001;
R(4,4) = .04;
R(2,2) = 5;
R(1,1) = 25*7;
R = R./7;


for j = 1:numTrials
    u = u_(:,:,j);
    y_data = y_(:,:,j);
    delta_u = [0;0];
    
    for k = 1:499
        % Prediction Step
        delta_x_hat_minus = F_tilde*delta_x_hat(:,k)+G_tilde*delta_u(:,k);
        P_minus = F_tilde*P(:,:,k)*F_tilde'+Q_KF;
        delta_u(:,k+1) = u(:,k+1) - u_nom_traj(:,k+1);
        
        % Correction Step
        K = P_minus*H_tilde'*inv(H_tilde*P_minus*H_tilde'+R);
        delta_y(:,k+1) = y_data(:,k+1)-y_nom_traj(:,k+1)+M_tilde*delta_u(:,k+1);
        %delta_y(:,k+1) = delta_y(:,k+1)';
        
        delta_x_hat(:,k+1) = delta_x_hat_minus+K*(delta_y(:,k+1)-H_tilde*delta_x_hat_minus);
        P(:,:,k+1) = (eye(6,6)-K*H_tilde)*P_minus;
        NEES(k) = (((NL_data(:,k,j)-x_Nom)-delta_x_hat(:,k+1))'*inv(P(:,:,k))*((NL_data(:,k,j)-x_Nom)-delta_x_hat(:,k+1)));
        NIS(k)  = delta_y(:,k+1)'*inv(H_tilde*P_minus*H_tilde'+R)*delta_y(:,k+1);
    end
    allNEES(j,:) = NEES;
    NEES = 0;
    allNIS(j,:) = NIS;
    NIS = 0;
    
    
end

%%DO NEES TEST:
epsNEESbar = mean(allNEES,1);
alphaNEES = 0.05; %%significance level
Nnx = numTrials*length(F_tilde);
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ numTrials;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ numTrials;

figure

plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, $\bar{\epsilon_x}$','FontSize',14,'Interpreter','latex')
xlabel('time step, k','FontSize',14)
%string = sprintf('NEES Estimation Results Q = %f',m);
%title(string,'FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')

figure

%%DO NIS TEST:
epsNISbar = mean(allNIS,1);
alphaNIS = 0.05; 
Nny = numTrials*size(H_tilde,1);
%%compute intervals:
r1y = chi2inv(alphaNIS/2, Nny )./ numTrials;
r2y = chi2inv(1-alphaNIS/2, Nny )./ numTrials;
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, $\bar{\epsilon}_y$','FontSize',14,'Interpreter','latex')
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
 figure
x_hat = delta_x_hat+ x_nom_traj(:,1);

string = {'xi','xi_dot','z','z dot','theta','theta dot'};
for k = 1:6
    subplot(3,2,k)
    
    plot(0.1:.1:50,x_hat(k,:))
    ylabel(string(k))
    hold on
    plot(0:.1:50,NL_data(k,:,end),'k')
    if k == 3
        plot(0.1:.1:50,2*(sqrt(reshape(P(k,k,:),1,[])))+20)
        plot(0.1:.1:50,-2*(sqrt(reshape(P(k,k,:),1,[])))+20)
    else
        plot(0.1:.1:50,2*(sqrt(reshape(P(k,k,:),1,[]))))
        plot(0.1:.1:50,-2*(sqrt(reshape(P(k,k,:),1,[]))))
    end
%     if k == 1
%         plot(0.1:.1:50,y_data(1,:,end),'r.')
%     elseif k == 3
%         plot(0.1:.1:50,y_data(2,:,end),'r.')
%     elseif k == 6
%         plot(0.1:.1:50,y_data(3,:,end),'r.')
%     end
    
end


