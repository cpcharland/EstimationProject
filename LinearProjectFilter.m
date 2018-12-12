function [P,delta_x_hat] = LinearProjectFilter(t,u,x0,R,P,Q)

delta_x_hat = x0;
delta_u = [0;0];

for k = t
    % Prediction Step
    delta_x_hat_minus = F_tilde*delta_x_hat(:,k)+G_tilde*delta_u(:,k);
    P_minus = F_tilde*P(:,:,k)*F_tilde'+Omega_tilde*Q*Omega_tilde';
    delta_u(:,k+1) = u(:,k+1) - u_nom_traj(:,k+1);
    
    % Correction Step
    K = P_minus*H_tilde'*inv(H_tilde*P_minus*H_tilde'+R);
    delta_y(:,k+1) = y_data(:,k+1)-y_nom_traj(:,k+1);
    
    delta_x_hat(:,k+1) = delta_x_hat_minus+K*(delta_y(:,k+1)-H_tilde*delta_x_hat_minus);
    P(:,:,k+1) = (eye(6,6)-K*H_tilde)*P_minus;
   
end

end

