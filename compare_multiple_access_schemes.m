function compare_multiple_access_schemes(M_x, M_y, user1_theta_deg, user2_theta_deg_vec)
% Description: 
% Input(s):
%           -
% Output(s):
%           - 
%

arguments
   M_x (1,1) double {mustBeNumeric} = 128
   M_y (1,1) double {mustBeNumeric} = 1
   user1_theta_deg (1,1) double {mustBeNumeric} = 10
   user2_theta_deg_vec (1,:) double {mustBeNumeric} = 10:0.001:130
end

% Define problem parameters
%  power
Pe_dBm = 30;
sigmasq_dBm = -101;
rho = 10^(Pe_dBm/10)/10^(sigmasq_dBm/10);

%  channel coefficients
alpha_vec = [8.3e-6, 7.4e-6];

%  power splitting (MB-NOMA)
gamma_vec = [0.2, 0.8];


% Define psi_f for MB-NOMA
%psi_f = @(theta_deg) 1i*2*pi*(M_x/2)*2*pi*cosd(theta_deg);
psi_f = @(theta_deg) 1i*2*pi*(M_x/2)*cosd(theta_deg);

M = M_x * M_y;

beta_1_2 = nan(size(user2_theta_deg_vec));
beta_no_abs_1_2 = nan(size(user2_theta_deg_vec));
beta_1_g = nan(size(user2_theta_deg_vec));
beta_2_g = nan(size(user2_theta_deg_vec));

sum_rate_SDMA = nan(size(user2_theta_deg_vec));
sum_rate_SB_NOMA = nan(size(user2_theta_deg_vec));
sum_rate_MB_NOMA = nan(size(user2_theta_deg_vec));

for k = 1:length(user2_theta_deg_vec)
    user2_theta_deg = user2_theta_deg_vec(k);
    
    % SDMA
    beta_1_2(k) = compute_beta(M_x, user1_theta_deg, user2_theta_deg);
    sum_rate_SDMA(k) = log2(1 + ( rho * alpha_vec(1)^2 * M )/( rho * alpha_vec(1)^2 * beta_1_2(k)^2 * M + 2)) + ...
        log2(1 + ( rho * alpha_vec(2)^2 * M )/( rho * alpha_vec(2)^2 * beta_1_2(k)^2 * M + 2));
    
    % SB-NOMA
    beta_1_g(k) = compute_beta(M_x, user1_theta_deg, (user1_theta_deg + user2_theta_deg)/2);
    beta_2_g(k) = compute_beta(M_x, user2_theta_deg, (user1_theta_deg + user2_theta_deg)/2);
    sum_rate_SB_NOMA(k) = log2(1 + rho * gamma_vec(1) * alpha_vec(1)^2 * beta_1_g(k)^2 * M) + ...
        log2(1 + ( rho * (1 - gamma_vec(1)) * alpha_vec(2)^2 * beta_2_g(k)^2 * M )/( rho * gamma_vec(1) * alpha_vec(2)^2 * beta_2_g(k)^2 * M + 1 ));

    % MB-NOMA
    dpsi = psi_f(user1_theta_deg) - psi_f(user2_theta_deg);
    [~, beta_no_abs_1_2(k)] = compute_beta(M_x/2, user1_theta_deg, user2_theta_deg);
    sum_rate_MB_NOMA(k) = log2(1 + rho * gamma_vec(1) * alpha_vec(1)^2 * abs(0.5 * (1 + beta_no_abs_1_2(k)))^2 * M) + ...
        log2(1 + ( rho * (1 - gamma_vec(1)) * alpha_vec(2)^2 * abs(0.5 * (beta_no_abs_1_2(k) + exp(dpsi)))^2 * M )/( rho * gamma_vec(1) * alpha_vec(2)^2 * abs(0.5 * (beta_no_abs_1_2(k) + exp(dpsi)))^2 * M + 1 ));
end

% Sort for plotting (connected lines)
[beta_sorted, sort_i] = sort(beta_1_2);

% Generate plots
figure(1);
plot(user2_theta_deg_vec - user1_theta_deg, beta_1_2);
xlabel('delta theta');
ylabel('beta');

delta_theta_deg = user2_theta_deg_vec - user1_theta_deg;
figure(2);
clf;
hold on;
grid on;
plot(delta_theta_deg, sum_rate_MB_NOMA);
plot(delta_theta_deg, sum_rate_SB_NOMA);
plot(delta_theta_deg, sum_rate_SDMA);
xlabel('${\Delta\theta}$ [deg]','Interpreter','Latex');
ylabel('Sum Rate');
legend('MB-NOMA', 'SB-NOMA', 'SDMA');
title('Sum Rate vs Angular Difference between UEs');
saveas(2, 'sumrate_deltatheta.png', 'png');

figure(3);
clf;
hold on;
grid on;
plot(beta_sorted, sum_rate_MB_NOMA(sort_i));
plot(beta_sorted, sum_rate_SB_NOMA(sort_i));
plot(beta_sorted, sum_rate_SDMA(sort_i));
xlabel('${\beta}$','Interpreter','Latex');
ylabel('Sum Rate');
legend('MB-NOMA', 'SB-NOMA', 'SDMA');
title('Sum Rate vs Spatial Interference between UEs');
saveas(3, 'sumrate_beta.png', 'png');

return