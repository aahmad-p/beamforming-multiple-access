function [beta, beta_no_abs] = compute_beta(M_x, user1_theta_deg, user2_theta_deg)
% Description: 
% Input(s):
%           -
% Output(s):
%           - 
%

arguments
   M_x (1,1) double {mustBeNumeric}
   user1_theta_deg (1,1) double {mustBeNumeric}
   user2_theta_deg (1,1) double {mustBeNumeric}
end

lambda_m = 1;
d_m = lambda_m/2;

% Function to generate steering vector (applicable to 1D antenna array)
a_f = @(theta_deg) transpose( exp(-1i*2*pi*(d_m/lambda_m)*cosd(theta_deg)*((1:M_x) - 1)) );

user1_a = a_f(user1_theta_deg);
user2_a = a_f(user2_theta_deg);

% MATLAB dot product: takes complex conjugate of first input
inner_product = dot(user1_a, user2_a);

% non-abs version used in MB-MIMO equation
beta_no_abs = inner_product/M_x;
beta = abs(beta_no_abs);

return