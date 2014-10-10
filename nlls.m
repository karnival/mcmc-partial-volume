% Constrained Gauss-Newton fitting for partial volumes

% S/S_0 = v1*exp(-TE/T2_1) + v2*exp(-TE/T2_2) + v3*exp(-TE/T2_3)
% v1+v2+v3 = 1

clear all;

% independent variable is time
TEs = [0:10E-3:160E-3];

TE_min = 50E-3;
TE_max = 150E-3;

% grey, white, csf
T2  = [110E-3; 80E-3; 400E-3];

% true volume fractions
v_t = [0.7; 0.2; 0.1];

% define func to create data and test against
func = @(v)( (v(1)*exp(-TEs/T2(1)) + v(2)*exp(-TEs/T2(2)) + v(3)*exp(-TEs/T2(3))) ...
              .* heaviside(TEs - TE_min) .* heaviside(TE_max - TEs) );

% generate noisy data -- don't bother saving the noise, though
sigma_n = 0.05;

% repeated measurements
repetitions = 10;

for i=1:repetitions
    y(i,:) = noise_generator(v_t, sigma_n, func, 1);

end

% minimise squared error over v_m
% constraints: sum(v_m)=1, v_m(i)>0
v_m = fmincon((@(v_m)(norm(y - (func(v_m)'*ones(1,repetitions))' ))), [0.6; 0.3; 0.1], [], [], [1 1 1], 1, [0; 0; 0], [1; 1; 1]);
