% Independent Metropolis-Hastings sampling to infer v_1, v_2 and v_3.

% TODO: infer T2_csf, sigma_n

clear all;

%% set up function and generate noisy data
% S/S_0 = v1*exp(-TE/T2_1) + v2*exp(-TE/T2_2) + v3*exp(-TE/T2_3)
% v1+v2+v3 = 1

clear all;

% independent variable is time

TE_min = 50E-3;
TE_max = 150E-3;

TEs = [TE_min:10E-3:TE_max];

% grey, white, csf
T2  = [110E-3; 80E-3; 400E-3];

% true volume fractions
v_t = [0.7; 0.2; 0.1];

% define func to create data and test against
func = @(v)( (v(1)*exp(-TEs/T2(1)) + v(2)*exp(-TEs/T2(2)) + v(3)*exp(-TEs/T2(3))) ...
              .* heaviside_asymm(TEs - TE_min) .* heaviside_asymm(TE_max - TEs) );

% generate noisy data -- don't bother saving the noise, though
% TODO: frame this in terms of SNR
sigma_n = 0.02;

% repeated measurements
repetitions = 10;

for i=1:repetitions
    y(i,:) = noise_generator(v_t, sigma_n, func, 1);

end


%% begin MCMC stuff

max_iterations = 100000;
burn_in = max_iterations*0.3;

accepted = 0;
rejected = 0;

samples = []; % accepted samples
draws   = []; % all proposed samples

% feasible region falls within one SD of sampling distrib => sets sigma_s
% sample two dimensions because equality constraint reduces to 2D problem
% this doesn't sample uniformly over the triangle, but hopefully that won't
% matter
sigma_s = sqrt(2/3)*eye(2,2);
mu_s    = [0 0];

centre_s= [1/3; 1/3; 1/3]; % centre of proposal distrib in xyz space

% initialise based on priors
% samples in x'y' space
samples_orig(1,:) = [0 0];
samples_xyz(1,:)  = [1/3 1/3 1/3];

for i=2:max_iterations
    % generate a feasible proposal
    feasible = 0;
    while feasible ~= 1
        % sample proposal in x'y' space
        prop = mvnrnd(mu_s, sigma_s);
        if(     prop(2) < sqrt(2/3) - prop(1)*sqrt(3) ...
                &&   prop(2) < sqrt(2/3) + prop(1)*sqrt(3) ...
                &&   prop(2) > -1/sqrt(6)  )
            feasible = 1;
        else
            feasible = 0;
        end
    end
    
    % translate prop (x' y') to xyz space
    prop_xyz = centre_s + (prop(1) * [-1/sqrt(2); 1/sqrt(2); 0]) + (prop(2) * [-1/sqrt(6); -1/sqrt(6); 2/sqrt(6)]);
    prop_xyz = prop_xyz';
    
    % test plot to show distribution on feasible surface
    % plot3(prop_xyz(1,:), prop_xyz(2,:), prop_xyz(3,:), '.')
    
    % normalising factor
    c = mvnpdf(samples_orig(i-1,:), mu_s, sigma_s) / ...
        mvnpdf(prop,                mu_s, sigma_s);
        
    % remember: log-likelihood here for stability's sake
    % will need to account for changing variance later
    A = -sum((y(1,:) - func(prop_xyz          )).^2)/(2*sigma_n^2);
    B = -sum((y(1,:) - func(samples_xyz(i-1,:))).^2)/(2*sigma_n^2);

    % priors are taken here as equal to transition distrib
    % prior_ratio = mvnpdf(prop, mu_s, sigma_s) / mvnpdf(samples_orig(i-1,:), mu_s, sigma_s);
    
    alpha = min(1, c*exp(A - B));
    
    u = rand();
    
    if u < alpha
        accepted = accepted + 1;
        samples_orig = [samples_orig; prop];
        samples_xyz  = [samples_xyz;  prop_xyz];
    else
        rejected = rejected + 1;
        samples_orig = [samples_orig; samples_orig(end,:)];
        samples_xyz  = [samples_xyz;  samples_xyz(end,:)];
    end
    
    draws = [draws; prop_xyz];
end

% trace plot
% plot(draws(:,1))
