% Independent Metropolis-Hastings sampling to infer v_1, v_2 and v_3.

% TODO: infer T2_csf, sigma_n

clear all;

max_iterations = 1000;
burn_in = max_iterations*0.3;

% feasible region falls within one SD of sampling distrib => sets sigma_s
% sample two dimensions because equality constraint reduces to 2D problem
% this doesn't sample uniformly over the triangle, but hopefully that won't
% matter
sigma_s = sqrt(2/3)*eye(2,2);
mu_s    = [0 0];

centre_s= [1/3; 1/3; 1/3]; % centre of proposal distrib in xyz space

% initialise based on priors
prop(1,:) = [0 0];

for i=2:max_iterations
    % generate a feasible proposal
    feasible = 0;
    while feasible ~= 1
        % sample proposal in x'y' space
        prop(i,:) = mvnrnd(mu_s, sigma_s);
        if(     prop(i,2) < sqrt(2/3) - prop(i,1)*sqrt(3) ...
                &&   prop(i,2) < sqrt(2/3) + prop(i,1)*sqrt(3) ...
                &&   prop(i,2) > -1/sqrt(6)  )
            feasible = 1;
        else
            feasible = 0;
        end
    end
    
    % translate prop (x' y') to xyz space
    prop_xyz(:,i) = centre_s + (prop(i,1) * [-1/sqrt(2); 1/sqrt(2); 0]) + (prop(i,2) * [-1/sqrt(6); -1/sqrt(6); 2/sqrt(6)]);

    % test plot to show distribution on feasible surface
    % plot3(prop_xyz(1,:), prop_xyz(2,:), prop_xyz(3,:), '.')
    
    % normalising factor
    c = mvnpdf(prop(i-1,:), mu_s, sigma_s) / ...
        mvnpdf(prop(i,  :), mu_s, sigma_s);
        
    % remember: log-likelihood here for stability's sake
    %log_r = lognpdf(prop_xyz(:,i))) - lognpdf(prop_xyz(:,i-1)))
end
