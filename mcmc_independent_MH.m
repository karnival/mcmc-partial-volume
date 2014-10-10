% Independent Metropolis-Hastings sampling to infer v_1, v_2 and v_3.

% TODO: infer T2_csf, sigma_n

clear all;

max_iterations = 1000;
iterations     = 0;
burn_in = max_iterations*0.3;

% feasible region falls within one SD of sampling distrib => sets sigma_s
% sample two dimensions because equality constraint reduces to 2D problem
% this doesn't sample uniformly over the triangle, but hopefully that won't
% matter
sigma_s = sqrt(2/3)*eye(2,2);
mu_s    = [0; 0];

centre_s= [1/3; 1/3; 1/3]; % centre of proposal distrib in xyz space


for i=1:max_iterations
    iterations = iterations+1;
    
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
    prop_xyz(:,i) = centre_s + (prop(1) * [-1/sqrt(2); 1/sqrt(2); 0]) + (prop(2) * [-1/sqrt(6); -1/sqrt(6); 2/sqrt(6)]);
end
