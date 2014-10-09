function y = noise_generator(m, sigma_n, func, new_noise)
%noise_generator() generates noise and returns noisy data. If new_noise==1 then it
%generates new noise and returns new y, with no data saved. If new_noise==0
%then it tries to load old noise -- generating and saving noise of the
%appropriate SD if necessary.
%   m - model parameters, as a vector
%   sigma_n - noise variance

%y = m(1) .* exp(-m(2)*x) + m(3) .* exp(-m(4)*x);

% noise-free y
%y = m(1) .* x.^2 .* (1-heaviside(x-m(3))) + m(1)*(m(2).*(x-m(3)) + 4).*(heaviside(x-m(3))-heaviside(x-5));

y = func(m);

% add noise
if sigma_n > 0 % TODO: add a provision to reuse pre-generated noise data -- lets us "compare apples with apples"...
%    noise = sigma_n * randn(1, length(y));
%    y = y + noise;

    if new_noise == 1
        fprintf('Generating new noise -- NOT USING any existing noise file.\n');
        noise = sigma_n * randn(1, length(y));
           
    else
        
        noise_folder = 'noise';
        cd(noise_folder);
        noise_filename = strcat('noise_',num2str(sigma_n),'.mat');
        noise_exists_for_sigma = exist(noise_filename,'file');
        
        if noise_exists_for_sigma == 2
            load(noise_filename,'noise');
        else
            fprintf('Noise does not yet exist for this variance. Generating noise.\n');
            noise = sigma_n * randn(1, length(y));
            fprintf('Saving new noise file.\n');
            save(noise_filename,'noise');
        end
        cd ..;
        
    end
    
    y = y + noise(1:length(y));
%    plot(x, y);
end

end
