function J=eval_jacobian_2(f,m)

epsilon = 1e-6; % delta
l_m = length(m); % length of m;

f0=feval(f,m); % caclulate f0


for j=1:l_m
        dm = [ zeros(j-1,1); epsilon; zeros(l_m-j,1)];
        J(:,j) = ((feval(f,m+dm) - feval(f,m))/epsilon);
end