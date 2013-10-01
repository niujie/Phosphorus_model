function X = cspoirnd(lam, n)
% This function will generate Poisson random variables with parameter
% lambda.

X = zeros(1, n);
for j = 1 : n
    b = rand(1);
    i = 0;
    while b >= exp(-lam)
        i = i + 1;
        b = b * rand(1);
    end
    X(j) = i;
end