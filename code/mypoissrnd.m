function x = mypoissrnd(lambda)
%MYPOISSRND  Poisson random numbers without toolboxes.
%   x = MYPOISSRND(lambda) returns Poisson draws with mean 'lambda'.
%   'lambda' can be scalar or an array; output matches its size.

    if any(lambda(:) < 0)
        error('mypoissrnd:NegativeLambda','Lambda must be >= 0.');
    end

    lambda = double(lambda);
    x = zeros(size(lambda));

    % Knuth's algorithm: fine for small vectors (like 181 days)
    for i = 1:numel(lambda)
        L = exp(-lambda(i));
        k = 0;
        p = 1;
        while p > L
            k = k + 1;
            p = p * rand;   % uses base MATLAB rand
        end
        x(i) = k - 1;
    end

    x = reshape(x, size(lambda));
end
