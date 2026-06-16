function [xi, w, gamma] = gaussQuadNodes(K, dist_type)
% gaussQuadNodes  Gauss quadrature nodes/weights for PCE coefficient projection.
%
% Returns K nodes and weights for the distribution:
%   'normal'  — Gauss-Hermite (probabilist), orthogonal w.r.t. N(0,1)
%   'uniform' — Gauss-Legendre, orthogonal w.r.t. Uniform[-1,1]
%
% Inputs:
%   K         — number of quadrature points (= PCE order + 1 for exact integration)
%   dist_type — 'normal' or 'uniform'
%
% Outputs:
%   xi    [K x 1] — quadrature nodes (standardized: N(0,1) or U[-1,1])
%   w     [K x 1] — quadrature weights (sum to 1 for normal, sum to 2 for uniform)
%   gamma [K x 1] — normalisation constants <Phi_n, Phi_n> for n=0..K-1
%                   Normal:  gamma(n) = n!   (probabilist Hermite)
%                   Uniform: gamma(n) = 2/(2n-1)  (Legendre, 1-indexed)

% Build symmetric tridiagonal Jacobi matrix — eigenvalues are nodes,
% eigenvectors give weights. Golub-Welsch algorithm.

switch lower(dist_type)
    case 'normal'
        % Probabilist Hermite He_n: recurrence b_n = sqrt(n)
        b = sqrt(1:K-1)';   % [K-1 x 1]
        J = diag(b, 1) + diag(b, -1);   % symmetric tridiagonal, diagonal=0
        gamma = factorial(0:K-1)';       % <He_n, He_n>_{N(0,1)} = n!

    case 'uniform'
        % Legendre P_n on [-1,1]: recurrence b_n = n/sqrt(4n^2-1)
        n_vec = (1:K-1)';
        b = n_vec ./ sqrt(4*n_vec.^2 - 1);
        J = diag(b, 1) + diag(b, -1);
        n_idx = (0:K-1)';
        gamma = 2 ./ (2*n_idx + 1);     % <P_n, P_n>_{U[-1,1]} = 2/(2n+1)

    otherwise
        error('gaussQuadNodes: dist_type must be ''normal'' or ''uniform''');
end

% Golub-Welsch: eigendecomposition of J
[V, D]   = eig(J, 'vector');
[xi, ix] = sort(D);
V        = V(:, ix);

switch lower(dist_type)
    case 'normal'
        % w_i = (V(1,i))^2  (first component squared, no extra sqrt(2*pi) needed
        % since He polynomials are orthogonal w.r.t. standard normal PDF)
        w = V(1,:)'.^2;    % weights sum to 1 (probability weights)
    case 'uniform'
        % w_i = 2*(V(1,i))^2   (Legendre: integral over [-1,1] = 2)
        w = 2 * V(1,:)'.^2;
end
end
