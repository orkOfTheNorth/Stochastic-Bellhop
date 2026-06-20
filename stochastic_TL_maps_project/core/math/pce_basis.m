function Phi = pce_basis(xi, K, dist_type)
% Evaluate PCE orthogonal basis polynomials Phi_0..Phi_K at nodes xi.
%
% Inputs:
%   xi        [N x 1]  standardised sample values
%                      Normal  → N(0,1) space
%                      Uniform → [-1,1] space
%   K         scalar   maximum polynomial order
%   dist_type string   'normal' (probabilist Hermite He_n)
%                      'uniform' (Legendre P_n)
%
% Output:
%   Phi [N x K+1]  columns are Phi_0, Phi_1, ..., Phi_K

N   = numel(xi);
xi  = xi(:);
Phi = zeros(N, K+1);

switch lower(dist_type)
    case 'normal'
        % Probabilist Hermite: He_0=1, He_1=x, He_{n+1}=x*He_n - n*He_{n-1}
        Phi(:,1) = ones(N,1);
        if K >= 1, Phi(:,2) = xi; end
        for n = 1:K-1
            Phi(:,n+2) = xi .* Phi(:,n+1) - n * Phi(:,n);
        end

    case 'uniform'
        % Legendre: P_0=1, P_1=x, P_{n+1}=((2n+1)*x*P_n - n*P_{n-1})/(n+1)
        Phi(:,1) = ones(N,1);
        if K >= 1, Phi(:,2) = xi; end
        for n = 1:K-1
            Phi(:,n+2) = ((2*n+1)*xi.*Phi(:,n+1) - n*Phi(:,n)) / (n+1);
        end

    otherwise
        error('pce_basis: dist_type must be ''normal'' or ''uniform''');
end
end
