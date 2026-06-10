function lb = chebyshevBound(EX, Var, FOM)
% Chebyshev lower bound on P(TL > FOM).
%
%   lb = max(0,  1 - Var / (Var + d^2))   where d = EX - FOM, d > 0
%
% Only pixels where EX > FOM get a non-zero bound (shadow candidates).
% All others return 0 (the point is already in the detection zone in expectation).

lb   = zeros(size(EX));
mask = (EX - FOM) > 0;
d    = EX(mask) - FOM;
V    = Var(mask);
lb(mask) = max(0, 1 - V ./ (V + d.^2));
end
