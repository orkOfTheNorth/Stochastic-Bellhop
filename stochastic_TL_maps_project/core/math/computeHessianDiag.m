function H = computeHessianDiag(TL_up, TL_nom, TL_dn, h)
% Diagonal Hessian d²TL/dθ² via 2nd-order central finite difference.
%
% Formula: H = (TL(θ+h) - 2*TL(θ) + TL(θ-h)) / h²
% Error is O(h²) — same order as the Jacobian central FD.
%
% Uses the same ±h Bellhop runs already computed for the Jacobian,
% so no additional Bellhop calls are needed.
H = (TL_up - 2*TL_nom + TL_dn) / h^2;
end
