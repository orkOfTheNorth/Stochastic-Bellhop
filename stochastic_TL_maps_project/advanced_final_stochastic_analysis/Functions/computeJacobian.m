function J = computeJacobian(TL_plus, TL_minus, h)
% Central finite-difference first-order derivative (Jacobian map).
%   J(z,r) = [TL(θ₀+h) − TL(θ₀−h)] / (2h)   [dB per unit of θ]
%
% TL_plus, TL_minus : Nz×Nr TL maps from +h and −h perturbation runs
% h                 : scalar step size in units of θ
%
% Swap this function to implement higher-order stencils or
% complex-step differentiation in the future.
    J = (TL_plus - TL_minus) / (2 * h);
end
