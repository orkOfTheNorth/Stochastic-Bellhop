function J = computeJacobian(TL_plus, TL_minus, h)
% Central finite-difference Jacobian of TL with respect to a parameter.
%
%   J = (TL_plus - TL_minus) / (2*h)     [dB / unit_of_parameter]
%
% Inputs:
%   TL_plus  (Nz x Nr) — TL at nominal + h
%   TL_minus (Nz x Nr) — TL at nominal - h
%   h        (scalar)  — step size in parameter units (Hz, m, or degC)

J = (TL_plus - TL_minus) / (2 * h);
end
