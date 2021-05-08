%% Introduction
% Uses the symbolic processing toolbox to analyse the quadrature curve.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Define wavelets and velocity

syms rho phi;

psi_e = (2.*rho.^2 - 1) ./ (rho.^2 + 1)^(5/2);
psi_o = (3.*rho) ./ (rho.^2 + 1)^(5/2);
psi_n = (2-rho.^2) ./ (rho.^2 + 1)^(5/2);

v_x = psi_e .* cos(phi) + psi_o .* sin(phi);
v_y = psi_o .* cos(phi) + psi_n .* sin(phi);


%% Check v_x square

v_x_square = v_x.^2;
v_x_square_norm = v_x_square * (1 + rho.^2).^5;
v_x_square_norm_hand = (4*rho^4 - 13*rho^2 + 1) * cos(phi)^2 +...
                       (6*rho^3 - 3*rho) * sin(2 * phi) +...
                       9*rho^2;
                   
% Check of ze hetzelfde zijn:
simplify(expand(v_x_square_norm_hand))
simplify(expand(v_x_square_norm))

%% Check v_y square

v_y_square = v_y.^2;
v_y_square_norm = v_y_square * (1 + rho.^2).^5;
v_y_square_norm_hand = (6*rho - 3*rho^3) * sin(2*phi) + ...
                       (13/2*rho^2 - 1/2*rho^4 -2) * cos(2*phi) + ...
                       5/2 * rho^2 + 1/2 * rho^4 + 2;
                   
% Check of ze hetzelfde zijn:
simplify(expand(v_y_square_norm_hand))
simplify(expand(v_y_square_norm))

%% Check quadrature curve

psi_quad_norm = (v_x_square_norm + 1 / 2 * v_y_square_norm);
psi_quad_norm_hand = (7/4*rho^4 - 13/4*rho^2 - 1/2) * cos(2*phi) + ...
                     (9/2*rho^3) * sin(2*phi) + ...
                     15/4*rho^2 + 9/4*rho^4 + 3/2;
simplify(expand(psi_quad_norm))
simplify(expand(psi_quad_norm_hand))

%% Check sym and skew component

psi_quad = sqrt(v_x.^2 + 1/2*v_y.^2);
psi_sym = ((7/4*rho^4 - 13/4*rho^2 - 1/2) * cos(2*phi) + 15/4*rho^2 + 9/4*rho^4 + 3/2) / (1 + rho^2)^5;
psi_skew = ((9/2*rho^3) * sin(2*phi)) / (1 + rho^2)^5;
psi_quad_hand = sqrt(psi_sym + psi_skew);

simplify(expand(psi_quad))
simplify(expand(psi_quad_hand))

%% Derivative psi_sym

syms C
psi_sym_norm = psi_sym / (1 + sin(phi)^2) * (1 + rho^2)^5 / C;
d_psi_sym_d_phi = diff(psi_sym_norm, phi);
d_psi_sym_d_phi = simplify(expand(d_psi_sym_d_phi))
solve(d_psi_sym_d_phi == 0)

%% Value of anchor points

phis = 0:0.1:2*pi;
psi_sym_norm = psi_sym / (1 + sin(phi)^2);
values = eval(subs(subs(psi_sym_norm, rho, 2/sqrt(5)), phi, phis));

plot(phis, values)
values(1)