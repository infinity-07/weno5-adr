%==========================================================================
% Approximate Dispersion Relation (ADR) of WENO5-JS and UW5
% Reproduces the quasi-linear ADR analysis of:
%   Pirozzoli, "On the spectral properties of shock-capturing schemes",
%   J. Comput. Phys. 219 (2006) 489-497.
%
% For a chosen reduced wavenumber phi, the procedure is:%   1) set v_j(0) = cos(j*phi)  on a periodic grid of N points
%   2) advance the linear advection eq.  v_t + a v_x = 0  with the scheme
%      under analysis for a very small time tau = sigma*h/a   (sigma << 1)
%   3) take the DFT and read off the complex amplitude at wavenumber phi
%   4) extract the modified wavenumber via
%         Phi(phi) = (i/sigma) * log( vhat(phi,tau) / vhat(phi,0) )
%      Re(Phi) -> approximate phase speed (dispersion)
%      Im(Phi) -> numerical dissipation (negative for stable schemes)
%
% Tested against the analytical UW5 modified wavenumber as a sanity check.
% Requires MATLAB R2016b or later (uses local functions in scripts).
%
%==========================================================================

clear; close all; clc;

%% -------- user parameters -----------------------------------------------
N        = 1000;        % number of grid points (periodic)
a        = 1.0;        % advection speed (a > 0 assumed)
h        = 1.0;        % grid spacing (acts as length unit; phi = w*h)
sigma    = 1e-2;       % a*tau/h  -- must be small (try 1e-3 .. 1e-2)
n_sub    = 1;          % RK3-TVD sub-steps inside [0,tau]; raise for accuracy
eps_weno = 1e-6;       % WENO5-JS regularization parameter (Jiang-Shu default)
%-------------------------------------------------------------------------

tau   = sigma*h/a;

% Discrete reduced wavenumbers supported by the grid
% (skip phi=0 and the Nyquist mode where the procedure is ill-conditioned)
j     = (0:N-1).';
modes = (1:N/2-1).';
phi   = 2*pi*modes/N;

Phi_WENO = zeros(numel(modes),1);
Phi_UW5  = zeros(numel(modes),1);

%% -------- sweep wavenumbers ---------------------------------------------
for k = 1:numel(modes)
    n  = modes(k);
    v0 = cos(j*phi(k));            % real, single-mode initial condition
    F0 = fft(v0);

    % WENO5-JS
    vW = advance(v0, h, tau, a, @(v,h,a) weno5js_rhs(v,h,a,eps_weno), n_sub);
    FW = fft(vW);
    Phi_WENO(k) = 1i/sigma * log( FW(n+1) / F0(n+1) );

    % UW5 (linear 5th-order upwind = WENO5 with frozen ideal weights)
    vU = advance(v0, h, tau, a, @uw5_rhs, n_sub);
    FU = fft(vU);
    Phi_UW5(k)  = 1i/sigma * log( FU(n+1) / F0(n+1) );
end

% Analytical UW5 modified wavenumber (closed form, for verification)
Phi_UW5_ana = uw5_analytical(phi);

fprintf('Sanity check: max |UW5_numerical - UW5_analytical| = %.2e\n', ...
        max(abs(Phi_UW5 - Phi_UW5_ana)));

%% -------- plot ----------------------------------------------------------
figure('Position',[100 100 780 720],'Color','w');

% ----- Real part (dispersion) -----
subplot(2,1,1); hold on; box on; grid on;
plot([0 pi],[0 pi],'k--','LineWidth',1.1);                     % spectral
plot(phi, real(Phi_UW5_ana),'b-' ,'LineWidth',1.6);            % UW5 analyt.
plot(phi, real(Phi_UW5)    ,'bo' ,'MarkerSize',5);             % UW5 numer.
plot(phi, real(Phi_WENO)   ,'r-^','MarkerSize',5,'LineWidth',1.2);
xlabel('\phi','FontSize',12);
ylabel('Re(\Phi)','FontSize',12);
title('Approximate dispersion relation - dispersion','FontSize',12);
legend('spectral','UW5 (analytical)','UW5 (ADR)','WENO5-JS (ADR)', ...
       'Location','NorthWest');
% xlim([0 pi]); ylim([0 pi+0.15]);
set(gca,'XTick',0:pi/4:pi,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'});

% ----- Imaginary part (dissipation) -----
subplot(2,1,2); hold on; box on; grid on;
plot([0 pi],[0 0],'k--','LineWidth',1.1);
plot(phi, imag(Phi_UW5_ana),'b-' ,'LineWidth',1.6);
plot(phi, imag(Phi_UW5)    ,'bo' ,'MarkerSize',5);
plot(phi, imag(Phi_WENO)   ,'r-^','MarkerSize',5,'LineWidth',1.2);
xlabel('\phi','FontSize',12);
ylabel('Im(\Phi)','FontSize',12);
title('Approximate dispersion relation - dissipation','FontSize',12);
legend('spectral','UW5 (analytical)','UW5 (ADR)','WENO5-JS (ADR)', ...
       'Location','SouthWest');
% xlim([0 pi]);
set(gca,'XTick',0:pi/4:pi,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'});

%==========================================================================
% Local functions
%==========================================================================

function v = advance(v0, h, T, a, rhs, n_sub)
% Integrate v_t = rhs(v) from 0 to T using n_sub steps of RK3-TVD.
    dt = T/n_sub;
    v  = v0;
    for k = 1:n_sub
        v = rk3_tvd_step(v, h, dt, a, rhs);
    end
end

function v = rk3_tvd_step(v, h, dt, a, rhs)
% One step of the 3rd-order TVD Runge-Kutta scheme of Shu & Osher.
    L0 = rhs(v ,h,a);   v1 = v  + dt*L0;
    L1 = rhs(v1,h,a);   v2 = 0.75*v + 0.25*v1 + 0.25*dt*L1;
    L2 = rhs(v2,h,a);   v  = (1/3)*v + (2/3)*v2 + (2/3)*dt*L2;
end

function dvdt = weno5js_rhs(v, h, a, eps_w)
% WENO5-JS (Jiang-Shu 1996) semi-discrete RHS for v_t + a v_x = 0, a>0.
% Periodic boundary conditions implemented via circshift.

    % Stencil values for f_{j+1/2} : v_{j-2}, v_{j-1}, v_j, v_{j+1}, v_{j+2}
    vm2 = circshift(v,  2);
    vm1 = circshift(v,  1);
    v0  = v;
    vp1 = circshift(v, -1);
    vp2 = circshift(v, -2);

    % Three candidate 3rd-order reconstructions at j+1/2
    q0 = (1/3)*vm2 - (7/6)*vm1 + (11/6)*v0;
    q1 = -(1/6)*vm1 + (5/6)*v0  + (1/3)*vp1;
    q2 = (1/3)*v0  + (5/6)*vp1 - (1/6)*vp2;

    % Jiang-Shu smoothness indicators
    b0 = (13/12)*(vm2 - 2*vm1 + v0 ).^2 + (1/4)*(vm2 - 4*vm1 + 3*v0).^2;
    b1 = (13/12)*(vm1 - 2*v0  + vp1).^2 + (1/4)*(vm1 - vp1).^2;
    b2 = (13/12)*(v0  - 2*vp1 + vp2).^2 + (1/4)*(3*v0 - 4*vp1 + vp2).^2;

    % Ideal (linear) weights
    d0 = 1/10;  d1 = 6/10;  d2 = 3/10;

    % Nonlinear weights
    a0 = d0 ./ (eps_w + b0).^2;
    a1 = d1 ./ (eps_w + b1).^2;
    a2 = d2 ./ (eps_w + b2).^2;
    asum = a0 + a1 + a2;
    w0 = a0./asum;  w1 = a1./asum;  w2 = a2./asum;

    fhat   = w0.*q0 + w1.*q1 + w2.*q2;    % f at j+1/2
    fhat_m = circshift(fhat, 1);           % f at j-1/2

    dvdt = -(a/h) * (fhat - fhat_m);
end

function dvdt = uw5_rhs(v, h, a)
% 5th-order linear upwind scheme (= WENO5 with ideal weights frozen in).
% Coefficients of f_{j+1/2}: (1/30, -13/60, 47/60, 27/60, -3/60)
    vm2 = circshift(v,  2);
    vm1 = circshift(v,  1);
    v0  = v;
    vp1 = circshift(v, -1);
    vp2 = circshift(v, -2);

    fhat = (1/30)*vm2 - (13/60)*vm1 + (47/60)*v0 + (27/60)*vp1 - (3/60)*vp2;
    dvdt = -(a/h) * (fhat - circshift(fhat,1));
end

function Phi = uw5_analytical(phi_arr)
% Closed-form modified wavenumber of UW5 (for verification).
%   Phi = -i * (1 - exp(-i phi)) * sum_l  c_l * exp(i l phi)
% where (c_{-2},...,c_{+2}) = (1/30, -13/60, 47/60, 27/60, -3/60).
    c   = [1/30, -13/60, 47/60, 27/60, -3/60];
    ell = -2:2;
    Phi = zeros(size(phi_arr));
    for k = 1:numel(phi_arr)
        S      = sum(c .* exp(1i*ell*phi_arr(k)));
        Phi(k) = -1i * (1 - exp(-1i*phi_arr(k))) * S;
    end
end
