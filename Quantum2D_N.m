%% Quantum2D_N ------------------------------------------------------------
%
%  Two-dimensional time-dependent Schroedinger equation solved with the
%  factored diagonal Pade propagator combined with high-order symmetric
%  splitting (chi/chi* composition of McLachlan, 2022).
%
%  Reproduces Fig. 3 of:
%    L. Freschi, C. Leveque, J. Caillat, R. Taieb, J. Dubois,
%    "High-order Crank-Nicolson method via diagonal Pade factorization",
%    Computer Physics Communications (submitted).
%
%  Problem:   u(T) = exp(-1i*(Kx + Ky + V)*T) * u(0)
%    Kx = (1/(2*dx^2)) * tridiag(-1, 2, -1)  along x
%    Ky = (1/(2*dy^2)) * tridiag(-1, 2, -1)  along y
%    V  = diag(V(x,y)),  V(x,y) = exp(-wx*x^2 - wy*y^2 + lambda*x*y)
%
%  One time step (operators listed in order of APPLICATION):
%    stage s odd  (chi* in the paper): Kx -> Ky -> V   with h = beta_s*dt
%    stage s even (chi  in the paper): V  -> Ky -> Kx  with h = beta_s*dt
%  for s = 1, ..., 2*Ns, with the palindromic sequence
%    beta_full = [beta_1, ..., beta_Ns, beta_Ns, ..., beta_1].
%
%  exp(-1i*h*Kx) and exp(-1i*h*Ky) are approximated by the factored Pade
%  propagator of order [q,q]
%    U_q(h,K) = prod_k [1 + 1i*alpha_k*h*K]^{-1} [1 - 1i*alpha_k*h*K],
%  applied to all grid lines at once (vectorized Thomas algorithm).
%  exp(-1i*h*V) is applied exactly (elementwise multiplication).
%
%  The reference solution is computed with expmv (MATLAB R2023b or later).
%
%  Storage: psi is an Ny x Nx array (rows = y, columns = x). MATLAB's
%  column-major vectorization psi(:) therefore has y varying fastest,
%  and the full 2D Hamiltonian is assembled accordingly:
%    H = kron(Kx_1d, Iy) + kron(Ix, Ky_1d) + diag(V(:)).
%

function Quantum2D_N

%% Parameters -------------------------------------------------------------
% Pade orders q to test
P.q                     =   1:5;

% Splitting scheme: 'Verlet', 'BM4' or 'BM6'
P.integrator            =   'BM4';

% Total propagation time and numbers of time steps
P.Tf                    =   4;
P.Nt                    =   2.^(0:15);

% Grid
P.dx                    =   0.1;
P.dy                    =   0.1;
P.Nx                    =   2^9;
P.Ny                    =   2^9;

% Initial Gaussian wavepacket
P.x0                    =   -5;
P.y0                    =   -5;
P.kx0                   =   2;
P.ky0                   =   2;
P.sx                    =   2;
P.sy                    =   2;

% Non-separable Gaussian potential
P.omegax                =   0.5;
P.omegay                =   0.5;
P.lambda                =   0.3;

%% Initialization ---------------------------------------------------------
P.h                     =   P.Tf ./ P.Nt;

x                       =   (-(P.Nx-1)/2 : (P.Nx-1)/2) * P.dx;
y                       =   (-(P.Ny-1)/2 : (P.Ny-1)/2) * P.dy;
[X, Y]                  =   meshgrid(x, y);              % Ny x Nx

psi0                    =   exp( -0.5*((X-P.x0)/P.sx).^2 ...
                                 -0.5*((Y-P.y0)/P.sy).^2 ) ...
                          .* exp( 1i*(P.kx0*X + P.ky0*Y) );
psi0                    =   psi0 ./ sqrt(sum(abs(psi0(:)).^2)*P.dx*P.dy);

P.Vmat                  =   exp(- P.omegax*X.^2 - P.omegay*Y.^2 + P.lambda*(X.*Y));

% Tridiagonal kinetic operators: sub-diagonal a, diagonal b, super-diagonal c
P.ax                    =   -0.5/P.dx^2 * ones(P.Nx-1,1);
P.bx                    =   +1.0/P.dx^2 * ones(P.Nx,1);
P.cx                    =   P.ax;
P.ay                    =   -0.5/P.dy^2 * ones(P.Ny-1,1);
P.by                    =   +1.0/P.dy^2 * ones(P.Ny,1);
P.cy                    =   P.ay;

% Palindromic splitting sequence
beta                    =   splitting_coefficients(P.integrator);
beta_full               =   [beta(:); flipud(beta(:))];

%% Reference solution -----------------------------------------------------
disp('> Computing reference solution with expmv...');
Kx_1d                   =   sparse( [1:P.Nx, 2:P.Nx,   1:P.Nx-1], ...
                                     [1:P.Nx, 1:P.Nx-1, 2:P.Nx  ], ...
                                     [P.bx',  P.ax',     P.cx'   ], P.Nx, P.Nx );
Ky_1d                   =   sparse( [1:P.Ny, 2:P.Ny,   1:P.Ny-1], ...
                                     [1:P.Ny, 1:P.Ny-1, 2:P.Ny  ], ...
                                     [P.by',  P.ay',     P.cy'   ], P.Ny, P.Ny );
N2D                     =   P.Nx * P.Ny;
H_2D                    =   kron(Kx_1d, speye(P.Ny)) + kron(speye(P.Nx), Ky_1d) ...
                          + spdiags(P.Vmat(:), 0, N2D, N2D);
tic;
u_exact                 =   reshape(expmv(-1i*H_2D, psi0(:), P.Tf), P.Ny, P.Nx);
fprintf('  done (%.2f s)\n', toc);

%% Propagation ------------------------------------------------------------
err                     =   zeros(length(P.q), length(P.h));
norm_dev                =   zeros(length(P.q), length(P.h));
comp_time               =   zeros(length(P.q), length(P.h));

for iq                  =   1:length(P.q)
    alpha               =   coefficients(P.q(iq));
    fprintf('> Pade [%d,%d]\n', P.q(iq), P.q(iq));

    for ih              =   1:length(P.h)
        h               =   P.h(ih);
        psi             =   psi0;

        tic;
        for it          =   1:P.Nt(ih)
            for s       =   1:length(beta_full)
                hs      =   beta_full(s) * h;
                if mod(s,2) == 1
                    psi =   apply_Kx(psi, P, alpha, hs);
                    psi =   apply_Ky(psi, P, alpha, hs);
                    psi =   apply_V (psi, P, hs);
                else
                    psi =   apply_V (psi, P, hs);
                    psi =   apply_Ky(psi, P, alpha, hs);
                    psi =   apply_Kx(psi, P, alpha, hs);
                end
            end
        end
        comp_time(iq,ih)=   toc;

        d               =   psi - u_exact;
        err(iq,ih)      =   sqrt(sum(abs(d(:)).^2)*P.dx*P.dy);
        norm_dev(iq,ih) =   abs(1 - sum(abs(psi(:)).^2)*P.dx*P.dy);
        fprintf('  dt = %-12g error = %.3e   |1-norm| = %.3e   (%.2f s)\n', ...
                h, err(iq,ih), norm_dev(iq,ih), comp_time(iq,ih));
    end
end

%% Display ----------------------------------------------------------------
cmap                    =   copper(length(P.q));
r                       =   splitting_order(P.integrator);

figure(1); clf;
subplot(2,1,1); hold on;
    for i               =   1:length(P.q)
        plot(P.h, err(i,:), 'o', 'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:));
        % slope r of the splitting, for every q (Proposition 2)
        plot(P.h, err(i,end)*(P.h/P.h(end)).^r, '-', 'Color', cmap(i,:));
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\delta t (a.u.)'); ylabel('error (a.u.)');
    title([P.integrator, ' + Pade [q,q]']);
subplot(2,1,2); hold on;
    for i               =   1:length(P.q)
        plot(P.h, norm_dev(i,:), 'o', 'Color', cmap(i,:), 'MarkerFaceColor', cmap(i,:));
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\delta t (a.u.)'); ylabel('|1 - norm|');

save(['data2D_' P.integrator '.mat'], 'P', 'err', 'norm_dev', 'comp_time', 'beta_full');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
function psi            =   apply_Kx(psi, P, alpha, h)
%APPLY_KX  psi <- U_q(h,Kx) psi, on all Ny rows at once.
psi                     =   pade_lines(psi.', P.ax, P.bx, P.cx, alpha, h).';

%--------------------------------------------------------------------------
function psi            =   apply_Ky(psi, P, alpha, h)
%APPLY_KY  psi <- U_q(h,Ky) psi, on all Nx columns at once.
psi                     =   pade_lines(psi, P.ay, P.by, P.cy, alpha, h);

%--------------------------------------------------------------------------
function psi            =   apply_V(psi, P, h)
%APPLY_V  psi <- exp(-1i*h*V) psi (exact, elementwise).
psi                     =   psi .* exp(-1i * h * P.Vmat);

%--------------------------------------------------------------------------
function D              =   pade_lines(D, a, b, c, alpha, h)
%PADE_LINES  Factored Pade propagator applied to every column of D.
%  The tridiagonal operator (a, b, c) acts along the first dimension of D;
%  each column is an independent one-dimensional problem.
n                       =   length(b);
for k                   =   1:length(alpha)
    fac                 =   -1i*h*alpha(k);
    D                   =   tridiag_mv(fac*a, ones(n,1)+fac*b, fac*c, D);
    fac                 =   +1i*h*alpha(k);
    D                   =   thomas(fac*a, ones(n,1)+fac*b, fac*c, D);
end

%--------------------------------------------------------------------------
function Y              =   tridiag_mv(a, b, c, V)
%TRIDIAG_MV  Y = A*V for tridiagonal A, applied to every column of V.
Y                       =   b .* V;
Y(1:end-1,:)            =   Y(1:end-1,:) + c .* V(2:end,:);
Y(2:end,:)              =   Y(2:end,:)   + a .* V(1:end-1,:);

%--------------------------------------------------------------------------
function D              =   thomas(a, b, c, D)
%THOMAS  Solve A*X = D for tridiagonal A, for every column of D.
%  Thomas algorithm, O(n) per column. All columns share the same matrix.
n                       =   length(b);
for i                   =   2:n
    w                   =   a(i-1) / b(i-1);
    b(i)                =   b(i) - w * c(i-1);
    D(i,:)              =   D(i,:) - w * D(i-1,:);
end
D(n,:)                  =   D(n,:) / b(n);
for i                   =   n-1:-1:1
    D(i,:)              =   (D(i,:) - c(i) * D(i+1,:)) / b(i);
end

%--------------------------------------------------------------------------
function alpha          =   coefficients(q)
%COEFFICIENTS  alpha_q^(k) = -1/z_q^(k), z_q^(k) roots of P_q (Table 1).
switch q
    case 1
        alpha           =   0.5;
    case 2
        alpha           =   [(3-1i*sqrt(3))/12, ...
                              (3+1i*sqrt(3))/12];
    case 3
        alpha           =   [0.215314423116112178, ...
                              0.142342788441943911 - 0.135799925708153803*1i, ...
                              0.142342788441943911 + 0.135799925708153803*1i];
    case 4
        alpha           =   [0.15843375973428236 - 0.04744101257110844*1i, ...
                              0.15843375973428236 + 0.04744101257110844*1i, ...
                              0.09156624026571764 - 0.11566261301312761*1i, ...
                              0.09156624026571764 + 0.11566261301312761*1i];
    case 5
        alpha           =   [0.13710881296519226, ...
                              0.11742725435969863 - 0.06104970381853820*1i, ...
                              0.11742725435969863 + 0.06104970381853820*1i, ...
                              0.06401833915770524 - 0.09834106917310919*1i, ...
                              0.06401833915770524 + 0.09834106917310919*1i];
    case 6
        alpha           =   [0.112981487353213770 - 0.023070678358901355*1i, ...
                              0.112981487353213770 + 0.023070678358901355*1i, ...
                              0.089573203698748516 - 0.062971624731703068*1i, ...
                              0.089573203698748516 + 0.062971624731703068*1i, ...
                              0.047445308948037714 - 0.084722574097167073*1i, ...
                              0.047445308948037714 + 0.084722574097167073*1i];
    case 7
        alpha           =   [0.100567464822504836, ...
                              0.092694772766323216 - 0.033882488943913687*1i, ...
                              0.092694772766323216 + 0.033882488943913687*1i, ...
                              0.070328696979243325 - 0.060773909166175561*1i, ...
                              0.070328696979243325 + 0.060773909166175561*1i, ...
                              0.036692797843181041 - 0.074059703707308632*1i, ...
                              0.036692797843181041 + 0.074059703707308632*1i];
    case 8
        alpha           =   [0.087372903654552227 - 0.0135661308688463049*1i, ...
                              0.087372903654552227 + 0.0135661308688463049*1i, ...
                              0.076688973859411999 - 0.038547151769477275*1i, ...
                              0.076688973859411999 + 0.038547151769477275*1i, ...
                              0.056629165024864069 - 0.057227482069531063*1i, ...
                              0.056629165024864069 + 0.057227482069531063*1i, ...
                              0.029308957461171706 - 0.065596184872823900*1i, ...
                              0.029308957461171706 + 0.065596184872823900*1i];
    case 9
        alpha           =   [0.0794026483914611415, ...
                              0.075504785130370233 - 0.0214077325481443664*1i, ...
                              0.075504785130370233 + 0.0214077325481443664*1i, ...
                              0.064202205264604698 - 0.040073605399966999*1i, ...
                              0.064202205264604698 + 0.040073605399966999*1i, ...
                              0.046581435869845285 - 0.053398589951849901*1i, ...
                              0.046581435869845285 + 0.053398589951849901*1i, ...
                              0.0240102495394492130 - 0.058762852440401573*1i, ...
                              0.0240102495394492130 + 0.058762852440401573*1i];
    otherwise
        error('Pade order q = %d not tabulated (q = 1, ..., 9).', q);
end

%--------------------------------------------------------------------------
function beta           =   splitting_coefficients(name)
%SPLITTING_COEFFICIENTS  Half-sequence {beta_s} of the chi/chi* composition
%  (Table 2). Valid for any number of operators S (McLachlan, 2022).
switch name
    case 'Verlet'
        beta            =   0.5;
    case 'BM4'
        beta            =   [ 0.0792036964311957; ...
                               0.1303114101821663; ...
                               0.2228614958676077; ...
                              -0.3667132690474257; ...
                               0.3246481886897062; ...
                               0.1096884778767498];
    case 'BM6'
        beta            =   [ 0.050262764400392; ...
                               0.098553683500650; ...
                               0.314960616927694; ...
                              -0.447346482695478; ...
                               0.492426372489876; ...
                              -0.425118767797691; ...
                               0.237063913978122; ...
                               0.195602488600053; ...
                               0.346358189850727; ...
                              -0.362762779254345];
    otherwise
        error('Unknown integrator %s (use Verlet, BM4 or BM6).', name);
end

%--------------------------------------------------------------------------
function r              =   splitting_order(name)
%SPLITTING_ORDER  Order r of the splitting scheme.
switch name
    case 'Verlet',  r   =   2;
    case 'BM4',     r   =   4;
    case 'BM6',     r   =   6;
end