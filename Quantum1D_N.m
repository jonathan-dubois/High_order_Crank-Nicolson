
%% exp_tridiagonal --------------------------------------------------------
%
%  Code in order to evaluate 
%
%       u = exp(-1i*A) v
%
%  with A a tridiagonal matrix and v a vector. We want to solve it using a
%  factorization
%
%       u = ( exp(-1i*h*A) )^nh v
%
%  where nh*h = 1. We do so by doing a truncation of the exponential
%  series
%
%       exp(B) = prod_n^N (1 + a_n B)
%
%  where the coefficients a_n are such that it converges to the exp.
%
%  The prototype code also compares the trick
%
%       ( exp(1i*h*A/2) )^nh u = ( exp(-1i*h*A/2) )^nh v
%
%  to enforce the norm conservation. The inverse is computed using the
%  Thomas algorithm.
%

function Quantum1D_N

%% Parameters -------------------------------------------------------------
% order of the approximation
P.N                     =   1:5;

% number of time steps for the highest h
P.Tf                    =   4;
P.Nt                    =   2.^(0:15);

% grid step
P.dx                    =   .1;
P.Nx                    =   2^9;

% initial conditions
P.x0                    =   -5;
P.E0                    =   2;
P.s0                    =   2;

%% Initialization ---------------------------------------------------------
% disp time step and integration time
P.h                     =   P.Tf ./ P.Nt;
disp(['> time step:     ', num2str(P.h)]);
disp(['> step number:   ', num2str(P.Nt)]);

% grid discretization
x                       =   (-(P.Nx-1)/2 : (P.Nx-1)/2) * P.dx;
disp(['> xmax:          ', num2str(max(x))]);

% initial condition
psi0                    =   transpose( exp(-.5*(x-P.x0).^2/P.s0^2) .* exp(+1i*x*sqrt(2*P.E0)) );
psi0                    =   psi0 ./ sqrt(sum(abs(psi0).^2)*P.dx);

% Hamilton operator
%P.V                     =   @(x)    - exp(-x.^2);
P.V                     =   @(x)    zeros(size(x));
P.a                     =   -.5/P.dx^2 * ones(P.Nx-1,1);
P.b                     =   +1/P.dx^2 * ones(P.Nx,1) + P.V(x');
P.c                     =   P.a;

% data to store
xnum                    =   zeros(length(P.N),length(P.h));
Enum                    =   zeros(length(P.N),length(P.h));
err                     =   zeros(length(P.N),length(P.h));
norm_vec                =   zeros(length(P.N),length(P.h));
psi_boundary            =   zeros(length(P.N),length(P.h));
comp_time               =   zeros(length(P.N),length(P.h));

% Sparse tridiagonal Hamiltonian & exact solution
i_sub                   =   (2:P.Nx)';
i_sup                   =   (1:P.Nx-1)';
i_main                  =   (1:P.Nx)';
H_sparse                =   sparse([i_main;i_sub;i_sup], ...
                                [i_main;i_sub-1;i_sup+1], ...
                                [P.b;P.a;P.c],P.Nx,P.Nx); 

% exact solution
disp(' ')
disp('> Exact solution via expmv');
tic;
u_exact                 =   expmv(-1i*H_sparse,psi0,P.Tf);
toc;

figure(2);
        plot(x,abs(psi0).^2,'k--','LineWidth',2);
    hold on;
        plot(x,abs(u_exact).^2,'k-','LineWidth',2);
drawnow;

%% Main -------------------------------------------------------------------
disp(' ');
disp('> Numerical solutions');
for iN                  =   1:length(P.N)
    % alpha_n coefficients
    disp([' >> N = ', num2str(P.N(iN))]);
    alpha               =   coefficients(P.N(iN));
    for ih              =   1:length(P.h)
        disp(['     >> h = ', num2str(P.h(ih))]);
        h               =   P.h(ih);
        u_app           =   psi0;
        tic;
        for i           =   1:P.Nt(ih)
            for j       =   1:P.N(iN)                                               
                fac     =   -1i*h*alpha(j);
                u_app   =   tridiag_mv(fac*P.a,ones(P.Nx,1)+fac*P.b,fac*P.c,u_app);
    
                fac     =   1i*h*alpha(j);
                u_app   =   thomas(fac*P.a,ones(P.Nx,1)+fac*P.b,fac*P.c,u_app); 
                % The inverse and multiplication commute => it can be made
                % sequentially
            end
        end
        comp_time(iN,ih)=   toc;
        err(iN,ih)      =   norm(u_app - u_exact) * sqrt(P.dx);
        Enum(iN,ih)     =   real(u_app'*tridiag_mv(P.a,P.b,P.c,u_app))*P.dx;
        xnum(iN,ih)     =   sum(abs(u_app).^2.*x'.*P.dx);
        norm_vec(iN,ih) =   sum(abs(u_app).^2.*P.dx);
        psi_boundary(iN,ih) ...
                        =   abs(u_app(1)) + abs(u_app(end));
    end
end

% display
cmap                    =   colormap(copper(length(P.N)));
figure(1);
for i                   =   1:length(P.N)
    set(gcf,Color='w');
    subplot(4,1,1);
        hold on;
        plot(P.h,err(i,:),'ko-','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:));
        set(gca,Color='w',XColor='k',YColor='k',ZColor='k',YScale='log',XScale='log');
    subplot(4,1,2);
        hold on;
        plot(P.h,abs(1-norm_vec(i,:)),'ko-','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:));
        set(gca,Color='w',XColor='k',YColor='k',ZColor='k',YScale='log',XScale='log');

    subplot(4,1,3);
        hold on;
        plot(P.h,psi_boundary(i,:),'ko-','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:));
        set(gca,Color='w',XColor='k',YColor='k',ZColor='k',YScale='log',XScale='log');

    subplot(4,1,4);
        hold on;
        plot(P.h,comp_time(i,:),'ko-','Color',cmap(i,:),'MarkerFaceColor',cmap(i,:));
        set(gca,Color='w',XColor='k',YColor='k',ZColor='k',YScale='log',XScale='log');
end
drawnow;

save('data.mat','P','xnum','Enum','err','norm_vec','psi_boundary','comp_time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function y              =   tridiag_mv(a,b,c,v)
%TRIDIAG_MV  Compute y = A*v for tridiagonal A defined by (a,b,c).
%
%   Input:
%     a : sub-diagonal (n-1 x 1)
%     b : main diagonal (n x 1)
%     c : super-diagonal (n-1 x 1)
%     v : vector (n x 1)
%
%   Output:
%     y : A*v  (n x 1)

n                       =   length(b);
y                       =   zeros(n,1);

% Main contribution
y                       =   b .* v;

% Add sub- and super-diagonal terms
y(1:end-1)              =   y(1:end-1) + c .* v(2:end);
y(2:end)                =   y(2:end)   + a .* v(1:end-1);

%
function x              =   thomas(a,b,c,d)
%THOMAS Solve a tridiagonal system A*x = d using the Thomas algorithm.
%
%   Input:
%     a : sub-diagonal (n-1 x 1)
%     b : main diagonal (n x 1)
%     c : super-diagonal (n-1 x 1)
%     d : right-hand side (n x 1)
%
%   Output:
%     x : solution vector (n x 1)
n                       =   length(b);
x                       =   zeros(n,1);

% Forward elimination
for i                   =   2:n
    w                   =   a(i-1) / b(i-1);
    b(i)                =   b(i) - w * c(i-1);
    d(i)                =   d(i) - w * d(i-1);
end

% Back substitution
x(n)                    =   d(n) / b(n);
for i                   =   n-1:-1:1
    x(i)                =   (d(i) - c(i) * x(i+1)) / b(i);
end

%
function alpha          =   coefficients(order)
    
switch order
    case 1
        alpha           =   .5;
    case 2
        alpha           =   [(3-1i*sqrt(3))/12, ...
                                (3+1i*sqrt(3))/12];
    case 3
        alpha           =   [0.21531442311611, ...
                                0.14234278844194 - 0.13579992570815*1i, ...
                                0.14234278844194 + 0.13579992570815*1i];
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
        alpha           =   [0.112981487353213770 - 0.023070678358901355.*1i, ...  
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
        error('Order not coded.');
end
