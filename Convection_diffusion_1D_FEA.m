% This program solves a one-dimensional convection-diffusion equation
%        a u_x - nu u_xx = f
% with Dirichlet boundary conditions using the finite element method
% with some stabilized formulations and with linear elements
% by Tariq Ridwan: tariq.ridwan@bsc.es && https://tariqridwan.github.io/
close all;
clear; clc
disp('This program solves a convection-diffusion equation: a u_x - nu u_xx = f')
disp('with 0<x<1 and essential boundary conditions on both ends')
disp('and linear or quadratic elements.')
disp('One of the following problems can be solved at once:')
disp('  [1]: boundary conditions: u(0)= 0, u(1) = 1. f = 0') % 0 source
disp('  [2]: boundary conditions: u(0)= 0, u(1) = 0. f = 1') % 1 source
disp('  [3]: boundary conditions: u(0)= 0, u(1) = 1. f = sin(pi x)');
disp('  [4]: boundary conditions: u(0)= 0, u(1) = 1. f = 10*exp(-5*x)-4*exp(-x)');
disp('  [5]: boundary conditions: u(0)= 0, u(1) = 1. f = 20*exp(-5*(x-1/8))-10*exp(-5*(x-1/4))');
problem = cinput('Select a problem to be solved:', 2);
disp('Choose element type: [1] Linear, [2] Quadratic');
ElemType = cinput('Select a problem to be solved:', 1);
%% Parameters
a = 1; % a constant, default = 1
% nu = 0.01; % viscosity, default = 0.01
nu = cinput('Diffusion coefficient nu', 0.01);
gamma = a/nu;
Ldom = [0,1]; % domain length
if ElemType == 1 % linear element
    % nElem = 10; % number of Element
    nElem = cinput('Number of elements',10);
    nPt = nElem + 1; % number of points, imax = 11
    h = ( Ldom(2)-Ldom(1) ) / nElem; % Distance between nodes / Each element's length
elseif ElemType == 2 % quadratic element
    % nElem = 5; % number of Element
    nElem = cinput('Number of elements',5);
    nPt = nElem*2 + 1; % number of points, imax = 11
    h = ( Ldom(2)-Ldom(1) ) / nElem/2; % Distance between nodes; Each element's length = 2*h
end
x = ( Ldom(1) : h : Ldom(2) ); % Location of x points
Pe = a*h/(2*nu); % Péclet number
disp(['Peclet number: ',num2str(Pe)]);
%% Dirichet boundary conditions
BC_D_loc = [1,nPt]; % Locations for Dirichet boundary conditions
if problem == 2
    BC_D = [0;0];
else
    BC_D = [0;1]; 
end
if ElemType == 1 % linear element
%% Method used for solving the problem
disp ('The problem can be solved using one of the following methods: ');
disp ('[0] Galerkin, [-1] FU, [1] SU, [2] SUPG, [3] GLS, [4] SGS');
method = input('Method = ');
%% 1. Element convection matrix
nGaussPt = 3; % number of gauss points
G_N1_dN1_dx = 0;
G_N2_dN1_dx = 0;
if nGaussPt == 2
    w_i = [1, 1];
    x_i = [sqrt(1/3), -sqrt(1/3)];
    for i = 1:nGaussPt
        G_N1_dN1_dx = G_N1_dN1_dx + w_i(i) * (1 - x_i(i));
        G_N2_dN1_dx = G_N2_dN1_dx + w_i(i) * (1 + x_i(i));
    end
    G_N1_dN2_dx = G_N1_dN1_dx;
    G_N2_dN2_dx = G_N2_dN1_dx;
elseif nGaussPt == 3
    w_i = [8/9, 5/9, 5/9];
    x_i = [0, sqrt(3/5), -sqrt(3/5)];
    for i = 1:nGaussPt
        G_N1_dN1_dx = G_N1_dN1_dx + w_i(i) * (1 - x_i(i));
        G_N2_dN1_dx = G_N2_dN1_dx + w_i(i) * (1 + x_i(i));
    end
    G_N1_dN2_dx = G_N1_dN1_dx;
    G_N2_dN2_dx = G_N2_dN1_dx;
end
N1_dN1_dx = -1/4 * G_N1_dN1_dx;
N1_dN2_dx = 1/4 * G_N1_dN2_dx;
N2_dN1_dx = -1/4 * G_N2_dN1_dx;
N2_dN2_dx = 1/4 * G_N2_dN2_dx;
Ce = a*[N1_dN1_dx, N1_dN2_dx;...
        N2_dN1_dx, N2_dN2_dx];
%% 2. Element diffusion matrix
dN1_dx_dN1_dx = 1/h;
dN1_dx_dN2_dx = -1/h;
dN2_dx_dN1_dx = -1/h;
dN2_dx_dN2_dx = 1/h;
if method == 0 % Galerkin
    beta = 0;
elseif method == -1 % Full-upwind
    beta = 1;
elseif method == 1 || method == 2 || method == 3 || method == 4 % SU / SUPG / GLS / SGS
    beta = coth(Pe) - 1/Pe;
end
Ke = (nu+beta*a*h/2)*[dN1_dx_dN1_dx, dN1_dx_dN2_dx; dN2_dx_dN1_dx, dN2_dx_dN2_dx];
CKe = Ce + Ke; % Element convection (Ce) + diffusion (Ke) matrix
%% 3. Element source matrix
if problem == 1 % source = 0, for all methods
    Fe = [0; 0];
elseif problem == 2 % source = 1
    if method == 0 || method == -1 || method == 1 % Galerkin / FU / SU
        Fe = [h/2; h/2];
    elseif method == 2 || method == 3 || method == 4 % SUPG / GLS / SGS
        Fe = [h/2 - beta*h/2; h/2 + beta*h/2]; % beta*h/2 = tau*a
    end
elseif problem == 3 || problem == 4 || problem == 5 % source = sin(pi*x) OR 10*exp(-5*x)-4*exp(-x) ...
    Fe = zeros(2,nElem);
    for n = 1:nElem
        x1 = x(n);
        x2 = x(n+1);
        if problem == 3
            s1 = sin(pi*x1);
            s2 = sin(pi*x2);
        elseif problem == 4
            s1 = 10*exp(-5*x1)-4*exp(-x1);
            s2 = 10*exp(-5*x2)-4*exp(-x2);
        elseif problem == 5
            s1 = 20*exp(-5*(x1-1/8))-10*exp(-5*(x1-1/4));
            s2 = 20*exp(-5*(x2-1/8))-10*exp(-5*(x2-1/4));
        end
        if method == 0 || method == -1 || method == 1 % Galerkin / FU / SU
            Fe(:,n) = [h/3*s1 + h/6*s2;... N1s
                       h/6*s1 + h/3*s2]; % N2s
        elseif method == 2 || method == 3 || method == 4 % SUPG / GLS / SGS
            tau = beta*h/(2*a);
            Fe(:,n) = [(h/3 - tau*a/2)*s1 + (h/6 - tau*a/2)*s2;... N1s
                       (h/6 + tau*a/2)*s1 + (h/3 + tau*a/2)*s2]; % N2s
        end
    end
end
%% Global convection, diffusion, & source matrix
CKg = zeros(nPt);
Fg = zeros(nPt,1);
if problem == 3 || problem == 4 || problem == 5 % source = sin(pi*x) OR 10*exp(-5*x)-4*exp(-x)
    for n = 1:nElem
        CKnew = zeros(nPt);
        CKnew(n:n+1,n:n+1) = CKe;
        CKg = CKg + CKnew; % Global convection, diffusion matrix
        Fnew = zeros(nPt,1);
        Fnew(n:n+1) = Fe(:,n);
        Fg = Fg + Fnew; % Global source matrix
    end
else % for constant source-term, problem == 1,2
    for n = 1:nElem
        CKnew = zeros(nPt);
        CKnew(n:n+1,n:n+1) = CKe;
        CKg = CKg + CKnew; % Global convection, diffusion matrix
        Fnew = zeros(nPt,1);
        Fnew(n:n+1) = Fe;
        Fg = Fg + Fnew; % Global source matrix
    end
end % Linear element calculation ends
elseif ElemType == 2 % quadratic element
%% Method used for solving the problem
disp ('The problem can be solved using one of the following methods: ');
disp ('[0] Galerkin, [1] SU, [2] SUPG');
method = input('Method = ');
%% 1. Element convection matrix
N1_dN1_dx = -1/2;
N1_dN2_dx = 2/3;
N1_dN3_dx = -1/6;
N2_dN1_dx = -2/3;
N2_dN2_dx = 0;
N2_dN3_dx = 2/3;
N3_dN1_dx = 1/6;
N3_dN2_dx = -2/3;
N3_dN3_dx = 1/2;
Ce = a*[N1_dN1_dx, N1_dN2_dx, N1_dN3_dx;...
        N2_dN1_dx, N2_dN2_dx, N2_dN3_dx;...
        N3_dN1_dx, N3_dN2_dx, N3_dN3_dx];
%% 2. Element diffusion matrix
dN1_dx_dN1_dx = 7/6/h;
dN1_dx_dN2_dx = -4/3/h;
dN1_dx_dN3_dx = 1/6/h;
dN2_dx_dN1_dx = -4/3/h;
dN2_dx_dN2_dx = 8/3/h;
dN2_dx_dN3_dx = -4/3/h;
dN3_dx_dN1_dx = 1/6/h;
dN3_dx_dN2_dx = -4/3/h;
dN3_dx_dN3_dx = 7/6/h;
if method == 0 % Galerkin
    beta = 0; beta_c = 0;
elseif method == 1 % SU
    beta = coth(Pe) - 1/Pe;
    beta_c = 2*((coth(Pe)-1/Pe)-(cosh(Pe))^2*(coth(2*Pe)-1/(2*Pe)))/(2-(cosh(Pe))^2);
elseif method == 2 % SUPG
    beta = coth(Pe) - 1/Pe;
    beta_c = ((2*Pe-1)*exp(3*Pe)+(-6*Pe+7)*exp(Pe)+(-6*Pe-7)*exp(-Pe)+(2*Pe+1)*exp(-3*Pe))...
            /( (Pe+3)*exp(3*Pe)+(-7*Pe-3)*exp(Pe)+ (7*Pe-3)*exp(-Pe)- (Pe+3)*exp(-3*Pe));
end
Ke = nu*[dN1_dx_dN1_dx, dN1_dx_dN2_dx, dN1_dx_dN3_dx;...
         dN2_dx_dN1_dx, dN2_dx_dN2_dx, dN2_dx_dN3_dx;...
         dN3_dx_dN1_dx, dN3_dx_dN2_dx, dN3_dx_dN3_dx]...
   +(a*h/2)*[beta_c*dN1_dx_dN1_dx, beta_c*dN1_dx_dN2_dx, beta_c*dN1_dx_dN3_dx;...
             beta  *dN2_dx_dN1_dx, beta  *dN2_dx_dN2_dx, beta  *dN2_dx_dN3_dx;...
             beta_c*dN3_dx_dN1_dx, beta_c*dN3_dx_dN2_dx, beta_c*dN3_dx_dN3_dx];
if method == 0 || method == 1 % Galerkin / SU
    CKe = Ce + Ke; % Element convection (Ce)+diffusion (Ke)
elseif method == 2 % SUPG
    dN1_dx_d2N1_dx2 = -1/h^2;
    dN1_dx_d2N2_dx2 = 2/h^2;
    dN1_dx_d2N3_dx2 = -1/h^2;
    dN2_dx_d2N1_dx2 = 0;
    dN2_dx_d2N2_dx2 = 0;
    dN2_dx_d2N3_dx2 = 0;
    dN3_dx_d2N1_dx2 = 1/h^2;
    dN3_dx_d2N2_dx2 = -2/h^2;
    dN3_dx_d2N3_dx2 = 1/h^2;
    nu_u_xx = nu*h/2*[beta_c*dN1_dx_d2N1_dx2, beta_c*dN1_dx_d2N2_dx2, beta_c*dN1_dx_d2N3_dx2;...
                      beta  *dN2_dx_d2N1_dx2, beta  *dN2_dx_d2N2_dx2, beta  *dN2_dx_d2N3_dx2;...
                      beta_c*dN3_dx_d2N1_dx2, beta_c*dN3_dx_d2N2_dx2, beta_c*dN3_dx_d2N3_dx2];
    CKe = Ce + Ke - nu_u_xx; % Element convection (Ce)+diffusion (Ke)+nu_u_xx matrix
end
%% 3. Element source matrix
inc = 0;
if problem == 1 % source = 0, for all methods
    Fe = [0; 0; 0];
elseif problem == 2 % source = 1
    if method == 0 || method == 1 % Galerkin / SU
        Fe = [h/3; 4*h/3; h/3];
    elseif method == 2 % SUPG
        Fe = [h/3 - beta*h/2; 4*h/3; h/3 + beta*h/2]; % beta*h/2 = tau*a
    end
elseif problem == 3 || problem == 4 || problem == 5 % source = sin(pi*x) OR 10*exp(-5*x)-4*exp(-x) ...
    Fe = zeros(3,nElem);
    for n = 1:nElem
        x1 = x(n+inc);
        x2 = x(n+1+inc);
        x3 = x(n+2+inc);
        if problem == 3
            s1 = sin(pi*x1);
            s2 = sin(pi*x2);
            s3 = sin(pi*x3);
        elseif problem == 4
            s1 = 10*exp(-5*x1)-4*exp(-x1);
            s2 = 10*exp(-5*x2)-4*exp(-x2);
            s3 = 10*exp(-5*x3)-4*exp(-x3);
        elseif problem == 5
            s1 = 20*exp(-5*(x1-1/8))-10*exp(-5*(x1-1/4));
            s2 = 20*exp(-5*(x2-1/8))-10*exp(-5*(x2-1/4));
            s3 = 20*exp(-5*(x3-1/8))-10*exp(-5*(x3-1/4));
        end
        if method == 0 || method == 1 % Galerkin / SU
            Fe(:,n) = [4*h/15*s1 +  2*h/15*s2 -   h/15*s3;... N1s
                       2*h/15*s1 + 16*h/15*s2 + 2*h/15*s3;... N1s
                        -h/15*s1 +  2*h/15*s2 + 4*h/15*s3]; % N3s
        elseif method == 2 % SUPG
            tau = beta*h/(2*a);
            Fe(:,n) = [(4*h/15-  tau*a/2)*s1 + (2*h/15-2*tau*a/3)*s2 + ( -h/15+  tau*a/6)*s3;... N1s
                       (2*h/15+2*tau*a/3)*s1 + 16*h/15           *s2 + (2*h/15-2*tau*a/3)*s3;... N1s
                       ( -h/15-  tau*a/6)*s1 + (2*h/15+2*tau*a/3)*s2 + (4*h/15+  tau*a/2)*s3]; % N3s
        end
        inc = inc + 1;
    end
end
%% Global convection, diffusion, & source matrix
CKg = zeros(nPt);
Fg = zeros(nPt,1);
inc = 0;
if problem == 3 || problem == 4 || problem == 5 % source = sin(pi*x) OR 10*exp(-5*x)-4*exp(-x)
    for n = 1:nElem
        CKnew = zeros(nPt);
        CKnew(n+inc:n+2+inc,n+inc:n+2+inc) = CKe;
        CKg = CKg + CKnew; % Global convection, diffusion matrix
        Fnew = zeros(nPt,1);
        Fnew(n+inc:n+2+inc) = Fe(:,n);
        Fg = Fg + Fnew; % Global source matrix
        inc = inc + 1;
    end
else % for constant source-term, problem == 1,2
    for n = 1:nElem
        CKnew = zeros(nPt);
        CKnew(n+inc:n+2+inc,n+inc:n+2+inc) = CKe;
        CKg = CKg + CKnew; % Global convection, diffusion matrix
        Fnew = zeros(nPt,1);
        Fnew(n+inc:n+2+inc) = Fe;
        Fg = Fg + Fnew; % Global source matrix
        inc = inc + 1;
    end
end % Quadratic element calculation ends
end
%% Calculation of u
Fg_Final = Fg - CKg(:,BC_D_loc)*BC_D;
CKeff = CKg(2:nPt-1,2:nPt-1);
Feff = Fg_Final(2:nPt-1);
ueff = CKeff\Feff;
u = zeros(nPt,1);
u(BC_D_loc) = BC_D;
u(2:nPt-1) = ueff;
%% Analytical solution
n_x_ex = 100; % number of divisions to divide L
dx_ex = Ldom(2)/n_x_ex; % length of the divisions
x_ex = 0:dx_ex:Ldom(2); % locations of the divisions
if problem == 1 % source = 0
    u_ex = (1-exp(x_ex*gamma))/(1-exp(gamma));
elseif problem == 2 % source = 1
    u_ex = 1/a * (x_ex - (1 - exp(gamma*x_ex)) / (1 - exp(gamma)) );
elseif problem == 3 % source = sin(pi*x)
    aux = pi*(a^2+nu^2*pi^2);
    e = exp(a/nu);
    c1 = (-aux+a*(e+1))/(aux*(e-1));
    c2 = (aux-2*a)/(aux*(e-1));
    u_ex = c1 + c2*exp(a*x_ex/nu) + nu*pi*(sin(pi*x_ex)-a*cos(pi*x_ex)/(nu*pi))/aux;
elseif problem == 4 % source = 10*exp(-5*x)-4*exp(-x)
    c1 = ( (2-2*exp(-5))/(a+5*nu) - (4-4*exp(-1))/(a+nu) - 1 ) / (nu/a * (1-exp(a/nu)));
    c2 = 2/(a+5*nu) - 4/(a+nu) - ( (2-2*exp(-5))/(a+5*nu) - (4-4*exp(-1))/(a+nu) - 1 ) / (1-exp(a/nu));
    u_ex = c1*nu/a*exp(x_ex*a/nu) + 4*exp(-x_ex)/(a+nu) - 2*exp(-5*x_ex)/(a+5*nu) + c2;
elseif problem == 5 % source = 20*exp(-5*(x-1/8))-10*exp(-5*(x-1/4))
    u_ex = 2*exp(-5*x_ex+5/8)*(-2+exp(5/8))/(5*nu+a) - ...     
        (5*nu*exp(5/8)^7-2*exp(5/8) - 4*exp(a/nu+5) + 2*exp(a/nu+45/8) + 4 + a*exp(5/8)^7)/...
        ((5*nu+a)*exp(5/8)^7*(exp(a/nu)-1)) + ...
        (4-2*exp(5/8)+5*nu*exp(5/8)^7+a*exp(5/8)^7-4*exp(5/8)^8+2*exp(5/8)^9)*exp(a*x_ex/nu)/...
        (exp(5/8)^7*(5*exp(a/nu)*nu+exp(a/nu)*a-5*nu-a));
end
%% Error calculation
Error = zeros(nPt,1);
for n = 2:nPt-1
    ana = u_ex( (n-1)*(h*n_x_ex) +1); % analytical solution at x location
%     Error(n) = abs((u(n) - ana)/ana);
    Error(n) = u(n) - ana;
end
%% Post-processing
yyaxis left
plot(x,u,'r-o',x_ex,u_ex,'k:','LineWidth',3,'MarkerSize',10)
% ylim ([-0 1])
xlabel('$L$','Interpreter','latex','fontsize',12)
ylabel('$u$','Interpreter','latex','fontsize',12)
grid on
yyaxis right
plot(x,Error,'r:','LineWidth',2)
ylabel('Error','Interpreter','latex','fontsize',12)
if method == 0 % Galerkin
    l = legend('Galerkin','exact','Error','Location','Best');
elseif method == -1 % Full-upwind
    l = legend('FU','exact','Error','Location','Best');
elseif method == 1 % SU
    l = legend('SU','exact','Error','Location','Best');
elseif method == 2 % SUPG
    l = legend('SUPG','exact','Error','Location','Best');
elseif method == 3 % GLS
    l = legend('GLS','exact','Error','Location','Best');
elseif method == 4 % SGS
    l = legend('SGS','exact','Error','Location','Best');
end
if problem == 1
    title('$au_x-\nu u_{xx}=0,u(0)=0,u(L)=1$','Interpreter','latex','fontsize',14)
elseif problem == 2
    title('$au_x-\nu u_{xx}=1,u(0)=u(L)=0$','Interpreter','latex','fontsize',14)
elseif problem == 3
    title('$au_x-\nu u_{xx}=sin(\pi x), u(0)=0, u(L)=1$','Interpreter','latex','fontsize',14)
elseif problem == 4
    title('$au_x-\nu u_{xx}=10e^{-5x}-4e^{-x},u(0)=0,u(L)=1$','Interpreter','latex','fontsize',14)
elseif problem == 5
    title('$au_x-\nu u_{xx}=20e^{-5(x-1)/8}-10e^{-5(x-1)/4},u(0,L)=0,1$','Interpreter','latex','fontsize',14)
end
set(l,'FontSize',12)
set(gca,'FontSize',12);
gtext (strcat('Pe =   ',num2str(Pe)), 'FontSize',16) % place the Pe at your desired location
% v = axis; % ==== % Matteo
% text (v(1)+0.05, v(4)-2.35, ['Pe = ',num2str(Pe)], 'FontSize',16, 'VerticalAlignment','Top')
print FEAfluids1D.png -dpng
print FEAfluids1D.eps -depsc
print FEAfluids1D.pdf -dpdf
savefig('FEAfluids1D.fig')
% hold on
