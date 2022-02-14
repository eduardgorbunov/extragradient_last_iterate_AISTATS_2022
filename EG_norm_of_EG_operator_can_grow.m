clear all; clc;
%
% In this example, we consider Extragradient method:
% w_{k+1} = w_k - gamma_2 * F(w_k - gamma_1 * F(w_k))
% with gamma_2 = gamma_1
%
% In particular, we show that ||F_{EG}(w_k)|| can be larger than
% ||F_{EG}(w_{k-1})|| even if F is cocoercive! We observe this phenomenon
% for gamma_1 = gamma_2 and for gamma_1 > gamma_2




% (0) Initialize an empty PEP
P=pep();

% (1) Set up the class of monotone inclusions
%param.L  =  1.0; param.mu = 0; % F is 1-Lipschitz and 0-strongly monotone
param.beta = 1.0;

%gamma1 = 1/(2*param.L);
gamma1 = param.beta/2;
gamma2 = gamma1/2;

%F = P.DeclareFunction('LipschitzStronglyMonotone',param);
F = P.DeclareFunction('Cocoercive',param);

% (2) Set up the starting points
w0=P.StartingPoint();
[ws, Fs] = F.OptimalPoint(); 

P.InitialCondition((ws-w0)^2<=1.0);  % Normalize the initial distance ||w0-ws||^2 <= 1

% (3) Algorithm

w12 = w0 - gamma1 * F.evaluate(w0);
w1 = w0 - gamma2 * F.evaluate(w12);
w22 = w1 - gamma1 * F.evaluate(w1);

% (4) Set up the performance measure: ||F_{EG}(w_k)||^2 - ||F_{EG}(w_{k-1})||^2
squared_norm_EG_diff = (F.evaluate(w22))^2 - (F.evaluate(w12))^2;
P.PerformanceMetric(squared_norm_EG_diff);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(squared_norm_EG_diff)   % worst-case squared norm