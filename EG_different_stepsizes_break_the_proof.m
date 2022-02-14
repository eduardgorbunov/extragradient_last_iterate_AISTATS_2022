clear all; clc;
%
% In this example, we consider Extragradient method:
% w_{k+1} = w_k - gamma_2 * F(w_k - gamma_1 * F(w_k))
% with gamma_2 = gamma_1
%
% In particular, we show that when gamma_1 and gamma_2 are not the same,
% then it is possible that ||F(w_k)|| < ||F(w_{k+1})|| even if F is
% cocoercive! This example illustrates the importance of having
% gamma_1 = gamma_2 for having ||F(w_{k+1})|| <= ||F(w_k)||




% (0) Initialize an empty PEP
P=pep();

% (1) Set up the class of monotone inclusions
%param.L  =  1.0; param.mu = 0; % F is 1-Lipschitz and 0-strongly monotone
param.beta = 1.0;

%gamma1 = 1/(2*param.L);
gamma1 = param.beta/2;
gamma2 = gamma1/4; %stepsizes are different

%F = P.DeclareFunction('LipschitzStronglyMonotone',param);
F = P.DeclareFunction('Cocoercive',param);

% (2) Set up the starting points
w0=P.StartingPoint();
[ws, Fs] = F.OptimalPoint(); 

P.InitialCondition((ws-w0)^2<=1.0);  % Normalize the initial distance ||w0-ws||^2 <= 1

% (3) Algorithm

w12 = w0 - gamma1 * F.evaluate(w0);
w1 = w0 - gamma2 * F.evaluate(w12);

% (4) Set up the performance measure: ||F(w_k)||^2 - ||F(w_{k-1})||^2
squared_norm_diff = (F.evaluate(w1))^2 - (F.evaluate(w0))^2;
P.PerformanceMetric(squared_norm_diff);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(squared_norm_diff)   % worst-case squared norm