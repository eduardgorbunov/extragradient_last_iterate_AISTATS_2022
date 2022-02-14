clear all; clc;
%
% In this example, we consider Extragradient method:
% w_{k+1} = w_k - gamma_2 * F(w_k - gamma_1 * F(w_k))
% 
% The goal is to find worst-case ratio ||x1 - y1||^2/||x0-y0||,
% where x1 is a result of 1 iteration of EG from x0 and
%       y1 is a result of 1 iteration of EG from y0
%
% The results are visualized in EG_showing_non_cocoercivity.ipynb



% (0) Initialize an empty PEP

beta = 1.0;

for gamma1 = linspace(beta/30, beta, 100)
   for gamma2 = [gamma1, gamma1/2, gamma1/4, gamma1/8, gamma1/10]
       P=pep();
       
       % (1) Set up the class of monotone inclusions
       %param.L  =  1; param.mu = 0; % F is 1-Lipschitz and 0-strongly monotone
       param.beta = beta;
       
       F = P.DeclareFunction('Cocoercive',param);
       
       % (2) Set up the starting points
       x0=P.StartingPoint();
       y0=P.StartingPoint();
       
       P.InitialCondition((x0-y0)^2==1);  % Normalize the initial distance ||w0-ws||^2 <= 1
       
       % (3) Algorithm
       x_F1 = F.evaluate(x0);
       y_F1 = F.evaluate(y0);
       x_F2 = F.evaluate(x0 - gamma1 * x_F1);
       y_F2 = F.evaluate(y0 - gamma1 * y_F1);
       
       % (4) Set up the performance measure: ||x1-y1||^2
       squared_norm = (x0 - gamma2 * x_F2 - y0 + gamma2 * y_F2)^2;
       P.PerformanceMetric(squared_norm);

       % (5) Solve the PEP
       P.solve()

       % (6) Evaluate the output
       double(squared_norm)   % worst-case squared norm

       res_norm = double(squared_norm);
       res_x0 = double(x0);
       res_y0 = double(y0);
       res_x_F1 = double(x_F1);
       res_y_F1 = double(y_F1);
       res_x_F2 = double(x_F2);
       res_y_F2 = double(y_F2);

       save(strcat('dump/EG_expansiveness_1e-1', sprintf('_%f_', gamma1), sprintf('%f', gamma2),'.mat'), 'res_norm', 'gamma1', 'gamma2', 'res_x0', 'res_y0', 'res_x_F1', 'res_y_F1', 'res_x_F2', 'res_y_F2');

       fprintf("======================================================\n")
       fprintf("gamma_1 = %20f, ", gamma1)
       fprintf("gamma_2 = %20f\n", gamma2)
       fprintf("||x - gamma * x_F2 - y + gamma * y_F2||^2 =  %20f\n", double(squared_norm));
   end
end

P=pep();

% (1) Set up the class of monotone inclusions
%param.L  =  1; param.mu = 0; % F is 1-Lipschitz and 0-strongly monotone
param.beta = 1.0;

%gamma1 = 1/(10*param.L);
gamma1 = param.beta/10;
gamma2 = gamma1;

%F = P.DeclareFunction('LipschitzStronglyMonotone',param);
F = P.DeclareFunction('Cocoercive',param);

% (2) Set up the starting points
x0=P.StartingPoint();
y0=P.StartingPoint();

P.InitialCondition((x0-y0)^2==1);  % Normalize the initial distance ||w0-ws||^2 <= 1

% (3) Algorithm

x_F1 = F.evaluate(x0);
y_F1 = F.evaluate(y0);
x_F2 = F.evaluate(x0 - gamma1 * x_F1);
y_F2 = F.evaluate(y0 - gamma1 * y_F1);

%P.AddConstraint((y0)^2 == 0);
%P.AddConstraint((x_F1 - y_F2)^2 == 0);
%P.AddConstraint((y_F1 - y0)^2 == 0);

%x12 = x0 - gamma1 * F.evaluate(x0);
%x1 = x0 - gamma2 * F.evaluate(x12);
%y12 = y0 - gamma1 * F.evaluate(y0);
%y1 = y0 - gamma2 * F.evaluate(y12);

% (4) Set up the performance measure: ||x1-y1||^2
squared_norm = (x0 - gamma2 * x_F2 - y0 + gamma2 * y_F2)^2;
P.PerformanceMetric(squared_norm);
%P.TraceHeuristic(1)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(squared_norm)   % worst-case squared norm

res_norm = double(squared_norm);
gamma = gamma1;
res_x0 = double(x0);
res_y0 = double(y0);
res_x_F1 = double(x_F1);
res_y_F1 = double(y_F1);
res_x_F2 = double(x_F2);
res_y_F2 = double(y_F2);

save('EG_expansiveness_1.mat', 'res_norm', 'gamma', 'res_x0', 'res_y0', 'res_x_F1', 'res_y_F1', 'res_x_F2', 'res_y_F2');

fprintf("||x - gamma * x_F2 - y + gamma * y_F2||^2 =  %20f\n", double(squared_norm));
%fprintf("x0 =  %10f\n", double(x0));
%fprintf("y0 =  %10f\n", double(y0));
%fprintf("x_F1 =  %10f\n", double(x_F1));
%fprintf("y_F1 =  %10f\n", double(y_F1));
%fprintf("x_F2 =  %10f\n", double(x_F2));
%fprintf("y_F2 =  %10f\n", double(y_F2));