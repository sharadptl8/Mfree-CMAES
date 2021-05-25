tic;
clc;
clear all;
close all;
% rng default;
lower=15;
upper=150; 
Maxit=800;
Parameter=10; % No. of decision variables (in this case this is 10 transmissivity values for 5 zones)
popmean1 = randsample(lower:upper,Parameter);    % objective variables initial point
popmean=popmean1';
sigma = 0.5;          % coordinate wise standard deviation (step size)
stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
stopeval = 1e3*Parameter^2;   % stop after stopeval number of function evaluations
  
% Strategy parameter setting: Selection  
lambda = 4+floor(3*log(Parameter));  % population size, offspring number
mu = lambda/2;               % number of parents/points for recombination
weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
mu = floor(mu);        
weights = weights/sum(weights);     % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4 + mueff/Parameter) / (Parameter+4 + 2*mueff/Parameter); % time constant for cumulation for C
cs = (mueff+2) / (Parameter+mueff+5);  % t-const for cumulation for sigma control
c1 = 2 / ((Parameter+1.3)^2+mueff);    % learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((Parameter+2)^2+mueff));  % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(Parameter+1))-1) + cs; % damping for sigma 
                                                      % usually close to 1
% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(Parameter,1); ps = zeros(Parameter,1);   % evolution paths for C and sigma
B = eye(Parameter,Parameter);                       % B defines the coordinate system
D = ones(Parameter,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2 
eigeneval = 0;                      % track update of B and D
chiN=Parameter^0.5*(1-1/(4*Parameter)+1/(21*Parameter^2));  % expectation of 
                                      %   ||Parameter(0,I)|| == norm(randn(Parameter,1))

iter=1; % initializing the iteration number
counteval = 0;
%% Loop for PSO
for it=1:Maxit
    % Generate and evaluate lambda offspring
    for k=1:lambda,
      pop(:,k) = popmean + sigma * B * (D .* randn(Parameter,1)); % m + sig * Normal(0,C) 
      counteval = counteval+1;
    end
    pop;
     [obj,object]=objfn(pop,iter); % objective function call
    
    % Sort by fitness and compute weighted mean into xmean
    [object, index] = sort(obj);  % minimization
    popold = popmean;
    popmean = pop(:,index(1:mu)) * weights;  % recombination, new mean value
    
    % Cumulation: Update evolution paths
    ps = (1-cs) * ps ... 
          + sqrt(cs*(2-cs)*mueff) * invsqrtC * (popmean-popold) / sigma; 
    hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/Parameter < 2 + 4/(Parameter+1);
    pc = (1-cc) * pc ...
          + hsig * sqrt(cc*(2-cc)*mueff) * (popmean-popold) / sigma; 

    % Adapt covariance matrix C
    artmp = (1/sigma) * (pop(:,index(1:mu)) - repmat(popold,1,mu));  % mu difference vectors
    C = (1-c1-cmu) * C ...                   % regard old matrix  
         + c1 * (pc * pc' ...                % plus rank one update
                 + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
         + cmu * artmp * diag(weights) * artmp'; % plus rank mu update 

    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
    
    % Update B and D from C
    if counteval - eigeneval > lambda/(c1+cmu)/Parameter/10  % to achieve O(Parameter^2)
      eigeneval = counteval;
      C = triu(C) + triu(C,1)'; % enforce symmetry
      [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
      D = sqrt(diag(D));        % D contains standard deviations now
      invsqrtC = B * diag(D.^-1) * B';
    end
[obj,object]=objfn(pop,iter); % Calculating the objective funtion value for
... the new target vector 

objectvalue(:,iter)=object(:,iter); % Saving obj. func. value for the 
... corresponding iteration number

iter=iter+1
end
toc
plot (iter,obj);
boxplot(pop') % "boxplot" is a command that plots the box plot; transpose 
... of pop (pop') is used because the box plot has to be made between its 
    ... own parameter values and not against different parameter values
xlswrite('fem_PSO_12.xlsx', objectvalue, 'Sheet1', 'A1')
xlswrite('fem_PSO_12.xlsx', pop, 'Sheet2', 'A1')
xlswrite('fem_PSO_12.xlsx', toc, 'Sheet3', 'A1')
