function [obj,object]=objfn(pop,iter)
% iterations=1;
Parameter=10; % no. of parameters to optimize
lambda= 4+floor(3*log(Parameter)); % population size (i.e no. of sets)

% observed head values after 365 days


% noise=0+2.5*(randn(49,1)); % to generate random numbers of mean zero and
... std. deviation 2.5
    
% hobs=noise+hobs; % observed head values with noise

for p=1:lambda
    % Using the target vector (variable name "pop") from "DEalgo", head 
    ... values are simulated using FEM simulator by calling function call 
        ... "femdiff"
[hsimu]=Mfreediff(pop,p);

hsim(:,p)=hsimu(:,p); % saving the simulated head values "hsimu" with 
... different variable name "hsim"
    
for ob=1:117
   
end
end

obj(p,1)=sum(objvalue(:,p)); % add objective function values of all the observation wells

end
object(:,iter)=obj(:,1); % saves the objective function values of all the iterations correspondingly
end
