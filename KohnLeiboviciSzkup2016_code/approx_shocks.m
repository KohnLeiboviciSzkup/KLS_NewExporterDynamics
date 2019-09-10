% 11-22-2010

% Function to approximate continuous process, iid(mu,sigma^2) or AR(1)(mu,sigma^2,rho<1), with normal distribution as a
% discrete markov process

% Takes as input:
% n (number of states of discrete process)
% rho (AR(1) coefficient), if =0, then iid process
% mu (unconditional mean of process y)
% sigma (variance of white noise)
% r (range of values for discrete process)

% IID process 
% y(t)~ N(0,sigma^2)
% Algorithm 1 in Nakajima (2007) -page 2-

% AR(1) process 
% y(t)= (1-rho)*mu + rho y(t-1) + epsilon(t)
% epsilon ~ N(0,sigma^2)
% rho < 1 (stationary process)
% Algorithm 3 in Nakajima (2007) -page 5-, based on Tauchen (1986)

% Returns:
% (z_1,...,z_n) Vector of values of states (increasing)
%  P(nxn) Matrix of transition probabilities



function [z,P,pi] = approx_shocks(n,c,rho,sigma,mu)

% Solve for sigma_z, SDev of stationaty process y (and discrete process z)

sigma_z=sigma/sqrt(1-rho^2);

% Solve for c, the number of S_Dev such that z_n-z_1=r
%c=r/(2*sigma_z);
r=(2*sigma_z)*c;

% Solve for values of discrete state space (equally distanced from each
% other)

z_1=mu-c*sigma_z;
e=(0:1:n-1)';
z=z_1+e*(r/(n-1));

% (The two lines above are a more efficient way to code what follows...)
% for i=1:n
%     z(i)=z_1+ (r/(n-1))*(i-1);
% end

% Constructing midpoints m(i), i=1...n-1

m=(z(2:n)+z(1:n-1))/2;
% (The two lines above are a more efficient way to code what follows...)
% for i=1:n-1
% m(i)=(z(i+1)+z(i))/2;
% end

% Define intervals Z(i) for i=1...n
% Z(1) = (-inf,m(1)]
% Z(i) = (m(i),m(i+1)], i=2,3,... n-1
% Z(n) = (m(n),inf)

% Approximate the transition probability P(i,j) as the probability that, conditional on z(i),
% z'=(1-rho)*mu+rho*z(i)+epsilon' falls into interval j (as defined above)

P=zeros(n,n);

for i=1:n
    P(i,1)=normcdf((m(1)-(1-rho)*mu-rho*z(i))/sigma);
       for j=2:n-1
            P(i,j)=normcdf((m(j)-(1-rho)*mu-rho*z(i))/sigma)-normcdf((m(j-1)-(1-rho)*mu-rho*z(i))/sigma);
        
       end
    P(i,n)=1-normcdf((m(n-1)-(1-rho)*mu-rho*z(i))/sigma);
end


% Obtaining Limiting distribution

% pi'*P=pi' => P'*pi=pi => (P'-I)*pi=0

evalc('[V,D]=eigs(P'',1)');
pi=V/sum(V); % limiting distribution (eigenvector associated to eigenvalue = 1 -the largest -)




