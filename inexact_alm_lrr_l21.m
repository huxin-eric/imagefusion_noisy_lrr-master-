function [Z,E] = inexact_alm_lrr_l21(X,A,lambda,display)
% This routine uses Inexact ALM algorithm to solve the following nuclear-norm optimization problem:
% min |Z|_*+lambda*|E|_2,1
% s.t., X = AZ+E
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        A -- D*M matrix of a dictionary, M is the size of the dictionary


if nargin<4
    display = false;
end
% display=fase,lambda=4.5
tol = 1e-8;
maxIter = 1e6;
[d, n] = size(X); % 16*16
m = size(A,2); % 16
rho = 1.1;
max_mu = 1e10;
mu = 1e-6;
atx = A'*X;
inv_a = inv(A'*A+eye(m));
%% Initializing optimization variables
% intialize
J = zeros(m,n);
Z = zeros(m,n);
E = sparse(d,n);  %全零稀疏矩阵: d×n
% S=sparse(X)—将矩阵X转化为稀疏矩阵的形式，
% 即矩阵X中任何零元素去除，非零元素及其下标（索引）组成矩阵S。

Y1 = zeros(d,n);
Y2 = zeros(m,n);
%% Start main loop
iter = 0;
if display
    disp(['initial,rank=' num2str(rank(Z))]);
end
while iter<maxIter
    iter = iter + 1;
    %update J
    temp = Z + Y2/mu;
    [U,sigma,V] = svd(temp,'econ');%奇异值分解[U, S, V] = svd(M);  则 U(m,m), S(m, n),  V(n, n).
    % 对于矩阵A(m*n)，存在U(m*m)，V(n*n)，S(m*n)，满足A = U*S*V’。
    % U和V中分别是A的奇异向量，而S是A的奇异值。
    % AA'的正交单位特征向量组成U，特征值组成S'S，
    % A'A的正交单位特征向量组成V，特征值（与AA'相同）组成SS'。
    % 核范数指的是A的奇异值之和
    sigma = diag(sigma); % 16*1
    svp = length(find(sigma>1/mu)); 
    %  find(X)找出矩阵X中的所有非零元素，并将这些元素的线性索引值（linear indices：按列）返回到向量ind中。
    if svp>=1
        sigma = sigma(1:svp)-1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)'; 
    %udpate Z
    Z = inv_a*(atx-A'*E+J+(A'*Y1-Y2)/mu);
    %update E
    xmaz = X-A*Z;
    temp = xmaz+Y1/mu;
    E = solve_l1l2(temp,lambda/mu);
    
    leq1 = xmaz-E;
    leq2 = Z-J;
    stopC = max(max(max(abs(leq1))),max(max(abs(leq2)))); % 两个参数误差最大值
    if display && (iter==1 || mod(iter,50)==0 || stopC<tol)
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
            ',rank=' num2str(rank(Z,1e-3*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
    end
    if stopC<tol  % 误差最大小于阈值，则循环终止
        break;
    else
        % 参数更新
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        mu = min(max_mu,mu*rho);
    end
end
