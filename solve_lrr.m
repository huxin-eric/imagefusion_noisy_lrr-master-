function [Z,E] = solve_lrr(X,A,lambda,reg,alm_type,display)
% 
% Aug 2013
% This routine solves the following nuclear-norm optimization problem,
% min |Z|_*+lambda*|E|_L
% s.t., X = AZ+E
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        A -- D*M matrix of a dictionary, M is the size of the dictionary
%        lambda -- parameter
%        reg -- the norm chosen for characterizing E,
%            -- reg=0 (default),                  use the l21-norm
%            -- reg=1 (or ther values except 0),  use the l1-norm
%        alm_type -- 0 (default)   use the exact ALM strategy
%                 -- 1             use the inexact ALM strategy
%
if nargin < 6 || isempty(display)  %这里输入了五个参数
    display = false;
end
if nargin<5 || isempty(alm_type)
    alm_type = 0 ;
end

if nargin<4 || isempty(reg)
    reg = 0;
end
% solve_lrr(X,A,lambda,reg,alm_type,display)
% solve_lrr(X,A,4.5,0,1,false)
norm_f_X = norm(X,'fro');
% norm函数可计算几种不同类型的矩阵范数
% 'fro’  A和A‘的积的对角线和的平方根，即sqrt(sum(diag(A'*A)))
norm_f_A = norm(A,'fro');

if norm_f_X>0 && norm_f_A>0
    Q = orth(A');  %返回矩阵A'正交基。
    B = A*Q;% 矩阵A乘A的转置的正交基

    if reg==0
        if alm_type ==0
            [Z,E] = exact_alm_lrr_l21v2(X,B,lambda,[],[],display);
        else
            [Z,E] = inexact_alm_lrr_l21(X,B,lambda,display); %当前采用这种算法
            % X=A，B=A*Q，Q是矩阵A的转置的正交基
        end
    else
        if alm_type == 0
            [Z,E] = exact_alm_lrr_l1v2(X,B,lambda,[],[],display);
        else
            [Z,E] = inexact_alm_lrr_l1(X,B,lambda,display);
        end
    end
    Z = Q*Z;
else
    Z = 0;
    E = 0;
end