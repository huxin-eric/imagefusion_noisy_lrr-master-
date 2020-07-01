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
if nargin < 6 || isempty(display)  %�����������������
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
% norm�����ɼ��㼸�ֲ�ͬ���͵ľ�����
% 'fro��  A��A���Ļ��ĶԽ��ߺ͵�ƽ��������sqrt(sum(diag(A'*A)))
norm_f_A = norm(A,'fro');

if norm_f_X>0 && norm_f_A>0
    Q = orth(A');  %���ؾ���A'��������
    B = A*Q;% ����A��A��ת�õ�������

    if reg==0
        if alm_type ==0
            [Z,E] = exact_alm_lrr_l21v2(X,B,lambda,[],[],display);
        else
            [Z,E] = inexact_alm_lrr_l21(X,B,lambda,display); %��ǰ���������㷨
            % X=A��B=A*Q��Q�Ǿ���A��ת�õ�������
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