% choose strategy for High freauency - LRR
% for wavelet ����С���任�ĸ�Ƶϵ���ں�
% lambda=4.5000,unit=16
function coe_f = fusionHighCoe(coe1,coe2,lambda, unit)

[m,n]=size(coe1);
coe_f = zeros(m,n);
% reshape(coe2(:,ii), [unit unit])
for ii=1:n
    % reshape��16*16
    coe1_t = reshape(coe1(:,ii), [unit unit]); 
    coe2_t = reshape(coe2(:,ii), [unit unit]);
    [Z1,~] = solve_lrr(coe1_t, coe1_t, lambda,0,1); % �������
    [Z2,~] = solve_lrr(coe2_t, coe2_t, lambda,0,1);
    %[Z3,~] = 
    LRR1 = sum(svd(Z1));
    LRR2 = sum(svd(Z2));  
    % choose-max ѡ��˷������ͼ��飨patch��
    % �˷�����������ֵ���
    % �����Z1��Z2��ͨ��IALM�õ���ͼ�����ֵ
    if LRR1>LRR2
        coe_t = reshape(coe1(:,ii), [unit unit])*Z1;
    else
        coe_t = reshape(coe2(:,ii), [unit unit])*Z2;
    end
    coe_f(:,ii) = coe_t(:);

end

end