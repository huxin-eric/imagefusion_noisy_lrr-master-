% choose strategy for High freauency - LRR
% for wavelet 用于小波变换的高频系数融合
% lambda=4.5000,unit=16
function coe_f = fusionHighCoe(coe1,coe2,lambda, unit)

[m,n]=size(coe1);
coe_f = zeros(m,n);
% reshape(coe2(:,ii), [unit unit])
for ii=1:n
    % reshape成16*16
    coe1_t = reshape(coe1(:,ii), [unit unit]); 
    coe2_t = reshape(coe2(:,ii), [unit unit]);
    [Z1,~] = solve_lrr(coe1_t, coe1_t, lambda,0,1); % 这个输入
    [Z2,~] = solve_lrr(coe2_t, coe2_t, lambda,0,1);
    %[Z3,~] = 
    LRR1 = sum(svd(Z1));
    LRR2 = sum(svd(Z2));  
    % choose-max 选择核范数大的图像块（patch）
    % 核范数就是奇异值求和
    % 这里的Z1和Z2是通过IALM得到的图像估计值
    if LRR1>LRR2
        coe_t = reshape(coe1(:,ii), [unit unit])*Z1;
    else
        coe_t = reshape(coe2(:,ii), [unit unit])*Z2;
    end
    coe_f(:,ii) = coe_t(:);

end

end