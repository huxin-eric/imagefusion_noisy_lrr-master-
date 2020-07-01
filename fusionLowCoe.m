% choose strategy for Low freauency - SF
function coe_f = fusionLowCoe(coe1,coe2, type, unit)
% type = 0
% unit = 16

[m,n]=size(coe1);% 256*56
if type==1
    coe_f = ones(unit)/2;
else
    coe_f = zeros(m,n); % 输出初始化
    for i=1:n
        
        var_A2_x1 = variance_block(reshape(coe1(:,i), [unit unit]));
        % 把每一列256个元素reshape为16*16
        % 分别计算16*16的框的空间频率
        var_A2_x2 = variance_block(reshape(coe2(:,i), [unit unit]));
        
        % 取空间频率大的那个16*16的低频图像块得到融合的低频系数
        if var_A2_x1>var_A2_x2
            coe_f(:,i) = coe1(:,i);
        else
            coe_f(:,i) = coe2(:,i);
        end
    end
    
end
end