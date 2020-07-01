% choose strategy for Low freauency - SF
function coe_f = fusionLowCoe(coe1,coe2, type, unit)
% type = 0
% unit = 16

[m,n]=size(coe1);% 256*56
if type==1
    coe_f = ones(unit)/2;
else
    coe_f = zeros(m,n); % �����ʼ��
    for i=1:n
        
        var_A2_x1 = variance_block(reshape(coe1(:,i), [unit unit]));
        % ��ÿһ��256��Ԫ��reshapeΪ16*16
        % �ֱ����16*16�Ŀ�Ŀռ�Ƶ��
        var_A2_x2 = variance_block(reshape(coe2(:,i), [unit unit]));
        
        % ȡ�ռ�Ƶ�ʴ���Ǹ�16*16�ĵ�Ƶͼ���õ��ںϵĵ�Ƶϵ��
        if var_A2_x1>var_A2_x2
            coe_f(:,i) = coe1(:,i);
        else
            coe_f(:,i) = coe2(:,i);
        end
    end
    
end
end