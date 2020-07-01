% Multi-focus noisy image fusion based on LRR
addpath('./lrr');
addpath(genpath('./methods'));

% multi-scale transform
methods = {'lrr_wave2', ...
        'shearlet', 'lrr_shearlet', 'nsst', 'lrr_nsst', ...
        'contourlet', 'lrr_contourlet', 'nsct', 'lrr_nsct'};
    
% noise type
noise_lambda={{'_gau_0005', '4.5'}, {'_gau_001', '3'}, {'_gau_005', '1'}, {'_gau_01', '1'}, ...
        {'_sp_01', '1.5'}, {'_sp_02', '1'},...
        {'_poi','2'}};

unit = 16; % windwos size
index = 1; % 通过index选择使用的融合方法和噪声类型
fusion_method = methods{index}; % lrr_wave2
noise_label = noise_lambda{index}{1}; % '_gau_0005'
lam = str2double(noise_lambda{index}{2}); % \lambda of LRR for different noises 4.500

disp([fusion_method,'-',noise_label]);
%for i=1:10
for i=1

    image_left = ['./mf_noise_images/mf_noise_images/image',num2str(i),noise_label,'_left.png'];
    image_right = ['./mf_noise_images/mf_noise_images/image',num2str(i),noise_label,'_right.png'];
    % 选择源图像
    sourceTestImage1 = imread(image_left);
    sourceTestImage2 = imread(image_right);
    %sourceTestImage3= imread

    tic
    eval(['fusion_image_LRR = ', fusion_method, '(sourceTestImage1,sourceTestImage2,lam, unit, 0);']);
    % 这里相当于执行“fusion_image_LRR = lrr_wave2(sourceTestImage1,sourceTestImage2,lam, unit, 0);”
    toc
    fused_path = ['./fused_images/fused',num2str(i),noise_label, '_', fusion_method, '_unit',num2str(unit), '.png'];
    imwrite(fusion_image_LRR,fused_path,'png');

end