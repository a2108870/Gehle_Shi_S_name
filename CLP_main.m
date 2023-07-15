load('.\utils\REC_groundtruth.mat');
real_rgb=REC_groundtruth;
LIST_FILES = load('.\utils\all_files_name_shi.mat');
%%
tic
now_ours_errors = [];
LSRS_errors = [];
for ii=1:size(real_rgb,1)
    current_illuminant = real_rgb(ii,:);
    current_illuminant = current_illuminant ./ norm(current_illuminant);
    name=['.\Gehle_Shi\',LIST_FILES{ii},'.tiff'];
    
    input_data = [];
    I=im2double(imread(name));
    
    I=I./max(I(:));
    mask_im2=(max(I,[],3)<=0.95);
    I=I.*double(mask_im2);
    
    I_xyz = [];
    I_xyz(:,:,1) = I(:,:,1)*0.4124564 + I(:,:,2)*0.3575761 + I(:,:,3)*0.1804375;
    I_xyz(:,:,2) = I(:,:,1)*0.2126729 + I(:,:,2)*0.7151522 + I(:,:,3)*0.0721750;
    I_xyz(:,:,3) = I(:,:,1)*0.0193339 + I(:,:,2)*0.1191920 + I(:,:,3)*0.9503041;
    
    current_illuminant_xyz = [current_illuminant(1)*0.4124564 + current_illuminant(2)*0.3575761 + current_illuminant(3)*0.1804375 ,
        current_illuminant(1)*0.2126729 + current_illuminant(2)*0.7151522 + current_illuminant(3)*0.0721750,
        current_illuminant(1)*0.0193339 + current_illuminant(2)*0.1191920 + current_illuminant(3)*0.9503041];
    current_illuminant_xyz_chorm_x = current_illuminant_xyz(1)./sum(current_illuminant_xyz);
    current_illuminant_xyz_chorm_y = current_illuminant_xyz(2)./sum(current_illuminant_xyz);
    
    I_x_chorm_tep = I_xyz(:,:,1) ./ sum(I_xyz,3);
    I_y_chorm_tep = I_xyz(:,:,2) ./ sum(I_xyz,3);
    I_x_chorm_tep(isnan(I_x_chorm_tep)) = 0;
    I_y_chorm_tep(isnan(I_y_chorm_tep)) = 0;
    % ori_params
    %     N_superpixels_num = 500;
    N_superpixels_num = 500;
    [I_xyz_superpixels,N] = superpixels(I_xyz,N_superpixels_num);
    
    I_x_chorm = [];
    I_y_chorm = [];
    for num_sp = 1:N_superpixels_num
        if ~isempty(I_x_chorm_tep(I_xyz_superpixels==num_sp)) ==1
            I_x_chorm = [I_x_chorm, mean(I_x_chorm_tep(I_xyz_superpixels==num_sp))];
            I_y_chorm = [I_y_chorm, mean(I_y_chorm_tep(I_xyz_superpixels==num_sp))];
        end
    end
    
    data = [I_x_chorm'; I_y_chorm'];
    %%
    % 输入数据：I_x_chorm, I_y_chorm为点集
    
    % 参数设置
    maxIterations = 1000; % 最大迭代次数
    thresholdDistance = 0.002; % 点到拟合直线的距离阈值
    inlierRatioThreshold = 0.05; % 内点比例阈值
    
    numPoints = numel(I_x_chorm); % 点的数量
    bestInlierCount = 0; % 最佳内点数量
    bestModel = struct(); % 最佳模型（直线）
    
    % RANSAC算法迭代
    for iteration = 1:maxIterations
        % 随机选择两个点作为直线的两个端点
        indices = randperm(numPoints, 2);
        x1 = I_x_chorm(indices(1));
        y1 = I_y_chorm(indices(1));
        x2 = I_x_chorm(indices(2));
        y2 = I_y_chorm(indices(2));
        
        % 根据两个点计算直线参数（斜率和截距）
        slope = (y2 - y1) / (x2 - x1);
        intercept = y1 - slope * x1;
        
        % 计算所有点到直线的距离
        distances = abs(I_y_chorm - (slope * I_x_chorm + intercept));
        
        % 统计落在距离阈值内的内点数量
        inlierIndices = find(distances <= thresholdDistance);
        inlierCount = numel(inlierIndices);
        
        % 判断内点数量是否超过阈值，并更新最佳模型和内点数量
        if inlierCount > inlierRatioThreshold * numPoints && inlierCount > bestInlierCount
            bestModel.slope = slope;
            bestModel.intercept = intercept;
            bestInlierCount = inlierCount;
        end
    end
    
    % 绘制拟合直线
    x = min(I_x_chorm):0.005:max(I_x_chorm);
    y = bestModel.slope * x + bestModel.intercept;
    
    %%
    % LSRS
    [estimated_ill_LSRS]=LSRS(I);
    %%
    estimated_ill_LSRS = estimated_ill_LSRS./norm(estimated_ill_LSRS);
    estimated_ill_LSRS_xyz = [estimated_ill_LSRS(1)*0.4124564 + estimated_ill_LSRS(2)*0.3575761 + estimated_ill_LSRS(3)*0.1804375 ,
        estimated_ill_LSRS(1)*0.2126729 + estimated_ill_LSRS(2)*0.7151522 + estimated_ill_LSRS(3)*0.0721750,
        estimated_ill_LSRS(1)*0.0193339 + estimated_ill_LSRS(2)*0.1191920 + estimated_ill_LSRS(3)*0.9503041];
    estimated_ill_LSRS_xyz_chorm_x = estimated_ill_LSRS_xyz(1)./sum(estimated_ill_LSRS_xyz);
    estimated_ill_LSRS_xyz_chorm_y = estimated_ill_LSRS_xyz(2)./sum(estimated_ill_LSRS_xyz);
    
    x0=x(1); x1=x(end);
    y0=y(1); y1=y(end);
    
    % 计算直线斜率和截距
    slope = (y1 - y0) / (x1 - x0);
    intercept = y0 - slope * x0;
    
    % 计算投影点的坐标
    projection_x = (estimated_ill_LSRS_xyz_chorm_x + slope * estimated_ill_LSRS_xyz_chorm_y - slope * intercept) / (slope^2 + 1);
    projection_y = slope * projection_x + intercept;
    
    %%
    est_ill_xyz = [projection_x, projection_y, 1-projection_x-projection_y];
    estimated_ill_rgb = [est_ill_xyz(1)*3.2404542 + est_ill_xyz(2)*-1.5371385 + est_ill_xyz(3)*-0.4985314 ,
        est_ill_xyz(1)*-0.9692660 + est_ill_xyz(2)*1.8760108 + est_ill_xyz(3)*0.0415560 ,
        est_ill_xyz(1)*0.0556434 + est_ill_xyz(2)*-0.2040259 + est_ill_xyz(3)*1.0572252]';
    ours_estimated_illuminant = estimated_ill_rgb./norm(estimated_ill_rgb);
    
    if ~(isnan(acos(ours_estimated_illuminant * current_illuminant') * (180/pi)))
        now_ours_errors = [now_ours_errors, acos(ours_estimated_illuminant * current_illuminant') * (180/pi)];
        LSRS_errors = [LSRS_errors, acos(estimated_ill_LSRS' * current_illuminant') * (180/pi)];
    end
    disp(['当前第',num2str(ii),'幅图像，LSRS误差为', num2str(acos(estimated_ill_LSRS' * current_illuminant') * (180/pi)),...
        ',修正后误差为', num2str(acos(ours_estimated_illuminant * current_illuminant') * (180/pi)),...
        ',LSRS平均误差为', num2str(mean(LSRS_errors)),...
        ',平均修正后误差为',num2str(mean(now_ours_errors))]);
end
[mn,md,tr,bst,wst,Q95,mx]=reportError_right(now_ours_errors);
toc


function [mn,md,tr,bst,wst,Q95,mx]=reportError_right(ae)
mn = mean(ae);
md = median(ae);
mx = max(ae);
ae_sort = sort(ae);
bst = mean(ae_sort(ceil(1:0.25*length(ae_sort))));
wst = mean(ae_sort(ceil(0.75*length(ae_sort)):end));
tr = (md+(ae_sort(ceil(0.25*length(ae_sort)))+ae_sort(ceil(0.75*length(ae_sort))))/2)/2;
Q95 = prctile(ae,95);
fprintf('Mean: %.2f, Median: %.2f, Trimean: %.2f, Best-25 %2.2f, Worst-25 %2.2f, Quantile-95: %2.2f, Max: %.2f\n', mn,md,tr,bst,wst,Q95,mx);
end