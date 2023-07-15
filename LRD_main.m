load('F:\image_dataset\picture\REC_groundtruth.mat');
real_rgb=REC_groundtruth;
LIST_FILES = load('F:\image_dataset\picture\all_files_name_shi.mat');
LIST_FILES=LIST_FILES.name;%加载每幅图像的名称
%%
tic
now_ours_errors = [];
for ii=1:size(real_rgb,1)
    %for ii=105
    current_illuminant = real_rgb(ii,:);
    current_illuminant = current_illuminant ./ norm(current_illuminant);
    name=['F:\image_dataset\Cube+PREprocess\',num2str(ii),'.png'];

    input_data = [];
    I=im2double(imread(name));    
    I=I./max(I(:));
    mask_im2=(max(I,[],3)<=0.95);
    I=I.*double(mask_im2);
    
    [Rx,Gx,Bx]=norm_derivative(I, 6, 2);
    input_data(:,:,1)=Rx;
    input_data(:,:,2)=Gx;
    input_data(:,:,3)=Bx;
    
    % 将图像均匀地分割成 4x4 块    
    % ori_params
%     num_rows = 4;
%     num_cols = 4;
    num_rows = 200;
    num_cols = 200;
    
    block_size_rows = floor(size(input_data, 1) / num_rows);
    block_size_cols = floor(size(input_data, 2) / num_cols);
    blocks = cell(num_rows, num_cols);
    for r = 1:num_rows
        row_start = (r-1) * block_size_rows + 1;
        row_end = min(r * block_size_rows, size(input_data, 1));
        for c = 1:num_cols
            col_start = (c-1) * block_size_cols + 1;
            col_end = min(c * block_size_cols, size(input_data, 2));
            blocks{r, c} = input_data(row_start:row_end, col_start:col_end, :);
        end
    end
    
    % 对每个图像块进行归一化操作
    for r = 1:num_rows
        for c = 1:num_cols
            block = blocks{r, c};
            max_vals = max(max(block, [], 1), [], 2); % 计算 RGB 通道中的最大值
            if mean(max_vals) ==0
                blocks{r, c} = double(block);
            else
                blocks{r, c} = double(block) ./ repmat(max_vals, size(block,1), size(block,2)); % 归一化操作
            end
        end
    end
    
    % 计算图像全局梯度的范数值，p=10
    p=10;
    gobal = cal_norm(input_data, p);
    
    % 计算所有图像块梯度的范数值之和locals_rgb，p=10
    locals_rgb = [0,0,0];
    for r = 1:num_rows
        for c = 1:num_cols
            block = blocks{r, c};
            if sum(block(:)) == 0
                locals_rgb = locals_rgb  + 0;
            else
                local_rgb = cal_norm(block, p);
                locals_rgb = locals_rgb  + local_rgb;
            end
        end
    end
    
    ours_estimated_illuminant = gobal./locals_rgb;
    ours_estimated_illuminant = ours_estimated_illuminant./norm(ours_estimated_illuminant);
    if ~(isnan(acos(ours_estimated_illuminant * current_illuminant') * (180/pi)))
        now_ours_errors =[now_ours_errors, acos(ours_estimated_illuminant * current_illuminant') * (180/pi)];
    end
    disp(['当前第',num2str(ii),'幅图像，误差为', num2str(acos(ours_estimated_illuminant * current_illuminant') * (180/pi)),'平均误差为',num2str(mean(now_ours_errors))]);
end
[mn,md,tr,bst,wst,Q95,mx]=reportError_right(now_ours_errors);

toc

function use_rgb = cal_norm(input_data, mink_norm)
kleur = power(input_data,mink_norm);
white_R = power(sum(sum(kleur(:,:,1))),1/mink_norm);
white_G = power(sum(sum(kleur(:,:,2))),1/mink_norm);
white_B = power(sum(sum(kleur(:,:,3))),1/mink_norm);
som=sqrt(white_R^2+white_G^2+white_B^2);
white_R=white_R/som;
white_G=white_G/som;
white_B=white_B/som;
use_rgb = [white_R, white_G, white_B];
% use_rgb = use_rgb./norm(use_rgb);
end

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