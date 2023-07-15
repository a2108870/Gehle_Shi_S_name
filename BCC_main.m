load('.\utils\REC_groundtruth.mat');
real_rgb=REC_groundtruth;
LIST_FILES = load('.\utils\all_files_name_shi.mat');
%%
tic
now_ours_errors = [];
for ii=1:size(real_rgb,1)
    %for ii=105
    current_illuminant = real_rgb(ii,:);
    current_illuminant = current_illuminant ./ norm(current_illuminant);
    name=['.\Gehle_Shi\',LIST_FILES{ii},'.tiff'];

    I=im2double(imread(name));

    I=I./max(I(:));
    mask_im2=(max(I,[],3)<=0.95&max(I,[],3)>=0.05);
    input_data=I.*double(mask_im2);
    
    GW_RGB = [mean(mean(input_data(:,:,1))),mean(mean(input_data(:,:,2))),mean(mean(input_data(:,:,3)))];
    GW_RGB = GW_RGB./norm(GW_RGB);
    input_data_correct = [];
    input_data_correct(:,:,1) = input_data(:,:,1)./GW_RGB(1);
    input_data_correct(:,:,2) = input_data(:,:,2)./GW_RGB(2);
    input_data_correct(:,:,3) = input_data(:,:,3)./GW_RGB(3);

    num_rows = 15;
    num_cols = 22;
    % ori_params
%     num_rows = 4;
%     num_cols = 7;

    block_size_rows = floor(size(input_data, 1) / num_rows);
    block_size_cols = floor(size(input_data, 2) / num_cols);
    blocks = cell(num_rows, num_cols);
    blocks_correct = cell(num_rows, num_cols);
    for r = 1:num_rows
        row_start = (r-1) * block_size_rows + 1;
        row_end = min(r * block_size_rows, size(input_data, 1));
        for c = 1:num_cols
            col_start = (c-1) * block_size_cols + 1;
            col_end = min(c * block_size_cols, size(input_data, 2));
            blocks{r, c} = input_data(row_start:row_end, col_start:col_end, :);
            blocks_correct{r, c} = input_data_correct(row_start:row_end, col_start:col_end, :);
        end
    end

    C_p = [0,0,0];
    C_p_back_up = [0,0,0];
    for r = 1:num_rows
        for c = 1:num_cols
            block = blocks{r, c};
            block_correct = blocks_correct{r, c};
            max_vals = reshape(max(max(block, [], 1), [], 2),1,3); % 计算 RGB 通道中的最大值
            max_vals = max_vals./norm(max_vals);
            u_p = reshape(mean(block_correct,[1,2]),1,3);
            u_p = u_p./norm(u_p);
            ill_patch = max_vals ./ u_p;
            ill_patch = ill_patch ./ norm(ill_patch);
            mean_vals_reshape = reshape(mean(block,[1,2]),1,3);
            mean_vals_reshape = mean_vals_reshape./norm(mean_vals_reshape);

            if (mean(max_vals)~=0 && mean(u_p)~=0) && acos(mean_vals_reshape * GW_RGB') * (180/pi)<=5
                C_p = C_p + ill_patch;
            end
            C_p_back_up = reshape(max_vals,1,3);
        end
    end
    
    if mean(C_p) == 0
        ours_estimated_illuminant = C_p_back_up;
    else
        ours_estimated_illuminant = C_p;
    end
    ours_estimated_illuminant = ours_estimated_illuminant./norm(ours_estimated_illuminant);
    if ~(isnan(acos(ours_estimated_illuminant * current_illuminant') * (180/pi)))
        now_ours_errors =[now_ours_errors, acos(ours_estimated_illuminant * current_illuminant') * (180/pi)];
    end
    disp(['当前第',num2str(ii),'幅图像，误差为', num2str(acos(ours_estimated_illuminant * current_illuminant') * (180/pi)),'平均误差为',num2str(mean(now_ours_errors))]);
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