% CEST processing using magnetization transfer ratio asymmetry (MTRasm)
% Jianpan Huang - jphuang@hku.hk, 20230818
clear all, close all, clc
addpath(genpath(pwd));

%% Parameter setting
filepath = 'data';
filename = 'simu_phant_5pools_amidediff.mat';
offs_intst = [3.5, 2.0]; % The interested offsets for amide and guan, in ppm
if_denoi = 1; % Set to 1 for CEST denoising using MLSVD, 0 for no denoising
rot_ang = 0; % Rotate the image when necessary (in angle degree)

%% Load data
load([filepath,filesep,filename]); % The mat file includes CEST images (img), frequency offsets (offs) and ROI (roi)
if ~exist("offs","var")
    error("Please define the variable 'offs' (list of saturatiion frequency offsets) before this if sentence")
end
if length(size(img)) == 3
    img_temp = img;
    clear img;
    img(:,:,1,:) = img_temp;
end
if ~exist("roi","var")
    roi = draw_load_roi(filepath, img(:,:,:,1), 'roi', 'polygon');
end

%% Seperate M0 and CEST images
img_m0 = img(:,:,:,1);
img(:,:,:,1) = [];
offs(1) = [];

%% Denoising
img = imrotate(img, rot_ang); % Rotate images when necessary
% img_all = fliplr(imrotate(img_all, rot_ang)); % Flip and rotate images when necessary
[xn, yn, sn, on] = size(img);
for s = 1:sn
    if if_denoi == 1
	    img(:,:,s,:) = mlsvd_denoi(squeeze(img(:,:,s,:)), round([xn/3, yn/3, on/4]));
    end
end

%% Z-spectrum normalization
img_norm = img./repmat(img_m0, [1,1,1,on]);

%% Generate deltaB0 and Correct B0
db0 = generate_load_db0_3d(filepath, 'b0map', img_norm, offs, roi, 'spline');
img_corr = correct_b0_3d(img_norm, db0, offs, roi, 'spline');
% figure, imagesc(db0),caxis([-1,1]),colorbar,colormap(jet);

%% Magnetization transfer ratio asymmetry (MTRasm)
if size(offs,1) == 1
    offs = offs';
end
offs_pos = offs(offs>=0);
if ~exist([filepath, filesep, 'mtrasym.mat'],'file') == 1
    %          Zi   A     G     dw   
    offs_max = max(abs(offs));
    offs_pos_interp = 0:0.01:offs_max;
    offs_neg_interp = 0:-0.01:-offs_max;
    img_mtrasym = zeros(xn, yn, sn, length(offs_pos));
    img_mtrasym_inv = zeros(xn, yn, sn, length(offs_pos));
    roi_vec = roi(:);
    roi_vec(roi_vec==0) = [];
    roi_len = length(roi_vec);
    % figure(1);
    h = waitbar(0, 'Doing MTRasym voxel by voxel, please wait >>>>>>'); 
    set(h, 'Units', 'normalized', 'Position', [0.4, 0.2, 0.3, 0.08])
    count = 1;
    tic
    for s = 1:sn
        for m = 1:xn
            for n = 1:yn
                if roi(m,n,s) == 1
                    count = count+1;
                    zspec = squeeze(img_corr(m,n,s,:));
                    zspec_pos_interp  = interp1(offs,zspec,offs_pos_interp,'spline');
                    zspec_neg_interp  = interp1(offs,zspec,offs_neg_interp,'spline');
                    % MTRasym
                    mtrasym_interp = zspec_neg_interp-zspec_pos_interp; 
                    mtrasym = interp1(offs_pos_interp,mtrasym_interp,offs_pos,'spline');
                    img_mtrasym(m,n,s,:) = mtrasym;
                    % MTRasym inverse
                    mtrasym_inv_temp = 1./(zspec_pos_interp+eps)-1./(zspec_neg_interp+eps);
                    mtrasym_inv = interp1(offs_pos_interp,mtrasym_inv_temp,offs_pos,'spline');
                    img_mtrasym_inv(m,n,s,:) = mtrasym_inv;

                    if mod(count,200)==0
                        figure(100),plot(offs,zspec,'-o',offs_pos,mtrasym,'-','LineWidth',1.5)
                        axis([min(offs(:)),max(offs(:)),-0.2,1.02]);
                        xlabel('Offset (ppm)'); ylabel('Z'); title('Z spectrum'); 
                        legend('Z', 'MTRasym','Location', 'east');
                        set(gca, 'Xdir', 'reverse', 'FontWeight','bold','FontSize',14,'LineWidth',3)
                        set(gcf,'color','w');
                    end
                    waitbar(count/roi_len,h)
                end
            end
        end   
    end
    toc
    close gcf;
    delete(h);
    save([filepath, filesep, 'mtrasym.mat'],'img_mtrasym','img_mtrasym_inv','offs','img','img_norm','img_m0','roi');
else
    load([filepath, filesep, 'mtrasym.mat']);
end

%% Find cest maps
[~,ind] = min(abs(offs_pos-offs_intst));
for m = 1:length(offs_intst)
    cest(:,:,:,m) = img_mtrasym(:,:,:,ind(m));
    cest_inv(:,:,:,m) = img_mtrasym_inv(:,:,:,ind(m));
end

%% Save mean values
roi_vec = roi(:);
for m = 1:length(offs_intst)
    cest_temp = cest(:,:,:,m);
    cest_vec = cest_temp(:);
    cest_inv_temp = cest_inv(:,:,:,m);
    cest_inv_vec = cest_inv_temp(:);
    cest_mean(1,m) = mean(cest_vec(roi_vec==1));
    cest_inv_mean(1,m) = mean(cest_inv_vec(roi_vec==1));
end
all_mean = [cest_mean;cest_inv_mean];
save_txt([filepath,filesep,'cestMeanValues_amide_guan_mtrasym.txt'], all_mean);

%% Display
cest_splice = zeros(xn,yn*sn,length(offs_intst));
cest_inv_splice = zeros(xn,yn*sn,length(offs_intst));
for m = 1:length(offs_intst)
    for s = 1:sn
        cest_splice(:,(s-1)*yn+1:s*yn,m) = cest(:,:,s,m);
        cest_inv_splice(:,(s-1)*yn+1:s*yn,m) = cest_inv(:,:,s,m);
    end
end
set(0,'defaultfigurecolor','w')
% CEST
figure, set(gcf,'unit','normalized','position',[0.3,0.2,0.2,0.6]);
roi_all = repmat(roi,[1,1,2]);
cest_splice(roi_all==0)=-1e8;
subplot(2,1,1), imagesc(cest_splice(:,:,1)),colorbar('fontsize',12,'fontweight','bold'),caxis([-0.1,0.1]),axis off,colormap(jet);title([num2str(offs_intst(1)),' ppm']);
subplot(2,1,2), imagesc(cest_splice(:,:,2)),colorbar('fontsize',12,'fontweight','bold'),caxis([-0.1,0.1]),axis off,colormap(jet);title([num2str(offs_intst(2)),' ppm']);
% CEST inv
cest_inv_splice(roi_all==0)=-1e8;
figure, set(gcf,'unit','normalized','position',[0.5,0.2,0.2,0.6]);
subplot(2,1,1), imagesc(cest_inv_splice(:,:,1)),colorbar('fontsize',12,'fontweight','bold'),caxis([-0.15,0.15]),axis off,colormap(jet);title([num2str(offs_intst(1)),' ppm (inv)']);
subplot(2,1,2), imagesc(cest_inv_splice(:,:,2)),colorbar('fontsize',12,'fontweight','bold'),caxis([-0.15,0.15]),axis off,colormap(jet);title([num2str(offs_intst(2)),' ppm (inv)']);
