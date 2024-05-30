% CEST processing using Lorentzian difference analysis (LDA)
% Jianpan Huang - jianpanhuang@outlook.com, 20230818
clear all, close all, clc
addpath(genpath(pwd));

%% Parameter setting
filepath = 'data';
filename = 'simu_phant_5pools_amidediff.mat';
offs_intst = [3.5, -3.5, 2.0]; % The interested offsets for amide, rNOE and guan, in ppm
fit_rgc = [-0.8, 0.8]; % Center range of Z-spectrum used for Lorentzian fitting, in ppm (0.5~1.5 ppm)
fit_rge = [-19, 19]; % Edge range of Z-spectrum used for Lorentzian fitting, in ppm
if_denoi = 0; % Set to 1 for CEST denoising using MLSVD, 0 for no denoising
rot_ang = 0; % Rotate the image when necessary (in angle degree)

%% Load data
load([filepath,filesep,filename]); % The mat file includes CEST images (img), frequency offsets (offs) and ROI (roi)
if length(size(img)) == 3
    img_temp = img;
    clear img;
    img(:,:,1,:) = img_temp;
end
% Seperate M0 and CEST images
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

%% Generate deltaB0 and correct B0
db0 = generate_load_db0_3d(filepath, 'b0map', img_norm, offs, roi, 'spline');
img_corr = correct_b0_3d(img_norm, db0, offs, roi, 'spline');
% figure, imagesc(db0),caxis([-1,1]),colorbar,colormap(jet);

%% Lorentzian difference analysis (LDA)
if size(offs,1) == 1
    offs = offs';
end
if ~exist([filepath, filesep, 'lda.mat'],'file') == 1
    %          Zi   A     G     dw   
    iv = [ 1    0.5   2     0 ];
    lb = [ 0.5  0.02  0.3  -0.1 ];
    ub = [ 1    1     100  +0.1 ]; % Here dw is set to a small range since B0 correction is done
    img_water_cur = zeros(xn, yn, sn, on);
    img_lda = zeros(xn, yn, sn, on);
    img_lda_inv = zeros(xn, yn, sn, on);
    fit_para = zeros(xn, yn, sn, length(iv));
    roi_vec = roi(:);
    roi_vec(roi_vec==0) = [];
    roi_len = length(roi_vec);
    % figure(1);
    h = waitbar(0, 'Doing Lorentzian fitting voxel by voxel, please wait >>>>>>'); 
    set(h, 'Units', 'normalized', 'Position', [0.4, 0.2, 0.3, 0.08])
    count = 1;
    tic
    for s = 1:sn
        for m = 1:xn
            for n = 1:yn
                if roi(m,n,s) == 1
                    count = count+1;
                    zspec = squeeze(img_corr(m,n,s,:));
                    fit_ind = get_fit_ind(zspec, offs, fit_rgc, fit_rge);
                    offs_trunc = offs(fit_ind);
                    zspec_trunc = zspec(fit_ind);
                    fit_para = lsqcurvefit(@lorentzian1pool,iv,offs_trunc,zspec_trunc,lb,ub);
                    zspec_fit = lorentzian1pool(fit_para,offs);
                    img_water_cur(m,n,s,:) = zspec_fit;
                    ld = zspec_fit-zspec; %lorentzian difference
                    img_lda(m,n,s,:) = ld;
                    img_lda_inv(m,n,s,:) = 1./(zspec+eps)-1./(zspec_fit+eps);
                    if mod(count,200)==0
                        figure(100),plot(offs,zspec,'-o',offs,zspec_fit,'-',offs,ld,'--','LineWidth',1.5)
                        axis([min(offs(:)),max(offs(:)),0,1.02]);
                        xlabel('Offset (ppm)'); ylabel('Z'); title('Z spectrum'); 
                        legend('Z', 'Zref', 'LD','Location', 'east');
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
    save([filepath, filesep, 'lda.mat'], 'img_water_cur','img_lda','img_lda_inv','offs','img','img_norm','img_m0','roi');
else
    load([filepath, filesep, 'lda.mat']);
end

%% Find CEST maps
[~,ind] = min(abs(offs-offs_intst));
for m = 1:length(offs_intst)
    cest(:,:,:,m) = img_lda(:,:,:,ind(m));
    cest_inv(:,:,:,m) = img_lda_inv(:,:,:,ind(m));
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
save_txt([filepath,filesep,'cestMeanValues_amide_rnoe_guan_lda.txt'], all_mean);

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
figure, set(gcf,'unit','normalized','position',[0.3,0.2,0.2,0.8]);
subplot(3,1,1), imagesc(cest_splice(:,:,1)),colorbar('fontsize',12,'fontweight','bold'),caxis([0,0.12]),axis off,colormap(jet);title([num2str(offs_intst(1)),' ppm']);
subplot(3,1,2), imagesc(cest_splice(:,:,2)),colorbar('fontsize',12,'fontweight','bold'),caxis([0,0.10]),axis off,colormap(jet);title([num2str(offs_intst(2)),' ppm']);
subplot(3,1,3), imagesc(cest_splice(:,:,3)),colorbar('fontsize',12,'fontweight','bold'),caxis([0,0.06]),axis off,colormap(jet);title([num2str(offs_intst(3)),' ppm']);
% CEST inv
figure, set(gcf,'unit','normalized','position',[0.5,0.2,0.2,0.8]);
subplot(3,1,1), imagesc(cest_inv_splice(:,:,1)),colorbar('fontsize',12,'fontweight','bold'),caxis([0,0.20]),axis off,colormap(jet);title([num2str(offs_intst(1)),' ppm inv']);
subplot(3,1,2), imagesc(cest_inv_splice(:,:,2)),colorbar('fontsize',12,'fontweight','bold'),caxis([0,0.18]),axis off,colormap(jet);title([num2str(offs_intst(2)),' ppm inv']);
subplot(3,1,3), imagesc(cest_inv_splice(:,:,3)),colorbar('fontsize',12,'fontweight','bold'),caxis([0,0.15]),axis off,colormap(jet);title([num2str(offs_intst(3)),' ppm inv']);