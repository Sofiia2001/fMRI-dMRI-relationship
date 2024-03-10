clear all
close all
clc

%% FUNCTIONAL MRI ANALYSIS

%% Load data

% Get the current directory as the base path
base_path = pwd;

% Add necessary directories to the MATLAB path
addpath(genpath(fullfile(base_path, 'Matlab_tools')))  % Add Matlab_tools directory and its subdirectories to the path
addpath(genpath(fullfile(base_path, 'Structural')))    % Add Structural directory and its subdirectories to the path
addpath(genpath(fullfile(base_path, 'Matlab_tools', 'Nifti_tools')))  % Add Nifti_tools directory and its subdirectories to the path

% Add Results directory and its subdirectories to the path
addpath(genpath(fullfile(base_path, 'Results')))

% Load Gray Matter (GM), White Matter (WM), and Cerebro-spinal Fluid (CSF) data
GM_data = load_untouch_nii('c1T1_2mm.nii');  % Load GM probability map
GM_data = double(GM_data.img);               % Convert GM image data to double for further processing
WM_data = load_untouch_nii('c2T1_2mm.nii');  % Load WM probability map
WM_data = double(WM_data.img);               % Convert WM image data to double for further processing
CSF_data = load_untouch_nii('c3T1_2mm.nii'); % Load CSF probability map
CSF_data = double(CSF_data.img);             % Convert CSF image data to double for further processing

% Load functional MRI (fMRI) data
addpath(genpath(fullfile(base_path, 'FMRI')))   % Add FMRI directory and its subdirectories to the path
fMRI = load_untouch_nii('FMRI_2_T1_2mm.nii.gz'); % Load fMRI data
fMRI_data = double(fMRI.img);                    % Convert fMRI image data to double for further processing

% Load T1-weighted structural MRI data
MRI = load_untouch_nii('T1_2mm.nii');   % Load T1-weighted structural MRI data
MRI_data = double(MRI.img);             % Convert MRI image data to double for further processing

% Load brain atlas data
addpath(genpath(fullfile(base_path, 'Atlas')))   % Add Atlas directory and its subdirectories to the path
atlas = load_untouch_nii('Hammers_2_T1_2mm_int.nii.gz'); % Load brain atlas data
atlas_data = double(atlas.img);                        % Convert atlas image data to double for further processing

% Clear unnecessary variables from the workspace to conserve memory
clear fMRI MRI atlas


%% Point 1: Data Preprocessing

%%%% a. Use SPM software to segment the T1w structural image into GM, WM and CSF tissues, obtaining the tissues probability maps. (This step may take a few minutes)

% The following instructions guide the user on how to perform tissue segmentation using SPM software:
% 1. Add the SPM Toolbox to the MATLAB search path, which can be done from the command line.
% 2. Open SPM software by typing 'spm' in the command line.
% 3. In the SPM interface, click on 'fMRI' or 'MRI', then select 'Segment'. 
% 4. Load the T1-weighted structural image file ('T1_2mm.nii.gz').
% 5. Use the default SPM settings for segmentation.
% 6. After segmentation, the resulting probability maps for GM, WM, and CSF tissues are obtained:
%    - c1T1_2mm.nii : Probability map for Gray Matter (GM)
%    - c2T1_2mm.nii : Probability map for White Matter (WM)
%    - c3T1_2mm.nii : Probability map for Cerebro-spinal Fluid (CSF)
% 7. Two additional images are discarded as they represent scalp tissues.
% 8. The correctness of the segmentation can be verified using the 'display' function in SPM.

% Now we can load the data, considering that file names might change if coregistration is applied.

% Define the repetition time (TR) for fMRI data
TR = 2.6;

% Get dimensions of the fMRI data
[nrows, ncols, nslices, nvolumes] = size(fMRI_data);

% Create a time vector for fMRI data
t = linspace(0, TR * nvolumes, nvolumes);

%%%% b. Threshold the three tissue probability maps and create for each tissue a binary mask.
%%%% Erode the masks (imerode matlab function) and compute the mean fMRI signal of WM and CSF.
%%%% Provide a justification of the employed thresholds and of the erosion parameters.

% Plot histograms of tissue probability maps
figure, histogram(GM_data(:), 80), title('Histogram Gray Matter')
ylim([0, 20000])
figure, histogram(WM_data(:), 80), title('Histogram White Matter')
ylim([0, 20000])
figure, histogram(CSF_data(:), 80), title('Histogram Cerebrospinal Fluid')
ylim([0, 20000])

% Provides good contrast between gray matter (dark gray) and white matter (lighter gray) tissues, 
% while CSF is void of signal (black). Water, such as CSF, as well as dense bone and air appear dark. 
%     Fat, such as lipids in the myelinated white matter, appears bright.
% Set thresholds for each tissue mask
th1 = 255 * 0.75;
th2 = 255 * 0.9;
th3 = 255 * 0.8;

% Justification for chosen thresholds:
% These thresholds were chosen after visually inspecting the tissue probability maps.
% The goal was to find thresholds that would accurately delineate the distributions of the three tissue types.

%Certainly! Let's provide reasoning for the chosen thresholds:

% The thresholds used for binary mask creation in the provided code are as follows:
% 
% - `th1 = 255 * 0.75`: Threshold for Gray Matter (GM) probability map.
% - `th2 = 255 * 0.9`: Threshold for White Matter (WM) probability map.
% - `th3 = 255 * 0.8`: Threshold for Cerebro-spinal Fluid (CSF) probability map.
% 
% Justification for these thresholds:
% 
% 1. **Gray Matter (GM)**:
%    - Gray matter typically exhibits higher signal intensity in T1-weighted MRI images compared to other tissues.
%    - By setting the threshold at 75% of the maximum intensity value (255 in a typical 8-bit image), 
% it's likely to capture most of the voxels representing gray matter while 
% minimizing inclusion of noise or non-gray matter regions with lower intensities.
% 
% 2. **White Matter (WM)**:
%    - White matter generally has lower signal intensity than gray matter in T1-weighted MRI images.
%    - Setting the threshold at 90% of the maximum intensity value ensures 
% that predominantly white matter voxels are included in the binary mask, 
% while potentially excluding regions of gray matter with higher intensity.
% 
% 3. **Cerebro-spinal Fluid (CSF)**:
%    - Cerebro-spinal fluid exhibits even lower signal intensity compared to gray and white matter.
%    - The threshold is set at 80% of the maximum intensity value, aiming 
% to capture voxels corresponding to CSF while minimizing inclusion of noise or misclassified tissue voxels.
% 
% These thresholds are chosen based on empirical observations and visual inspection 
% of the tissue probability maps. Adjustments to the thresholds might be necessary based on 
% the specific characteristics of the MRI data and the segmentation results.

% Erosion is a morphological operation commonly used in image processing to reduce the size 
% of objects in a binary image. It works by "eroding away" the boundaries of regions of 
% foreground pixels (typically represented as white) in a binary image. 
% 
% Here's how erosion works:
% 
% 1. **Kernel or structuring element**: Erosion requires a small matrix called a 
% kernel or structuring element. This kernel defines the neighborhood around each pixel 
% that is considered during the erosion operation. Common kernel shapes include squares, circles, or crosses.
% 
% 2. **Operation**: For each pixel in the input image, the corresponding pixel in the output 
% (eroded) image is set to white only if all the pixels in the neighborhood defined by the kernel 
% in the input image are white. Otherwise, the pixel is set to black (background). 
% This effectively "shrinks" the foreground regions in the image.
% 
% 3. **Boundary effects**: At the boundaries of objects, erosion can cause thinning or 
% complete removal of the object if the kernel extends beyond the object boundaries.
% 
% In MATLAB, erosion can be performed using the `imerode` function, which takes the 
% input binary image and the structuring element as input arguments.
% 
% Erosion is often used in conjunction with other morphological operations, such as dilation, 
% opening, and closing, to achieve desired effects such as noise reduction, object separation, or 
% feature extraction in binary images. In the context of the provided code, erosion is used after 
% thresholding the tissue probability maps to create binary masks, likely to refine the masks and ensure 
% they accurately represent the tissue regions while reducing any noise or artifacts.

% Create binary masks for each tissue using the defined thresholds
GM_mask = GM_data > th1;
WM_mask = WM_data > th2;
CSF_mask = CSF_data > th3;

% For each voxel in the WM probability map (WM_data), 
% if the intensity value exceeds 90% of the maximum intensity value (255), 
%     % set the corresponding pixel in the WM_mask to 1 (true); otherwise, 
%     % set it to 0 (false).
% WM_mask will contain 1s where the probability of being white matter is 
% relatively high, and 0s elsewhere.


% WM_mask_neg=ones(size(GM_data));
% WM_mask_neg=WM_mask_neg-GM_mask;
% 
% implay(GM_mask)
% implay(WM_mask)
% implay(WM_mask_neg)
% implay(CSF_mask)

% For each pixel in the binary mask, the corresponding pixel in the output (eroded) 
% mask is set to 1 (true) only if all the pixels in the neighborhood defined by the 
% structuring element in the input mask are also 1. Otherwise, the pixel is set to 0 (false).
% This operation effectively reduces the size of the foreground regions in the binary mask, 
% making them shrink towards the center of mass.
% Only immediate horizontal pixels are checked

se = strel('square',2);
GM_mask = imerode(GM_mask,se);
WM_mask = imerode(WM_mask,se);
CSF_mask = imerode(CSF_mask,se);

% We chose a structural element for the erosion such that it would remove
% most of the outliers without removing too many informative pixels. We found 
% that the square structural element of size 2-by-2 is the best compromise
% between the two.

for i = 1: nvolumes
    tmp = squeeze(fMRI_data(:,:,:,i)); 
    data2D_WM(i,:)  = tmp(WM_mask); 
    data2D_CSF(i,:) = tmp(CSF_mask); 
end

mean_WM_signal = nanmean(data2D_WM,2); % mean WM signal
mean_CSF_signal = nanmean(data2D_CSF,2); % mean CSF signal

%%%% c. Create the sumEPI image by summing the EPI volumes in the 4th dimension. 
%%%% Create a binary mask for the sumEPI image and provide a justification 
%%%% of the employed threshold.

sumEPI = sum(fMRI_data, 4);

figure, histogram(sumEPI(:),100), title('histogram sumEPI')
ylim([0,20000])
th4 = 100000;

% We chose this threshold because it's the point that divides the two gaussian
% distributions, since it's the lowest point between the two curves in the 
% histogram (expecially clear when we increase the number of bins).

sumEPI_mask = sumEPI > th4;
% implay(CSF_mask)

%%%% d. Mask the Hammers atlas with two masks: the GM binary mask and the sumEPI mask.
atlas_masked = atlas_data .* GM_mask;
atlas_masked = atlas_masked .* sumEPI_mask;

%%%% e. ROI time activity curve extraction:for each masked ROI of the atlas, 
%%%% extract the mean fMRI signal. 
%%%% Using the Hammers_labels.pdf file, discard for the following analyses 
%%%% the masked ROIs with less than 15 voxels and those belonging from: amygdala (3,4),
%%%% cerebellum (17,18), brainstem (19), corpus callosum (44), substantia nigra (74,75), ventricles (45,46,47,48,49).

%Discarding the ROI that belong to regions of no interest OR have less then
%15 voxels

disp('Removing regions of non interest from atlas')

% amyg (3 4); cer (17 18); bstem (19); cc (44); snigra (74 75); ventr (45 46 47 48 49)
% Initialize indices of ROIs to be excluded
idx_exl = [3 4 17 18 19 44 45 46 47 48 49 74 75];

% Loop through each possible index value present in the atlas data
for i = 1:max(atlas_data(:))
    % Check if the current index is present in the exclusion list
    if sum(idx_exl == i) ~= 0
        % If index is in the exclusion list, zero out corresponding values in the atlas_masked array
        % those not equal to i will be set to 1 (so we leave them)
        % those equal to i will be get rid of 
        atlas_masked = atlas_masked .* (atlas_masked ~= i);
    else
        % If index is not in the exclusion list
        % Check if the number of voxels belonging to the current ROI is less than 15
        if sum(atlas_masked == i, 'all') < 15
            % If fewer than 15 voxels, add index to exclusion list and zero out corresponding values
            idx_exl = [idx_exl i];
            atlas_masked = atlas_masked .* (atlas_masked ~= i);
        end
    end
end


% Extract the mean fMRI signal for each remaining ROI

% Find indices of remaining ROIs
remainedROI_idx = nonzeros(unique(atlas_masked));
nROI_diff = length(remainedROI_idx);

% Display the number of discarded ROIs
disp(['Number of discarded ROIs: ' num2str(length(unique(idx_exl)))]);

% Initialize matrix to store mean fMRI signal for each remaining ROI
mean_ROI_TAC_original = nan(nvolumes, length(remainedROI_idx));

% Loop through each remaining ROI
for i = 1:length(remainedROI_idx)
    % Get the index of the current ROI
    name_ROI = remainedROI_idx(i);
    disp(['Working on ROI #' num2str(name_ROI)])

    % Extract the binary mask for the current ROI
    mask_ROI = (atlas_masked == remainedROI_idx(i));

    % Initialize matrix to store fMRI signal for the current ROI
    data2D_GM = zeros(size(fMRI_data,4), sum(mask_ROI(:)));

    % Extract fMRI signal for the current ROI for each time point
    for tt = 1:size(fMRI_data,4)
        tmp = squeeze(fMRI_data(:,:,:,tt));
        data2D_GM(tt,:) = tmp(logical(mask_ROI));
    end

    % Calculate the mean fMRI signal for the current ROI
    mean_ROI_TAC_original(:,i) = nanmean(data2D_GM,2);

    % Plot the Time Activity Curve (TAC) for the current ROI (optional)
%     plot(t, mean_ROI_TAC_original(:,i))
%     xlabel('time')
%     ylabel('u.a')
%     title(['Time Activity Curve for ROI: ', num2str(name_ROI)])
%     pause(1.5)
end


clear th1 th2 th3 th4 GM_mask WM_mask se tmp sumEPI_mask sumEPI i CSF_data
clear tt mask_ROI name_ROI data2D_GM data2D_WM data2D_CSF atlas_data 
clear fMRI_data MRI_data nrows ncols nslices idx_exl GM_data WM_data


%% PART 2: DATA DENOISING

%% a. Noise regression
load('FMRI/MOCOparams.mat');

% Calculating temporal derivatives; considering the time domain as 2.6 sec
% difference between time stamps
[~, tempDerivatives] = gradient(newMOCOparams, TR);

% Linear Regression 
% size(X) = (225, 14); size(Y) = (225, 64)
X = [newMOCOparams tempDerivatives mean_WM_signal mean_CSF_signal];
Y = mean_ROI_TAC_original;

% Transforming the data into Z-scores
GC = zscore(X);

% beta - coefficients from Linear Regression; size(beta) = (14, 64)
beta_hat = (GC' * GC) \ (GC' * Y);

% Regressing out the noise
Y_regressed = Y - GC * beta_hat;

% Visualizing the regression matrix
figure
imagesc(GC);  % For a simple intensity plot
title('Regression matrix');
colormap gray;  % Optional: Set the colormap
colorbar;  % Optional: Display the colorbar

clear newMOCOparams mean_WM_signal mean_CSF_signal beta_hat

%% b. Temporal filtering

% Visualizing the frequency spectrum of one ROI signal

hippocampus_idx = find(remainedROI_idx==1);
thalamus_idx = find(remainedROI_idx==40);
% FFT
ROI_fft_original = Y_regressed(:, thalamus_idx) - mean(Y_regressed(:, thalamus_idx));
fft_signal = fft(ROI_fft_original);

% Frequency axis calculation
L = length(ROI_fft_original);  % Length of the signal
f = (0:(L/2))/L*(1/TR);  % Frequency axis

% Take the magnitude and normalize
fft_signal_mag = abs(fft_signal/L);

% Plot the frequency spectrum
figure
plot(f, fft_signal_mag(1:L/2+1));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Thalamus ROI');

%% Optionally, you can also plot the power spectrum (squared magnitude)
figure
fft_signal_pow = (abs(fft_signal)/L).^2;
plot(f, fft_signal_pow(1:L/2+1));
grid on
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Thalamus ROI');

%% FFT of the average ROI

% FFT on the average 
ROI_fft_original = mean(Y_regressed - mean(Y_regressed, 1), 2);
fft_signal = fft(ROI_fft_original);

% Frequency axis calculation
L = length(ROI_fft_original);  % Length of the signal
f = (0:(L/2))/L*(1/TR);  % Frequency axis

% Take the magnitude and normalize
fft_signal_mag = abs(fft_signal/L);

% Plot the frequency spectrum
figure
plot(f, fft_signal_mag(1:L/2+1));
grid on
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Average ROI');

% Optionally, you can also plot the power spectrum (squared magnitude)
figure
fft_signal_pow = (abs(fft_signal)/L).^2;
plot(f, fft_signal_pow(1:L/2+1));
grid on
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Averge ROI');

% We selected thalamus region because we use it in the task 4 but no
% visibile drift can be seen, so by only that is not possible to provide 
% any further comment about the denoising step. So we also checked for 
% hippocampus region, where is present an actual drift that, after the 
% process of denoising is noticible reduced.

%% Cut-off frequencies 

% From the class [0.008, 0.08] Hz is the range related to neural activity in
% BOLD signal of GM areas. We choose a band pass from 0.008 to 0.15 Hz in
% order to mantain also higher frequencies that could be also related to
% hemodynamic response

% Looking for the order of the filter(tuning parameters)
Fs = 1/TR;
Wp = [0.008 0.15]/(Fs/2);
Ws = [0.003 0.19]/(Fs/2);
Rp = 1;
Rs = 40;

%Cut-off frerquencies and filter order
[n, Wn] = buttord(Wp, Ws, Rp, Rs); % Wn =[0.064 Hz 0.1581] Hz

% Designing a Butterworth filter
[b, a] = butter(n,Wn);

figure
freqz(b,a,[],(1/TR))

filtered_signal = filtfilt(b, a, Y_regressed);

%% FFT on filtered signal (average of ROIs)

% FFT on the average 
ROI_fft_original = mean(filtered_signal, 2);
fft_signal = fft(ROI_fft_original);

% Frequency axis calculation
L = length(ROI_fft_original);  % Length of the signal
f = (0:(L/2))/L*(1/TR);  % Frequency axis

% Take the magnitude and normalize
fft_signal_mag = abs(fft_signal/L);

% Plot the frequency spectrum
figure
semilogy(f, fft_signal_mag(1:L/2+1));
grid on
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of filtered ROI');

% Optionally, you can also plot the power spectrum (squared magnitude)
figure
fft_signal_pow = (abs(fft_signal)/L).^2;
semilogy(f, fft_signal_pow(1:L/2+1));
grid on
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Filtered ROI');

%% Plot the signal in time domain(average of ROIs)
figure

plot(mean(Y_regressed - mean(Y_regressed, 1), 2));
hold on
plot(mean(filtered_signal, 2));
xlabel('Time [samples]');
ylabel('Amplitude');
title('Regressed signal vs. Filtered regressed signal');
legend('regressed','regressed + filtered')

clear f Fs Wn Wp Ws Rp Rs b a n L fft_signal thalamus_idx tempDerivatives
clear X Y ROI_fft_original 

clear  fft_signal_pow fft_signal_mag 
% clear Y_regressed

%% 3. VOLUME CENSORING

% Discard the volumes that are affected by motion artefacts

load('FMRI/FDparams.mat');
FD_values = FD(:,1); % FD values in mm are reported in the first column
affectedIndices = find(FD_values > 0.2); % framewise displacement 0.2 mm  

% For each identified volume to discard, remove also one before and two after volumes
artifacts_volumes = unique([affectedIndices-1, affectedIndices, affectedIndices+1, affectedIndices+2]);
final_volumes = setdiff(1:nvolumes,artifacts_volumes);

% Update the number of volumes of data matrix
nvolumes_reduced = length(final_volumes);

% I keep only volumes without motion artifacts 
mean_ROI_signal = filtered_signal(final_volumes,:); % or Y_regressed
mean_ROI_TAC_original = mean_ROI_TAC_original(final_volumes,:);
Y_regressed = Y_regressed(final_volumes,:);


clear FD_values FD affectedIndices final_volumes artifacts_volumes nvolumes filtered_signal

%% 4. Visual Check of the preprocessing step
% Plot the original time-course of the left thalamus region and what 
% was obtained after each denoising step. Do you see a drift in the original signal? If so, is the denoising able to remove it?

n_idx_thalamus = 40; %insert ROI index of the LEFT thalamus 
idx_thalamus=find(remainedROI_idx==n_idx_thalamus);

figure
hold on
plot([1:TR:nvolumes_reduced*TR]/60 , mean_ROI_TAC_original(:,idx_thalamus)-mean(mean_ROI_TAC_original(:,idx_thalamus)),'LineWidth',2)
plot([1:TR:nvolumes_reduced*TR]/60 , Y_regressed(:,idx_thalamus)-mean(Y_regressed(:,idx_thalamus)),'LineWidth',2)
plot([1:TR:nvolumes_reduced*TR]/60 , mean_ROI_signal(:,idx_thalamus)-mean(mean_ROI_signal(:,idx_thalamus)),'LineWidth',2)
plot([t(1),t(end)]./60,[0,0])

title('Left Thalamus Time Course after performing the denoising steps')
legend('Original epi signal','Epi signal after noise regression','Preprocessed epi signal')
xlabel('time [min]')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'fontsize',16,'fontweight','bold')
hold off

% In this region we did not observe any drift in the signal. Nevertheless
% checking other regions of the brain (for example the one of the
% hippocampus) it is possible to observe such drift, and that it is corrected by
% the temporal filtering.

clear n_idx_thalamus idx_thalamus nvolumes_reduced mean_ROI_TAC_original t TR Y_regressed

%% 5. Static FC Matrix Computation
[FC, p] = corr(mean_ROI_signal);
zFC     = atanh(FC);

figure
subplot(121)
imagesc(zFC)
axis square
colormap jet
caxis([-1 2])
colorbar
title('z-Fisher transformed FC ')
set(gca,'fontsize',16,'fontweight','bold')
subplot(122)
imagesc(p)
colormap jet
axis square
caxis([0 1])
colorbar
title('p-values matrix')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'fontsize',16,'fontweight','bold')
% saveas(4,fullfile(output_file_prefx, ['zFC_matrix.png']))

clear mean_ROI_signal FC

%% 6. Multiple Comparison Correction
alpha = 0.05; % Set the significance level

% Trying False Discovery Rate (FDR) correction
[h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p,alpha,'dep','yes');
% 'dep' indicates that the data are not all positive correlated or independent
N_fdr = sum(h_fdr(:)); % Count the number of significant results
kept_fdr = N_fdr/(size(zFC,1)*size(zFC,2)); % Calculate the proportion of significant results

% Sparsification based on FDR correction
zFC_corr = zFC; % Initialize the corrected matrix
zFC_corr(~h_fdr) = nan; % Set non-significant elements to NaN

% Binarization of the corrected matrix
zFC_binary = ~isnan(zFC_corr); % Binarize the matrix: 1 for significant, 0 for non-significant

% Justification of the chosen correction method
% We chose the False Discovery Rate (FDR) approach 
% because it provides better control over the reduction of 
% Type II error (false negatives) compared to Type I error (false positives). 
% This is particularly important in fMRI data analysis where it's more critical t
% o avoid missing true activations than to risk identifying false ones.
% Additionally, Bonferroni correction is overly conservative because it assumes 
% independence among all tests, which is not typically the case with fMRI data due to 
% spatial correlation. FDR correction is more flexible and accounts for dependencies 
% among tests, making it a preferable choice in this context.


clear alpha h_fdr N_fdr kept_fdr h_bonferroni N_bonferroni kept_bonferroni number_tests
clear adj_p adj_ci_cvrg crit_p p

%% 7. GRAPH MEASURES

% To summarize the functional connectivity in terms of node centrality, for each ROI compute
% the node degree, the eigenvector centrality (with eigenvector_centrality_und provided
% function) and the normalized betweenness centrality (with betweenness_wei.m provided
% function). In the metrics computation, consider only the statistically significative functional
% connections obtained after the multiple comparison correction at point 6. Plot the node
% degree, the eigenvector centrality and the normalized betweenness centrality of the ROIs
% using the stem matlab function. Which are for each metric the 10 ROIs with the higher
% metrics values? Provide the indices of these regions.

% node degree
DEG=sum(zFC_binary)';

% eigenvector centrality
EC = eigenvector_centrality_und(double(zFC_binary))';

% normalized betweenness centrality
% Node betweenness centrality is the fraction of all shortest paths in 
% the network that contain a given node. Nodes with high values of 
% betweenness centrality participate in a large number of shortest paths.

% since higher correlations are more naturally interpreted as shorter distances
% the input matrix should be some inverse of the connectivity matrix:
G = 1./zFC_corr;

% in conclusion nodes with high values of betweenness centrality
% participate in a large number of connection with high correlation values.
BTW_NORM = betweenness_wei(G);

% Betweenness centrality may be normalised to the range [0,1] as
% BC/[(N-1)(N-2)], where N is the number of nodes in the network.
nROI = length(remainedROI_idx);
BTW_NORM = BTW_NORM/((nROI - 1)*(nROI - 2));

figure
subplot(3,1,1)
stem(remainedROI_idx,DEG, 'LineWidth', 1.2)
title('node degree')
xlabel('ROIs')
subplot(3,1,2)
stem(remainedROI_idx,EC, 'LineWidth', 1.2)
title('eigenvectors centrality')
xlabel('ROIs')
subplot(3,1,3)
stem(remainedROI_idx,BTW_NORM, 'LineWidth', 1.2)
title('normalized betweenness centrality')
xlabel('ROIs')

[ ~, degree_idx] = sort(DEG, 'descend');
degree_idx =  degree_idx(1:10);
node_degree = remainedROI_idx(degree_idx); %  8  30  11  52  61  66  82  14  53  60

[ ~, eigenvector_idx] = sort(EC, 'descend');
eigenvector_idx =  eigenvector_idx(1:10);
eigenvector_centrality = remainedROI_idx(eigenvector_idx); %  8  30  11  82  66  61  67  13  52  64

[ ~, BC_idx] = sort(BTW_NORM, 'descend');
BC_idx =  BC_idx(1:10);
BC = remainedROI_idx  (BC_idx); %  32  65  33  35  41  1  38  40  39  27

% the numbers above are the indexes for each metric of the 10 ROIs with
% the higher values

clear BC node_degree eigenvector_centrality BC_idx eigenvector_idx degree_idx
clear zFC_binary G nROI


%% DIFFUSION MRI ANALYSIS
%% 1. Diffusion signal visualization & understanding
%% a.

% Load the diffusion volumes, the bvals file and the bvecs file
addpath(genpath(fullfile(base_path, 'DMRI')))

load('bvals')
load('bvecs')
dMRI_volumes=load_untouch_nii('diffusion_volumes.nii'); %[120×120×90×103 single]
dMRI_volumes=double(dMRI_volumes.img); %4D double

% How many different DWIs have been acquired?
nVols=size(dMRI_volumes,4); % this is the number of DWIs: the fourth dimension is 103
nSlices=size(dMRI_volumes,3);
nVox=size(dMRI_volumes,1);

% Excluding b=0, how many diffusion shells does this acquisition feature? (consider a
% small tolerance α=±20 s/mm2 in the shell definition).

tolerance = 20;
% Exclude b=0 value
nonZeroBvals = bvals(bvals > 0);
% Sort the non-zero b-values
sortedBvals = sort(nonZeroBvals);
% Determine the boundaries of each shell based on the tolerance
shell_values = [sortedBvals(1)];
for i = 2:length(sortedBvals)
    if abs(sortedBvals(i) - shell_values(end)) > tolerance
        shell_values = [shell_values, sortedBvals(i)];
    end
end

numShells = length(shell_values);
disp("Number of diffusion shells (excluding b=0): " + numShells);

clear i nonZeroBvals numShells sortedBvals tolerance base_path

%% b.
% Plot the diffusion signal of a voxel populated principally with cerebrospinal fluid
% (hint: you can skip steps 1a-1b and get back to them when you have computed the
% DTI metrics, which can aid you in the selection); is the diffusion signal ordered by
% its b-value? If not, sort it so that the signal points corresponding to the same shell
% are shown consequently (and shells are ordered in an ascending fashion).
% Provide both the plot of the unsorted signal and the sorted one.

voxel_CSF_selected=[83,70,38];  %FA(83,70,38)=0.0368, atlas_mask(83,70,38)=46, CSF_mask(83,70,38)=1

% Explanation of the process of selecting a voxel that is likely populated with cerebrospinal fluid (CSF) based on several aspects. Let's go through each aspect:
% DTI Metrics (FA): The fractional anisotropy (FA) is a metric derived from diffusion tensor imaging (DTI) that provides information about the degree of diffusion 
% anisotropy in tissues. In this case, the FA value at the coordinates (83,70,38) is reported to be low (FA(83,70,38)=0.0368). A low FA value is expected 
% in a voxel populated principally with CSF because CSF exhibits isotropic diffusion (FA<0.18 for isotropic grey matter).
% Anatomical Atlas: The anatomical atlas provides information about the anatomical regions corresponding to specific coordinates. In this case, the anatomical 
% atlas indicates that at the coordinates (83,70,38), the value is 46 (atlas_mask(83,70,38)=46), which corresponds to the lateral ventricle. The ventricles are 
% regions in the brain where the cerebrospinal fluid flows, indicating a higher likelihood of CSF presence in this voxel.
% By considering these different aspects, the voxel at the coordinates (83,70,38) is selected as a voxel populated principally with CSF.

signal_CSF_selected=squeeze(dMRI_volumes(voxel_CSF_selected(1),voxel_CSF_selected(2),voxel_CSF_selected(3),:));
figure
plot(signal_CSF_selected,'o-')
title(['plot of unsorted signal for voxel  ',num2str(voxel_CSF_selected)])
ylabel('diffusion value')

[~,idx_sort]=sort(bvals);
signal_CSF_selected_sorted=signal_CSF_selected(idx_sort);
figure
plot(signal_CSF_selected_sorted,'o-')
title(['plot of sorted signal for voxel  ',num2str(voxel_CSF_selected)])
ylabel('diffusion value')

clear voxel_CSF_selected signal_CSF_selected idx_sort CSF_mask signal_CSF_selected_sorted

%% c. By visually inspecting the sorted signal, provide a brief comment both on the inter
%%b-value and on the intra b-value variabilities. Why do these signal variations occur?

% By visually inspecting the sorted diffusion signal, we can observe both 
% inter-b-value and intra-b-value variabilities. Considering that CSF has 
% free and isotropic diffusion, intra and inter b-value variations will be
% more influenced by technical and acquisition factors rather than biological
% tissue characteristics. 
% The inter-b-value variability is due to the relationship between the diffusion
% MRI signal and the b-value. The signal attenuation is described by the equation
% S = S0 * exp(-b * D), where S is the acquired signal, S0 is the signal
% without diffusion weighting, b is the b-value, and D is the diffusion coefficient. 
% As the b-value increases, the diffusion weighting intensifies, resulting in
% a lower signal intensity. This is true for both isotropic and anisotropic diffusion.
% Therefore, the inter-b-value variability reflects the sensitivity of the diffusion signal
% to changes in the b-value. On the other hand, intra-b-value variability can 
% be influenced by various factors such as tissue properties, partial volume 
% effects, imaging artifacts, and noise. Since we are considering the diffusion 
% signal of a voxel populated principally with CSF (which exhibits isotropic diffusion),
% any signal fluctuations will primarily be attributed to noise and artifacts. 
% Indeed, noise and imaging artifacts in the diffusion MRI acquisition process 
% can introduce random fluctuations in the signal, causing intra-b-value variability.

%% 2. Diffusion tensor computation
%% a.
% From the entirety of the diffusion volumes data, create a new 4D matrix containing
% only the volumes corresponding to b=0 s/mm2 and to the shell closest to b=1000
% s/mm^2 identified at the point 1a.
diff_values=abs(shell_values-1000);
i=find(diff_values==min(diff_values));
shell_closest=shell_values(i);
idx_valid=zeros(1,nVols); % =0 if the volume does not correspond neither to b=0 nor to b=shell_closest
for i=1:nVols
    if bvals(i)==0 || (bvals(i)>shell_closest && bvals(i)<=shell_closest+20)
        idx_valid(i)=1; % =1 if the volume correspond to b=0 or to b=shell_closest
    end
end
dMRI_volumes_reduced=dMRI_volumes(:,:,:,idx_valid==1); % 4D matrix reduced
bvals_reduced=bvals(idx_valid==1); % there is the need of reducing also bvals and bvecs for next steps
bvecs_reduced=bvecs(:,idx_valid==1);
nVols_reduced=sum(idx_valid); 

clear diff_values shell_closest idx_valid i dMRI_volumes bvals bvecs shell_values nVols

%% b.
% Fit the voxel-wise diffusion tensor (using the linear least square approach seen in
% class) on the whole brain diffusion data created at point 2a. When performing the
% log(S/S0) transformation of the signal, use as S0 the voxel-wise mean value of all
% b=0 volumes of the available dataset. Use the eigenvalue/eigenvector
% decomposition to recover the FA / MD / RD indices.
S0=mean(dMRI_volumes_reduced(:,:,:,bvals_reduced==0),4);
Slog=zeros(nVox,nVox,nSlices,nVols_reduced);
for i=1:nVols_reduced
    Slog(:,:,:,i)=log((dMRI_volumes_reduced(:,:,:,i)./S0)+eps); % log of the normalized signal for the fit
end

% loading of the brain mask
brain_mask=load_untouch_nii('diffusion_brain_mask.nii'); % it is a binary mask
brain_mask=brain_mask.img > 0.5;
% histogram(brain_mask.img,10) <-- to choose the threshold

brain_mask=brain_mask & S0>0; % removing voxels for which S0=0 due to an error in the acquisition step
% There are some voxel for which the signal is constant and equal to 0, probably due to an error
% in the acquisition step, so we decided to remove them.

brain_mask(:,[37,38],:)=false; 
% removing voxels for which there was an error in the acquisition step
% By visual inspection of the results at a first stage we noticed that theese voxels (that are
% peripheral voxels in the brain mask) assume unusually high values for MD index, probably due to an error
% in the acquisition step, so we decided to remove them to avoid that they influence the results.

% build the B design matrix for the linear least squares approach
B=zeros(nVols_reduced,6); % Initialize the B matrix with zeros, where nVols_reduced is the number of reduced volumes
B(:,1:3)=bvecs_reduced'.^2; % Assign the squares of the x, y, and z components of the reduced b-vectors to the first three columns of B
B(:,4)=bvecs_reduced(1,:).*bvecs_reduced(2,:); % Assign the product of the x and y components of the reduced b-vectors to the fourth column of B
B(:,5)=bvecs_reduced(1,:).*bvecs_reduced(3,:); % Assign the product of the x and z components of the reduced b-vectors to the fifth column of B
B(:,6)=bvecs_reduced(2,:).*bvecs_reduced(3,:); % Assign the product of the y and z components of the reduced b-vectors to the sixth column of B
B=B.*bvals_reduced'; % Multiply each column of B by the corresponding reduced b-values


% initialize the structures which will be used to contain DTI indexes
FA=zeros(nVox,nVox,nSlices);
MD=zeros(nVox,nVox,nSlices);
RD=zeros(nVox,nVox,nSlices);

FirstX=zeros(nVox, nVox, nSlices);
FirstY=zeros(nVox, nVox, nSlices);
FirstZ=zeros(nVox, nVox, nSlices);

% start the cycle to fit the voxel-wise diffusion tensor
for k=1:nSlices
    % print fitting progress
    disp([' Fitting Slice ',num2str(k)])
    for i=1:1:nVox
        for j=1:1:nVox
            % check if current voxel belongs to the mask
            if (brain_mask(i,j,k))
                
                % extract the signal from each voxel
                VoxSignal=squeeze(Slog(i,j,k,:));
                % fit the DTI using LLS approach
                D=-(B'*B)\B'*VoxSignal;
                % reconstruct the diffusion tensor from the fitted parameters
                T=[D(1) D(4)/2 D(5)/2;
                   D(4)/2 D(2) D(6)/2;
                   D(5)/2 D(6)/2 D(3)];
               
                % compute eigenvalues and eigenvectors
                [eigenvects,eigenvals]=eig(T); 
                eigenvals=diag(eigenvals);
                % manage negative eigenvals as shown in laboratory:
                if eigenvals(1)<0 && eigenvals(2)<0 && eigenvals(3)<0
                    eigenvals=abs(eigenvals); % if all <0 -> take the absolute value
                end
                % otherwise -> put negatives to zero
                if eigenvals(1)<0, eigenvals(1)=0; end
                if eigenvals(2)<0, eigenvals(2)=0; end
                if eigenvals(3)<0, eigenvals(3)=0; end
                
                % compute Fractional Anisotropy index
                FA(i,j,k)=(1/sqrt(2))*sqrt((eigenvals(1)-eigenvals(2)).^2+(eigenvals(2)-eigenvals(3)).^2 + (eigenvals(1)-eigenvals(3)).^2) ...
                    ./sqrt(eigenvals(1).^2+eigenvals(2).^2+eigenvals(3).^2);
                % compute Mean Diffusivity index
                MD(i,j,k)=mean(eigenvals);

                % Diffusion direction color encoded maps
                % sort eigenvalues, eigenvectors
                [sorted, idx_sort]=sort(eigenvals);
                eigenvals=eigenvals(idx_sort);
                eigenvects=eigenvects(:,idx_sort);

                % RD index
                RD(i,j,k)=mean(eigenvals([1,2]));

                % I take the principal eigenvector and I decompose it
                First=eigenvects(:,3); % it's the eigenvector corresponding to the higher eigenvalue
                FirstX(i,j,k)=abs(First(1))*FA(i,j,k);
                FirstY(i,j,k)=abs(First(2))*FA(i,j,k);
                FirstZ(i,j,k)=abs(First(3))*FA(i,j,k);
            end
        end
    end
end

% Explanation of row: brain_mask(:,[37,38],:)=false; 
% for ii=1:size(MD,2)
%     imagesc(squeeze(MD(:,ii,:)))
%     disp(['iteration n. ',num2str(ii)])
%     pause
%     close
% end 

clear B i j k VoxSignal D T eigenvals eigenvects idx_sort principal_eigenvector brain_mask
clear bvals_reduced bvecs_reduced dMRI_volumes_reduced First ii nVols_reduced sorted Slog S0

%% c.
% Provide the visualization of the FA and MD maps for a central slice.
% FA for the central slice
figure
imagesc(imrotate(squeeze(FA(:,:,45)),90)) 
colorbar
title('FA index for slice 45')
% MD for the central slice
figure
imagesc(imrotate(squeeze(MD(:,:,45)),90))
colorbar
title('MD index for slice 45')
% Its colour encoded MD map weighted for the FA
% Creating the colour encoded directional map
color(:,:,:,1)=FirstX; 
color(:,:,:,2)=FirstY;
color(:,:,:,3)=FirstZ;
figure
for ii=1:1:90 % Number of slices
    slice=reshape(color(:,:,ii,:),nVox, nVox, 3);
    image(imrotate(slice,90))
    pause(0.1)
end

% Visualizing the central slice (3rd dimension)
central_slice=reshape(color(:,:,45,:),nVox, nVox, 3);
figure
image(imrotate(central_slice, 90))

clear nSlice_FA nSlice_MD color central_slice central_slice2 FirstX FirstY FirstZ
clear ii slice

%% d.
% Mask the FA and MD maps (as done for the Hammers atlas in fMRI: point 1.d),
% extract their mean values in each ROI (fMRI: point 1.e.). 

% initialize the structures which will be used to contain ROI-mean DTI indexes
mean_ROI_FA=zeros(nROI_diff,1); 
mean_ROI_MD=zeros(nROI_diff,1);
for i=1:nROI_diff
    name_ROI=remainedROI_idx(i);
    disp(['Working on ROI #',num2str(name_ROI)])
    ROI_mask=atlas_masked==name_ROI;
    FA_values=[];
    for v=1:nVox
        for vv=1:nVox
            for vvv=1:nSlices
                if ROI_mask(v,vv,vvv)==1
                    FA_values=[FA_values, FA(v,vv,vvv)];
                end
            end
        end
    end
    mean_ROI_FA(i)=mean(FA_values); 
    MD_ROI=MD(ROI_mask==1); 
    mean_ROI_MD(i)=mean(MD_ROI); 
end

% visualization of the results
figure
subplot(211), stem(remainedROI_idx,mean_ROI_FA, 'LineWidth',1.2), xlabel('ROIs'), title('mean FA')
subplot(212), stem(remainedROI_idx,mean_ROI_MD, 'LineWidth',1.2), xlabel('ROIs'), title('mean MD')

clear atlas_masked i MD_ROI FA_values name_ROI v vv vvv ROI_mask remainedROI_idx nROI_diff
clear nSlices nVox



%% dMRI/fMRI INTEGRATION

%% 2. Quantitative results

% Compute and provide the Spearman correlation between the six pairs of variables of point
% 1. Discuwss the results: is there a statistically significant relationship between any pair of
% variables?

ROI_total = [DEG, EC', BTW_NORM];

v1 = zeros(1,6);
p1 = zeros(1,6);
xl = {'DEG vs FA', 'EC vs FA', 'BTW NORM vs FA','DEG vs MD', 'EC vs MD', 'BTW NORM vs MD'};

for i = 1:3
    [h,p] = corr(ROI_total(:,i),mean_ROI_FA, 'type','Spearman');
    v1(1,i) = h;
    p1(1,i) = p;
end

for i = 1:3
    [h,p] = corr(ROI_total(:,i),mean_ROI_MD, 'type', 'Spearman');
    v1(1,i+3) = h;
    p1(1,i+3) = p;
end

figure
bar(v1');
ylim([-1 1])
ylabel('correlation values')
title('correlations between fMRI and dMRI metrics')
set(gca,'xticklabel',xl)

save Results\results.mat

%% Visual Inspection

ROI_total = [DEG, EC', BTW_NORM];
xl = {'Node degree', 'Eigenvector centrality', 'Betweenness centrality'};

figure('Name', 'Visual Inspection', 'NumberTitle', 'off');
sgtitle('fMRI metrics vs. dMRI metrics');

for i = 1:3
    subplot(2, 3, i);
    scatter(ROI_total(:,i), mean_ROI_FA);
    hold on;
    xlabel(xl(i));
    ylabel('ROIs FA');
    lsline;
    grid on;
    h = lsline;
    h.Color = 'red';  % Set regression line color to red
    h.LineWidth = 2;  % Set regression line thickness

    % Use pre-computed p-value and Spearman correlation value
    pval = p1(i);
    spearman_corr = v1(i);
    text(0.5, 0.9, ['SP = ' num2str(spearman_corr, '%.2f')], 'Units', 'normalized');
    text(0.5, 0.8, ['p-value = ' num2str(pval, '%.4f')], 'Units', 'normalized');
    
    hold off;
end

for i = 1:3
    subplot(2, 3, i + 3);
    scatter(ROI_total(:,i), mean_ROI_MD);
    hold on;
    xlabel(xl(i));
    ylabel('ROIs MD');
    lsline;
    grid on;
    h = lsline;
    h.Color = 'red';  % Set regression line color to red
    h.LineWidth = 2;  % Set regression line thickness

    % Use pre-computed p-value and Spearman correlation value
    pval = p1(i+3);
    spearman_corr = v1(i+3);
    text(0.5, 0.9, ['SP = ' num2str(spearman_corr, '%.2f')], 'Units', 'normalized');
    text(0.5, 0.8, ['p-value = ' num2str(pval, '%.4f')], 'Units', 'normalized');
    
    hold off;
end

clear xl;

%% Results
% Fractional Anisotropy (FA) is negatively correlated with the Node Degree (DEG) at a
% statistically significant level (p-value=0.0028). 

% Fractional Anisotropy is negatively correlated with the Eigenvector Centrality (EC)
% at a statistically significant level (p-value=0.0029). 

% Fractional Anisotropy is not correlated with the Normalized Betweenness Centrality 
% (BTW_NORM) at a statistically significant level (p-value=0.086).

% Mean Diffusivity (MD) is correlated with the Node Degree at a statistically significant
% level (p-value=0.0170). 

% Mean Diffusivity (MD) is correlated with the Eigenvector Centrality at a statistically 
% significant level (p-value=0.0231).

% Mean Diffusivity is not correlated with the Normalized Betweenness Centrality
% at a statistically significant level (p-value=0.3542).