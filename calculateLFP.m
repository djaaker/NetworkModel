function [LFP] = calculateLFP(SNN)
    %% Calculating LFP from the spikes,ge and gi alpha = 1.65;
    dt = SNN.dt;
    tau = 0.006;
    tau_shift = tau/dt;
    
    E_e = SNN.neurons(1).E_e;
    E_i = SNN.neurons(1).E_i;
    
    LFP.lfp = abs(circshift(SNN.ge_neurons.*(E_e-SNN.v_neurons),[0 tau_shift]) - 1.65*(SNN.gi_neurons.*(E_i-SNN.v_neurons)));
    LFP.lfp = LFP.lfp(1:SNN.N_E,:); % only excitatory neurons
    
    LFP.lfpz = ( LFP.lfp - repmat(nanmean(LFP.lfp,2),[1 size(LFP.lfp,2)]) )./repmat(nanstd(LFP.lfp,[],2),[1 size(LFP.lfp,2)]);
    
    LFP.lfp = reshape(LFP.lfp,SNN.grid_size_e,SNN.grid_size_e,[]);
    LFP.lfpz = reshape(LFP.lfpz,SNN.grid_size_e,SNN.grid_size_e,[]);
    %% Geting LFP by averaging over a patch of 10x10 neurons 
    LFP.lfp_patch_size = 10;
    T = size(LFP.lfpz,3);
    
    N_blocks = floor(SNN.grid_size_e / LFP.lfp_patch_size); % Number of full blocks
    row_remainder = mod(SNN.grid_size_e, LFP.lfp_patch_size); % Remaining rows
    col_remainder = mod(SNN.grid_size_e, LFP.lfp_patch_size); % Remaining columns
    
    % Process full lfp_patch_size x lfp_patch_size blocks
    lfpz_main = LFP.lfpz(1:N_blocks*LFP.lfp_patch_size, 1:N_blocks*LFP.lfp_patch_size, :);
    lfpz_reshaped_main = reshape(lfpz_main, LFP.lfp_patch_size, N_blocks, LFP.lfp_patch_size, N_blocks, T);
    LFP.LFP_main = squeeze(mean(mean(lfpz_reshaped_main, 1), 3)); % Full blocks
    
    % Handle remaining rows and columns if necessary
    if row_remainder > 0
        lfpz_row = LFP.lfpz(N_blocks*LFP.lfp_patch_size+1:end,1:N_blocks*LFP.lfp_patch_size, :);
        lfpz_reshaped_row = reshape(lfpz_row, row_remainder, 1, LFP.lfp_patch_size, N_blocks, T);
        LFP.LFP_row = squeeze(mean(mean(lfpz_reshaped_row, 1), 3)); % Remaining rows
    end
    if col_remainder > 0
        lfpz_col = LFP.lfpz(1:N_blocks*LFP.lfp_patch_size,N_blocks*LFP.lfp_patch_size+1:end, :);
        lfpz_reshaped_col = reshape(lfpz_col, LFP.lfp_patch_size, N_blocks,col_remainder,1,T);
        LFP.LFP_col = squeeze(mean(mean(lfpz_reshaped_col, 1), 3)); % Remaining cols
    end
    if row_remainder > 0 && col_remainder > 0
        LFP.LFP_corner = mean(LFP.lfpz(N_blocks*LFP.lfp_patch_size+1:end, N_blocks*LFP.lfp_patch_size+1:end, :), [1,2]); % Bottom-right corner
    end
    
    % Combining all
    N_blocks_combined = ceil(SNN.grid_size_e / LFP.lfp_patch_size);
    LFP.N_blocks_combined = N_blocks_combined;
    LFP.N_blocks = N_blocks;
    LFP.LFP_combined = zeros(N_blocks_combined,N_blocks_combined,T);
    
    LFP.LFP_combined(1:N_blocks,1:N_blocks,:) = LFP.LFP_main;
    LFP.LFP_combined(end,1:N_blocks,:) = LFP.LFP_row;
    LFP.LFP_combined(1:N_blocks,end,:) = LFP.LFP_col;
    LFP.LFP_combined(end,end,:) = LFP.LFP_corner;

    %% Downsampling the LFP to 1000 Hz 
    sampling_freq = 1/dt;
    target_sampling_freq = 1000;
    downsampling_ratio = sampling_freq/target_sampling_freq;

    LFP.sampling_freq = sampling_freq;
    LFP.target_sampling_freq = target_sampling_freq;
    % I am using downsample. Dont know if I should be worried about
    % applying an anti aliasing filter. May need to use decimate 

    LFP.LFP_downsampled = LFP.LFP_combined(:,:,1:downsampling_ratio:end);

    LFP.LFP_decimated = zeros(N_blocks_combined,N_blocks_combined,T/downsampling_ratio);

    for i = 1:N_blocks_combined % rows
        for j = 1:N_blocks_combined % cols
            LFP.LFP_decimated(i,j,:) = decimate(squeeze(LFP.LFP_combined(i,j,:)),downsampling_ratio);
        end
    end

    %% Applying a Gaussian kernel to the LFP data
    % Create the Gaussian filter
    gauss_sigma = 2;
    fprintf('The sigma for the gaussian filter for LFP is : %d% \n', gauss_sigma*(SNN.grid_length/SNN.grid_size_e)*LFP.lfp_patch_size);
    % Apply convolution with the kernel (preserving the original size)
    LFP.LFP_gaussian_downsampled = zeros(size(LFP.LFP_downsampled));
    LFP.LFP_gaussian = zeros(size(LFP.LFP_combined));
    
    for i=1:size(LFP.LFP_gaussian_downsampled,3)
        LFP.LFP_gaussian_downsampled(:,:,i) = imgaussfilt(squeeze(LFP.LFP_downsampled(:,:,i)),gauss_sigma);
    end

    for i=1:size(LFP.LFP_gaussian,3)
        LFP.LFP_gaussian(:,:,i) = imgaussfilt(squeeze(LFP.LFP_combined(:,:,i)),gauss_sigma);
    end

end