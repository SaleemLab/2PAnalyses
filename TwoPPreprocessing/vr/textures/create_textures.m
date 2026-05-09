%%---------------------------PARAMETERS---------------------------------
% Corridor's parameters
corridorL = 200; % corridor length in cm
corridorH = 12; % VR corridor height in cm

% Landmark parameters
texwidth = 0.04; % fraction of the corridor's length covered by a landmark
sf_base = 4; % no.of vertical bars visible
tex_contrast = 0.5;

% BG periodocity and chunk parameters
BG_contrast = 0.5; % geneated with 0.35 on 22/02/26 BG1_35C; was set to 0.25 in the new version. 
BGperiodicity = 0.12; %fraction of the corridor's visible length
BGchunkWidth = 0.02; %fraction of the corridor's visible length 0.02
BGdensity = 20; %# dots per chunk
BG_endGray = 0.5;

% Seed for random background
random_seed = 0;
s = RandStream('mt19937ar','Seed',random_seed);

%%------------------------------Landmarks--------------------------------
texsize = 512;
% Making horizontal and vertical axes asymmetric because the texture aspect raio is 1:1.5 (W=8cm, H=12cm)
sf_H = sf_base / texwidth * (corridorH / corridorL); % no.of horizonal bars visible
sf_V = sf_base; % no.of vertical bars visible

% Gray
textures(1).matrix = BG_endGray*ones(64,64);

% Unfiltered Whitenoise
textures(2).matrix = rand(s, 16, 512);

% Vertical grating
textures(6).matrix = 0.5*(1 + tex_contrast*repmat(sin(0:((2*sf_V*pi)/texsize):(2*sf_V)*pi-(((2*sf_V)*pi)/texsize)),texsize,1));

% Horizontal grating
textures(7).matrix = 0.5*(1 + tex_contrast*repmat(sin(0:(((2*sf_H)*pi)/texsize):(2*sf_H)*pi-(((2*sf_H)*pi)/texsize))',1,texsize));

% Plaid
textures(8).matrix = (textures(6).matrix+textures(7).matrix)/2;

% c = [0:0.1:1.0];

% % Horixontal grating at 0, 10 to 100% contrast
% for n = 6:16
%     textures(n).matrix = 0.5+c(n-5)*0.5*repmat(sin(0:(((2*sf)*pi)/texsize):(2*sf)*pi-(((2*sf)*pi)/texsize)),texsize,1);
% end

% % Anti-phase Vertical, horizontal gratings
% textures(17).matrix = 1-textures(4).matrix;
% textures(20).matrix = 1-textures(3).matrix;
% textures(21).matrix = 1-textures(17).matrix;


%%---------------------Background textures-------------------------------
% Calculates how many "columns" of dots fit in one 20% repeat (e.g., 0.20/0.02 = 10 chunks).
BGchunkNb = BGperiodicity / BGchunkWidth;
% Horizontal resolution of the final image 
finalBGlength = 2048; % 2048; Size of the BG textures in pixels. Should be power of 2
% Calculates vertical resolution based on the corridor's aspect ratio (12:200).
finalBGheight = 2^round(log2(finalBGlength * corridorH / (corridorL)));

%Size of BG pattern in pixels.
% That's just for construction. The BG texture will be interpolated 
% to finalBGlength anyway
BGchunklength = 80;
BGlength = BGchunklength * BGchunkNb;
% Sets the vertical height of the noise 
BGheight = round(BGlength * 1 / BGperiodicity * corridorH / corridorL);

% Smoothing Gaussian filter
sigma = BGchunklength / 4; % lowpass cutoff at half chunk size 
sigma1 = sigma;
filtSize = sigma * 10; %  x10 to avoid edge artifacts
BGheight_pad = BGheight + filtSize*2;

% Building the smoothing filter once: this is what smooths the dots into
% blobs 
x = 1:filtSize;
y = exp(-((x-round(filtSize/2)).^2)/(2*sigma^2));
y1 = exp(-((x-round(filtSize/2)).^2)/(2*sigma1^2));
y = y./sum(y);
y1 = y1./sum(y1);
filt2 = y'*y1;

% Generates random dots: bBecause randi(s, ...) and rand(s) are called inside the texID loop, Texture 2, 3, 4, and 5 all get different random dots.
for texID = 2:5
    % unfiltered BG pattern. BGdensity non-zero values every BGchunklength
    % columns. The rest stays at 0.5.
    Im_BG = 0.5*ones(BGheight_pad,BGlength);
    rand_idx = BGchunklength/2:BGchunklength:BGlength;
    for j = rand_idx
        for k = 1:BGdensity
            if mod(k, 2) == 0
                Im_BG(filtSize+randi(s, BGheight), j) = 0.5 + 0.5*rand(s); 
            else
                Im_BG(filtSize+randi(s, BGheight), j) = 0.5 - 0.5*rand(s);
            end
        end
    end
    
    % padded full BG texture with repeats
    pad_repeats = ceil(filtSize / BGlength);
    Im = repmat(Im_BG, [1, ceil(1 / BGperiodicity) + 2*pad_repeats]);
    
    % Convolving the BG texture with filter
    Imf = conv2(Im, filt2,'same');
    
    % Cropping padding repeats
    Imf = Imf(filtSize+1:end-filtSize,pad_repeats*BGlength + 1:end-pad_repeats*BGlength);
    
    % Cropping to actual number of repeats along the actual corridor's length
    Imf = Imf(:, 1:round(1 / BGperiodicity * BGlength));
    
    % Interpolating to get final texture of length finalTexLength pixels
    H = size(Imf, 1);
    L = size(Imf, 2);
    [u, v] = meshgrid(1:L, 1:H);
    [uq, vq] = meshgrid(linspace(1,L, finalBGlength), linspace(1,H, finalBGheight));
    Imf = interp2(u, v, Imf, uq, vq);
    
    %Normalizing to [-1, 1]
    Imf = (Imf - 0.5);
    Imf = Imf./max(abs(Imf), [], 'all');
    %Scaling by BG_contrast and scaling to [0, 1]
    Im_new = 0.5*(1 + BG_contrast*Imf);
    
    %Saving to textures
    textures(texID).matrix = Im_new;
end

% pause(1);
% close;

% clear Im ImF filt2 texID x y


%% plot textures

% cd C:\Users\m.morimoto\Documents\GitHub\SaleemLab-VR\VRCentral\data
% load('textures_hf_Mik4.mat')

figure; 

subplot(7,1,1)
tex=textures(2).matrix;    
imagesc(tex, [0 1]);
title({['background 1 ', num2str(size(tex,1)),'x', num2str(size(tex,2))],...
    ['contrast=', num2str(BG_contrast), ' filtsize=', num2str(filtSize),' sigma=', num2str(sigma), ' sigma1=', num2str(sigma1)]})
colormap gray; axis equal; box off; axis off

subplot(7,1,2)
tex=textures(3).matrix;    
imagesc(tex, [0 1]);
title({['background 2 ', num2str(size(tex,1)),'x', num2str(size(tex,2))],...
    ['contrast=', num2str(BG_contrast), ' filtsize=', num2str(filtSize),' sigma=', num2str(sigma), ' sigma1=', num2str(sigma1)]})
colormap gray; axis equal; box off; axis off

subplot(7,1,3)
tex=textures(4).matrix;    
imagesc(tex, [0 1]);
title({['background 3 ', num2str(size(tex,1)),'x', num2str(size(tex,2))],...
    ['contrast=', num2str(BG_contrast), ' filtsize=', num2str(filtSize),' sigma=', num2str(sigma), ' sigma1=', num2str(sigma1)]})
colormap gray; axis equal; box off; axis off

subplot(7,1,4)
tex=textures(5).matrix;    
imagesc(tex, [0 1]);
title({['background 4 ', num2str(size(tex,1)),'x', num2str(size(tex,2))],...
    ['contrast=', num2str(BG_contrast), ' filtsize=', num2str(filtSize),' sigma=', num2str(sigma), ' sigma1=', num2str(sigma1)]})
colormap gray; axis equal; box off; axis off

subplot(7,1,5); 
tex=textures(6).matrix;    
imagesc(tex, [0 1]);
title({['Horizontal grating', num2str(size(tex,1)),'x', num2str(size(tex,2))], ['sf=', num2str(sf_H)]})
colormap gray; axis equal; box off; axis off

subplot(7,1,6); 
tex=textures(7).matrix;    
imagesc(tex, [0 1]);
title({['Vertical grating', num2str(size(tex,1)),'x', num2str(size(tex,2))], ['sf=', num2str(sf_V)]})
colormap gray; axis equal; box off; axis off

subplot(7,1,7); 
tex=textures(8).matrix;    
imagesc(tex, [0 1]);
title({['plaid ', num2str(size(tex,1)),'x', num2str(size(tex,2))], ['sf V=', num2str(sf_V),' sf H=', num2str(sf_H)]})
colormap gray; axis equal; box off; axis off


figure; 

grating = textures(6).matrix;
plaid = textures(8).matrix;
final_width = round(finalBGlength * texwidth);
final_height = finalBGheight;
H = size(plaid, 1);
L = size(plaid, 2);
[u, v] = meshgrid(1:L, 1:H);
[uq, vq] = meshgrid(linspace(1,L, final_width), linspace(1,H, final_height));
grating = interp2(u, v, grating, uq, vq);
plaid = interp2(u, v, plaid, uq, vq);

i = 0;
for k = 2:5
    i = i + 1;
    subplot(4,1,i)
    tex=textures(k).matrix(:,1:finalBGlength);
    grating_start1 = round(finalBGlength * 0.20 - final_width/2);
    grating_start2 = round(finalBGlength * 0.60 - final_width/2);
    plaid_start1 = round(finalBGlength * 0.40 - final_width/2);
    plaid_start2 = round(finalBGlength * 0.80 - final_width/2);
    tex(:, grating_start1:grating_start1+final_width-1) = grating;
    tex(:, grating_start2:grating_start2+final_width-1) = grating;
    tex(:, plaid_start1:plaid_start1+final_width-1) = plaid;
    tex(:, plaid_start2:plaid_start2+final_width-1) = plaid;

    imagesc([0, 100], [0, 100 * corridorH / corridorL], tex(:,1:size(tex, 2)), [0 1]);
    title({['Visible background ', num2str(k-1), ' ', num2str(size(tex,1)),'x', num2str(size(tex,2))],...
        ['BG periodocity =', num2str(100*BGperiodicity), '%', ' contrast = ', num2str(BG_contrast)],...
        ['contrast=', num2str(BG_contrast), ' sigma=', num2str(100 * sigma / BGchunklength), '% of chunk width']})
    colormap gray; axis equal; box off; axis off
end

%% save
% 
% tex1 = textures(2).matrix;
% tex2 = textures(3).matrix;
% tex3 = textures(4).matrix;
% tex4 = textures(5).matrix;
% 
% grating_v = textures(6).matrix;
% plaid = textures(8).matrix;
% grey = textures(1).matrix;
% figure;
% imagesc(tex, [0 1]);
% axis equal; axis off; colormap(gray);
% 
% imwrite(tex1, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-20PercPeriodicity\BG1_50C20P.jpg');
% imwrite(tex2, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-20PercPeriodicity\BG2_50C20P.jpg');
% imwrite(tex3, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-20PercPeriodicity\BG3_50C20P.jpg');
% imwrite(tex4,  '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-20PercPeriodicity\BG4_50C20P.jpg');
% imwrite(grating_v, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-20PercPeriodicity\grating_vertical_newRun.jpg');
% imwrite(plaid, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-20PercPeriodicity\plaid_newRun.jpg');
% imwrite(grey, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-20PercPeriodicity\grey_newRun.jpg')
% close all