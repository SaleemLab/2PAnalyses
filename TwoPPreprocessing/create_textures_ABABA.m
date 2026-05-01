%%---------------------------PARAMETERS---------------------------------
corridorL = 200; 
corridorH = 12; 
texwidth = 0.04; 
sf_base = 4; 
tex_contrast = 0.5;

BG_contrast = 0.5; 
BGperiodicity = 0.20; 
BGchunkWidth = 0.02; 
BGdensity = 20; 
BG_endGray = 0.5;

random_seed = 0;
s = RandStream('mt19937ar','Seed',random_seed);

%%------------------------------Landmarks--------------------------------
texsize = 512;
sf_H = sf_base / texwidth * (corridorH / corridorL); 
sf_V = sf_base; 

textures(1).matrix = BG_endGray*ones(64,64);
textures(2).matrix = rand(s, 16, 512);

% Vertical grating
textures(6).matrix = 0.5*(1 + tex_contrast*repmat(sin(0:((2*sf_V*pi)/texsize):(2*sf_V)*pi-(((2*sf_V)*pi)/texsize)),texsize,1));
% Horizontal grating
textures(7).matrix = 0.5*(1 + tex_contrast*repmat(sin(0:(((2*sf_H)*pi)/texsize):(2*sf_H)*pi-(((2*sf_H)*pi)/texsize))',1,texsize));
% Plaid
textures(8).matrix = (textures(6).matrix+textures(7).matrix)/2;

%%---------------------Background textures-------------------------------
BGchunkNb = BGperiodicity / BGchunkWidth;
finalBGlength = 2048; 
finalBGheight = 2^round(log2(finalBGlength * corridorH / (corridorL)));

% Using your original 80-pixel chunk length to match smoothing grain
BGchunklength = 80;
BGlength_master = BGchunklength * BGchunkNb; 
BGheight = round(BGlength_master * 1 / BGperiodicity * corridorH / corridorL);

sigma = BGchunklength / 4; 
sigma1 = sigma;
filtSize = sigma * 10; 
BGheight_pad = BGheight + filtSize*2;

x = 1:filtSize;
y = exp(-((x-round(filtSize/2)).^2)/(2*sigma^2));
y1 = exp(-((x-round(filtSize/2)).^2)/(2*sigma1^2));
y = y./sum(y);
y1 = y1./sum(y1);
filt2 = y'*y1;

for texID = 2:5
    % create noise A leading to the grating
    Im_A = 0.5*ones(BGheight_pad, BGlength_master);
    for j = BGchunklength/2:BGchunklength:BGlength_master
        for k = 1:BGdensity
            row = filtSize+randi(s, BGheight);
            if mod(k, 2) == 0
                Im_A(row, j) = 0.5 + 0.5*rand(s); 
            else
                Im_A(row, j) = 0.5 - 0.5*rand(s);
            end
        end
    end
    
    % noise B leading to the plaid
    Im_B = 0.5*ones(BGheight_pad, BGlength_master);
    for j = BGchunklength/2:BGchunklength:BGlength_master
        for k = 1:BGdensity
            row = filtSize+randi(s, BGheight);
            if mod(k, 2) == 0
                Im_B(row, j) = 0.5 + 0.5*rand(s); 
            else
                Im_B(row, j) = 0.5 - 0.5*rand(s);
            end
        end
    end
    
    % A B A B A 
    Im_Combined = [Im_A, Im_B, Im_A, Im_B, Im_A];
    
    % 4. Filter with padding logic to match original sharpness
    pad = Im_A; 
    Im_Full = [pad, Im_Combined, pad];
    Imf = conv2(Im_Full, filt2, 'same');
    Imf = Imf(filtSize+1:end-filtSize, BGlength_master+1 : 6*BGlength_master);
    
    % 5. Interpolate to 2048px (Only one interpolation to avoid double-blurring)
    [u, v] = meshgrid(1:size(Imf,2), 1:size(Imf,1));
    [uq, vq] = meshgrid(linspace(1,size(Imf,2), finalBGlength), linspace(1,size(Imf,1), finalBGheight));
    Imf = interp2(u, v, Imf, uq, vq);
    
    Imf = (Imf - 0.5);
    Imf = Imf./max(abs(Imf), [], 'all');
    textures(texID).matrix = 0.5*(1 + BG_contrast*Imf);
end

%%
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

figure; 
grating = textures(6).matrix;
plaid = textures(8).matrix;

% Landmark width in pixels
final_width = round(finalBGlength * texwidth); 
final_height = finalBGheight;

% Rescale landmarks to the background height
[u, v] = meshgrid(1:size(plaid,2), 1:size(plaid,1));
[uq, vq] = meshgrid(linspace(1, size(plaid,2), final_width), linspace(1, size(plaid,1), final_height));
grating_res = interp2(u, v, grating, uq, vq);
plaid_res = interp2(u, v, plaid, uq, vq);

for k = 2:5
    subplot(4,1,k-1)
    % Extract the full 200cm background (2000px if using the 10px/cm fix)
    tex = textures(k).matrix(:, 1:finalBGlength);
    
    % Landmark centers in cm: 40, 80, 120, 160
    % Landmark centers in pixels: (cm / 200) * finalBGlength
    centers_px = [0.20, 0.40, 0.60, 0.80] * finalBGlength;
    
    for idx = 1:4
        c = round(centers_px(idx));
        start_col = c - round(final_width/2) + 1;
        
        % Place gratings at 40 and 120, plaids at 80 and 160
        if mod(idx,2) == 1
            tex(:, start_col:start_col+final_width-1) = grating_res;
        else
            tex(:, start_col:start_col+final_width-1) = plaid_res;
        end
    end
    
    % DISPLAY IN CM: [0 to 200] on X, [0 to 12] on Y
    imagesc([0, 200], [0, 12], tex, [0 1]);
    
    title({['Visible background ', num2str(k-1)], ...
           ['X-Axis in cm (Landmarks at 40, 80, 120, 160)']});
       
    xlabel('Distance (cm)');
    ylabel('Height (cm)');
    colormap gray; 
    axis tight; % Ensure no extra white space
    box off;
end


% save

tex1 = textures(2).matrix;
tex2 = textures(3).matrix;
tex3 = textures(4).matrix;
tex4 = textures(5).matrix;

grating_v = textures(6).matrix;
plaid = textures(8).matrix;
grey = textures(1).matrix;
figure;
imagesc(tex, [0 1]);
axis equal; axis off; colormap(gray);

imwrite(tex1, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-ABABA\BG1_50C20P.jpg');
imwrite(tex2, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-ABABA\BG2_50C20P.jpg');
imwrite(tex3, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-ABABA\BG3_50C20P.jpg');
imwrite(tex4,  '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-ABABA\BG4_50C20P.jpg');
imwrite(grating_v, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-ABABA\grating_vertical_newRun.jpg');
imwrite(plaid, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-ABABA\plaid_newRun.jpg');
imwrite(grey, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\BGTextures-ABABA\grey_newRun.jpg')
% close all
