function binData = import_BonVisionParamsSparseNoise_bin(filename,grid_size)
% import SparseNoise binary file
% Sam Solomon 
fileID=fopen(filename);
thisBinFile=fread(fileID);
fclose(fileID);


% Translate stimulus into -1:1 scale
stim_matrix = zeros(1,length(thisBinFile));
stim_matrix(thisBinFile==0)=-1;
stim_matrix(thisBinFile==255)=1;
stim_matrix(thisBinFile==128)=0;


% Make a NxM grid from the stimulus log
stim_matrix = reshape(stim_matrix, [grid_size(1), grid_size(2), length(thisBinFile)/grid_size(1)/grid_size(2)]);
stim_matrix = stim_matrix(:,:,1:end-1); % The last 'stimulus'


% Make it a table
binData = table;
for thisTrial = 1:size(stim_matrix,3)
    datatable.stim_matrix{thisTrial,1} = squeeze(stim_matrix(:,:,thisTrial));
end
end