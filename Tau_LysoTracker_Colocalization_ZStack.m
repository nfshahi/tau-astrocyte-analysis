clear; clc; close all

% Set input and output folders
mainFolder = fullfile(pwd, 'Rawdata'); 
outputFolder = fullfile(pwd, 'OUTPUT');

% Check that the input folder exists
if ~exist(mainFolder, 'dir')
    error('The folder "Rawdata" does not exist in the current directory (%s).', pwd);
end

% Create the output folder if needed
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Find subfolders containing the image data
dirContents = dir(mainFolder);
isSubdir = [dirContents.isdir];
subFolders = dirContents(isSubdir);

% Remove the current and parent directory entries
validFolderIdx = ~strcmp({subFolders.name}, '.') & ~strcmp({subFolders.name}, '..');
subFolders = subFolders(validFolderIdx);

if isempty(subFolders)
    error('No subfolders found inside "Rawdata". Please ensure your data is organized in subfolders.');
end

% Initialize storage for correlation results
compiledData = struct();
maxRows = 0;

% =========================================================================
% Analysis parameters
% =========================================================================
% Median filter settings
useMedianFilter = true;
medFilterSize   = 3;

% Polynomial background correction settings
usePolyCorrection = true;
polyOrder         = 2;

% Intensity thresholds applied after image correction
thresh_c001    = 0;
thresh_c002    = 0;
% =========================================================================

% Calculate correlations for each subfolder
for i = 1:length(subFolders)
    currentSubfolderName = subFolders(i).name;
    currentSubfolderPath = fullfile(mainFolder, currentSubfolderName);
    
    fprintf('Processing subfolder: %s...\n', currentSubfolderName);
    
    % Find c001 images
    c1Files = dir(fullfile(currentSubfolderPath, '*c001*.tif'));
    if isempty(c1Files)
        c1Files = dir(fullfile(currentSubfolderPath, '*C001*.tif')); 
    end
    
    % Store correlations for the current subfolder
    folderCorrelations = [];
    
    % Match each c001 image with its corresponding c002 image
    for j = 1:length(c1Files)
        c1Name = c1Files(j).name;
        
        c2Name = strrep(c1Name, 'c001', 'c002');
        c2Name = strrep(c2Name, 'C001', 'C002'); 
        
        c2Path = fullfile(currentSubfolderPath, c2Name);
        
        % Process the image pair if the matching c002 image exists
        if exist(c2Path, 'file')
            c1Path = fullfile(currentSubfolderPath, c1Name);
            
            % Read images as double precision
            img1 = double(imread(c1Path));
            img2 = double(imread(c2Path));
            
            % Apply median filtering
            if useMedianFilter
                img1 = medfilt2(img1, [medFilterSize medFilterSize]);
                img2 = medfilt2(img2, [medFilterSize medFilterSize]);
            end
            
            % Apply polynomial background correction
            if usePolyCorrection && polyOrder >= 1
                img1 = img1 - fit2DPolynomial(img1, polyOrder);
                img2 = img2 - fit2DPolynomial(img2, polyOrder);
                
                % Remove negative values after background subtraction
                img1(img1 < 0) = 0;
                img2(img2 < 0) = 0;
            end
            
            % Convert images to vectors for correlation analysis
            vec1 = img1(:);
            vec2 = img2(:);
            
            % Apply intensity thresholds
            validPixels = (vec1 > thresh_c001) & (vec2 > thresh_c002);
            
            % Calculate Pearson correlation for the remaining pixels
            if sum(validPixels) > 10
                corrMatrix = corrcoef(vec1(validPixels), vec2(validPixels));
                r_val = corrMatrix(1, 2);
            else
                r_val = NaN;
            end
            
            % Store the correlation value
            folderCorrelations = [folderCorrelations; r_val]; %#ok<AGROW>
        end
    end
    
    % Store results for the current subfolder
    validFieldName = matlab.lang.makeValidName(currentSubfolderName);
    compiledData.(validFieldName).originalName = currentSubfolderName;
    compiledData.(validFieldName).values = folderCorrelations;
    
    % Track the maximum number of image pairs
    if length(folderCorrelations) > maxRows
        maxRows = length(folderCorrelations);
    end
end

% Build the results table
totalRows = 1 + maxRows + 1 + 1; 
fields = fieldnames(compiledData);
finalCellArray = cell(totalRows, length(fields));

% Populate the table
for col = 1:length(fields)
    fieldName = fields{col};
    folderInfo = compiledData.(fieldName);
    
    % Set the column header
    finalCellArray{1, col} = folderInfo.originalName;
    
    % Add individual correlation values
    numVals = length(folderInfo.values);
    for row = 1:numVals
        finalCellArray{1 + row, col} = folderInfo.values(row);
    end
    
    finalCellArray{1 + maxRows + 1, col} = '';
    
    % Calculate the average correlation
    if numVals > 0
        avgVal = mean(folderInfo.values, 'omitnan');
        finalCellArray{totalRows, col} = avgVal;
    else
        finalCellArray{totalRows, col} = NaN;
    end
end

% Add row labels
rowLabels = cell(totalRows, 1);
rowLabels{1} = 'Pair Index';
for r = 1:maxRows
    rowLabels{1 + r} = ['Pair ', num2str(r)];
end
rowLabels{1 + maxRows + 1} = '';
rowLabels{totalRows} = 'AVERAGE';

finalCellArray = [rowLabels, finalCellArray];

% Export results to Excel
outputTable = cell2table(finalCellArray(2:end, :), 'VariableNames', finalCellArray(1, :));
excelFilePath = fullfile(outputFolder, 'All_Subfolders_Correlations.xlsx');
writetable(outputTable, excelFilePath);

fprintf('-----------------------------------------\n');
fprintf('All processing complete!\nSaved compiled results to: OUTPUT/All_Subfolders_Correlations.xlsx\n');

% =========================================================================
% Helper function for polynomial background correction
% =========================================================================
function bg = fit2DPolynomial(img, order)
    [h, w] = size(img);
    [X, Y] = meshgrid(1:w, 1:h);
    X_vec = X(:);
    Y_vec = Y(:);
    Z_vec = img(:);
    
    % Build the polynomial design matrix
    A = [];
    for i = 0:order
        for j = 0:(order-i)
            A = [A, (X_vec.^i) .* (Y_vec.^j)]; %#ok<AGROW>
        end
    end
    
    % Calculate polynomial coefficients
    coeffs = A \ Z_vec;
    
    % Reconstruct the background surface
    bg = reshape(A * coeffs, h, w);
end
