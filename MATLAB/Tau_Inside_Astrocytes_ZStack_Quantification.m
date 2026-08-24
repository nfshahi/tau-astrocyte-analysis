% =========================================================================
% Quantification of Tau Within Astrocytes from Live-Cell Z-Stack Images
% =========================================================================
clear; clc; close all;

%% 1. Set Input, Output, and Analysis Parameters
rootInputFolder = fullfile('RawData');
ctrlFolder = fullfile('RawData', 'Ctrl1'); 
baseOutputFolder = 'AdvancedSheeeeet';

% File naming patterns used to identify the cell, tau, and control images
cellSearchPattern = 'c002';
tauSearchPattern = 'c001';
ctrlSearchPattern = 'c001';

% Parameters used for cell segmentation and tau detection
intensityThreshold = 0.15;
minCellArea = 300;
shrinkPercent = 10;
tauAmplification = 3.0;

%% 2. Create Output Folders
if ~isfolder(rootInputFolder), error('Root input folder does not exist.'); end
if ~isfolder(ctrlFolder), warning('Control folder not found at: %s', ctrlFolder); end

dirInverted = fullfile(baseOutputFolder, 'inverted_dark_masks');
dirComplement = fullfile(baseOutputFolder, 'complement_shrunk_cells');
dirTauMasks = fullfile(baseOutputFolder, 'tau_isolated_masks');
dirDualOverlay = fullfile(baseOutputFolder, 'dual_color_overlays');
dirStats = fullfile(baseOutputFolder, 'quantification_results');

subDirs = {dirInverted, dirComplement, dirTauMasks, dirDualOverlay, dirStats};
for s = 1:length(subDirs)
    if ~isfolder(subDirs{s}), mkdir(subDirs{s}); end
end

%% 3. Identify Input Subfolders and Control Images
subFolders = dir(rootInputFolder);
subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name}, {'.', '..'}));

if isempty(subFolders)
    subFolders = struct('name', {''}, 'folder', {rootInputFolder});
end

% Load the control images used to determine the tau intensity threshold
ctrlFiles = dir(fullfile(ctrlFolder, ['*', ctrlSearchPattern, '*.tif*']));
if isempty(ctrlFiles)
    warning('No control files found matching pattern "*%s*.tif" in %s.', ctrlSearchPattern, ctrlFolder);
end

summaryData = table();

%% 4. Process Each Input Subfolder
for f = 1:length(subFolders)
    if isempty(subFolders(f).folder)
        currentInFolder = subFolders(f).name;
        subFolderName = '';
    else
        currentInFolder = fullfile(subFolders(f).folder, subFolders(f).name);
        subFolderName = subFolders(f).name;
    end
    
    cellFiles = dir(fullfile(currentInFolder, ['*', cellSearchPattern, '*.tif*']));
    if isempty(cellFiles), continue; end
    
    fprintf('--> Processing subfolder: %s (%d images found)\n', subFolderName, length(cellFiles));
    
    if ~isempty(subFolderName)
        outInvSub = fullfile(dirInverted, subFolderName);
        outCompSub = fullfile(dirComplement, subFolderName);
        outTauSub = fullfile(dirTauMasks, subFolderName);
        outDualSub = fullfile(dirDualOverlay, subFolderName);
        if ~isfolder(outInvSub), mkdir(outInvSub); end
        if ~isfolder(outCompSub), mkdir(outCompSub); end
        if ~isfolder(outTauSub), mkdir(outTauSub); end
        if ~isfolder(outDualSub), mkdir(outDualSub); end
    else
        outInvSub = dirInverted;
        outCompSub = dirComplement;
        outTauSub = dirTauMasks;
        outDualSub = dirDualOverlay;
    end
    
    %% 5. Process Each Cell and Tau Image Pair
    for i = 1:length(cellFiles)
        cellFileName = cellFiles(i).name;
        [~, baseName, ~] = fileparts(cellFileName);
        
        tauFileName = strrep(cellFileName, cellSearchPattern, tauSearchPattern);
        cellFullPath = fullfile(currentInFolder, cellFileName);
        tauFullPath = fullfile(currentInFolder, tauFileName);
        
        if ~isfile(tauFullPath)
            fprintf('  Skipping %s: Matching Tau file not found.\n', cellFileName);
            continue;
        end
        
        if ~isempty(ctrlFiles)
            ctrlIdx = mod(i-1, length(ctrlFiles)) + 1;
            ctrlFileName = ctrlFiles(ctrlIdx).name;
            ctrlFullPath = fullfile(ctrlFolder, ctrlFileName);
        else
            ctrlFileName = 'None';
            ctrlFullPath = '';
        end
        
        try
            % Calculate the tau detection threshold from the control image
            if ~isempty(ctrlFullPath) && isfile(ctrlFullPath)
                ctrlImg = imread(ctrlFullPath);
                if size(ctrlImg, 3) == 3, ctrlImg = rgb2gray(ctrlImg); end
                meanCtrlIntensity = mean(im2double(ctrlImg(:)));
            else
                meanCtrlIntensity = 0.05; 
            end
            tauThreshold = meanCtrlIntensity * tauAmplification;
            
            % Generate a mask identifying the cell regions
            imgCell = imread(cellFullPath);
            if size(imgCell, 3) == 3, imgCellGray = rgb2gray(imgCell); else, imgCellGray = imgCell; end
            imgCellDouble = im2double(imgCellGray);
            
            binMaskInverted = imgCellDouble < intensityThreshold;
            binMaskInverted = bwareaopen(binMaskInverted, minCellArea);
            binMaskInverted = imfill(binMaskInverted, 'holes');
            
            % Invert the mask to obtain the cell regions
            binMaskComplement = ~binMaskInverted;
            
            % Shrink the cell regions before tau quantification
            if shrinkPercent > 0
                cc = bwconncomp(binMaskComplement);
                refinedComplement = false(size(binMaskComplement));
                for c = 1:cc.NumObjects
                    singleCellIdx = cc.PixelIdxList{c};
                    singleCellMask = false(size(binMaskComplement));
                    singleCellMask(singleCellIdx) = true;
                    stats = regionprops(singleCellMask, 'EquivDiameter');
                    if ~isempty(stats)
                        cellDiam = stats.EquivDiameter;
                        retractionRadius = max(1, round((cellDiam / 2) * (shrinkPercent / 100)));
                        se = strel('disk', retractionRadius);
                        erodedCell = imerode(singleCellMask, se);
                        if any(erodedCell(:))
                            refinedComplement = refinedComplement | erodedCell;
                        else
                            refinedComplement = refinedComplement | singleCellMask;
                        end
                    end
                end
                binMaskComplement = refinedComplement;
            end
            
            % Create a mask of tau-positive pixels
            imgTau = imread(tauFullPath);
            if size(imgTau, 3) == 3, imgTauGray = rgb2gray(imgTau); else, imgTauGray = imgTau; end
            tauMask = im2double(imgTauGray) > tauThreshold;
            tauMask = bwareaopen(tauMask, 10);
            
            % Quantify the amount of tau-positive area within the cell mask
            totalCellPixels = sum(binMaskComplement(:));
            tauInsideCells = tauMask & binMaskComplement;
            totalTauInsideCells = sum(tauInsideCells(:));
            
            if totalCellPixels > 0
                pctTauInCells = (totalTauInsideCells / totalCellPixels) * 100;
            else
                pctTauInCells = 0;
            end
            
            % Create images showing the masks and tau/cell overlap
            if isa(imgCellGray, 'uint8'), maxVal = 255;
            elseif isa(imgCellGray, 'uint16'), maxVal = 65535;
            else, maxVal = 1; end
            
            imgRGB = cat(3, imgCellGray, imgCellGray, imgCellGray);
            
            % Inverted cell mask
            R1 = imgRGB(:,:,1); G1 = imgRGB(:,:,2); B1 = imgRGB(:,:,3);
            R1(binMaskInverted) = maxVal; G1(binMaskInverted) = 0; B1(binMaskInverted) = 0;
            imwrite(cat(3, R1, G1, B1), fullfile(outInvSub, [baseName, '_inverted.tif']), 'Compression', 'none');
            
            % Shrunk cell mask
            R2 = imgRGB(:,:,1); G2 = imgRGB(:,:,2); B2 = imgRGB(:,:,3);
            R2(binMaskComplement) = maxVal; G2(binMaskComplement) = 0; B2(binMaskComplement) = 0;
            imwrite(cat(3, R2, G2, B2), fullfile(outCompSub, [baseName, '_complement_shrunk.tif']), 'Compression', 'none');
            
            % Tau-positive mask
            tauRGB = cat(3, imgTauGray, imgTauGray, imgTauGray);
            RTau = tauRGB(:,:,1); GTau = tauRGB(:,:,2); BTau = tauRGB(:,:,3);
            RTau(tauMask) = maxVal; GTau(tauMask) = 0; BTau(tauMask) = 0;
            imwrite(cat(3, RTau, GTau, BTau), fullfile(outTauSub, [baseName, '_tau_isolated.tif']), 'Compression', 'none');
            
            % Overlay showing cell regions and tau-positive regions
            R_dual = imgCellGray; G_dual = imgCellGray; B_dual = imgCellGray;
            R_dual(binMaskComplement) = 0; G_dual(binMaskComplement) = maxVal; B_dual(binMaskComplement) = 0;
            R_dual(tauMask) = maxVal; G_dual(tauMask) = 0; BTau(tauMask) = 0;
            dualOverlay = cat(3, R_dual, G_dual, B_dual);
            imwrite(dualOverlay, fullfile(outDualSub, [baseName, '_dual_tau_cells.tif']), 'Compression', 'none');
            
            % Save the measurements for each image pair
            newRow = table(string(subFolderName), string(cellFileName), string(ctrlFileName), meanCtrlIntensity, tauThreshold, totalCellPixels, totalTauInsideCells, pctTauInCells, ...
                'VariableNames', {'SubFolder', 'CellImage', 'ControlUsed', 'CtrlMeanIntensity', 'ComputedTauThreshold', 'CellPixels', 'TauInsidePixels', 'PercentTauInCells'});
            summaryData = [summaryData; newRow];
            
            fprintf('    Done: %s | Tau in Cells: %.2f%%\n', cellFileName, pctTauInCells);
            
        catch ME
            fprintf('    Error on %s: %s\n', cellFileName, ME.message);
        end
    end
end

%% 6. Export Quantification Results
csvPath = fullfile(dirStats, 'hierarchical_tau_cell_quantification_summary.csv');
writetable(summaryData, csvPath);

excelPath = fullfile(dirStats, 'hierarchical_tau_cell_quantification_summary.xlsx');
writetable(summaryData, excelPath, 'Sheet', 'All_Files');

% Calculate average measurements for each subfolder
if ~isempty(summaryData)
    try
        folderNames = unique(summaryData.SubFolder);
        avgSubfolderTable = table();
        
        for k = 1:length(folderNames)
            fName = folderNames(k);
            matches = summaryData.SubFolder == fName;
            
            meanPct = mean(summaryData.PercentTauInCells(matches));
            meanCellPix = mean(summaryData.CellPixels(matches));
            meanTauPix = mean(summaryData.TauInsidePixels(matches));
            fileCount = sum(matches);
            
            avgRow = table(fName, fileCount, meanCellPix, meanTauPix, meanPct, ...
                'VariableNames', {'SubFolder', 'TotalImages', 'AvgCellPixels', 'AvgTauInsidePixels', 'AvgPercentTauInCells'});
            avgSubfolderTable = [avgSubfolderTable; avgRow];
        end
        
        writetable(avgSubfolderTable, excelPath, 'Sheet', 'Subfolder_Averages');
        fprintf('Successfully generated subfolder averages summary sheet in Excel.\n');
    catch
        fprintf('Note: Could not create Subfolder_Averages sheet, but main Excel file was saved.\n');
    end
end

fprintf('\nBatch processing complete across all subfolders! Summary files saved to %s\n', dirStats);
