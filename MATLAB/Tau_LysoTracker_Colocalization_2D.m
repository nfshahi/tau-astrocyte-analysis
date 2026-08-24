%% ---------------------------------------------------------
% Batch correlation across subfolders
% Reference-based thresholding:
% CY5 threshold from C*, GFP threshold from B*
% Correlation only for selected parameters (e.g., A* and D*)
% Export results to Excel and plot vs time
% ---------------------------------------------------------

clear; clc; close all

%% -------------------------------
% USER PARAMETERS
%% -------------------------------
mainFolder = 'RAWDATA';  % Main folder containing subfolders 1,2,3,...
CY5_refLetter = 'C';              % Reference group for CY5 threshold
GFP_refLetter = 'B';              % Reference group for GFP threshold
analyzeLetters = {'A','B','C','D'};       % Only analyze these parameters

% Output Excel file
excelFile = 'CorrelationResults.xlsx';

%% -------------------------------
% List all subfolders (time points)
%% -------------------------------
subFolders = dir(mainFolder);
subFolders = subFolders([subFolders.isdir]);
subFolders = subFolders(~ismember({subFolders.name},{'.','..'}));

% Convert folder names to numbers and sort
folderNums = str2double({subFolders.name});   % convert names to numeric
[~, sortIdx] = sort(folderNums);              % get sorting index
folderNames = {subFolders(sortIdx).name};    % sorted folder names as strings

%% -------------------------------
% Preallocate result storage
%% -------------------------------
allParams = {};  % store all A*, D* parameter names encountered
timePoints = folderNums(sortIdx);  % sorted times
allResults = []; % will store numeric correlation results

%% -------------------------------
% Iterate over subfolders (time points)
%% -------------------------------
for f = 1:length(folderNames)
    subFolder = folderNames{f};
    subPath = fullfile(mainFolder, subFolder);

    fprintf('\nProcessing folder: %s\n', subFolder);

    % -------------------------------
    % Load images
    % -------------------------------
    files = dir(fullfile(subPath,'*.tif'));
    if isempty(files)
        warning('No TIFF files in folder %s', subFolder);
        continue;
    end

    data = struct();
    for k = 1:length(files)
        fname = files(k).name;
        fpath = fullfile(subPath, fname);

        token = regexp(fname, '^([A-D]\d+)_.*_(\d+)_(CY5|GFP)_', 'tokens');
        if isempty(token)
            continue;
        end

        param   = token{1}{1};
        tileNum = str2double(token{1}{2});
        channel = token{1}{3};

        img = im2double(imread(fpath));

        if ~isfield(data,param)
            data.(param).CY5 = cell(1,16);
            data.(param).GFP = cell(1,16);
        end
        data.(param).(channel){tileNum} = img;
    end

    % -------------------------------
    % Compute reference thresholds
    % -------------------------------
    CY5_vals = [];
    GFP_vals = [];
    params = fieldnames(data);
    for p = 1:length(params)
        param = params{p};

        % CY5 reference
        if startsWith(param,CY5_refLetter)
            CY5_tiles = data.(param).CY5;
            for t = 1:numel(CY5_tiles)
                if ~isempty(CY5_tiles{t})
                    CY5_vals = [CY5_vals; CY5_tiles{t}(:)]; %#ok<SAGROW>
                end
            end
        end

        % GFP reference
        if startsWith(param,GFP_refLetter)
            GFP_tiles = data.(param).GFP;
            for t = 1:numel(GFP_tiles)
                if ~isempty(GFP_tiles{t})
                    GFP_vals = [GFP_vals; GFP_tiles{t}(:)]; %#ok<SAGROW>
                end
            end
        end
    end

    CY5_thresh = mean(CY5_vals);
    GFP_thresh = mean(GFP_vals);

    %% -------------------------------
    % Correlation calculation for selected parameters
    %% -------------------------------
    corrResults = struct();
    for p = 1:length(params)
        param = params{p};

        if ~any(startsWith(param, analyzeLetters))
            continue;
        end

        % Keep track of all parameters
        if ~ismember(param, allParams)
            allParams{end+1} = param; %#ok<SAGROW>
        end

        CY5_tiles = {};
        GFP_tiles = {};

        for t = 1:length(data.(param).CY5)
            hasCY5 = ~isempty(data.(param).CY5{t});
            hasGFP = ~isempty(data.(param).GFP{t});
        
            if hasCY5 && hasGFP
                CY5_tiles{end+1} = data.(param).CY5{t}; %#ok<SAGROW>
                GFP_tiles{end+1} = data.(param).GFP{t}; %#ok<SAGROW>
            else
                % Print ignored tiles
                if ~hasCY5 && hasGFP
                    fprintf('%s: Tile %d missing CY5 → ignored in folder %s\n', param, t, subFolder);
                elseif hasCY5 && ~hasGFP
                    fprintf('%s: Tile %d missing GFP → ignored in folder %s\n', param, t, subFolder);
% % % %                 elseif ~hasCY5 && ~hasGFP
% % % %                     fprintf('%s: Tile %d missing CY5 and GFP → ignored in folder %s\n', param, t, subFolder);
                end
            end
        end

        if isempty(CY5_tiles)
            corrVal = NaN;
        else
            CY5_stitched = cat(2, CY5_tiles{:});
            GFP_stitched = cat(2, GFP_tiles{:});

            v1 = CY5_stitched(:);
            v2 = GFP_stitched(:);

            % Apply reference thresholds
            mask = (v1 > CY5_thresh) & (v2 > GFP_thresh);
            v1_use = v1(mask);
            v2_use = v2(mask);

            if isempty(v1_use)
                corrVal = NaN;
            else
%%%%%%----------%%%%%%  matlab corrcoef function (pearson correlation)
                R = corrcoef(v1_use,v2_use);
                corrVal = R(1,2);
            end
        end

        corrResults.(param) = corrVal;
    end

    % -------------------------------
    % Store results in numeric array row
    % -------------------------------
    row = NaN(1,length(allParams));  % numeric row
    for a = 1:length(allParams)
        paramName = allParams{a};
        if isfield(corrResults,paramName)
            row(a) = corrResults.(paramName);
        end
    end

    allResults(f,:) = row;
end

%% -------------------------------
% Create final table
%% -------------------------------
T = table(timePoints', 'VariableNames', {'Time'});
for a = 1:length(allParams)
    T.(allParams{a}) = allResults(:,a);
end

% -------------------------------
% Export to Excel
% -------------------------------
writetable(T, excelFile);
fprintf('\nResults exported to %s\n', excelFile);

%% -------------------------------
% Plot correlation vs time
%% -------------------------------
figure('Name','Correlation vs Time','NumberTitle','off'); hold on;
colors = lines(length(allParams));

for a = 1:length(allParams)
    plot(T.Time, T.(allParams{a}), '-o','Color',colors(a,:),'DisplayName',allParams{a});
end

xlabel('Time (hours)');
ylabel('Pearson correlation');
title('Correlation vs Time');
legend('show','Location','best');
grid on;

