clear;
eeglab;
close all;

% Cargar FieldTrip
addpath('C:\Users\mdomi\Documents\Toolbox\eeglab2023.0');
eeglab nogui;
addpath(genpath('C:\Users\mdomi\Documents\Toolbox\fieldtrip-20230118'));
ft_defaults;

% Cargar máscaras y matrices grupales
load('electrodePairsByBand.mat');  
load('avgMatricesByBand.mat');     

% Configuración
basePath = 'E:\EEGO\datos\EEG_data\TF data\EmotionRegulationforanalysis3';
conditions = {'REAPPRAISE', 'NEUTRAL', 'NEGATIVE', 'SUPPRESS'};
neutralIdx = find(strcmp(conditions, 'NEUTRAL'));
timeWindow = [0.3, 4.0];
bands = fieldnames(electrodePairsByBand);

% Crear máscara combinada por banda
maskByBand = struct();
for b = 1:length(bands)
    bandName = bands{b};
    combinedMask = [];

    for c = 1:length(conditions)
        condName = conditions{c};
        if ~strcmp(condName, 'NEUTRAL') && isfield(electrodePairsByBand.(bandName), condName)
            currentMask = electrodePairsByBand.(bandName).(condName).maskedMatrix;
            if isempty(combinedMask)
                combinedMask = ~isnan(currentMask) & currentMask > 0;
            else
                combinedMask = combinedMask | (~isnan(currentMask) & currentMask > 0);
            end
        end
    end
    maskByBand.(bandName) = combinedMask;
end

% Inicializar tabla de resultados
results = {};

% Iterar por banda
for b = 1:length(bands)
    bandName = bands{b};
    fprintf('Procesando banda: %s\n', bandName);

    % Límites de frecuencia
    switch bandName
        case 'delta', bandLimits = [1 3];
        case 'theta', bandLimits = [4 8];
        case 'alpha', bandLimits = [9 12];
        case 'beta',  bandLimits = [15 30];
    end

    maskMatrix = maskByBand.(bandName);

    % Iterar por condiciones ≠ NEUTRAL
    for c = 1:length(conditions)
        conditionName = conditions{c};
        if strcmp(conditionName, 'NEUTRAL')
            continue;
        end

        fprintf('  Condición: %s\n', conditionName);

        condPath = fullfile(basePath, conditionName);
        neutralPath = fullfile(basePath, 'NEUTRAL');
        files = dir(fullfile(condPath, 'S*.set'));

        for subjIdx = 1:numel(files)
            subjectName = files(subjIdx).name;
            subjectID = erase(subjectName, '.set');
            fprintf('    Sujeto: %s\n', subjectID);

            % --- Cargar condición experimental ---
            EEG = pop_loadset('filename', subjectName, 'filepath', condPath);
            ft_data = eeglab2fieldtrip(EEG, 'raw');

            cfg = [];
            cfg.latency = timeWindow;
            ft_data = ft_selectdata(cfg, ft_data);

            cfg = [];
            cfg.method = 'mtmfft';
            cfg.taper = 'dpss';
            cfg.output = 'fourier';
            cfg.tapsmofrq = 2;
            cfg.foilim = bandLimits;
            freq = ft_freqanalysis(cfg, ft_data);

            cfg = [];
            cfg.method = 'wpli_debiased';
            connExp = ft_connectivityanalysis(cfg, freq);
            connMatrixExp = squeeze(mean(connExp.wpli_debiasedspctrm, 3));

            % --- Cargar condición NEUTRAL del mismo sujeto ---
            % Reemplazar sufijo de condición por _NEUTRAL
            neutralFile = regexprep(subjectName, ['_' upper(conditionName)], '_NEUTRAL');

            try
                EEG = pop_loadset('filename', neutralFile, 'filepath', neutralPath);
            catch
                warning('No se encontró el archivo NEUTRAL para %s. Saltando...', subjectName);
                continue;
            end

            ft_data = eeglab2fieldtrip(EEG, 'raw');
            cfg.latency = timeWindow;
            ft_data = ft_selectdata(cfg, ft_data);

            cfg = [];
            cfg.method = 'mtmfft';
            cfg.taper = 'dpss';
            cfg.output = 'fourier';
            cfg.tapsmofrq = 2;
            cfg.foilim = bandLimits;
            freq = ft_freqanalysis(cfg, ft_data);

            cfg = [];
            cfg.method = 'wpli_debiased';
            connNeu = ft_connectivityanalysis(cfg, freq);
            connMatrixNeu = squeeze(mean(connNeu.wpli_debiasedspctrm, 3));

            % --- Diferencia: condición – neutral ---
            diffMatrix = connMatrixExp - connMatrixNeu;

            % --- Extraer valores significativos según máscara grupal ---
            valsDiff = diffMatrix(maskMatrix);
            meanDiff = mean(valsDiff, 'omitnan');

            % --- Guardar ---
            results = [results; {subjectID, bandName, conditionName, meanDiff}];
        end
    end
end

% Exportar CSV
T = cell2table(results, 'VariableNames', ...
    {'Subject', 'Band', 'Condition', 'Connectivity_Difference'});
writetable(T, 'subject_connectivity_differences.csv');
fprintf('\n✅ Archivo exportado: subject_connectivity_differences.csv\n');
