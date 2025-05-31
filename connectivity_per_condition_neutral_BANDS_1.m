clear;
eeglab;
close all;

% Cargar FieldTrip
addpath('C:\Users\mdomi\Documents\Toolbox\eeglab2023.0');
eeglab nogui;
addpath(genpath('C:\Users\mdomi\Documents\Toolbox\fieldtrip-20230118'));
ft_defaults;

% Configuración inicial
basePath = 'E:\EEGO\datos\EEG_data\TF data\EmotionRegulationforanalysis';
conditions = {'REAPPRAISE', 'NEUTRAL', 'NEGATIVE', 'SUPPRESS'};
timeWindow = [0.3, 4.0]; % Ventana de tiempo en segundos
bands = struct(...
    'delta', [1 3], ...
    'theta', [4 8], ...
    'alpha', [9 12], ...
    'beta', [15 30] ...
);
numConditions = length(conditions);

% Inicializar para almacenar matrices promedio por banda
avgMatricesByBand = struct();

% Procesar por banda y condición
for bandName = fieldnames(bands)'
    bandName = bandName{1}; % Nombre de la banda (e.g., 'delta')
    bandLimits = bands.(bandName); % Límites de frecuencia para la banda
    fprintf('Procesando banda: %s (%d-%d Hz)\n', bandName, bandLimits(1), bandLimits(2));

    avgMatrices = cell(numConditions, 1); % Inicializar matrices promedio para esta banda

    for condIdx = 1:numConditions
        condPath = fullfile(basePath, conditions{condIdx});
        files = dir(fullfile(condPath, 'S*.set')); % Buscar archivos por sujeto

        fprintf('  Procesando condición: %s\n', conditions{condIdx});
        allMatrices = [];

        for subjIdx = 1:numel(files)
            subjectName = files(subjIdx).name;
            fprintf('    Procesando sujeto: %s\n', subjectName);

            % Cargar datos
            EEG = pop_loadset('filename', subjectName, 'filepath', condPath);

            % Convertir a FieldTrip
            ft_data = eeglab2fieldtrip(EEG, 'raw');

            % Filtrar por ventana de tiempo
            cfg = [];
            cfg.latency = timeWindow;
            ft_data = ft_selectdata(cfg, ft_data);

            % Análisis de frecuencia
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.taper = 'dpss';
            cfg.output = 'fourier';
            cfg.tapsmofrq = 2;
            cfg.foilim = bandLimits; % Banda de frecuencia actual
            freq = ft_freqanalysis(cfg, ft_data);

            % Análisis de conectividad
            cfg = [];
            cfg.method = 'wpli_debiased';
            connectivity = ft_connectivityanalysis(cfg, freq);

            % Matriz de conectividad promedio
            connMatrix = squeeze(mean(connectivity.wpli_debiasedspctrm, 3));
            allMatrices = cat(3, allMatrices, connMatrix); % Acumular matrices
        end

        % Calcular la matriz promedio para la condición
        if ~isempty(allMatrices)
            avgMatrices{condIdx} = mean(allMatrices, 3, 'omitnan');
        else
            avgMatrices{condIdx} = NaN; % Si no hay datos, asignar NaN
        end
    end

    % Guardar las matrices promedio para esta banda
    avgMatricesByBand.(bandName) = avgMatrices;

    % Visualizar las matrices promedio por condición para esta banda
    for condIdx = 1:numConditions
        if ~isempty(avgMatrices{condIdx}) && ~all(isnan(avgMatrices{condIdx}(:)))
            figure;
            imagesc(avgMatrices{condIdx}, [0, 0.5]); % Ajusta el rango de color
            colorbar;
            title(sprintf('Average Matrix - Band: %s, Condition: %s', bandName, conditions{condIdx}));
            xlabel('Electrodes');
            ylabel('Electrodes');
            set(gca, 'FontSize', 14);                   % Aumenta tamaño de texto en los ejes
            set(colorbar, 'FontSize', 14);              % Aumenta tamaño de texto del colorbar

            axis square;
            colormap(jet);
        else
            fprintf('Sin datos válidos para la condición: %s en la banda: %s\n', conditions{condIdx}, bandName);
        end
    end
end

% Visualizar las diferencias entre condiciones experimentales y neutral por banda
neutralIdx = find(strcmp(conditions, 'NEUTRAL'));

for bandName = fieldnames(bands)'
    bandName = bandName{1};
    avgMatrices = avgMatricesByBand.(bandName);

    for condIdx = 1:numConditions
        if condIdx ~= neutralIdx && ~isempty(avgMatrices{neutralIdx}) && ~isempty(avgMatrices{condIdx}) ...
                && ~all(isnan(avgMatrices{neutralIdx}(:))) && ~all(isnan(avgMatrices{condIdx}(:)))

            % Calcular la diferencia
            diffMatrix = avgMatrices{condIdx} - avgMatrices{neutralIdx};

            % Visualización de la matriz de diferencias
            figure;
            set(gcf, 'Color', 'w');  % Establece fondo blanco para la figura

            imagesc(diffMatrix, [-0.1, 0.1]); % Ajusta el rango de color según tus datos
            cb = colorbar;
            ylabel(cb, 'Δ wPLI');
            title(sprintf('Difference - Band: %s, %s vs Neutral', bandName, conditions{condIdx}));
            
            xlabel('Electrodes');
            ylabel('Electrodes');
            xtickangle(90);
            set(gca, 'FontSize', 18);                   % Aumenta tamaño de texto en los ejes
            set(cb, 'FontSize', 18);              % Aumenta tamaño de texto del colorbar

            axis square;
            colormap(jet);

            % Incluir etiquetas de los electrodos
            if exist('EEG', 'var') && isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
                electrodeLabels = {EEG.chanlocs.labels};
                set(gca, 'XTick', 1:length(electrodeLabels), 'XTickLabel', electrodeLabels, ...
                         'YTick', 1:length(electrodeLabels), 'YTickLabel', electrodeLabels);
                 axis square;
            end
        else
            fprintf('Sin datos válidos para la diferencia entre %s y Neutral en la banda: %s\n', conditions{condIdx}, bandName);
        end
    end
end

%% ------------------------------------------- UMBRAL -------------------------------------------
% --------------------------------------------

% ------------------------------------------- UMBRAL POR BANDA -------------------------------------------
% --------------------------------------------

% Configuración del umbral
percentileThreshold = 95; % Percentil usado como umbral
electrodeLabels = {EEG.chanlocs.labels}; % Asegúrate de que esté definido previamente
neutralIdx = find(strcmp(conditions, 'NEUTRAL'));

% Inicializar estructura para almacenar resultados
electrodePairsByBand = struct();
figureCount = 0; % Contador de figuras generadas

% Iterar sobre bandas
for bandName = fieldnames(avgMatricesByBand)'
    bandName = bandName{1};
    fprintf('Procesando banda: %s\n', bandName);

    % Obtener matrices promedio para la banda actual
    avgMatrices = avgMatricesByBand.(bandName);

    % Inicializar lista para almacenar todos los valores positivos de diferencia
    allPositiveDiffs = [];

    % Recopilar diferencias positivas para calcular umbral por banda
    for condIdx = 1:numel(conditions)
        if condIdx ~= neutralIdx && ~isempty(avgMatrices{neutralIdx}) && ~isempty(avgMatrices{condIdx}) ...
                && ~all(isnan(avgMatrices{neutralIdx}(:))) && ~all(isnan(avgMatrices{condIdx}(:)))

            diffMatrix = avgMatrices{condIdx} - avgMatrices{neutralIdx};
            diffValues = diffMatrix(:);
            diffValues = diffValues(~isnan(diffValues)); % Eliminar NaN
            allPositiveDiffs = [allPositiveDiffs; diffValues(diffValues > 0)];
        end
    end

    % Calcular umbral único por banda
    if isempty(allPositiveDiffs)
        fprintf('No hay diferencias positivas para calcular umbral en la banda %s.\n', bandName);
        continue;
    end
    threshold = prctile(allPositiveDiffs, percentileThreshold);

    % Almacenar el umbral por banda
    electrodePairsByBand.(bandName).threshold = threshold;

    % Aplicar el umbral a cada comparación dentro de la banda
    for condIdx = 1:numel(conditions)
        if condIdx ~= neutralIdx && ~isempty(avgMatrices{neutralIdx}) && ~isempty(avgMatrices{condIdx}) ...
                && ~all(isnan(avgMatrices{neutralIdx}(:))) && ~all(isnan(avgMatrices{condIdx}(:)))

            diffMatrix = avgMatrices{condIdx} - avgMatrices{neutralIdx};

            % Atenuar valores por debajo del umbral
            attenuatedMatrix = diffMatrix;
            attenuatedMatrix(attenuatedMatrix < threshold) = attenuatedMatrix(attenuatedMatrix < threshold) * 0.3;
            attenuatedMatrix(attenuatedMatrix < 0) = NaN; % Excluir decrementos

            % Guardar resultados
            electrodePairsByBand.(bandName).(conditions{condIdx}) = struct(...
                'maskedMatrix', attenuatedMatrix ...
            );

            % Mostrar matriz enmascarada
            if any(~isnan(attenuatedMatrix(:)))
                figureCount = figureCount + 1;
                figure;
                set(gcf, 'Color', 'w');  % Establece fondo blanco para la figura

                imagesc(attenuatedMatrix, [0, max(allPositiveDiffs)]);
                cb = colorbar;
                ylabel(cb, 'Δ wPLI');
                set(cb, 'FontSize', 18);

                title(sprintf('Band: %s, %s vs Neutral (Threshold: %.4f)', ...
                    bandName, conditions{condIdx}, threshold));
              

                xlabel('Electrodes');
                ylabel('Electrodes');
                axis square;
                colormap(jet);
                set(gca, 'XTick', 1:length(electrodeLabels), 'XTickLabel', electrodeLabels, ...
                         'YTick', 1:length(electrodeLabels), 'YTickLabel', electrodeLabels);
                xtickangle(90);
                set(gca, 'FontSize', 18);                   % Aumenta tamaño de texto en los ejes
                set(cb, 'FontSize', 18);              % Aumenta tamaño de texto del colorbar

            else
                fprintf('La matriz atenuada para la banda %s y condición %s no tiene valores válidos.\n', ...
                        bandName, conditions{condIdx});
            end
        end
    end
end

% Mostrar umbrales por banda
for bandName = fieldnames(electrodePairsByBand)'
    bandName = bandName{1};
    fprintf('Banda: %s - Umbral aplicado: %.4f\n', ...
        bandName, electrodePairsByBand.(bandName).threshold);
end

fprintf('Total de figuras generadas: %d\n', figureCount);
%%
% Mostrar pares de electrodos con conectividad aumentada por banda y condición
for bandName = fieldnames(electrodePairsByBand)'
    bandName = bandName{1};
    fprintf('\nPares de electrodos con conectividad aumentada en la banda: %s\n', bandName);

    threshold = electrodePairsByBand.(bandName).threshold;
    avgMatrices = avgMatricesByBand.(bandName);

    for condIdx = 1:numel(conditions)
        condName = conditions{condIdx};

        if condIdx ~= neutralIdx && isfield(electrodePairsByBand.(bandName), condName)
            % Calcular matriz de diferencia original
            diffMatrix = avgMatrices{condIdx} - avgMatrices{neutralIdx};

            % Encontrar pares que superan el umbral
            [rowIdx, colIdx] = find(diffMatrix > threshold);

            fprintf('  Condición: %s\n', condName);
            if isempty(rowIdx)
                fprintf('    Sin pares de electrodos por encima del umbral (%.4f).\n', threshold);
            else
                for i = 1:length(rowIdx)
                    el1 = electrodeLabels{rowIdx(i)};
                    el2 = electrodeLabels{colIdx(i)};
                    val = diffMatrix(rowIdx(i), colIdx(i));
                    fprintf('    %s–%s: %.4f\n', el1, el2, val);
                end
            end
        end
    end
end
%%
% Guardar resultados para usar en el segundo script
save('electrodePairsByBand.mat', 'electrodePairsByBand');
save('avgMatricesByBand.mat', 'avgMatricesByBand');
