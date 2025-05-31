clear;
clc;

% Carpeta con los archivos Excel
dataFolder = 'C:\Users\mdomi\Documents\Submissions\Emotion Regulation - Connectivity\frequencies\';  % <-- cambia esta ruta según corresponda
bands = {'delta', 'theta', 'alpha', 'beta'};

% Inicializar
models = struct();

% Recorrer cada banda
for i = 1:length(bands)
    band = bands{i};
    filePath = fullfile(dataFolder, [band '.xlsx']);

    % Leer datos
    opts = detectImportOptions(filePath, 'NumHeaderLines', 0);
    opts = setvartype(opts, {'AVD', 'ANX', 'wPLI'}, 'double');
    T = readtable(filePath, opts);

    % Asegurar que Condition sea categórica
    T.Condition = categorical(T.Condition);

    % Ajustar modelo LME
    lme = fitlme(T, ...
        'wPLI ~ Condition * AVD + Condition * ANX + (1 | SubjectID)');

    % Mostrar resumen
    fprintf('\n===== %s band =====\n', upper(band));
    disp(lme);

    % Guardar el modelo
    models.(band) = lme;
end
