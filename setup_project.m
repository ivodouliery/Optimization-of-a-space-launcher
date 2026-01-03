
% Script to add all project folders to MATLAB path
% Run this once per session or add to your startup.m

currentPath = fileparts(mfilename('fullpath'));
addpath(fullfile(currentPath, 'SQP'));
addpath(fullfile(currentPath, 'Cas', 'MHW4D'));
addpath(fullfile(currentPath, 'Cas', 'Ariane1_Masse'));
addpath(fullfile(currentPath, 'Cas', 'PE'));
addpath(fullfile(currentPath, 'Cas', 'Lanceur'));

fprintf('Le code est prêt à être utilisé\n');