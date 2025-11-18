clear all; close all; clc;
%% source folders containing scripts not in this folder
addpath(genpath('FE_routines'));addpath(genpath('functions'));
addpath(genpath('mesh_utilities'));
addpath(genpath('optimization'));
addpath(genpath('utilities'));
addpath(genpath('plotting'));

global OPT FE

path_data = 'data/val/';

OPT.beta_min = 0;
OPT.beta_max = 30;
OPT.iter = 1;
OPT.eta_H = 0.5;
%% Initialization
% ============================
no_input = 50; % number of samples
% Latin Hypercube Sampling in [0,1]
lhs = lhsdesign(no_input, 1); % no_input x 2
% Scale to your parameter ranges
H_max = 10;
TR_min_List = lhs*0.9*H_max;        % scale to [0,30]
% Display results
OPT.parameters.slimit = 12;
for index = 1:no_input
    OPT.TR_min = TR_min_List(index);
    OPT.TR_max = OPT.TR_min  + 0.5;
    get_inputs();
    OPT.options.max_iter = 200;
    init_FE();
    init_optimization();
    %% Analysis
    perform_analysis();
    compute_stress_U();
    F = FE.elem_node'; % matrix of faces to be sent to patch function
    V = FE.coords'; % vertex list to be sent to patch function
    % ===================================================
    if index ==1
        X = zeros(no_input,FE.mesh_input.elements_per_side(2),FE.mesh_input.elements_per_side(1));
        Y = zeros(no_input,FE.mesh_input.elements_per_side(2),FE.mesh_input.elements_per_side(1));
        Y_Mises = zeros(no_input,FE.mesh_input.elements_per_side(2),FE.mesh_input.elements_per_side(1));

        A = zeros(FE.mesh_input.elements_per_side(2), FE.mesh_input.elements_per_side(1));
        x_c = FE.centroids(1, :) * FE.dim_scale ;
        y_c = FE.centroids(2, :) * FE.dim_scale ;
        rows = floor(y_c) + 1;   % y → row index
        cols = floor(x_c) + 1;   % x → column index
        idx = sub2ind(size(A), rows, cols);

    end
    X(index,idx) = FE.svm;

    % ------------------------------------
    %% Optimization
    runmma(OPT.dv, @(x)obj(x), @(x)nonlcon(x));
    % Combine them into one 3-channel array
    Y(index,idx) = OPT.pen_rho_e;
    Y_Mises(index,idx) = FE.svm;



end
%
% Save to .mat file for later Python training
save(fullfile(path_data,  'X.mat'), 'X');
save(fullfile(path_data, 'Y.mat'), 'Y');
save(fullfile(path_data,  'Y_Mises.mat'), 'Y_Mises');

% A(idx) = FE.svm;
% A = flipud(A);
% figure;
% imagesc(A);              % automatically scales colors
% axis equal tight;
% colormap(jet);           % use jet colormap
% colorbar;                % show color scale
% caxis([0, 12]);   % set color range from 0 to max_value
% axis off;
