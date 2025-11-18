clear all; close all; clc;
%% source folders containing scripts not in this folder
addpath(genpath('FE_routines'));addpath(genpath('functions'));
addpath(genpath('mesh_utilities'));
addpath(genpath('optimization'));
addpath(genpath('utilities'));
addpath(genpath('plotting'));

global OPT FE

path_data = 'data/train/';

OPT.beta_min = 0;
OPT.beta_max = 30;
OPT.iter = 1;
OPT.eta_H = 0.5;
%% Initialization
% ============================
no_input = 200; % number of samples
% Latin Hypercube Sampling in [0,1]
lhs = lhsdesign(no_input, 1); % no_input x 2
% Scale to your parameter ranges
H_max = 10;
TR_min_List = lhs*0.9*H_max;        % scale to [0,30]
% Display results
OPT.parameters.slimit = 12;
OPT.TR_min = TR_min_List(1);
OPT.TR_max = OPT.TR_min  + 0.5;
get_inputs();
OPT.options.max_iter = 200;
init_FE();
init_optimization();
%% Analysis
perform_analysis();


% ------------------------------------
%% Optimization
t1 = tic;   % start timer

runmma(OPT.dv, @(x)obj(x), @(x)nonlcon(x));

time = toc(t1)   % end timer and store elapsed time
