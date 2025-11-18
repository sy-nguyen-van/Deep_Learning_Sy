function [] = perform_analysis()
%
% Filter and penalize densities, solve the finite
% element problem for the displacements and reaction forces, and then
% evaluate the relevant functions.
%
global OPT FE

% Filter densities
dv_macro = OPT.H * OPT.dv;

    coeff = OPT.beta_min + OPT.iter*(OPT.beta_max-OPT.beta_min)/OPT.options.max_iter;
    [OPT.filt_rho_e, OPT.diff_Hea] = Heaviside(dv_macro,coeff,OPT.eta_H);


[OPT.pen_rho_e, OPT.dpen_rho_e] = penalize(OPT.filt_rho_e, ...
    OPT.parameters.penalization_param, ...
    OPT.parameters.penalization_scheme, ...
    FE.material.rho_min);


% Perform FE analysis
FE_analysis();

% Evaluate objective and constraint functions
evaluate_relevant_functions();

end

