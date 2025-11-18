function compute_stress_U()
%
% This function computes the mean compliance and its sensitivities
% based on the last finite element analysis
%
global FE OPT

% compute the element (centroidal) von Mises stress from the displacement.
% Note that this is not the FEA stress because Me is
% computed with the fully-solid modulus, not with the ersatz modulus.

se = zeros(FE.n_elem, FE.nloads);
% Initialize displacement at element centroids
strain = zeros(3, FE.n_elem, FE.nloads);  % for 2D problems: [εxx; εyy; γxy]

for iload=1:FE.nloads
    Ul = FE.U(:,iload);
    for j=1:FE.n_elem
        Ue = Ul(FE.edofMat(j,:));
        % While element matrices Me should be positive semidefinite,
        % numerically they can have very small but negative eigenvalues, hence the
        % von Mises stress can end up being a complex number (with a near-zero
        % complex part). Therefore, we no longer compute se as
        % sqrt(Ue'*Me*Ue), but we compute first the FE stress sigma, and
        % then the von Mises stress as se = sqrt(sigma'*V*sigma).
        B0e = FE.B0e(:,:,j);
        sigma = (FE.material.C*B0e*Ue);
        se(j,iload) = sqrt(sigma'*FE.V*sigma);
        % ==============
         % Compute strain at element centroid
        eps_e = B0e * Ue;                 % [εxx; εyy; γxy]
        
        strain(:, j, iload) = eps_e;      % store strain
    end
end

% Relaxed stress
p = OPT.parameters.penalization_param;
q = OPT.parameters.relaxation_param;
rhomin = FE.material.rho_min;
[re,dredrhof] = relaxdens(OPT.filt_rho_e,p,q,rhomin, 'stdpq');
FE.svm = re.*se;
FE.strain = strain;

end