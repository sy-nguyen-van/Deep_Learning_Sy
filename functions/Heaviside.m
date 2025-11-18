
function [rho_H, d_rho_H] = Heaviside(rho,beta,eta)
rho_H = (tanh(beta*eta)+tanh(beta*(rho-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta)));
d_rho_H =(beta*(sech(beta*(rho-eta))).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));
end


