function [N_xi, N_eta] = Quad_grad(aa, xi, eta)
    % Gradients of shape functions for triangular element
    if aa == 1
        N_xi = -1; N_eta = -1;
    elseif aa == 2
        N_xi = 1; N_eta = 0;
    elseif aa == 3
        N_xi = 0; N_eta = 1;
    end
end