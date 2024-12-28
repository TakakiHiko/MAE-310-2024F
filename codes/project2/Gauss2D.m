function [xi, eta, weight] = Gauss2D(n_int)
    % Gaussian quadrature points for a triangle
    if n_int == 1
        % Single integration point
        xi = 1/3; eta = 1/3; weight = 1/2;
    elseif n_int == 3
        % Three-point integration
        xi = [1/6, 2/3, 1/6];
        eta = [1/6, 1/6, 2/3];
        weight = [1/6, 1/6, 1/6];
    else
        error('Unsupported number of integration points for triangles.');
    end
end