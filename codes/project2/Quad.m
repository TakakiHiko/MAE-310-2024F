function N = Quad(aa, xi, eta)
    % Shape functions for triangular element
    if aa == 1
        N = 1 - xi - eta;
    elseif aa == 2
        N = xi;
    elseif aa == 3
        N = eta;
    end
end