function [N, dN_dxi] = shape(xi, eta)
    % 输入:
    %   xi  - 自然坐标 \xi
    %   eta - 自然坐标 \eta
    %
    % 输出:
    %   N       - 形函数值 (1x4)
    %   dN_dxi  - 形函数梯度 (2x4), 第一行是对xi的导数，第二行是对eta的导数
    
    % 形函数
    N = [
        0.25 * (1 - xi) * (1 - eta);  % N1
        0.25 * (1 + xi) * (1 - eta);  % N2
        0.25 * (1 + xi) * (1 + eta);  % N3
        0.25 * (1 - xi) * (1 + eta)   % N4
    ];
    
    % 形函数对自然坐标的导数
    dN_dxi = [
        -0.25 * (1 - eta),  0.25 * (1 - eta),  0.25 * (1 + eta), -0.25 * (1 + eta);  % dN/dxi
        -0.25 * (1 - xi),  -0.25 * (1 + xi),   0.25 * (1 + xi),  0.25 * (1 - xi)    % dN/deta
    ];
end