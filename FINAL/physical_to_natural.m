function [xi, eta] = physical_to_natural(x, y, coords)
    % 输入：
    % x, y - 物理坐标
    % coords - 节点坐标 (4x2 矩阵)，每行是 [x_i, y_i]
    % 输出：
    % xi, eta - 自然坐标

    % 初始猜测
    xi = 0; eta = 0;
    tol = 1e-6; % 收敛阈值
    max_iter = 50; % 最大迭代次数
    
    for iter = 1:max_iter
        % 计算形函数及其导数
        [N, dN_dxi] = quad4_shape_functions(xi, eta);
        
        % 计算物理坐标和雅可比矩阵
        x_phys = sum(N .* coords(:, 1));
        y_phys = sum(N .* coords(:, 2));
        
        J = dN_dxi * coords; % 雅可比矩阵
        detJ = det(J);
        if abs(detJ) < 1e-12
            error('雅可比矩阵行列式接近零，可能无法收敛。');
        end
        
        % 计算 f 和修正
        f = [x - x_phys; y - y_phys];
        delta = J \ f;
        
        % 更新自然坐标
        xi = xi + delta(1);
        eta = eta + delta(2);
        
        % 检查收敛
        if norm(delta) < tol
            return;
        end
    end
    error('牛顿迭代未收敛，请检查初始猜测或单元形状。');
end

function [N, dN_dxi] = quad4_shape_functions(xi, eta)
    % 形函数
    N = [
        0.25 * (1 - xi) * (1 - eta);  % N1
        0.25 * (1 + xi) * (1 - eta);  % N2
        0.25 * (1 + xi) * (1 + eta);  % N3
        0.25 * (1 - xi) * (1 + eta)   % N4
    ].';

    % 形函数导数
    dN_dxi = [
        -0.25 * (1 - eta),  0.25 * (1 - eta),  0.25 * (1 + eta), -0.25 * (1 + eta); % dN/dxi
        -0.25 * (1 - xi),  -0.25 * (1 + xi),   0.25 * (1 + xi),  0.25 * (1 - xi)   % dN/deta
    ];
end