function [J, dN_dx] = compute_jacobian(dN_dxi, node_coords)
    % 输入:
    %   dN_dxi     - 形函数对自然坐标的导数 (2x4)
    %   node_coords - 单元节点坐标 (4x2), 每行是 (x, y)
    %
    % 输出:
    %   J       - 雅可比矩阵 (2x2)
    %   dN_dx   - 形函数对物理坐标的导数 (2x4)
    
    % 计算雅可比矩阵
    J = dN_dxi * node_coords;
    
    % 雅可比矩阵的逆
    J_inv = inv(J);
    
    % 形函数对物理坐标的导数
    dN_dx = J_inv * dN_dxi;
end