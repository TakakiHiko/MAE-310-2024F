function ID = generate_ID_from_IEN(IEN, n_dof_per_node)
    % 输入:
    %   IEN            - 节点连接矩阵 (n_element x n_node_per_element)
    %   n_dof_per_node - 每个节点的自由度数量
 
    % 输出:
    %   ID - 自由度编号矩阵 (n_element x (n_node_per_element * n_dof_per_node))

    [n_element, n_node_per_element] = size(IEN);
    ID = zeros(n_element, n_node_per_element * n_dof_per_node);

    for i = 1:n_element
        for j = 1:n_node_per_element
            % 当前单元的第 j 个节点
            global_node = IEN(i, j);
            for k = 1:n_dof_per_node
                % 对应的全局自由度
                global_dof = (global_node - 1) * n_dof_per_node + k;
                % 填入 ID 矩阵
                ID(i, (j - 1) * n_dof_per_node + k) = global_dof;
            end
        end
    end
end