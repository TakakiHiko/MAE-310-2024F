filename = 'mesh.txt';
fileID = fopen(filename, 'r');

% 初始化存储
nodes = [];
elements = [];

% 逐行读取文件
while ~feof(fileID)
    line = fgetl(fileID);
    if contains(line, '$ParametricNodes')
        % 读取节点部分
        numNodes = str2double(fgetl(fileID));
        nodes = zeros(numNodes, 4); % 假设节点有4列数据
        for i = 1:numNodes
            data = sscanf(fgetl(fileID), '%f')';
            nodes(i, :) = data(1:4); % 取前4列
        end
    elseif contains(line, '$Elements')
        % 读取单元部分
        numElements = str2double(fgetl(fileID));

        elements = cell(numElements, 1);
        for i = 1:numElements
            data = sscanf(fgetl(fileID), '%f')';
            elements{i} = data; % 保存整个单元定义
        end
    end
end

fclose(fileID);
IEN=[];
lineflag=1;
for i=1:numElements
vec=cell2mat(elements(i,1));
if vec(1,2)==3
tips=vec(1,3);
nodes=vec(1,4+tips:length(vec));
IEN(lineflag,1)=lineflag;
IEN(lineflag,2:1+length(nodes))= nodes;
lineflag=lineflag+1;
end
end
IEN=IEN(1:length(IEN),2:5);%generate IEN matrix


freedom=2;%节点自由度

ID=generate_ID_from_IEN(IEN,freedom);

function ID = generate_ID_from_IEN(IEN, n_dof_per_node)
    % 输入:
    %   IEN            - 节点连接矩阵 (n_element x n_node_per_element)
    %   n_dof_per_node - 每个节点的自由度数量
    %
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


