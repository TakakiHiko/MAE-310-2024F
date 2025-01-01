function [nodes, elements] = readmsh(filename)
    % 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件');
    end

    % 初始化变量
    nodes = [];
    elements = [];

    while ~feof(fid)
        line = fgetl(fid);
        
        if contains(line, '$Nodes')
            % 读取节点数据
            numNodes = str2double(fgetl(fid));  % 节点数目
            nodes = zeros(numNodes, 4);
            for i = 1:numNodes
                nodes(i, :) = fscanf(fid, '%f', 4);
            end
        elseif contains(line, '$Elements')
            % 读取元素数据
            numElements = str2double(fgetl(fid));  % 元素数目
            elements = zeros(numElements, 5);  % 假设每个元素有5个数据项
            for i = 1:numElements
                elements(i, :) = fscanf(fid, '%f', 5);
            end
        end
    end
    
    % 关闭文件
    fclose(fid);
end

