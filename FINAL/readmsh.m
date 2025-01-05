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

