filename = 'mesh.txt'; % 替换为你的文件名
fid = fopen(filename, 'r');

% 初始化变量
elements = [];
readingElements = false;

% 遍历文件逐行读取
while ~feof(fid)
    line = strtrim(fgetl(fid)); % 读取一行并去掉空格
    if strcmp(line, '$Elements')
        readingElements = true; % 进入 Elements 部分
        numElements = str2double(fgetl(fid)); % 下一行是元素数量
        continue;
    elseif strcmp(line, '$EndElements')
        readingElements = false; % 退出 Elements 部分
        break;
    end
    
    if readingElements
        elementData = str2double(line); 
        elements = [elements; elementData]; 
    end
end

% 关闭文件
fclose(fid);

% 提取有用信息
elementIDs = elements(:, 1); % 单元编号
elementTypes = elements(:, 2); % 单元类型
elementTags = elements(:, 4); % 标签信息
elementNodes = elements(:, 5:end); % 节点编号（IEN矩阵核心部分）

triangleElements = elementTypes == 3; % 只保留三角形单元
IEN_matrix= elementNodes(triangleElements, :); % 对应节点关系

% 打印结果
disp('IEN 矩阵:');
disp(IEN_matrix);