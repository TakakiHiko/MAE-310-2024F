
filename = 'mesh.txt';
fileID = fopen(filename, 'r');
po=0.3;%泊松比
ym=1e9;%杨氏模量
f=1e4;%边界条件10kpa
% 初始化存储
nodes = [];
elements = [];
elementnodes=4;%rectangle elements
upperid=[];
leftid=[];
rightid=[];
lowerid=[];%准备储存边界节点编号
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
%坐标识别边界且记录其编号
 for i=1:numNodes
 cood=nodes(i,:);
  if cood(1,2)==-1
   lowerid=[lowerid,i];
  end
 if cood(1,2)==1
   upperid=[upperid,i];
 end
 if cood(1,3)==1
  rightid=[rightid,i];
 end
 if cood(1,3)==-1
  leftid=[leftid,i];
  end
 end
      


IEN=[];
lineflag=1;
for i=1:numElements
vec=cell2mat(elements(i,1));
if vec(1,2)==3
tips=vec(1,3);
nodess=vec(1,4+tips:length(vec));
IEN(lineflag,1)=lineflag;
IEN(lineflag,2:1+length(nodess))= nodess;
lineflag=lineflag+1;
end
end
IEN=IEN(1:length(IEN),2:5);%generate IEN matrix
%找出边界单元编号
force_element=zeros(length(rightid)-1,3);

numElements=lineflag-1;
ks=1;
for i=1:numElements
ar=IEN(i,:);
   flag=0;
 for j=1:length(rightid)
   if ismember(rightid(1,j),ar)
       flag=flag+1;
   end
 end
 if flag==2
force_element(ks,1)=i;
flag=0;
ks=ks+1;
 elseif flag==1
     flag=0;
 end
end
for i=1:length(force_element)
    k=2;
for j=1:4
if ismember(IEN(force_element(i,1),j),rightid)
force_element(i,k)=IEN(force_element(i,1),j);
k=k+1;
end
end
end
freedom=2;%节点自由度

ID=generate_ID_from_IEN(IEN,freedom);

fclose(fileID);

%local stiffness matrix
sizee=size(IEN);
Kee1=zeros(8,8,lineflag-1);
%Local force matrix
F=zeros(numNodes*freedom,1); 
elemntflag=1;
%无体积力，仅有表面力
F_local=zeros(8,numElements);
Gausspoint=6;
[xi_gauss,w_gauss]=Gauss(Gausspoint,-1,1);
example=0;
for i=1:length(force_element)
f_local=zeros(4,1);
load=[f;0];
elementindex=force_element(i,1);
coo1=[nodes(force_element(i,2),2),nodes(force_element(i,2),3)];
coo2=[nodes(force_element(i,3),2),nodes(force_element(i,3),3)];
L=norm(coo1-coo2);
Jaco=L/2;
for j = 1:length(xi_gauss)
    % 当前高斯点局部坐标
    xi = xi_gauss(j);
    
    % 形状函数 (线性单元)
    N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
    
    % 将形状函数扩展到自由度 (二维问题)
    N_ext = kron(N, eye(2)); % 生成 [N1 0 N2 0; 0 N1 0 N2]
    
    % 累加表面力贡献
    f_local = f_local + N_ext' * load * w_gauss(j) * Jaco;
end
F_local(1:4,elementindex)=f_local;
end


for i=1:lineflag-1
coo=zeros(4,2);
for j=1:4
    index=IEN(i,j);
coo(j,1)=nodes(index,2);
coo(j,2)=nodes(index,3);
end
Kee1(:,:,i)=quad4_stiffness(ym,po,coo);

end
k_global=zeros(numNodes*freedom,numNodes*freedom);
for e=1:numElements
elementnodes=IEN(e,:);
freedommap=zeros(1,numel(elementnodes)*freedom);
 for i = 1:numel(elementnodes)
        freedommap(2*i-1:2*i) = [(elementnodes(i)-1)*freedom + 1, ...
                             (elementnodes(i)-1)*freedom + 2];
 end

 Ke=Kee1(:,:,e);
 Fe=F_local(:,e);
 for i=1:length(freedommap)
    for j=1:length(freedommap)
     k_global(freedommap(i), freedommap(j)) = k_global(freedommap(i), freedommap(j)) + Ke(i, j);
    end
    F(freedommap(i),1)=F(freedommap(i),1)+Fe(i,1);
 end
end
%有对称性可得，对于左侧存在x位移位移为0，下侧y位移为0
%将对应的自由度以及刚度矩阵调0

% 
% p1=leftid.*2-1;
% p2=lowerid.*2-1;
% p=[p1 p2];
% k_global(p, :) = [];
% k_global(:, p) = [];
% F(p,:)=[];
U=F\k_global;

function K = quad4_stiffness(E, nu, coords)
    % 输入参数:
    % E: 弹性模量
    % nu: 泊松比
    % coords: 4x2 矩阵，节点坐标 [x1, y1; x2, y2; x3, y3; x4, y4]

    % 构建弹性矩阵 D (平面应力假设)
    D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    
    % 高斯积分点及权重 
    Gausspoints=6;
    [gauss_pts,weights] = Gauss(Gausspoints,-1,1);
    
    
    % 初始化刚度矩阵
    K = zeros(8, 8);  % 4节点，每个节点2个自由度
    
    % 循环计算高斯积分
    for i = 1:Gausspoints
        for j = 1:Gausspoints
            % 当前积分点 (ξ, η)
            xi = gauss_pts(i);
            eta = gauss_pts(j);
            
            % 计算形函数及其导数
            [N, dN_dxi] = shape_function_quad4(xi, eta);
            
            % Jacobian 矩阵
            J = dN_dxi * coords;
            
            % Jacobian 行列式和逆
            detJ = det(J);
            invJ = inv(J);
            
            % 计算形函数对 x, y 的导数
            dN_dxy = invJ * dN_dxi;
            
            % 构造应变-位移矩阵 B
            B = zeros(3, 8);
            for k = 1:4
                B(1, 2*k-1) = dN_dxy(1, k);
                B(2, 2*k) = dN_dxy(2, k);
                B(3, 2*k-1) = dN_dxy(2, k);
                B(3, 2*k) = dN_dxy(1, k);
            end
            
            % 高斯点的权重
            weight = weights(i) * weights(j);
            
            % 累积刚度矩阵
            K = K + B' * D * B * detJ * weight;
        end
    end
end

function [N, dN_dxi] = shape_function_quad4(xi, eta)
    % 四节点四边形单元形函数及其导数
    N = 0.25 * [(1 - xi)*(1 - eta);
                (1 + xi)*(1 - eta);
                (1 + xi)*(1 + eta);
                (1 - xi)*(1 + eta)];
    
    dN_dxi = 0.25 * [-1 + eta, -1 + xi;
                      1 - eta, -1 - xi;
                      1 + eta,  1 + xi;
                     -1 - eta,  1 - xi]';
end

