

ux1=0;
u0=1;
u=@(x) x.^5;%exact function
calu=1/11; %u平方在0到1的积分
calux=16/9;%u对x的偏导的平方在0到1的积分值

ux=@(x) 5*x^4;%x的一次导数
uxx=@(x) 20*x^3;%对x的二次导数

nelmax=16;%giving the biggest number of elements
n_node=3;%number of nodes in each elements
n_int=16;%number of point for gaussquadrature
degree=2;%number of polyshape oreder
pp=degree;
mistake=zeros(2,nelmax/2-1);
[x,w]=Gauss(n_int,-1,1);
for n_el=2:2:nelmax
%     f = @(x) -20*x.^3; % f(x) is the source
% g = 1.0;           % u    = g  at x = 1
% h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
n_en = pp + 1;         % number of element or local nodes
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations


hh = 1.0 / (n_np - 1); % space between two adjacent nodes
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

IEN = zeros(n_el, n_en);

for ee = 1 : n_el
  for flag1 = 1 : n_en
    IEN(ee, flag1) = (ee - 1) * pp + flag1;
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% allocate the stiffness matrix
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
F = zeros(n_eq, 1);

% Assembly of the stiffness matrix and load vector
for ee = 1 : n_el
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
  f_ele = zeros(n_en, 1);    % allocate a zero element load vector

  x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

  % quadrature loop
  for qua = 1 : n_int    
    dx_dxi = 0.0;
    x_l = 0.0;
    for flag1 = 1 : n_en
      x_l    = x_l    + x_ele(flag1) * PolyShape(pp, flag1, xi(qua), 0);
      dx_dxi = dx_dxi + x_ele(flag1) * PolyShape(pp, flag1, xi(qua), 1);
    end
    dxi_dx = 1.0 / dx_dxi;

    for flag1 = 1 : n_en
      f_ele(flag1) = f_ele(flag1) + weight(qua) * PolyShape(pp, flag1, xi(qua), 0) * u(x_l) * dx_dxi;
      for bb = 1 : n_en
        k_ele(flag1, bb) = k_ele(flag1, bb) + weight(qua) * PolyShape(pp, flag1, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
      end
    end
  end
 
  % Assembly of the matrix and vector based on the ID or LM data
  for flag1 = 1 : n_en
    P = ID(IEN(ee,flag1));
    if(P > 0)
      F(P) = F(P) + f_ele(flag1);
      for bb = 1 : n_en
        Q = ID(IEN(ee,bb));
        if(Q > 0)
          K(P, Q) = K(P, Q) + k_ele(flag1, bb);
        else
          F(P) = F(P) - k_ele(flag1, bb) * u0; % handles the Dirichlet boundary data
        end
      end
    end
  end
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) +ux1;

% Solve Kd = F equation
d_temp = K \ F;

uh=d_temp;
d=[uh;u0];

err1=0;
err2=0;

for flag_el=1:n_el
    x_ele=zeros(n_en,1);
    for flag=1:n_en
        x_ele(flag)=x_coor(IEN(flag_el,flag));
    end

    L1=0;
    H2=0;

   
        for l = 1 : n_int
            u_h = 0.0;
            uhx = 0.0;
            x_l = 0.0;
            changer = 0.0;
            
            for flag1 = 1 : n_en
                changer = changer + x_ele(flag1) * PolyShape(degree,flag1, xi(l), 1);
            end
            
            for flag1 = 1 : n_en
               u_h = u_h +u(x_ele(flag1)) * PolyShape(degree,flag1, xi(l), 0); 
               
                uhx = uhx + u(x_ele(flag1)) * PolyShape(degree,flag1, xi(l), 1) / changer;
      
                x_l = x_l + x_ele(flag1) * PolyShape(degree,flag1, xi(l), 0);
                 
            end
            ureal=u(x_l);
            uxreal=ux(x_l);
            L1 = L1 + w(l) * (u_h - ureal)^2 * changer;
             H2 = H2 + w(l) * (uhx - uxreal)^2 * changer;

        end
        err1=err1+L1;
        err2=err2+H2;

end
f_L1=(err1)^0.5/calu;
f_H2=(err2)^0.5/calux;
mistake(1,(n_el/2))=f_L1;
mistake(2,(n_el/2))=f_H2;
end
mistake=log(mistake);
mesh_number=4:4:n_eq;
mesh_number=log(mesh_number);
figure;
subplot(2,1,1);
plot(mesh_number,mistake(1,:),"-*");
 xlabel('log(mesh_size)');
    ylabel('log(Error u(x))');
    title('Plot of Error  vs. mesh size');

 subplot(2,1,2);
plot(mesh_number,mistake(2,:),"-o");
 xlabel('log(mesh_size)');
    ylabel('log(Error ux(x))');
    title('Plot of Error  vs. mesh size');
k1=polyfit(mesh_number,mistake(1,:),1)
k2=polyfit(mesh_number,mistake(2,:),1)
%展示线性相关系数
r1=corr(mesh_number',mistake(1,:)')
r2=corr(mesh_number',mistake(2,:)')
%可见对于u（x）计算相对误差之时，其误差随着mesh—size的减小而呈3阶逼近。而ux（x）2阶。
%而且其线性拟合的关系较好。可以认为分析准确


