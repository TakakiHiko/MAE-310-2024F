clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 3;               % number of nodes in an element
n_el_x = 60;               % number of elements in x-dir
n_el_y = 60;               % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array
% IEN = zeros(n_el, n_en);
% for ex = 1 : n_el_x
%   for ey = 1 : n_el_y
%     ee = (ey-1) * n_el_x + ex; % element index
%     IEN(ee, 1) = (ey-1) * n_np_x + ex;
%     IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
%     IEN(ee, 3) =  ey    * n_np_x + ex + 1;
%     IEN(ee, 4) =  ey    * n_np_x + ex;
%   end
% end
% IEN array for triangular elements
IEN = zeros(2 * n_el, 3); % two triangles per quadrilateral element
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee_quad = (ey-1) * n_el_x + ex; % quadrilateral element index
    n1 = (ey-1) * n_np_x + ex;
    n2 = (ey-1) * n_np_x + ex + 1;
    n3 = ey * n_np_x + ex + 1;
    n4 = ey * n_np_x + ex;

    ee_tri1 = 2 * (ee_quad - 1) + 1; % first triangle
    ee_tri2 = 2 * (ee_quad - 1) + 2; % second triangle

    IEN(ee_tri1, :) = [n1, n2, n3];
    IEN(ee_tri2, :) = [n1, n3, n4];
  end
end
% ID array
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter;

LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
% 组装局部刚度矩阵和载荷向量
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en, n_en);
  f_ele = zeros(n_en, 1);
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
    end

    % Jacobian determinant for triangle
    detJ = abs((x_ele(2) - x_ele(1)) * (y_ele(3) - y_ele(1)) - ...
               (x_ele(3) - x_ele(1)) * (y_ele(2) - y_ele(1))) / 2;
    
    % Assemble local stiffness and load contributions
    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        k_ele(aa, bb) = k_ele(aa, bb) + weight(ll) * detJ * kappa * (Na_xi * Nb_xi + Na_eta * Nb_eta);
      end
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
    end
  end

  % Assemble into global system
  % (similar to your current assembly code)
end
 
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end


% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");

% EOF