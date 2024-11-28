clear;
clc;
hold on;
ux1=0;
u0=1;
u=@(x) x.^5;%exact function
calu=1/11; %u平方在0到1的积分
calux=16/9;%u对x的偏导的平方在0到1的积分值

ux=@(x) 5*x^4;%x的一次导数
uxx=@(x) 20*x^3;%对x的二次导数

nelmax=16;%giving the biggest number of elements
n_node=3;%number of nodes in each elements
n_int=15;%number of point for gaussquadrature
for n_el=2:1:32% giving the mesh size
hh = 1 / n_el / n_node;
    x_coor = 0 : hh : 1;

 %generate Id
     n_point=n_el*n_node+1;
     ID = 1 : n_point;
    ID(end) = 0;
    n_eq = n_point- 1;

    %generate IEN
    IEN = zeros(4, n_el);
     for ee = 1 : n_el
        IEN(1,ee) = 3*ee - 2;
        IEN(2,ee) = 3*ee - 1;
        IEN(3,ee) = 3*ee;
        IEN(4,ee) = 3*ee + 1;
     end

    %GAUSS Quadrature
[x,w]=Gauss(n_int,-1,1);

%generate K and F
K=zeros(n_eq,n_eq);
F=zeros(n_eq,1);
   for ee = 1 : n_el

        k_e = zeros(4,4); 
        f_e = zeros(4,1);

        x_ele = zeros(4,1);
        for aa = 1 : 4
            x_ele(aa) = x_coor(IEN(aa,ee)); % A = IEN(a,e)
        end

    for l = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        for aa = 1 : 4
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(aa, xi(l), 1);
            x_l = x_l + x_ele(aa) * PolyShape(aa, xi(l), 0);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : 4
            for bb = 1 : 4
                k_e(aa,bb) = k_e(aa,bb) + weight(l) * PolyShape(aa, xi(l), 1) * PolyShape(bb, xi(l), 1) * dxi_dx;
            end
        end

        for aa = 1 : 4
            f_e(aa) = f_e(aa) + weight(l) * PolyShape(aa, xi(l), 0) * f(x_l) * dx_dxi;
        end

    end

    % Now we need to put element k and f into global K and F
    for aa = 1 : 4
        AA = IEN(aa,ee);
        PP = ID(AA);
        if PP > 0
            F(PP) = F(PP) + f_e(aa);
            for bb = 1 : 4
                BB = IEN(bb,ee);
                QQ = ID(BB);
                if QQ > 0
                    K(PP,QQ) = K(PP,QQ) + k_e(aa,bb);
                else
                    F(PP) = F(PP) - k_e(aa,bb) * g;
                end
            end
        end
    end

        if ee == 1
        F(ID(IEN(1,ee))) = F(ID(IEN(1,ee))) + h;
        end
    end
% Now we have K and F
% Solve Kd = F
uh = K \ F;
d = [uh; g];
%tips: all the code depend on course material,below part calculates the
%error

end
