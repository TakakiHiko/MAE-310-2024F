R = 0.3;// 给定了圆的半径
L = 1.0;// 给定正方形的边长

Point(1) = {L, -L, 0};
Point(2) = {L, L, 0};
Point(3) = {-L, L, 0};
Point(4) = {-L, -L, 0};// 这四行给定矩形的四个顶点
Point(5) = {-L + R, -L, 0};
Point(6) = {-L, -L + R, 0};//定义了两个四分之一圆的起点
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};//利用了三角函数给出四分之一圆弧的中点

Circle(1) = {5, 4, 7};//以点4为圆心，连接点5和7的圆弧
Circle(2) = {7, 4, 6};//上同

Line(3) = {6, 3};//定义线，连接6和3，下同
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};
Line(7) = {2, 7};

Curve Loop(1) = {4, 7, 2, 3};//定义闭合曲线
Plane Surface(1) = {1};//根据闭合曲线来定义面

Curve Loop(2) = {7, -1, -6, -5};
Plane Surface(2) = {2};//上同

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3;、=//进行有限元映射，准备生成网络，在每条线段上生成3个节点

Transfinite Surface{1};
Transfinite Surface{2};//使用结构化网格

Recombine Surface{1};
Recombine Surface{2};//转化为四边形网格


Mesh.ElementOrder = 1;//线性插值
Mesh.Algorithm = 8;//算法编号

// EOF
