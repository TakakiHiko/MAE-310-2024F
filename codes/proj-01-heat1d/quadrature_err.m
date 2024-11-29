
lmax=5;
mistake=zeros(1,lmax);
for n_int =1:1:lmax;

[xi, weight] = Gauss(n_int, -1, 1);

n = 6;%u积分为x^6/6;
exact = (1 - (-1)^(n+1)) / (n+1);
approxmation = 0;
for l = 1 : n_int
    approxmation = approxmation + weight(l) * xi(l)^n;
end

mistake(1,n_int) = abs(exact - approxmation);
end
plot(1:lmax,log(mistake));
xlabel("number of quadrature points");
ylabel("log(error)")
%可以见到在qudrature
%points的数量小于4时，误差随着点的数量增加而缓慢地减缓，直到点的数量大于4时，误差跳跃，说明逼近积分到达6阶精度。