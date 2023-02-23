clear;
close all;
load('./di_wsrdp3.mat');
load('./di2.mat');
load('./dixy2.mat');

n = 200;
%获取样本的中心点
Allcenters_data = [];

%wsrdp参数
a=1.03043553225395;
c=0.41365838905871;
m=-0 ;
theta=pi/(7);
dc=0.01;

data = [sample_cell{1, 1}(:,1), sample_cell{1, 1}(:,2), sample_cell{1, 1}(:,3)];
data = double(data)/255;
[idx,~]= kmeans(data,3);

sample_x = sample_xy{1, 1}(:,1);
sample_y = sample_xy{1, 1}(:,2);
%figure,gscatter(sample_x,sample_y,idx);
flag = 0;
Allcenters = [];
cluster_idx = cell(1, n);
for j=1:n
    idx = sample_cell_wsrdp{1,j};
    cluster_idx{1, j}=idx;
    data = [sample_cell{1, j}(:,1), sample_cell{1, j}(:,2), sample_cell{1, j}(:,3)];
    sample_x = sample_xy{1, j}(:,1);
    sample_y = sample_xy{1, j}(:,2);
    Coordinate = [sample_xy{1, j}(:,1),sample_xy{1, j}(:,2)];
    
    %figure,gscatter(sample_x, sample_y, idx);
    label = unique(cluster_idx{1, j});
    for i=1:length(unique(label))
        idx_data = data(find(cluster_idx{1, j} == i),:);                       % idx_data          不同类别的所有点的值
        Coordinate2 = Coordinate(find(cluster_idx{1, j} == i),:);              % Coordinate2       不同类别里所有点的位置
        Allcenters = [Allcenters;mean(Coordinate2)];             % Allcenters        中心点的位置
        Allcenters_data = [Allcenters_data;mean(idx_data)];      % Allcenters_data   中心点的值
        %{
        [len,~]= size(idxdata);
        if len == 1
            continue
            Allcenters_data = [Allcenters_data;idxdata];
            flag = flag + 1;
        else
            Allcenters_data = [Allcenters_data;mean(idxdata)];
        end
        %}
    end
end
%{
B = Allcenters_data(~isnan(Allcenters_data));
[len,~] = size(B);
di_centers = reshape(B,len/3,3);
%}
%di_centers = Allcenters;
%di_centers = Allcenters_data;
%temp = di_centers/10;
%[cl,Allcenters_idx] = WSRDP_Rotation_DPeaks(temp,dc,a,c,m,theta);
%[Allcenters_idx,~] = kmeans(di_centers,3);

min_v = 0;
max_v = 10;
accuracy = 0.01;
[x, y] = deal(min_v:accuracy:max_v, min_v:accuracy:max_v);
norm_Allcenters=mapminmax(Allcenters',min_v,max_v)';
% norm_Allcenters = Allcenters_data;

[Allcenters_idx,~] = kmeans(Allcenters,3);
figure()
gscatter(norm_Allcenters(:,1),norm_Allcenters(:,2), Allcenters_idx);
hold off

meankds = [];
COVS = [];
label = unique(Allcenters_idx);
len = length(unique(label));
for i=1:len
    labelname = label(i);
    Allcenters2 = norm_Allcenters(find(Allcenters_idx == labelname));
    meankd = mean(Allcenters2);
    meankds = [meankds;meankd];
    a =  Allcenters2;
    %a = 0;
    covkd = cov(a);
    COVS = [COVS;covkd];
end
[sp,product,result] = Sp(meankds,COVS);
fprintf('样本结果\n%d\n',result);

meankds = [];
%求高斯密度函数
label = unique(Allcenters_idx);
figure()
for i = 1:length(unique(label))
    labelname = label(i);
    if labelname == -1
        continue
    else
        kd = norm_Allcenters(find(Allcenters_idx == labelname),:);
        meankd = mean(kd);     %二维均值
        meankds = [meankds;meankd];
        % 绘制双变量正态分布概率密度曲线
%         x = 0 : 1 : 10;
%         y = 0 : 1 : 10;
%         x = 0 : 1 : 480;     % 140:1:280,40:1:220
%         y = 0 : 1 : 320 ;   
        [X,Y] = meshgrid(x, y);
        %y均值 方差 均方差
        U1 = meankd(1);
        DX = var(kd(:,1));
        dx = sqrt(DX);
        %y均值 方差 均方差
        U2 = meankd(2);
        DY = var(kd(:,2));
        dy = sqrt(DY);
        
        %协方差矩阵 2*2
        covkd = cov(kd(:,1), kd(:,2));
        r = covkd(2) / (dx * dy);
        % X Y的协方差的行列式 
        COV = det(covkd);
        %COVS = [COVS;covkd];
        %part1 1/(2*pi*协方差的行列式)
        part1 = 1 / ( 2 * pi * dx * dy * sqrt( 1 - r^ 2 ));
        %part1 = 1 / ( 2 * pi * COV);%e前面部分
        
        %part2 e上面部分 -1/2(x-均值)的转置（协方差矩阵的逆）(x-均值)
        %invkd = inv(covkd);%求逆 弃之采用建议的除法
        x_u = kd-meankd;
        %part2 =  -1 / (2*x_u.'/covkd*x_u);
        p1 = - 1 / ( 2 * ( 1 - r^ 2 ));
        px = (X - U1).^ 2. / DX;
        py = (Y - U2).^ 2. / DY;
        pxy = 2 * r.* (X - U1).* (Y - U2)./ (dx * dy);
        
        Z{1,i} = part1 * exp(p1 * (px - pxy + py));
        zlim([0 10]);
        mesh(X, Y, Z{1,i})
        hold on
       % xlabel('X1')
       % ylabel('X2')
        zlabel('Probability Density')
        shading interp 
        set(gca,'FontSize',20);
    end
end



figure()
AX = norm_Allcenters(:,1);
AY = norm_Allcenters(:,2);
XX = 2*min(AX)-max(AX):1:2*max(AX)-min(AX);
YY = 2*min(AY)-max(AY):1:2*max(AY)-min(AY);
figure()
set(gcf,'unit','normalized','position',[0.2,0.2,0.10,0.13])

[~, ~, XMesh, YMesh, ZMesh, colorList]=density2C(AX, AY, XX,YY);  %[~, ~, XMesh, YMesh, ZMesh, colorList]=density2C(AX, AY,XX,YY);
colormap(colorList)
contourf(XMesh, YMesh, ZMesh, 10)
xlim([0 10])
ylim([0 10])
%  xlim([150 300])
%  ylim([0 260])
set(gca,'FontSize',15);


