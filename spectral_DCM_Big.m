clear;
close all;
%加载总体数据
dataall=load('ConfLongDemo.txt');
% dataall=csvread('data(1).csv'); 
% ori_label=csvread('ori_label.csv');

%布局
% min_v = 0;
% max_v = 10;
% accuracy = 0.01;
% %数据规范化
% data=mapminmax(dataall', min_v, max_v)';
%坐标规范化
%[x, y] = deal(min_v:accuracy:max_v, min_v:accuracy:max_v);

%集中式聚类
%参数，需要调参
% sigma=3;
% k=5;
% X = dataall;
% [row,col] = size(X);
%计算临接矩阵W,采用高斯核函数RBF
% W = zeros(row,row);
% index_all = zeros(row,row);
% for i = 1:row                   % 全连接构造相似度矩阵W
%     for j=1:row
%         if i ~= j
%             W(i,j) = exp((-sum((X(i,:)-X(j,:)).^2))/(2*sigma.^2));
%         end
%     end
% end
% %tic
% cluster_idx = spectral_clustering(W,k);
% %toc
% %作图
% figure(),gscatter(data(:,1), data(:,2), cluster_idx);


n = 500;
sample = 500;
%布局
min_v = 0;
max_v = 10;
accuracy = 0.01;
data=mapminmax(dataall',min_v,max_v)';
[x, y] = deal(min_v:accuracy:max_v, min_v:accuracy:max_v);

%样本空间di
ori_sample_cell = cell(1, n);
norm_sample_cell = cell(1, n);
%存放每个点的分类情况
cluster_idx = cell(1, n);     % cluster_idx  存放聚类标签
spectral_idx = cell(1, n);

for i = 1:n
    ori_sample_cell{1, i} = datasample(dataall,sample); %数据随机抽取样本
    norm_sample_cell{1, i} = mapminmax(ori_sample_cell{1, i}',min_v,max_v)';
end
%中心点矩阵
ori_Allcenters =[];
norm_Allcenters =[];

%对样本的spectral聚类
for i = 1:n
    m=i;
    ori_sample = ori_sample_cell{1, i};
    norm_sample = norm_sample_cell{1, i};
    %调参
    sigma=4;   %sigma=3;   12
    k=7;
    X = ori_sample;
    [row,col] = size(X);
    %计算临接矩阵W,采用高斯核函数RBF
    W = zeros(row,row);
    index_all = zeros(row,row);
    for i = 1:row                   % 全连接构造相似度矩阵W
        for j=1:row
            if i ~= j
                W(i,j) = exp((-sum((X(i,:)-X(j,:)).^2))/(2*sigma.^2));
            end
        end
    end
    
    %tic
    cluster_idx{1, i} = spectral_clustering(W,k);
    
    %sample的聚类标签 ：spectral_idx  (original sample)
    spectral_idx1 = cluster_idx{1,sample};
    [spectral_idx{1, m},~] = spectral_idx1;
    %toc
   
    %打印
%      figure(),drawshapes3d(ori_sample,cluster_idx{1,i},sample);
%      hold on
    
    label = unique(cluster_idx{1, i});
    length(label);
    %取每个类的点
    norm_dikind = [];
    ori_dikind = [];
    
    for j = 1:length(label)
        labelname = label(j);
        %去噪
        if labelname == -1
            continue
        else
            norm_dikind = [norm_sample(find(cluster_idx{1, i} == labelname),:)];
            ori_dikind =  [ori_sample(find(cluster_idx{1, i} == labelname),:)];
            C = mean(norm_dikind, 1);
            C1 = mean(ori_dikind, 1);
            %打印中心点
            %plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
            norm_Allcenters = [norm_Allcenters; C];
            ori_Allcenters = [ori_Allcenters; C1];
        end
    end
end

%对中心点聚类
%cluster_idx = dbscan(norm_Allcenters, epsilonCenter, minptsCenter);%%采用dbscan聚类(,范围,点数)
%调参
[norm_Allcenters_idx,~] =  kmeans(norm_Allcenters,k);
[ori_Allcenters_idx,~] = kmeans (ori_Allcenters,k);

figure(),drawshapes3d(norm_Allcenters,norm_Allcenters_idx,3500);
figure(),drawshapes3d(ori_Allcenters,ori_Allcenters_idx,3500);
%[clustCent,cluster_idx,clustMembsCell] = MeanShiftCluster(norm_Allcenters',1.5);
%打印中心点聚类结果
% figure(),gscatter(norm_Allcenters(:,1),norm_Allcenters(:,2), norm_Allcenters_idx);
% figure(),gscatter(ori_Allcenters(:,1),ori_Allcenters(:,2), ori_Allcenters_idx);
%hold off

label = unique(norm_Allcenters_idx);
meankds = [];
COVS = [];
for i = 1:length(unique(label))
    labelname = label(i);
    if labelname == -1
        continue
    end
    
    kd = norm_Allcenters(find(norm_Allcenters_idx == labelname),:);
    check = size(kd);  
    %meankd = mean(kd);    %二维均值
    
    if check ~= 1
        meankd = mean(kd);    %二维均值
    else
        continue
    end
    meankds = [meankds;meankd];
    covkd = cov(kd(:,1), kd(:,2));
    COVS = [COVS;covkd];
end
fprintf('计算\n');
[sp,product,result] = Sp(meankds,COVS);
fprintf('样本结果\n %d\n',result);

%求高斯密度函数
label = unique(norm_Allcenters_idx);
figure()
for i = 1:length(unique(label))
    labelname = label(i);
    if labelname == -1
        continue
    else
        kd = norm_Allcenters(find(norm_Allcenters_idx == labelname),:);
        meankd = mean(kd);     %二维均值
        meankds = [meankds;meankd];
        % 绘制双变量正态分布概率密度曲线
        %x = -0.5 : 0.001 : 0.5 ;  
        %y = -0.5 : 0.001 : 0.5 ;      
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
        %zlim([0 5]);
        mesh(X, Y, Z{1,i})
        hold on
       % xlabel('X1')
       % ylabel('X2')
        zlabel('Probability Density')
        shading interp 
        set(gca,'FontSize',20);
    end
end

%打印中心点的等高线图
figure()
[~, ~, XMesh, YMesh, ZMesh, colorList]=density2C(norm_Allcenters(:,1), norm_Allcenters(:,2), x, y);
colormap(colorList)
contourf(XMesh, YMesh, ZMesh, 10)
set(gca,'FontSize',20);