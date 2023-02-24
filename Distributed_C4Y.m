clear;
close all;

% loading overall data
% load('7.mat','points');
% dataall = points;
dataall=load('synthesis_3.txt');
%dataall=load('synthesis_2.txt');
%dataall=load('unbalance2.txt');

min_v = 0;
max_v = 10;
accuracy = 0.01;
% data Normalization
data=mapminmax(dataall',min_v,max_v)';
% coordinate Normalization
%[x, y] = deal(min_v:accuracy:max_v, min_v:accuracy:max_v);

% parameters
n = 10;                                        % number of samples
sample = 900;                                  % the number of data points contained in each sample
k=5;                                           % number of clustering categories

[epsilonD,minptsD]=deal( 0.040 , 15 );         
[epsilonDi,minptsDi]=deal( 0.045 , 14 );       % DBSCAN parameters
%d=0.1                                         % DPeak parameters
%bandwidth1=0.15;                                 
%bandwidth2=0.178;                             % meanshift parameters, bandwidth1:overall data, bandwidth2: samples

% spectral clustering parameters
%sigma=3;
%X = dataall;
%[row,col] = size(X);
%计算临接矩阵W,采用高斯核函数RBF
%W = zeros(row,row);
%index_all = zeros(row,row);
%for i = 1:row                   % 全连接构造相似度矩阵W
%    for j=1:row
%        if i ~= j
%            W(i,j) = exp((-sum((X(i,:)-X(j,:)).^2))/(2*sigma.^2));
%        end
%    end
%end                                        

figure(),gscatter(dataall(:,1), dataall(:,2));

% overall data clustering results, cluster_idx：clustering labels
dione = dataall;
cluster_idx{1, 1} = dbscan(dione, epsilonD, minptsD);           %DBSCAN
%[cluster_idx{1, i} ,~]= kmeans(dataall, k);                    %kmeans
%[cluster_idx,p,s,K,C,Klist] = Dpeak(dataall, d, k);            %DPeak
%[clustCent,cluster_idx,clustMembsCell] = MeanShiftCluster(dataall',bandwidth1);   %meanshift
%cluster_idx = spectral_clustering(W,k);                        %spectral clustering
figure(),gscatter(data(:,1), data(:,2), cluster_idx);


ori_sample_cell = cell(1, n);                                      %original sample sets that stors all original samples   
norm_sample_cell = cell(1, n);                                     %normalized samle sets, where all elements of samples \in [0,10]  
for i = 1:n
    ori_sample_cell{1, i} = datasample(dataall, sample);           %random sample of data
    norm_sample_cell{1, i} = mapminmax(ori_sample_cell{1, i}',min_v,max_v)';
end

cluster_idx = cell(1, n);                  %clustering labels
norm_Allcenters = [];                      %stors all clustering centroids
ori_Allcenters = [];                       %all centeroids in [0,10]

%clustering of samples
for i = 1:n
    ori_sample = ori_sample_cell{1, i};                                       %fetch the ith original sample  
    norm_sample = norm_sample_cell{1, i};                                     %fetch the ith normalized sample
   
   
    cluster_idx{1,i}=dbscan(ori_sample,epsilonDi,minptsDi);                  % DBSCAN
    %[cluster_idx{1, i} ,~]= kmeans(ori_sample, k);                          % kmeans  
    %[cluster_idx{1, i},p,s,K,C,Klist] = Dpeak(ori_sample, d, k);            % Dpeak
    %[clustCent,cluster_idx{1,i},clustMembsCell] = MeanShiftCluster(ori_sample',bandwidth2);   %meanshift
    
    %X = ori_sample;
    %[row,col] = size(X);
    %计算临接矩阵W,采用高斯核函数RBF
    %W = zeros(row,row);
    %index_all = zeros(row,row);
    %for i = 1:row                   % 全连接构造相似度矩阵W
    %   for j=1:row
    %        if i ~= j
    %            W(i,j) = exp((-sum((X(i,:)-X(j,:)).^2))/(2*sigma.^2));
    %       end
    %    end
    %end
    $cluster_idx{1, i} = spectral_clustering(W,k);
    %spectral_idx  (original sample)
    %spectral_idx1 = cluster_idx{1,sample};
    %[spectral_idx{1, m},~] = spectral_idx1;                           %spectral clustering
    
    %draw clustering result on original sample
    %figure(),gscatter(ori_sample(:,1),ori_sample(:,2),cluster_idx{1,i});
    %draw clutering result on normalized sample
    %figure(),gscatter(norm_sample(:,1),norm_sample(:,2),cluster_idx{1,i});
    %hold on
   
    label = unique(cluster_idx{1, i});                                    %number of categories
    length(label);
    norm_dikind = [];                                                     %points in each category
    ori_dikind = [];
    
    for j = 1:length(label)
        labelname = label(j);
        if labelname == -1
            continue                                 %Denoising
        else
            norm_dikind = [norm_sample(find(cluster_idx{1, i} == labelname),:)];
            ori_dikind =  [ori_sample(find(cluster_idx{1, i} == labelname),:)];
            C = mean(norm_dikind, 1);               %calculation center point
            C1 = mean(ori_dikind, 1);               %calculation center point in [0,10]
            %plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
            norm_Allcenters = [norm_Allcenters; C];
            ori_Allcenters = [ori_Allcenters; C1];
        end
    end
end



%clustering of centroids
[norm_Allcenters_idx,~] =  kmeans(norm_Allcenters,k);
[ori_Allcenters_idx,~] = kmeans (ori_Allcenters,k);

%polt centroid clustering results
figure(),gscatter(norm_Allcenters(:,1),norm_Allcenters(:,2), norm_Allcenters_idx);
figure(),gscatter(ori_Allcenters(:,1),ori_Allcenters(:,2), ori_Allcenters_idx);
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
    %meankd = mean(kd);    
    
    if check ~= 1
        meankd = mean(kd);    %mean value
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

% the Gaussian density function
label = unique(norm_Allcenters_idx);
figure()
for i = 1:length(unique(label))
    labelname = label(i);
    if labelname == -1
        continue
    else
        kd = norm_Allcenters(find(norm_Allcenters_idx == labelname),:);
        meankd = mean(kd);    
        meankds = [meankds;meankd];
        % plot bivariate normal distribution probability density curves
        %x = -0.5 : 0.001 : 0.5 ;  
        %y = -0.5 : 0.001 : 0.5 ;      
        [X,Y] = meshgrid(x, y);
        %x Mean    Variance    Mean Variance
        U1 = meankd(1);
        DX = var(kd(:,1));
        dx = sqrt(DX);
        %y Mean    Variance    Mean Variance
        U2 = meankd(2);
        DY = var(kd(:,2));
        dy = sqrt(DY);
        
        % covariance matrix 2*2
        covkd = cov(kd(:,1), kd(:,2));
        r = covkd(2) / (dx * dy);
        COV = det(covkd);
        part1 = 1 / ( 2 * pi * dx * dy * sqrt( 1 - r^ 2 ));
        x_u = kd-meankd;
        %part2 =  -1 / (2*x_u.'/covkd*x_u);
        p1 = - 1 / ( 2 * ( 1 - r^ 2 ));
        px = (X - U1).^ 2. / DX;
        py = (Y - U2).^ 2. / DY;
        pxy = 2 * r.* (X - U1).* (Y - U2)./ (dx * dy);
        
        Z{1,i} = part1 * exp(p1 * (px - pxy + py));
        zlim([0 80]);
        mesh(X, Y, Z{1,i})
        hold on
       % xlabel('X1')
       % ylabel('X2')
        zlabel('Probability Density')
        shading interp 
        set(gca,'FontSize',20);
    end
end

% polt contour map of the center point
figure()
[~, ~, XMesh, YMesh, ZMesh, colorList]=density2C(norm_Allcenters(:,1), norm_Allcenters(:,2), x, y);
colormap(colorList)
contourf(XMesh, YMesh, ZMesh, 10)
set(gca,'FontSize',20);
