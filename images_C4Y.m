close all;
%I = imread('1.jpg');
%I = imread('2.jpg');
%I = imread('6.png');
%I = imread('3.png');% k = 3 4 5
I = imread('4.png');% k = 3 4 5 
%I = imread('5.png');% k = 3 4 5
figure, imshow(I);

k = 3;
%报错可能要改投影图像画图精度
%Center_accuracy = 0.01;

[originalRAW,originalCOL,R_double, G_double, B_double, A_h, B_h, C_h, row_gray_re, A_l, B_l, C_l,A_row_gradient,B_row_gradient,C_row_gradient,A_line_gradient,B_line_gradient,C_line_gradient] = Characteristics(I);
[x,y]=meshgrid( (1 : originalRAW) ,fliplr(1:originalCOL));

%ori_data =[R_double, G_double, B_double, A_h, B_h, C_h, row_gray_re, A_l, B_l, C_l];
ori_data =[R_double, G_double, B_double];
ori_data=mapminmax(ori_data',1,10)';
%ori_data =[R_double, G_double, B_double, A_h, B_h, C_h, row_gray_re, A_l, B_l, C_l, A_row_gradient, B_row_gradient, C_row_gradient, A_line_gradient, B_line_gradient, C_line_gradient];
[idx,Centers] = kmeans(ori_data,k);
ori_data;

label = unique(idx);
length(label);
reshape(idx,originalRAW,originalCOL);

%idx_image = arrayfun(@(x,y)sprintf( '(%i,%i)' ,[x,y]),x,y, 'UniformOutput' , false );
%打印中心点
figure,gscatter(reshape(x, 1, originalRAW*originalCOL), reshape(y, 1, originalRAW*originalCOL), idx);
for j = 1:length(label)
    labelname = label(j);
    C = [x(idx == labelname), y(idx == labelname)];
    C = abs(mean(C));
    hold on
    plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
end

%聚类结果

%采样
Sample_value = 15440;   %15440
data = [reshape(x, originalRAW*originalCOL,1),reshape(y, originalRAW*originalCOL,1), R_double, G_double, B_double,A_h, B_h, C_h, row_gray_re, A_l, B_l, C_l];

n = 5;
%样本空间di
sample_cell = cell(1, n);
Allcenters = [];
Allcentersi = cell(1, length(label));
%fprintf('开始采样\n');
%sample_centers = cell(1, 3);
sample_centers = [];

for i = 1:n
    %fprintf('第%d次采样\n',i);
    %数据随机抽取样本
    sample_cell{1, i} = datasample(data, Sample_value);
    sample_x =  sample_cell{1, i}(:,1);
    sample_y =  sample_cell{1, i}(:,2);
    
    %ori_data = [di{1, i}(:,3), di{1, i}(:,4), di{1, i}(:,5), di{1, i}(:,6), di{1, i}(:,7), di{1, i}(:,8), di{1, i}(:,9), di{1, i}(:,10), di{1, i}(:,11), di{1, i}(:,12)];
    ori_data = [sample_cell{1, i}(:,3), sample_cell{1, i}(:,4), sample_cell{1, i}(:,5)];
    [idx,Centers]=kmeans(ori_data, k);
    sample_centers = [sample_centers;Centers];


    label = unique(idx);
    length(label);
    %reshape(idx,originalRAW,originalCOL);
    %画图

    xlim([0 1460])
    figure,gscatter(sample_x,sample_y,idx);
    set(gcf,'unit','normalized','position',[.2 .2 .12 .15])
    set(gca,'FontSize',15);  
    xlabel('')
    ylabel('')


    %中心点
    for j = 1:length(label)
        labelname = label(j);
        C = [sample_x(idx == labelname),sample_y(idx == labelname)];
        C = mean(C);
        Allcenters = [Allcenters;C];
        %hold on
       %plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
    end
    
end
%fprintf('采样结束\n');
%一 二维坐标
%对中心点聚类
Allcenters_idx = kmeans(Allcenters, k);%%采用kmeans聚类
%打印中心点聚类结果
figure()
gscatter(Allcenters(:,1),Allcenters(:,2), Allcenters_idx);
hold off

label = unique(Allcenters_idx);
meankds = [];
COVS = [];
for i = 1:length(unique(label))
    labelname = label(i);
    kd = Allcenters(find(Allcenters_idx == labelname),:);
    meankd = mean(kd);     %二维均值
    meankds = [meankds;meankd];
    covkd = cov(kd(:,1), kd(:,2));
    COVS = [COVS;covkd];
end
%fprintf('二维坐标下计算\n');
[sp,product,result] = Sp(meankds,COVS);
fprintf('样本结果\n%d\n',result);


%%按类打印中心点的等高线图
%{
Allcenters2 = [];
for i = 1:length(unique(label))
    figure()    
    labelname = label(i);
    Allcenters2 = Allcenters(Allcenters_idx == labelname,:);
    AX = Allcenters2(:,1);
    AY = Allcenters2(:,2);
    XX = 2*min(AX)-max(AX):1:2*max(AX)-min(AX);
    YY = 2*min(AY)-max(AY):1:2*max(AY)-min(AY);
    [~, ~, XMesh, YMesh, ZMesh, colorList]=density2C(AX, AY, XX, YY);
    colormap(colorList)
    contourf(XMesh, YMesh, ZMesh, 10)
end
%}
%打印中心点的等高线图

%fprintf('打印中心点的等高线图\n');
figure()

AX = Allcenters(:,1);
AY = Allcenters(:,2);
%XX = 2*min(AX)-max(AX):1:2*max(AX)-min(AX);
%YY = 2*min(AY)-max(AY):1:2*max(AY)-min(AY);

%[~, ~, XMesh, YMesh, ZMesh, colorList]=density2C(AX, AY,XX,YY);
[~, ~, XMesh, YMesh, ZMesh, colorList]=density2C(AX, AY,450:1:750,100:1:600);
colormap(colorList)
contourf(XMesh, YMesh, ZMesh, 10)




