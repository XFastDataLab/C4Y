function [idx,p,s,K,C,Klist] = Dpeak(dataPts,d,kinds)
%d聚类距离，k聚类个数 dataPts n*2
D=pdist2(dataPts,dataPts);
[row,col] = size(dataPts);
n = row;
idx = zeros(n,1);
%局部密度 从大到小排序n*1
p = [];
%记录点集在原数据的位置n*1
xy = [];
%局部密度对应聚类中心距离n*1
s = [];
%p*s n*1;
K = [];

%计算p,记录xy
for i = 1:n
    %xy = [xy;i];
    sum = 0;
    for j=1:n
        if i == j
            continue
        else
            %di = DS(dataPts(i,:),dataPts(j,:));
            di = D(i,j);
            if di < d
                sum = sum + 1;
            end
        end
    end
    p = [p;sum];
end
%对p排序
[p,xy]=sort(p,'descend');
%计算s
%p
%p(1)
%计算最大密度点到所有点中的最大距离
maxd = 0;
for i = 1:n
    Di  = D(xy(1),xy(i));
    if maxd < Di
        maxd = Di;
    end
end
s = [s;maxd/maxd];
%计算密度大于本身点的点集中的最小距离
for i = 2:n
    mind = 999999999;
    for j = 1:i-1 %密度大的点集索引
        Di  = D(xy(i),xy(j));
        if mind >= Di
            mind = Di;
        end
    end
    s = [s;(mind)/(maxd)];
end

%s
%计算K
%figure(),gscatter(p, s);
for i = 1:n
    ki = s(i)*p(i);
    K = [K;ki];
end
%K,xy
%idx;
%选取前k个K
[K2,xy2] = sort(K,'descend');
%figure(),gscatter(dataPts(:,1), dataPts(:,2));
for i = 1:kinds
    %聚类中心索引
    ik = xy(xy2(i));
    hold on;
    C =  dataPts(ik,:);
    %plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
end
Centre = [];
%{
for i = 1:kinds
    %聚类中心索引
    ik = xy(xy2(i));
    Centre = [Centre;dataPts(ik,:)];
    for j = 1:n
        di = DS(dataPts(ik,:),dataPts(xy(j),:));
        if di < d.^2
            idx(xy(j)) = i;
        end
    end
end
%}
%聚类
X = dataPts;
Klist = cell(1,kinds);         %记录类别的索引，以便于分类
IDX=zeros(n,1);            %定义一个n行1列的矩阵
visited=false(n,1);        %创建一维的标记数组，全部初始化为false，代表还未被访问
for i=1:n %遍历1~n个所有的点
    %是否访问
    ki = xy2(n-i+1);%从最小k值点开始往前遍历
    if visited(ki) %被访问跳过
        continue
    end
    [root,list]= Findroot(ki);%root = 1~k
    Klist{1,root} = [Klist{1,root};list];
end
    function [root,list] = Findroot(ki)
        flag = 0;
        for j=1:kinds
            if ki == xy2(j)
                flag = 1;
                break
            end
        end
        if flag == 0 %递归调用
            mind = 9999999;
            for i=1:ki-1
                Di  = D(xy(i),xy(ki));
                if mind >= Di
                    root = i;
                    mind = Di;
                end
            end
            [root,listtmp] = Findroot(root);
        else
            root = j;
            listtmp = [];
        end 
        list = [listtmp;ki];        
    end
for i=1:kinds
    idxxy = Klist{1,i};
    for j=1:size(idxxy)
        idx = xy(idxxy(j,:));
        IDX(idx,1) =  i;
    end
end
idx = IDX';
   
end