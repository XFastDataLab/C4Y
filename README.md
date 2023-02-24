# C4Y
C4Y is a novel distributed clustering index, based on the sampling theories: 
(1) The data of each distributed node is considered as a sample of the entire data set. Each sample should have similar data 
distribution, as well as category distribution, to the entire data set; 
(2) the centers of each category in all samples conform to Gaussian
distribution. Thus, the proposed metric works by determining the extent to which the centers of each category in different samples
conform to Gaussian distribution, and the degree of difference between those categories.

***********************************************************************************
The C4Y program was compiled under Windows using matlab R2016b.
***********************************************************************************

Files
===================================================================================
These program mainly containing:

-startup code named `Distrubted_C4Y.m` and `images_C4Y.m`.

-two main functions of C4Y named `SP.m` and `density2C.m`.

-some data sets.

Dataset Format
===================================================================================
The dataset should be given in a text file of the following format:
-First, prepare the data set, no column number and row number are required.When using malab to read the given mat format data in the project, you can use the `load` function. By default, the variable name of the data matrix is named `points`.
For instance, the first 10 lines of the sample dataset "4.mat"(whose data number is 3412 and dimension is 2) are shown as below:

-14.1897478152872	157.776071507798
-13.1823191128589	156.783555810430
-15.1822635126559	156.768642805368
-12.1748904104294	155.791040113062
-14.1748348102266	155.776127108000
-16.1747792100238	155.761214102939
-26.1745012090096	155.686649077631
-28.1744456088068	155.671736072569
-11.1674617080001	154.798524415694
-13.1674061077973	154.783611410632

```matlab
load('4.mat','points');
dataall = points;
```

An example of quick start
===================================================================================
Step1:
Open the startup code `Distrubted_C4Y.m`, and load the data sets we prepared already.
```matlab
load('4.mat','points');
dataall = points;
```

Step2：
Choose to use DBSCAN, kmeans, meanshift, DPeak or spectral clustering to perform distributed clustering, and adjust the parameters. (the following is an example of DBSCAN.)
```matlab
% parameters
n = 10;                                        % number of samples
sample = 900;                                  % the number of data points contained in each sample
k=5;                                           % number of clustering categories
[epsilonD,minptsD]=deal( 0.040 , 15 );         
[epsilonDi,minptsDi]=deal( 0.045 , 14 );       % DBSCAN parameters
dione = dataall;
cluster_idx{1, 1} = dbscan(dione, epsilonD, minptsD);           %DBSCAN

%random sample
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
    figure(),gscatter(ori_sample(:,1),ori_sample(:,2),cluster_idx{1,i});
   
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
```

Step2：
Use kmeans clustering for all clustering centroids and call `SP.m` to calculate the value of C4Y and draw its corresponding probability density function plot and contour plot.
```matlab
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
```
Step4: -Press the “Run” button in Matlab and the code will be run.

Output Format
===================================================================================
The first 3 columns are randomly selected results from 500 samples; the fourth 
column represents the distribution of all cluster centers; the last two columns 
are its corresponding PDF and contour, respectively.As shown below:
![](https://github.com/XFastDataLab/C4Y/blob/main/result.png)

DBSCAN with ε = 30, MinPts = 7, which has C4Y = 5.74.
