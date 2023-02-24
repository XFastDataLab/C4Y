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
Choose to use DBSCAN, kmeans, meanshift, DPeak or spectral clustering to perform distributed clustering, and adjust the parameters
```matlab

