function [sp,product,result] = Sp(X,Covariance)
%求分子高低维都可以
%分母用各个协方差矩阵的行列式的绝对值开平方
%X是k行n维的数据 类中心点间的距离
%Covariance是k行1列的协方差行列式
[row,col] = size(X);
sp = 2 / ((row^2)-row);
sum = 0;
for i=1:row
    for j=i+1:row
        Dist = sqrt(X(i,:).^2*ones(size(X(j,:)'))+ones(size(X(i,:)))*(X(j,:)').^2-2*X(i,:)*X(j,:)');
        sum = sum + Dist;
    end
end
sp = sp * sum;
[row,col] = size(Covariance);
product = 1;
for i=1:col:row
    %product = product * log(sqrt(abs(det(Covariance(i:i+col-1,:))))+exp(1));
    %product = product * sqrt(abs(det(Covariance(i:i+col-1,:)))); %初始
    product = product * sqrt(abs(det(Covariance(i:i+col-1,:)))); 
end
%product = log(product+exp(1));
result = sp / product;
result =  log10((result)+1);
end