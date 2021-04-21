
clc
clear all
close(figure(gcf))
X=cell2mat(struct2cell(load('Xlineplain.mat')));
% main function
[Q,Qd,obj,b]=LRR(X,2,0.0001,2,1,1);

obj(50)
NUM1=find(Q(:,1)==1);
NUM2=find(Q(:,2)==1);
[l1,l]=size(NUM1);
[l2,l]=size(NUM2);
X1=zeros(3,l1);
X2=zeros(3,l2);

for i=1:l1
    X1(:,i)=X(:,NUM1(i));
end
for i=1:l2
    X2(:,i)=X(:,NUM2(i));
end

scatter3(X1(1,:),X1(2,:),X1(3,:),'o','filled','r');
hold on;
scatter3(X2(1,:),X2(2,:),X2(3,:),'*','b');

