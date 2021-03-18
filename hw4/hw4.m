clear all; close all; clc;
load fisheriris

[train_image, train_label] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
train_image = im2double(reshape(train_image, size(train_image,1)*size(train_image,2), []).');
train_label = im2double(train_label);
train_image = train_image';

[test_image,  test_label] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');
test_image = im2double(reshape(test_image, size(test_image,1)*size(test_image,2), []).');
test_label = im2double(test_label);
test_image = test_image';


meanFort = mean(train_image,2);
train_image = double(train_image)-repmat(mn,1,length(train_image)); 

[U, S, V ] = svd(train_image, 'econ');

en = 0;
result = sum(diag(S));
thres = 0.8; 
r = 0;

while en < thres
    r = r + 1;
    en = en + S(r,r)/result;
end

rank = r;

train_image = (U(:, 1:rank))'*train_image; 
test_image = (U(:, 1:rank))'*test_image;

lamd = diag(S).^2;

for i = [5,8]
    Projection_idx = train_image(:, find(train_label == i));
    plot3(Projection_idx(1,:), Projection_idx(2,:), Projection_idx(3,:),...
          'o', 'DisplayName', num2str(i)); 
    hold on
end
legend show


sc = get(groot, 'ScreenSize'); 


X = train_image;
T = test_image;


dimen = size(train_image,1);
s1 = zeros(dimen);
s2 = zeros(dimen); 
trainsize = size(train_image, 2);
testsize = size(test_image, 2); 
Mu = mean(train_image, 2);  

for i = [0,1]
    
    mask = (train_label ==  i);
    x = X(:, mask);
    ni = size(x, 2);
    pi = ni / trainsize;
    mu_i = mean(x, 2);

    Si = (x - repmat(mu_i, [1,ni]))*(x - repmat(mu_i, [1,ni]))';

    s1 = s1 + Si ;
    s2 = s2 + (mu_i - Mu) * (mu_i - Mu)';
end

M = pinv(s1) * s2; 
[U, D, V] = svd(M);

G2 = U(:,1:rank);
Y2 = G2' * X;

2dfig = figure('Name', '2-D Plot');
set(2dfig,[60 60 sc(3)-120 sc(4) - 140]);

for number = [0,1]
    
    mask = (train_label ==  number);
    a = Y2(1,mask);
    b = Y2(2,mask);
    
    d = [a'; b'];
    plot(d, 1*number*ones(size(d)),'o',...
        'DisplayName', num2str(number)); hold on 
    
end
legend show

%% accuracy

function [accu] = classifyNN(test_data, train_data, test_label, train_label)

train_size = size(train_data, 2);
test_size = size(test_data, 2);
c = zeros(test_size, 1);

parfor test_digit = 1:1:test_size

    test_mat = repmat(test_data(:, test_digit), [1,train_size]);
    distance = sum(abs(test_mat - train_data).^2);
    [M,I] = min(distance);
    if train_label(I) == test_label(test_digit)
        c(test_digit) = c(test_digit) + 1;
    end
end

accu = double(sum(c)) / test_size;
end

n1 = 0;
n2 = 1;

Y = G2' * X;
Y_t = G2'* T;

train_n= Y(:, find(train_label == n1|train_label ==n2));
test_n = Y_t(:, find(test_label == n1 |test_label ==n2)); 

accuracy = classifyNN(test_n, train_n,...
    test_label(find(test_label == n1 |test_label ==n2)), ...
    train_label(find(train_label == n1 |train_label ==n2)));

%% SVM
feat = 50;
z0 = [];
o1 = [];
t2 = [];
for i = 1:10000
    if train_labels(i) == 2
        t2(:,end+1) = PCA(1:feature,i);
    elseif train_labels(i) == 1
        o1(:,end+1) = PCA(1:feature,i);
    elseif train_labels(i) == 0
        z0(:,end+1) = PCA(1:feature,i);
end end
size = min(size(z0,2),min(size(o1,2),size(t2,2)));
zeroone  = [z0 o1]
zerotwo  = [t2 o1]
svm1 = fitcsvm(zeroone',zero1one);
cv1 = crossval(svm1);
svma1 = kfoldLoss(CV_1_0)
svm2 = fitcsvm(zerotwo',zero2two);
cv2 = crossval(svm2);
svma2 = kfoldLoss(cv2)

%% decision tree
twave = dcwavelet(train_images);
tree=fitctree(twave',train_labels,'CrossVal','on');
terror = kfoldLoss(tree)
label1  = [zeros([1 991]) ones([1 991])]
label2  = [zeros([1 991]) 2*ones([1 991])]
tree1=fitctree(zero_one', label1,'CrossVal','on');
tree2=fitctree(zero_two', label2,'CrossVal','on');
terror1 = kfoldLoss(tree1)
terror2 = kfoldLoss(tree2)
  
