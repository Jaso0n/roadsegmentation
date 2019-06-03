clear  all; close all; clc
[fileName,pathName] = uigetfile('*.*','Please select an image');%uigetfile打开一个模拟对话框，返回文件的名称和路径。点×两个参数都返回0

if(fileName)
    fileName = strcat(pathName,fileName);%strcat：水平串联字符串
    fileName = lower(fileName);%一致的小写字母形式
else
    msgbox('Please select an image');
    return; %退出程序
end
 
I = imread(fileName);
I = imresize(I,0.4);
subplot(2,2,1),imshow(I);title('原图');
if( ~( size(I,3)-3 ))%如果是RGB图像
    I = rgb2gray(I);%转化为单通道灰度图
    [M, N] = size(I);
end
I = im2double(I);%图像灰度值归一化到[0,1]之间

%种子点的交互式选择
if( exist('x','var') == 0 && exist('y','var') == 0)
    subplot(2,2,2),imshow(I,[]);title('灰度图');
    hold on;
    [y,x] = getpts;%鼠标取点  回车确定
    x = round(x);%四舍五入为最近的整数
    y = round(y);
end

J = zeros(size(I)); %记录区域生长所得到的区域
reg_mean = I(x,y);%表示分割好的区域的平均值，初始化为种子的灰度值
reg_size = 1; %分割到的区域，初始化只有种子点一个
neg_free = 10000; %动态分配内存的时候每次申请的连续空间大小
neg_list = zeros(neg_free,3);
%定义邻域列表，并预分配用于存储待分析的像素点的坐标值和灰度值的空间
neg_pos = 0;%用于记录neg_list中的待分析像素点个数
pixdist = 0;%记录最新像素点增加到分割区域后的距离测度
%如果当前坐标为（x, y),那么通过neigb可以得到四邻域像素的位置
neigb = [-1 0;
          1  0;
          0 -1;
          0  1];

%开始进行区域生长,当所有待分析的邻域像素点和已经分割好的区域像素点的灰度值距离
%大于reg_max,区域生长结束

while pixdist<0.15 && reg_size<numel(I)
    %增加新的邻域像素到neg_list中
    for j = 1:4
        xn = x + neigb(j,1);
        yn = y + neigb(j,2);
        %检查邻域像素是否超过了图像的边界
        ins = (xn>=1)&&(yn>=1)&&(xn<=M)&&(yn<=N);
        %如果邻域像素在图像内部且未分割，添加到邻域列表中
        if ins && J(xn,yn)==0
            neg_pos = neg_pos+1;%待分析像素点+1
            neg_list(neg_pos,:) = [xn, yn, I(xn, yn)];%存储对应点的位置和灰度值
            J(xn, yn) = 1;%标注该位置像素点已经被访问，但不表示已经被分割
        end
    end
    
%如果分配的存储空间不够，申请新的内存空间
    if neg_pos +10 > neg_free
        neg_free = neg_free +10000;
        neg_list((neg_pos+1):neg_free,:) = 0;
    end
    %从所有待分析的像素点中选择一个像素点，该点的灰度值和已经分割好区域灰度平均值的
    %差的绝对值是待分析像素中最小的
    dist = abs(neg_list(1:neg_pos,3) - reg_mean);
    [pixdist, index] = min(dist);
    %计算区域的新的均值
    reg_mean = (reg_mean*reg_size +neg_list(index,3))/(reg_size + 1);
    reg_size = reg_size + 1;
    %将旧的种子点标记为已经分割好的区域像素点
    J(x,y)=2;%标志该像素点已经是分割好的像素点
    x = neg_list(index, 1);
    y = neg_list(index, 2);

    %将新的种子点从待分析的邻域像素列表中移除
    neg_list(index,:) = neg_list(neg_pos,:);%新的种子点元素值由最后一个待分析邻域值取代
    neg_pos = neg_pos - 1;
end

J = (J==2);%我们之前将分割好的像素点转化为逻辑数组
hold off;
subplot(2,2,3),imshow(J);title('道路分割结果');
subplot(2,2,4),imshow(I+J);title('灰度图覆盖');
