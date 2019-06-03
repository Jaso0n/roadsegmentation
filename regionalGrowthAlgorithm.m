clear  all; close all; clc
[fileName,pathName] = uigetfile('*.*','Please select an image');%uigetfile��һ��ģ��Ի��򣬷����ļ������ƺ�·���������������������0

if(fileName)
    fileName = strcat(pathName,fileName);%strcat��ˮƽ�����ַ���
    fileName = lower(fileName);%һ�µ�Сд��ĸ��ʽ
else
    msgbox('Please select an image');
    return; %�˳�����
end
 
I = imread(fileName);
I = imresize(I,0.4);
subplot(2,2,1),imshow(I);title('ԭͼ');
if( ~( size(I,3)-3 ))%�����RGBͼ��
    I = rgb2gray(I);%ת��Ϊ��ͨ���Ҷ�ͼ
    [M, N] = size(I);
end
I = im2double(I);%ͼ��Ҷ�ֵ��һ����[0,1]֮��

%���ӵ�Ľ���ʽѡ��
if( exist('x','var') == 0 && exist('y','var') == 0)
    subplot(2,2,2),imshow(I,[]);title('�Ҷ�ͼ');
    hold on;
    [y,x] = getpts;%���ȡ��  �س�ȷ��
    x = round(x);%��������Ϊ���������
    y = round(y);
end

J = zeros(size(I)); %��¼�����������õ�������
reg_mean = I(x,y);%��ʾ�ָ�õ������ƽ��ֵ����ʼ��Ϊ���ӵĻҶ�ֵ
reg_size = 1; %�ָ�����򣬳�ʼ��ֻ�����ӵ�һ��
neg_free = 10000; %��̬�����ڴ��ʱ��ÿ������������ռ��С
neg_list = zeros(neg_free,3);
%���������б���Ԥ�������ڴ洢�����������ص������ֵ�ͻҶ�ֵ�Ŀռ�
neg_pos = 0;%���ڼ�¼neg_list�еĴ��������ص����
pixdist = 0;%��¼�������ص����ӵ��ָ������ľ�����
%�����ǰ����Ϊ��x, y),��ôͨ��neigb���Եõ����������ص�λ��
neigb = [-1 0;
          1  0;
          0 -1;
          0  1];

%��ʼ������������,�����д��������������ص���Ѿ��ָ�õ��������ص�ĻҶ�ֵ����
%����reg_max,������������

while pixdist<0.15 && reg_size<numel(I)
    %�����µ��������ص�neg_list��
    for j = 1:4
        xn = x + neigb(j,1);
        yn = y + neigb(j,2);
        %������������Ƿ񳬹���ͼ��ı߽�
        ins = (xn>=1)&&(yn>=1)&&(xn<=M)&&(yn<=N);
        %�������������ͼ���ڲ���δ�ָ��ӵ������б���
        if ins && J(xn,yn)==0
            neg_pos = neg_pos+1;%���������ص�+1
            neg_list(neg_pos,:) = [xn, yn, I(xn, yn)];%�洢��Ӧ���λ�úͻҶ�ֵ
            J(xn, yn) = 1;%��ע��λ�����ص��Ѿ������ʣ�������ʾ�Ѿ����ָ�
        end
    end
    
%�������Ĵ洢�ռ䲻���������µ��ڴ�ռ�
    if neg_pos +10 > neg_free
        neg_free = neg_free +10000;
        neg_list((neg_pos+1):neg_free,:) = 0;
    end
    %�����д����������ص���ѡ��һ�����ص㣬�õ�ĻҶ�ֵ���Ѿ��ָ������Ҷ�ƽ��ֵ��
    %��ľ���ֵ�Ǵ�������������С��
    dist = abs(neg_list(1:neg_pos,3) - reg_mean);
    [pixdist, index] = min(dist);
    %����������µľ�ֵ
    reg_mean = (reg_mean*reg_size +neg_list(index,3))/(reg_size + 1);
    reg_size = reg_size + 1;
    %���ɵ����ӵ���Ϊ�Ѿ��ָ�õ��������ص�
    J(x,y)=2;%��־�����ص��Ѿ��Ƿָ�õ����ص�
    x = neg_list(index, 1);
    y = neg_list(index, 2);

    %���µ����ӵ�Ӵ����������������б����Ƴ�
    neg_list(index,:) = neg_list(neg_pos,:);%�µ����ӵ�Ԫ��ֵ�����һ������������ֵȡ��
    neg_pos = neg_pos - 1;
end

J = (J==2);%����֮ǰ���ָ�õ����ص�ת��Ϊ�߼�����
hold off;
subplot(2,2,3),imshow(J);title('��·�ָ���');
subplot(2,2,4),imshow(I+J);title('�Ҷ�ͼ����');
