data = xlsread('IMU_GPS.csv');
X_ACC = data(:,4)*9.8;%X����ٶ� 
Y_ACC = data(:,5)*9.8;%Y����ٶ�
earth = 6378137;%����뾶
L = length(X_ACC);

lat_1 = zeros(L,1);
lon_1 = zeros(L,1);
lat = zeros(L,1);
lon = zeros(L,1);
d_x = zeros(L,1);
d_y = zeros(L,1);
d_lon = zeros(L,1);
d_lat = zeros(L,1);

%תΪ��γ��
for i = 1:L
    lon_1(i) = floor(data(i,2));
    lat_1(i) = floor(data(i,3));
    lon(i) = lon_1(i) + (mod(data(i,2),1)/60)*100;
    lat(i) = lat_1(i) + (mod(data(i,3),1)/60)*100;
end
%%%%%%%%%%%%%%%%%%%
%%��������ϵ����
for i = 1:L
   d_lon(i) = lon(i)-lon(1); 
   d_lat(i) = lat(i)-lat(1);
   d_x(i) = 2 * 3.14 * earth * (d_lon(i)/360);
   d_y(i) = 2 * 3.14 * (earth / cos(d_lat(i)/360 * 2 * 3.14))*(d_lat(i)/360);
end

pm_x = d_x; %GPS���Ȳ�
am_x = X_ACC;%X����ٶ�����
pm_y = d_y;  %GPSγ�Ȳ�
am_y = Y_ACC;%Y����ٶ�����
t = 0.02; %��������
X1kf = zeros(3,1);
X2kf = zeros(3,1);
%״̬ת�ƾ����Э�������
F = [1,t,0.5*t*t;0,1,t;0,0,1];%״̬ת�ƾ���
P1 = [10,0,0;0,2,0;0,0,0.8]; %���徭��γ��Э�������
P2 = [10,0,0;0,2,0;0,0,0.8];
Q = [10,0,0;0,0.8,0;0,0,0.8];%��������
H = [1,0,0;0,0,1];%�۲����
R = [100,0;0,100]; %�۲�����
I = [1,0,0;0,1,0;0,0,1];
XX = zeros(L,1);
YY = zeros(L,1);
%�Ծ�γ�ȷֱ��ں�
for i = 1:L
    %Ԥ��
    X1_pre = F * X1kf;
    P1_pre = F*P1*F'+Q;
    
    Z1 = [pm_x(i);5];
    e1 = Z1 - H * X1_pre;
    Kg1 = P1_pre * H' * inv(H*P1_pre*H' + R);
    %����
    X1kf = X1_pre + Kg1*e1;
    P1 = (I - Kg1*H)*P1_pre;
    XX(i) = X1kf(1);
end
%�Ծ�γ�ȷֱ��ں�
for i = 1:L
    %Ԥ��
    X2_pre = F * X2kf;
    P2_pre = F*P2*F'+Q;
    Z2 = [pm_y(i);5];
    e2 = Z2 - H * X2_pre;
    Kg2 = P2_pre * H' * inv(H*P2_pre*H' + R);
    %����
    X2kf = X2_pre + Kg2*e2;
    P2 = (I - Kg2*H)*P2_pre;
    YY(i) = X2kf(1);
end
%��ͼ
figure
hold on; box on;
plot(d_y,d_x,'r');
plot(YY,XX,'b');
legend('ԭʼ�켣','EKF�켣');
