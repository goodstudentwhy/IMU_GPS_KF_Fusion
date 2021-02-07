data = xlsread('IMU_GPS.csv');
X_ACC = data(:,4)*9.8;%X轴加速度 
Y_ACC = data(:,5)*9.8;%Y轴加速度
earth = 6378137;%地球半径
L = length(X_ACC);

lat_1 = zeros(L,1);
lon_1 = zeros(L,1);
lat = zeros(L,1);
lon = zeros(L,1);
d_x = zeros(L,1);
d_y = zeros(L,1);
d_lon = zeros(L,1);
d_lat = zeros(L,1);

%转为经纬度
for i = 1:L
    lon_1(i) = floor(data(i,2));
    lat_1(i) = floor(data(i,3));
    lon(i) = lon_1(i) + (mod(data(i,2),1)/60)*100;
    lat(i) = lat_1(i) + (mod(data(i,3),1)/60)*100;
end
%%%%%%%%%%%%%%%%%%%
%%计算距离上电距离
for i = 1:L
   d_lon(i) = lon(i)-lon(1); 
   d_lat(i) = lat(i)-lat(1);
   d_x(i) = 2 * 3.14 * earth * (d_lon(i)/360);
   d_y(i) = 2 * 3.14 * (earth / cos(d_lat(i)/360 * 2 * 3.14))*(d_lat(i)/360);
end

pm_x = d_x; %GPS经度差
am_x = X_ACC;%X轴加速度数据
pm_y = d_y;  %GPS纬度差
am_y = Y_ACC;%Y轴加速度数据
t = 0.02; %运行周期
X1kf = zeros(3,1);
X2kf = zeros(3,1);
%状态转移矩阵和协方差矩阵
F = [1,t,0.5*t*t;0,1,t;0,0,1];%状态转移矩阵
P1 = [10,0,0;0,2,0;0,0,0.8]; %定义经度纬度协方差矩阵
P2 = [10,0,0;0,2,0;0,0,0.8];
Q = [10,0,0;0,0.8,0;0,0,0.8];%过程噪声
H = [1,0,0;0,0,1];%观测矩阵
R = [100,0;0,100]; %观测噪声
I = [1,0,0;0,1,0;0,0,1];
XX = zeros(L,1);
YY = zeros(L,1);
%对经纬度分别融合
for i = 1:L
    %预测
    X1_pre = F * X1kf;
    P1_pre = F*P1*F'+Q;
    
    Z1 = [pm_x(i);5];
    e1 = Z1 - H * X1_pre;
    Kg1 = P1_pre * H' * inv(H*P1_pre*H' + R);
    %更新
    X1kf = X1_pre + Kg1*e1;
    P1 = (I - Kg1*H)*P1_pre;
    XX(i) = X1kf(1);
end
%对经纬度分别融合
for i = 1:L
    %预测
    X2_pre = F * X2kf;
    P2_pre = F*P2*F'+Q;
    Z2 = [pm_y(i);5];
    e2 = Z2 - H * X2_pre;
    Kg2 = P2_pre * H' * inv(H*P2_pre*H' + R);
    %更新
    X2kf = X2_pre + Kg2*e2;
    P2 = (I - Kg2*H)*P2_pre;
    YY(i) = X2kf(1);
end
%画图
figure
hold on; box on;
plot(d_y,d_x,'r');
plot(YY,XX,'b');
legend('原始轨迹','EKF轨迹');
