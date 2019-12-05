% Script file: test_GSO_GSOhj.m
%
% Purpose:
% 测试“单独使用GSO”和“GSO+hookejeeves”算法的搜索速度（评价指标：评价函数的指标大小和评估次数）
% 
% Record of revisions:
% Date         Programmer    Description of change
% ==========   ==========    =====================
% 2019-11-21   Chaofan Yu    Original code
%
%% ============================================ GSO =========================================
clear all,clc,close all
global PopSize
global fname
NUMBER1 = 1;
NUMBER2 = 30;
rangersPercent = 0.2;
pursuitAngleCoefficient = 2;
turingAngleCoefficient = 8;
initAngle = pi/4;
aCoefficient = 1;
bCoefficient = 1;
%                 1  2   3   4   5
maxIter=50;  % 迭代次数
% 1:SCH;    2:FON   3:POL   4:KUR   5:ZDT1
%       1 2 3 4 5
flag = [1 1 1 1 1];
% =====================f2================================
if flag(2)==1
    PopSize = 30;
    lmaxCoefficient = 1;
    initAngle = pi/4;
    aCoefficient = 1;
    bCoefficient = 1;
    fname = 'function2';  % 函数入口
    fprintf('===============================FON===============================\n');
    NDim = 24;
    numObjec=2;
    initProducer=[];
    direcDul = 6;  NumCoeffi = 1;  tempFlag = 1;   c1 = 1; c2 = 1.0;   NumCoeffi = 0;  numShift1 = 1;
    if NDim > 6
        flagDirec = 1;      %1则用修改后的，0用原有的坐标转换方程
    else
        flagDirec = 0;      %1则用修改后的，0用原有的坐标转换方程
    end
    global hvcounts;
    hviters=zeros(2,30);
    for i=1:30
        hvcounts=1;
        [ fbestvals, bestmembers, archiveNew, fvaluesNew, fvaluesAll, archiveAll,hv,hvcounts] = GSOMP_3(fname,NDim,maxIter,flagDirec,numObjec,initProducer,rangersPercent,pursuitAngleCoefficient,turingAngleCoefficient,lmaxCoefficient,initAngle,aCoefficient,bCoefficient,direcDul,c1,c2,NumCoeffi,numShift1,hvcounts);
        hviters(1,i)=hvcounts;
        hviters(2,i)=hv(end,1);
    end
    fprintf('\n');
    save('ycf_GSO.mat');
    %====================================
    figure(1)
    fvaluesNew(:,2)=fvaluesNew(:,2)*10000+1450;
    size(fvaluesNew);
    plot( fvaluesNew(:,1), fvaluesNew(:,2), '.r' );
    %fvaluesNew(:,3)=10^5-fvaluesNew(:,3);
    grid on
    title('FON--Pareto Front Obtained by FGSOMP');
end
hold on;
fprintf('hvcounts==%d\n',hvcounts);
figure(2)
plot(hv);

%% =============================================== GSO+hookejeeves =============================================
clear all,clc,close all
global PopSize
global fname
NUMBER1 = 1;
NUMBER2 = 30;
rangersPercent = 0.2;
pursuitAngleCoefficient = 2;
turingAngleCoefficient = 8;
initAngle = pi/4;
aCoefficient = 1;
bCoefficient = 1;
%                 1  2   3   4   5
maxIter=50;  % 迭代次数
% 1:SCH;    2:FON   3:POL   4:KUR   5:ZDT1
%       1 2 3 4 5
flag = [1 1 1 1 1];
% =====================f2================================
if flag(2)==1
    PopSize = 30;
    lmaxCoefficient = 1;
    initAngle = pi/4;
    aCoefficient = 1;
    bCoefficient = 1;
    tic
    fname = 'function2';  % 函数入口
    fprintf('===============================FON===============================\n');
    NDim = 24;
    numObjec=2;
    initProducer=[];
    direcDul = 6;  NumCoeffi = 1;  tempFlag = 1;   c1 = 1; c2 = 1.0;   NumCoeffi = 0;  numShift1 = 1;
    if NDim > 6
        flagDirec = 1;      %1则用修改后的，0用原有的坐标转换方程
    else
        flagDirec = 0;      %1则用修改后的，0用原有的坐标转换方程
    end
    global hvcounts;
    hviters=zeros(2,30);
    for i=1:30
        hvcounts=1;
        [ fbestvals, bestmembers, archiveNew, fvaluesNew, fvaluesAll, archiveAll,hv,hvcounts] = GSOMP_2hv(fname,NDim,maxIter,flagDirec,numObjec,initProducer,rangersPercent,pursuitAngleCoefficient,turingAngleCoefficient,lmaxCoefficient,initAngle,aCoefficient,bCoefficient,direcDul,c1,c2,NumCoeffi,numShift1,hvcounts);
        hviters(1,i)=hvcounts;
        hviters(2,i)=hv(end,1);
    end
    fprintf('\n');
    t= toc;
    save('ycf_GSOhj.mat');
    %====================================
    figure(1)
    fvaluesNew(:,2)=fvaluesNew(:,2)*10000+1450;
    size(fvaluesNew);
    plot( fvaluesNew(:,1), fvaluesNew(:,2), '.r' );
    %fvaluesNew(:,3)=10^5-fvaluesNew(:,3);
    grid on
    title('FON--Pareto Front Obtained by FGSOMP');
end
hold on;
fprintf('hvcounts==%d\n',hvcounts);
figure(2)
plot(hv);

%% ====================================== 作图分析两种算法的性能 ===============================
clear all;clc;close all;
% 初始化
iter=1:30;
hv_gso=zeros(2,30);
hv_gsohj=zeros(2,30);
% 载入"GSO"算法的hv指标
hv_gso=load('ycf_GSO.mat','hviters');
%  GSO算法---hv评价函数的评价次数
figure(1)
plot(iter,hv_gso.hviters(1,:),'b-',iter,hv_gso.hviters(2,:),'r--');
axis([0 30 0 280]);
title('hv indicator of GSO algorithm');
xlabel('number of iterations');
ylabel('hv indicators');
legend('hv函数调用次数','hv评价指标');
% 
% 载入"GSO+hookejeeves"算法的hv指标
hv_gsohj=load('ycf_GSOhj.mat','hviters');
%  GSO+hookejeeves算法---hv评价函数的评价次数
figure(2)
plot(iter,hv_gsohj.hviters(1,:),'b-',iter,hv_gsohj.hviters(2,:),'r--');
axis([0 30 0 280]);
title('hv indicator of GSOhj algorithm');
xlabel('number of iterations');
ylabel('hv indicators');
legend('hv函数调用次数','hv评价指标');
% 
% 对比两种算法的hv评价次数
figure(3)
plot(iter,hv_gso.hviters(2,:),'r-+',iter,hv_gsohj.hviters(2,:),'b-*');
title('Detail comparison of GSO and GSOhj');
xlabel('number of iterations');
ylabel('hv indicators');
legend('gso--hv评价指标','gsohj--hv评价指标');
%
% 显示两种算法的评价指标的平均值和均方差
hv_mu1=mean(hv_gso.hviters(2,:),2);
hv_std1=std(hv_gso.hviters(2,:));
fprintf('GSO algorithm: hv_mu==%d\t,hv_std==%d\n',hv_mu1,hv_std1);
%
hv_mu2=mean(hv_gsohj.hviters(2,:),2);
hv_std2=std(hv_gsohj.hviters(2,:));
fprintf('GSOhj algorithm: hv_mu==%d\t,hv_std==%d\n',hv_mu2,hv_std2);

