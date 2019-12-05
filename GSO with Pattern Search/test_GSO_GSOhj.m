% Script file: test_GSO_GSOhj.m
%
% Purpose:
% ���ԡ�����ʹ��GSO���͡�GSO+hookejeeves���㷨�������ٶȣ�����ָ�꣺���ۺ�����ָ���С������������
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
maxIter=50;  % ��������
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
    fname = 'function2';  % �������
    fprintf('===============================FON===============================\n');
    NDim = 24;
    numObjec=2;
    initProducer=[];
    direcDul = 6;  NumCoeffi = 1;  tempFlag = 1;   c1 = 1; c2 = 1.0;   NumCoeffi = 0;  numShift1 = 1;
    if NDim > 6
        flagDirec = 1;      %1�����޸ĺ�ģ�0��ԭ�е�����ת������
    else
        flagDirec = 0;      %1�����޸ĺ�ģ�0��ԭ�е�����ת������
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
maxIter=50;  % ��������
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
    fname = 'function2';  % �������
    fprintf('===============================FON===============================\n');
    NDim = 24;
    numObjec=2;
    initProducer=[];
    direcDul = 6;  NumCoeffi = 1;  tempFlag = 1;   c1 = 1; c2 = 1.0;   NumCoeffi = 0;  numShift1 = 1;
    if NDim > 6
        flagDirec = 1;      %1�����޸ĺ�ģ�0��ԭ�е�����ת������
    else
        flagDirec = 0;      %1�����޸ĺ�ģ�0��ԭ�е�����ת������
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

%% ====================================== ��ͼ���������㷨������ ===============================
clear all;clc;close all;
% ��ʼ��
iter=1:30;
hv_gso=zeros(2,30);
hv_gsohj=zeros(2,30);
% ����"GSO"�㷨��hvָ��
hv_gso=load('ycf_GSO.mat','hviters');
%  GSO�㷨---hv���ۺ��������۴���
figure(1)
plot(iter,hv_gso.hviters(1,:),'b-',iter,hv_gso.hviters(2,:),'r--');
axis([0 30 0 280]);
title('hv indicator of GSO algorithm');
xlabel('number of iterations');
ylabel('hv indicators');
legend('hv�������ô���','hv����ָ��');
% 
% ����"GSO+hookejeeves"�㷨��hvָ��
hv_gsohj=load('ycf_GSOhj.mat','hviters');
%  GSO+hookejeeves�㷨---hv���ۺ��������۴���
figure(2)
plot(iter,hv_gsohj.hviters(1,:),'b-',iter,hv_gsohj.hviters(2,:),'r--');
axis([0 30 0 280]);
title('hv indicator of GSOhj algorithm');
xlabel('number of iterations');
ylabel('hv indicators');
legend('hv�������ô���','hv����ָ��');
% 
% �Ա������㷨��hv���۴���
figure(3)
plot(iter,hv_gso.hviters(2,:),'r-+',iter,hv_gsohj.hviters(2,:),'b-*');
title('Detail comparison of GSO and GSOhj');
xlabel('number of iterations');
ylabel('hv indicators');
legend('gso--hv����ָ��','gsohj--hv����ָ��');
%
% ��ʾ�����㷨������ָ���ƽ��ֵ�;�����
hv_mu1=mean(hv_gso.hviters(2,:),2);
hv_std1=std(hv_gso.hviters(2,:));
fprintf('GSO algorithm: hv_mu==%d\t,hv_std==%d\n',hv_mu1,hv_std1);
%
hv_mu2=mean(hv_gsohj.hviters(2,:),2);
hv_std2=std(hv_gsohj.hviters(2,:));
fprintf('GSOhj algorithm: hv_mu==%d\t,hv_std==%d\n',hv_mu2,hv_std2);

