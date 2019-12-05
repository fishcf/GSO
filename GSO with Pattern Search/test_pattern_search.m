% Script file: test_pattern_search.m
%
% Purpose:
% 1.pattern search ���������ķ�ά���㣨��ͬ��ά�����ò�ͬ����������������
% 2.GSO �� pattern search����ϣ�GSO��������Ϊ����ʱpattern search����Ч����ѣ�
% Record of revisions:
% Date         Programmer    Description of change
% ==========   ==========    =====================
% 2019-11-23   Chaofan Yu    Original code
%
%% ===================================== ��ά��������ģʽ���� =======================================
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
% ��������
maxIter=50;  
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
    [ fbestvals, bestmembers, archiveNew, fvaluesNew, fvaluesAll, archiveAll,hv] = GSOMP_2_divide_dimention(fname,NDim,maxIter,flagDirec,numObjec,initProducer,rangersPercent,pursuitAngleCoefficient,turingAngleCoefficient,lmaxCoefficient,initAngle,aCoefficient,bCoefficient,direcDul,c1,c2,NumCoeffi,numShift1);
    fprintf('\n');
    % ====================================
    figure(1)
    fvaluesNew(:,2)=fvaluesNew(:,2)*10000+1450;
    size(fvaluesNew);
    plot( fvaluesNew(:,1), fvaluesNew(:,2), '.r' );
    %fvaluesNew(:,3)=10^5-fvaluesNew(:,3);
    grid on
    title('FON--Pareto Front Obtained by Divide-Step-GSOMP');
    xlabel('������'); ylabel('��·ԣ��');
end
hold on;
figure(2)
plot(hv);
title('Hypervolume����');
xlabel('��������'); ylabel('hvָ��');





