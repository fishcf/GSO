function [ fbestdirectionval, bestdirection ] = getBestOfSampleValue( SampleValue, SamplePosition, produceri )
% for SOB-GSOMP with constraint
% reference: K Deb: NSGA-II
%��SampleValue�л����Ѹ���-->producer
%   Input:
%       SampleValue     -	4*objectNum;   4�У�4��direction��������ΪĿ�����
%       SamplePosition  -   NDim*4      
%   Output:
%       fbestdirctionval-   1*objectNum;    ��Ѹ����Ŀ�꺯��ֵ
%       bestdirection   -   1*1;            ��Ѹ����λ�ã�����SampleValue�еĵڼ��У�ȡֵ��Χ[1,2,3,4]

% last modified by zjp in 2010/08/21

global epsConst
    front = paretofront(SampleValue);   %4*1    1:belong to the front;  0:not belong to the front
    %% ����һ����ȡ��һ��front
    bestdirection = find(front==1,1);     %��ȡ��һ��front ��non-dominated solution��
    fbestdirectionval = SampleValue(bestdirection, :);
    %% ��������
    % �����