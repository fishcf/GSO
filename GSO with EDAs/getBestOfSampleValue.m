function [ fbestdirectionval, bestdirection ] = getBestOfSampleValue( SampleValue, SamplePosition, produceri )
% for SOB-GSOMP with constraint
% reference: K Deb: NSGA-II
%从SampleValue中获得最佳个体-->producer
%   Input:
%       SampleValue     -	4*objectNum;   4行（4个direction），列数为目标个数
%       SamplePosition  -   NDim*4      
%   Output:
%       fbestdirctionval-   1*objectNum;    最佳个体的目标函数值
%       bestdirection   -   1*1;            最佳个体的位置，即在SampleValue中的第几行，取值范围[1,2,3,4]

% last modified by zjp in 2010/08/21

global epsConst
    front = paretofront(SampleValue);   %4*1    1:belong to the front;  0:not belong to the front
    %% 方案一：简单取第一个front
    bestdirection = find(front==1,1);     %简单取第一个front （non-dominated solution）
    fbestdirectionval = SampleValue(bestdirection, :);
    %% 方案二：
    % 待添加