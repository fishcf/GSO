function [ archiveNew, fvaluesNew, fvaluesAll, archiveAll, fbests, indexs ] = getArchiveProducer( ...
    fvalues, population, fvaluesOld, archiveOld, fvaluesAll, archiveAll, numObjec )
%GETARCHIVE_PRODUCER Summary of this function goes here
%��SampleValue�л����Ѹ���-->producer
%   Input:
%       fvalues   	-	n*objectNum;   ����ΪĿ�����,n is a certain number
%       population	-   NDim*?1;       ��fvalues���Ӧ�ĸ���
%   Output:
%       archiveNew  -   NDim*n          ��ռ�Ž⼯��n is a certain number
%       fvaluesNew  -   n*objectNum     n is a certain number
% 
%GETPRODUCERS for SOB-GSOMP
% get the best f1,f2 with the constraint satisfied
% population:   NDim*popSize
% f1, f2:       1*popSize
% fvalMulti:    numObjec*popSize
%   fbests:     1*2
%   indexs:     1*2
% last modified by zjp in 2011/01/31
global distPop
global fname
distPop = [];
if ~isempty(archiveOld)
    exexutefunction=strcat(fname,'(archiveOld)');
    [fvaluesOld1, fvaluesOld2]=eval(exexutefunction);
    fvaluesOld = [fvaluesOld1; fvaluesOld2];
%     fvaluesOld = costEmission(archiveOld);
    fvaluesOld = fvaluesOld';%????
end
fvaluesAll = [];
archiveAll = [];
% % fvaluesAll = [fvaluesAll; fvaluesOld];
% % archiveAll = [archiveAll archiveOld];
if size(archiveOld,2)<1     %��һ��
    
    %GETPRODUCERS===================================
    fvalMulti = fvalues';
    for i = 1:numObjec
        [fbests(1, i), indexs(1, i)] = min(fvalMulti(i, :));
    end    
    %GETPRODUCERS===================================
    front = paretofront(fvalues);   %������
    archiveNew = population(:, find(front==1));  
    fvaluesNew = fvalues(find(front==1),:);
else    %�ǵ�һ�Σ���ԭ�е�pareto�����������������©pareto��
    
    %GETPRODUCERS===================================
    fvalMulti = fvalues';
    for i = 1:numObjec
        [fbests(1, i), indexs(1, i)] = min(fvalMulti(i, :));
    end    
    %GETPRODUCERS===================================
    fvalues = [fvaluesOld;fvalues];
    front = paretofront(fvalues);    %������

    population = [archiveOld population];
    archiveNew = population(:, find(front==1)); 
    fvaluesNew = fvalues(find(front==1),:);
end