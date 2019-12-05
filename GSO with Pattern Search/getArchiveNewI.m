function archiveNewI = getArchiveNewI(archiveNew, fvaluesNew, tourSize)
% Input:
%   archiveNew    -   NDim*n          非占优解集，n is a certain number
%   fvaluesNew    -   n*objectNum     n is a certain number
% Outpurt:
%   archiveNewI   -   NDim*1
% last modified by zjp in 2010/08/21
%     global equFlag
    global conFlag
    global epsConst
    global distPop%     tourSize = 5;
    global xPre
        % 方案一：scrounger随机取一个pareto解作为跟随对象
        num = size(archiveNew, 2);
        r = rand(1);
        bd = [0 1/num];
        i = 1;
        for ite = 1:10
            if (r>=bd(1) & r<bd(2))
                break;            
            else
                bd = bd+1/num;
            end
            i = i+1;
        end
        archiveNewI = archiveNew(:,i);
    % 方案二：scrounger选取欧式距离最近的producer作为跟随对象
    % 待添加
end