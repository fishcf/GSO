function archiveNewI = getArchiveNewI(archiveNew, fvaluesNew, tourSize)
% Input:
%   archiveNew    -   NDim*n          ��ռ�Ž⼯��n is a certain number
%   fvaluesNew    -   n*objectNum     n is a certain number
% Outpurt:
%   archiveNewI   -   NDim*1
% last modified by zjp in 2010/08/21
%     global equFlag
    global conFlag
    global epsConst
    global distPop%     tourSize = 5;
    global xPre
        % ����һ��scrounger���ȡһ��pareto����Ϊ�������
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
    % ��������scroungerѡȡŷʽ���������producer��Ϊ�������
    % �����
end