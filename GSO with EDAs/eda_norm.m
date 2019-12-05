function [X1,X2,BestF1,BestF2,iters]=eda_norm(X,Max_Iter)
% EDA_NORM
%
% Purpose:
% 一种基于正态分布的分布估计算法，概率模型中将正态分布模型和余弦函数
% 结合在一起，并且采用遗传算法中交叉的思想，设计了随机交叉环节
%
% Record of revisions:
% Date Programmer Description of change
% ==== =========== ==========================
% 2019/11/28 Yu Chaofan Original code
%
% Define variables:
% N --- 种群个体总数
% m --- 优势种群个数
% mu --- 均值
% sigma --- 标准差
% pc -- 交叉概率
% X --- 种群
% Xbest --- 最优个体
% u1、u2 --- (0,1)中均匀分布的随机数
% a --- [0,1]中分布的随机数
%
% Initialize variables：
Population=X;
% Set Flags to indicate the update times of curtailment(弃风量) and branch margin（线路裕度）
Flag1=0;
Flag2=0;
% 种群数目:
N=30;
% 选择优势种群数目:
% m=N/4;
m=8;
% 交叉概率:
pc=0.1;
iters=0;
%
%% ========================================== BestF1 =====================================
while iters<Max_Iter
    iters=iters+1;
    % ==== 选择优势种群 ====
    % 1.计算每个个体的适应度
    [F1,~]=fun(X);
    % 2.对种群进行排序
    [~,idx]=sort(F1);
    % 3.选取前m个优势个体
    X=X(:,idx(1:m));
    %
    % 构建概率模型：统计前m个个体包含的信息，计算其均值、标准差
    mu=mean(X,2);
    sigma=std(X,0,2);
    %
    % 生成新种群：从构建的概率模型中随机产生 N 个新样本，构成新的种群
    u1=rand(24,30); u2=rand(24,30);
    X=repmat(mu,1,30)+repmat(sigma,1,30).*((-2*log(u1)).^(1/2)).*cos(2*pi*u2);
    [F1,~]=fun(X);
    % 随机交叉：根据交叉概率pc,随机挑选N*pc个个体与最好的个体Xbest交叉
    [BestF1,idx_best]=min(F1);
    Xbest=X(:,idx_best);
    a=unifrnd(0,1,24,30);
    % 随机挑选N*pc个个体
    p=randperm(N*pc);
    Xnew=a.*repmat(Xbest,1,30)+(ones(24,30)-a).*repmat(X(:,p),1,10);
    [Fnew,~]=fun(Xnew);
    % 更新
    if min(Fnew)<BestF1
        X=Xnew;
        Flag1=Flag1+1;
    end
end
if Flag1==0
    [F1,~]=fun(Population);
    [BestF1,idx_best]=min(F1);
    Xbest=X(:,idx_best);
    X1=Xbest;
else
    [F1,~]=fun(X);
    [BestF1,idx_best]=min(F1);
    Xbest=X(:,idx_best);
    X1=Xbest;
end
%
%% ===================================== BestF2 =============================
X=Population;
iters=0;
while iters<Max_Iter
    iters=iters+1;
    % ==== 选择优势种群 ====
    % 1.计算每个个体的适应度
    [~,F2]=fun(X);
    % 2.对种群进行排序
    [~,idx]=sort(F2);
    % 3.选取前m个优势个体
    X=X(:,idx(1:m));
    %
    % 构建概率模型：统计前m个个体包含的信息，计算其均值、标准差
    mu=mean(X,2);
    sigma=std(X,0,2);
    %
    % 生成新种群：从构建的概率模型中随机产生 N 个新样本，构成新的种群
    u1=rand(24,30); u2=rand(24,30);
    X=repmat(mu,1,30)+repmat(sigma,1,30).*((-2*log(u1)).^(1/2)).*cos(2*pi*u2);
    [~,F2]=fun(X);
    % 随机交叉：根据交叉概率pc,随机挑选N*pc个个体与最好的个体Xbest交叉
    [BestF2,idx_best]=min(F2);
    Xbest=X(:,idx_best);
    a=unifrnd(0,1,24,30);
    p=randperm(N*pc);
    Xnew=a.*repmat(Xbest,1,30)+(ones(24,30)-a).*repmat(X(:,p),1,10);
    [~,Fnew]=fun(Xnew);
    % 更新
    if min(Fnew)<BestF2
        X=Xnew;
        Flag2=Flag2+1;
    end
end
if Flag2==0
    [~,F2]=fun(Population);
    [BestF2,idx_best]=min(F2);
    Xbest=X(:,idx_best);
    X2=Xbest;
else
    [~,F2]=fun(X);
    [BestF2,idx_best]=min(F2);
    Xbest=X(:,idx_best);
    X2=Xbest;
end
%
fprintf('Flag1=%d ,Flag2=%d\n',Flag1,Flag2);
%% ===================================== Sub Function =======================
% 适应度函数
    function [f_1,f_2]=fun(x)
        %   Bus Data Format
        %       1   bus number (1 to 29997)
        %       2   bus type
        %               PQ bus          = 1
        %               PV bus          = 2
        %               reference bus   = 3
        %               isolated bus    = 4
        
        %               zip load bus    = 9
        
        %       3   Pd, real power demand (MW)
        %       4   Qd, reactive power demand (MVAR)
        %       5   Gs, shunt conductance (MW (demanded?) at V = 1.0 p.u.)
        %       6   Bs, shunt susceptance (MVAR (injected?) at V = 1.0 p.u.)
        %       7   area number, 1-100
        %       8   Vm, voltage magnitude (p.u.)
        %       9   Va, voltage angle (degrees)
        %   (-)     (bus name)
        %       10  baseKV, base voltage (kV)
        %       11  zone, loss zone (1-999)
        %   (+) 12  maxVm, maximum voltage magnitude (p.u.)
        %   (+) 13  minVm, minimum voltage magnitude (p.u.)
        %
        tmpbus = [
            1   3   0.0     0.0     0.0 0.0 1   1.0855  0.0000  345.0000    1   1.0500  0.9500;
            2   2   21.7    12.7    0.0 0.0 1   1.0653  0.0000  345.0000    1   1.1000  0.9500;
            3   1   2.40    1.20    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            4   1   7.60    1.60    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            5   2   94.2    19      0.0 0.0 1   1.0333  0.0000  345.0000    1   1.1000  0.9500;
            6   1   0.0     0.0     0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            7   1   22.8    10.9    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            8   2   30.0    30.0    0.0 0.0 1   1.0386  0.0000  345.0000    1   1.1000  0.9500;
            9   1   0.0     0.0     0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            10   1   5.80    2.00    0.0 19  1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            11   2   0.0     0.0     0.0 0.0 1   1.0848  0.0000  345.0000    1   1.1000  0.9500;
            12   1   11.2    7.50    0.0 0.1 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            13   2   0.0     0.0     0.0 0.0 1   1.0512  0.0000  345.0000    1   1.1000  0.9500;
            14   1   6.2     1.60    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            15   1   8.2     2.50    0.0 0.1 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            16   1   3.5     1.80    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            17   1   9.0     5.80    0.0 0.1 1   1.0000  0.0000  345.0000    1   1.0500  0.1000;
            18   1   3.2     0.90    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            19   1   9.5     3.40    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            20   1   2.2     0.70    0.0 0.1 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            21   1   17.5    11.2    0.0 0.1 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            22   1   0.0     0.0     0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            23   1   3.20    1.60    0.0 0.1 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            24   1   8.70    6.70    0.0 4.3 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            25   1   0.0     0.0     0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            26   1   3.50    2.30    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            27   1   0.0     0.0     0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            28   1   0.0     0.0     0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            29   1   2.40    0.90    0.0 0.1 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            30   1   10.6    1.90    0.0 0.0 1   1.0000  0.0000  345.0000    1   1.0500  0.9500;
            ];
        
        casename = 'case30test';
        global bus;
        global  windenvironment;   % Num*5的矩阵,一行一次抽样
        global  numberofsample;
        global  possibility;
        if nargin==1
            %注意：一列表示一个个体
            pop_size=size(x,2);
            cur=zeros(numberofsample,pop_size);
            Slim_max=zeros(numberofsample,pop_size);
            % 模拟400次 每次计算20个个体------->计算开销巨大,削减为20
            for i = 1:numberofsample
                windpositions=[2 7 10 16 24];
                bus(windpositions,3)=tmpbus(windpositions,3)-windenvironment(i,:)';
                pforscen=bus(windpositions,3);
                for j=1:pop_size
                    %对于每个个体,每次进行一个样本的计算时,首先要重置为完全消纳的状态
                    bus(windpositions,3)=pforscen;
                    [cur(i,j),Slim_max(i,j)]=CalFitnessnew(casename,x(:,j),i); % 第一个目标
                end
            end
            %===========================================================================
            f_1 = sum(repmat(possibility',1,pop_size).*cur)/numberofsample;
            f_2 = sum(repmat(possibility',1,pop_size).*Slim_max)/numberofsample;
        else
            error('Only one inputs');
        end
    end
end







