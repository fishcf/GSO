function [X1,X2,BestF1,BestF2,iters]=eda_norm(X,Max_Iter)
% EDA_NORM
%
% Purpose:
% һ�ֻ�����̬�ֲ��ķֲ������㷨������ģ���н���̬�ֲ�ģ�ͺ����Һ���
% �����һ�𣬲��Ҳ����Ŵ��㷨�н����˼�룬�����������滷��
%
% Record of revisions:
% Date Programmer Description of change
% ==== =========== ==========================
% 2019/11/28 Yu Chaofan Original code
%
% Define variables:
% N --- ��Ⱥ��������
% m --- ������Ⱥ����
% mu --- ��ֵ
% sigma --- ��׼��
% pc -- �������
% X --- ��Ⱥ
% Xbest --- ���Ÿ���
% u1��u2 --- (0,1)�о��ȷֲ��������
% a --- [0,1]�зֲ��������
%
% Initialize variables��
Population=X;
% Set Flags to indicate the update times of curtailment(������) and branch margin����·ԣ�ȣ�
Flag1=0;
Flag2=0;
% ��Ⱥ��Ŀ:
N=30;
% ѡ��������Ⱥ��Ŀ:
% m=N/4;
m=8;
% �������:
pc=0.1;
iters=0;
%
%% ========================================== BestF1 =====================================
while iters<Max_Iter
    iters=iters+1;
    % ==== ѡ��������Ⱥ ====
    % 1.����ÿ���������Ӧ��
    [F1,~]=fun(X);
    % 2.����Ⱥ��������
    [~,idx]=sort(F1);
    % 3.ѡȡǰm�����Ƹ���
    X=X(:,idx(1:m));
    %
    % ��������ģ�ͣ�ͳ��ǰm�������������Ϣ���������ֵ����׼��
    mu=mean(X,2);
    sigma=std(X,0,2);
    %
    % ��������Ⱥ���ӹ����ĸ���ģ����������� N ���������������µ���Ⱥ
    u1=rand(24,30); u2=rand(24,30);
    X=repmat(mu,1,30)+repmat(sigma,1,30).*((-2*log(u1)).^(1/2)).*cos(2*pi*u2);
    [F1,~]=fun(X);
    % ������棺���ݽ������pc,�����ѡN*pc����������õĸ���Xbest����
    [BestF1,idx_best]=min(F1);
    Xbest=X(:,idx_best);
    a=unifrnd(0,1,24,30);
    % �����ѡN*pc������
    p=randperm(N*pc);
    Xnew=a.*repmat(Xbest,1,30)+(ones(24,30)-a).*repmat(X(:,p),1,10);
    [Fnew,~]=fun(Xnew);
    % ����
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
    % ==== ѡ��������Ⱥ ====
    % 1.����ÿ���������Ӧ��
    [~,F2]=fun(X);
    % 2.����Ⱥ��������
    [~,idx]=sort(F2);
    % 3.ѡȡǰm�����Ƹ���
    X=X(:,idx(1:m));
    %
    % ��������ģ�ͣ�ͳ��ǰm�������������Ϣ���������ֵ����׼��
    mu=mean(X,2);
    sigma=std(X,0,2);
    %
    % ��������Ⱥ���ӹ����ĸ���ģ����������� N ���������������µ���Ⱥ
    u1=rand(24,30); u2=rand(24,30);
    X=repmat(mu,1,30)+repmat(sigma,1,30).*((-2*log(u1)).^(1/2)).*cos(2*pi*u2);
    [~,F2]=fun(X);
    % ������棺���ݽ������pc,�����ѡN*pc����������õĸ���Xbest����
    [BestF2,idx_best]=min(F2);
    Xbest=X(:,idx_best);
    a=unifrnd(0,1,24,30);
    p=randperm(N*pc);
    Xnew=a.*repmat(Xbest,1,30)+(ones(24,30)-a).*repmat(X(:,p),1,10);
    [~,Fnew]=fun(Xnew);
    % ����
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
% ��Ӧ�Ⱥ���
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
        global  windenvironment;   % Num*5�ľ���,һ��һ�γ���
        global  numberofsample;
        global  possibility;
        if nargin==1
            %ע�⣺һ�б�ʾһ������
            pop_size=size(x,2);
            cur=zeros(numberofsample,pop_size);
            Slim_max=zeros(numberofsample,pop_size);
            % ģ��400�� ÿ�μ���20������------->���㿪���޴�,����Ϊ20
            for i = 1:numberofsample
                windpositions=[2 7 10 16 24];
                bus(windpositions,3)=tmpbus(windpositions,3)-windenvironment(i,:)';
                pforscen=bus(windpositions,3);
                for j=1:pop_size
                    %����ÿ������,ÿ�ν���һ�������ļ���ʱ,����Ҫ����Ϊ��ȫ���ɵ�״̬
                    bus(windpositions,3)=pforscen;
                    [cur(i,j),Slim_max(i,j)]=CalFitnessnew(casename,x(:,j),i); % ��һ��Ŀ��
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







