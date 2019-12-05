function  [f_1,f_2] = function2(x)
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
    %run('C:\Users\YUANZHENG\Desktop\GSO_OPF\case30test.m');

    
    %注意：一列表示一个个体
    pop_size=size(x,2);
    cur=zeros(numberofsample,pop_size);
    Slim_max=zeros(numberofsample,pop_size);
    %模拟400次 每次计算30个个体------->计算开销巨大,削减为30
    for i = 1:numberofsample       
       windpositions=[2 7 10 16 24];
       bus(windpositions,3)=tmpbus(windpositions,3)-windenvironment(i,:)'; 
       pforscen=bus(windpositions,3);
       %加入了无功的处理
%        tanfai=0.3287;
%        bus(windpositions,4)=bus(windpositions,4)-windenvironment(i,:)'*tanfai;

       for j=1:pop_size
           %对于每个个体,每次进行一个样本的计算时,首先要重置为完全消纳的状态
           bus(windpositions,3)=pforscen;
%            disp(bus(windpositions,3))
           [cur(i,j),Slim_max(i,j)]=CalFitnessnew(casename,x(:,j),i); % 第一个目标
           
       
       
       end
       
    end
    %===========================================================================

    f_1 = sum(repmat(possibility',1,pop_size).*cur)/numberofsample;
    f_2 = sum(repmat(possibility',1,pop_size).*Slim_max)/numberofsample;
    
else
    error('Only one inputs');
end
% save('x');
end
   
