function [ fbestvals, bestmembers, archiveNew, fvaluesNew, fvaluesAll, archiveAll,hv,hvcounts] = GSOMP_2(fname, ...
    NDim,MaxIter,flagDirec,numObjec,initProducer,rangersPercent,pursuitAngleCoefficient, ...
    turingAngleCoefficient,lmaxCoefficient,initAngle,aCoefficient,bCoefficient,direcDul,c1,c2,NumCoeffi,numShift1,hvcounts)
% Group Search Optimizer with Multiple Producers for multiobject optimization in Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified in 2019-11-21 by Chaofan YU, Guangxi University
% 1. only use the GSO algorithm to search the Pareto front
% 2. add HVCOUNTS to show the number of calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified in 2013-05-20 by Junpeng Zhan, Zhejiang University
% 1. the equations' mistakes have been corrected by zjp
% 2. add a function: getPopInspace() to keep particles in variable space
% 3. add inequality constraint handling (2010-10-23) ---- search the words (3 places): considering the inequality constraints
% 4. circshift(distance.*direction, -unidrnd(NumCoeffi0+floor(iteration/25)) (2010-10-26)---- search the words: floor(iteration numShift)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
%   fname       - the name of the evaluation .m function
%   NDim        - dimension of the evalation function
%   MaxIter     - maximum iteration
%   numObjec - number of producers
% Output Arguments:
%       archiveNew  -   NDim*n          非占优解集，n is a certain number
%       fvaluesNew  -   n*objectNum     n is a certain number

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

global bus;
%保留的典型场景数,初始实质抽了400个样本
bus = [
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
global numberofsample;
numberofsample=20;

global  windenvironment;   % wind environment in each generartion for MOGSO
global  possibility;
%% secnario reduction产生风速误差场景--------->详见函数scenarioReduction函数,返回5*Num的矩阵Wind_P和Wind_Q以及1*Num的概率
% windspeed_mean=[4.3 7.5 5.6 8.9 6.5];------->当做高斯分布的参考值
% windspeed_sd=0.08*windspeed_mean; % 风速预测均方差（标准差）
% mu=windspeed_mean;
% sigma=windspeed_sd;
% windpositions=[2 7 10 16 24];

load('wind_data'); %内含 Posibility,Wind_P,Wind_Q,final,initial
clear final initial
windenvironment= Wind_P';
possibility=Posibility;

    global PopSize
    global finalFlag
    global timeHorizon
% flagDirec = 0;      %1则用修改后的，0用原有的坐标转换方程
% direcDul = 50;      %常用zjp  100则前10几代收敛比较快
% direcDul = 10;      %常用zjp  100则前10几代收敛比较快
numShift = Inf;
% NumCoeffi = 1;
circUpDown = 1;    %-1则向下，1则向上循环移动
    fvaluesNewIter = struct('fvaluesNewIter',{1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10, ...
        1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10, ...
        1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,101,102});
    fiter = 1;
    lamda = zeros(1, timeHorizon);
    finalFlag = 0;
    outFlagCount = 0;
    flag=0;
    iteration = 0;
    tourSize = 2;
    NumCoeffi0 = 0;
% PopSize=46;     % population of members
% PopSize=98;     % population of members

    %=====initAngle======================================================
    angle=initAngle.*ones(NDim-1,PopSize);   % Initialize head angle
    %====================================================================
    leftangle=angle; 
    rightangle=angle;

% Bound=eval(fname);
    Bound=[20,80;15,50;10,35;10,30;12,40;0.950000000000000,1.05000000000000;0.950000000000000,1.10000000000000;0.950000000000000,1.10000000000000;0.950000000000000,1.10000000000000;0.950000000000000,1.10000000000000;0.950000000000000,1.10000000000000;1,17;1,17;1,17;1,17;0,6;0,6;0,6;0,6;0,6;0,6;0,6;0,6;0,6;];
% Defined lower bound and upper bound.
    LowerBound = zeros(NDim,PopSize);
    UpperBound = zeros(NDim,PopSize);   
    for i=1:PopSize
        LowerBound(:,i)=Bound(:,1);
        UpperBound(:,i)=Bound(:,2);
    end

    DResult = 1e-1;    % Desired results (never used)
    population =  rand(NDim, PopSize).*(UpperBound-LowerBound) + LowerBound;     % Initialize swarm population
    % initial the population according to each generation unit's maximal generation capacity
    if ~isempty(initProducer)
        if size(initProducer, 2)>PopSize
            population = initProducer(:, 1:PopSize);
        else
            population(:,1:size(initProducer,2)) = initProducer;   %getInitProducer('ZDT1',30);
        end
    end
    % population(1,:) = 0:0.01:0.97;%===================
    vmax = ones(NDim,PopSize);
    for i=1:NDim
        vmax(i,:)=(UpperBound(i,:)-LowerBound(i,:));
    end
    r = sqrt( sum(vmax(:,1).^2) );
    clear vmax
    %r=norm(vmax(:,1));%distance l
    %=====lmaxCoefficient================================================
    r = r*lmaxCoefficient; 
    r=r/2;
    %====================================================================
    distance=r*repmat(ones(1,PopSize),NDim,1);

    %a= round(((NDim+1)^.5))*ones(1,PopSize);
    a = aCoefficient * round( ((NDim+1)^.5) );
    b = bCoefficient * round( ((NDim+1)^.5) );
%% Polar2Cartesian_coordinate_transformation 极坐标-笛卡尔坐标变换
if flagDirec == 1                
    i = NDim - direcDul;
    j = 1;
    angleRep = repmat(angle(:,j), 1, direcDul);
    angleTri = angleRep(end-direcDul+1:end, :);
    % L = tril(X,k) returns the elements on and below the kth diagonal of X. k = 0 is 
    % the main diagonal, k > 0 is above the main diagonal, and k < 0 is below the main diagonal.
    angleTril = angleTri.*(tril(ones(direcDul, direcDul), 0)); %转为下三角型矩阵
    % If A is a nonempty matrix, then prod(A) treats the columns of A as vectors 
    % and returns a row vector of the products of each column.
    direction(i:NDim-1, j) = sin(angle(i-1:NDim-2, j)) .* prod(cos(angleTril))';
    direction(NDim, j) = sin(angle(NDim-1, j));
    direction = repmat(direction, 1, PopSize);
else
    for j=1:PopSize
        direction(1,j)=prod(cos(angle(1:NDim-1,j)));
        for i=2:NDim-1
            direction(i,j)=sin(angle(i-1,j)).*prod(cos(angle(i:NDim-1,j)));
        end
        direction(NDim,j)=sin(angle(NDim-1,j));
    end %equation(1)
end

    exexutefunction=strcat(fname,'(population)');
    
    % Prevent members from flying outside search space
    OutFlag = population<LowerBound | population>UpperBound;
    % population = population - OutFlag.*distance.*direction;  
    population = population - OutFlag.*circshift(distance.*direction,0);  
    [ population, outFlagCount ] = getPopInspace( population, LowerBound, UpperBound, outFlagCount);
    
    %Evaluate initial population
    [fvalue1, fvalue2] = eval(exexutefunction);
    fvalMulti = [fvalue1; fvalue2];
% ===================== getArchive =================================================
    archiveOld = [];
    fvaluesOld = [];
    fvaluesAll = [];
    archiveAll = [];
    % [ archiveNew, fvaluesNew, fvaluesAll ] = getArchive( [fvalue1' fvalue2'], population, fvaluesOld, archiveOld, fvaluesAll );
    % [ archiveNew, fvaluesNew, fvaluesAll, archiveAll ] = getArchive( fvalMulti', population, fvaluesOld, archiveOld, fvaluesAll, archiveAll );
    % zjp 2011-01-31
    [ archiveNew, fvaluesNew, fvaluesAll, archiveAll, fbestvals(1,1:numObjec) indexs(1,1:numObjec) ] = ...
        getArchiveProducer( fvalMulti', population, fvaluesOld, archiveOld, fvaluesAll, archiveAll, numObjec );
    archiveOld = archiveNew;
    fvaluesOld = fvaluesNew;
    % Finding best member in initial population
    bestmembers = population(:,indexs);
    oldangle=angle;
    oldindexs=indexs;    % oldindex=index;
    badcounters = zeros(1,PopSize);    % badcounter=0;

    while(flag == 0) && (iteration < MaxIter)
        % rand('state',sum(100*clock))
        iteration = iteration +1;
        NumCoeffi0 = NumCoeffi0 + NumCoeffi;
        numShift = 1;  
        produceri = 0;
         for j=1:PopSize
            R1 = randn(1); % randn(n) returns an n-by-n matrix of normally distributed random numbers.
            R2 = rand(NDim-1,1);
            R3 = rand(NDim, 1);
            R4 = rand(NDim, 1);
    % ============ numObjec ==========================================
            if ~isempty(find(indexs==j))
    %       if j==index % Stop and search around
    % ===========================================================
                produceri = produceri + 1;
                SamplePosition=[];
                SampleAngle=[];
                SampleValue=[];
                %========= pursuitAngleCoefficient =======================================================
                leftangle=pursuitAngleCoefficient*(pi/(a^2)/2).*R2+angle(:,j);  %trun left, the angle range between new and old direction is \in [0,pi/a^2]
                rightangle=-pursuitAngleCoefficient*(pi/(a^2)/2).*R2+angle(:,j); %trun right, the angle range between new and old direction is \in [0,pi/a^2]
                %===========================================================================================
            % distance(:,j)=r*R1;
            distance(:,j)=r*randn(size(distance, 1), 1);
            if flagDirec == 1 % 坐标变换
                i = NDim - direcDul;
                angleRep = repmat(angle(:,j), 1, direcDul);
                angleTri = angleRep(end-direcDul+1:end, :);
                angleTril = angleTri.*(tril(ones(direcDul, direcDul), 0));
                direction(i:NDim-1, j) = sin(angle(i-1:NDim-2, j)) .* prod(cos(angleTril))';
                direction(NDim, j) = sin(angle(NDim-1, j));
            else
                direction(1,j)=prod(cos(angle(1:NDim-1,j)));
                for i=2:NDim-1
                    direction(i,j)=sin(angle(i-1,j)).*prod(cos(angle(i:NDim-1,j)));
                end
                direction(NDim,j)=sin(angle(NDim-1,j));
            end
    % ==========================numObjec====================================
                % if badcounter>=a        % If the producer can not find a better area after $a$ iterations 
                if badcounters(1,j)>=a    % If the producer can not find a better area after $a$ iterations 
                    angle(:,j)=oldangle(:,j); % it will turn its head back to zero degree
                end
                SamplePosition = [SamplePosition,population(:,j)];
                SampleAngle = [SampleAngle,angle(:,j)];
                SampleValue = [SampleValue; fvalMulti(:, j)'];
                %SampleValue = [SampleValue; [fvalue1(j),fvalue2(j)] ];  %SampleValue = [SampleValue,fvalue(j)];

           %% Look Straight
            % StraightPosition=population(:,j)+distance(:,j).*direction(:,j);
            % Unidrnd()generates random numbers for the discrete uniform distribution with maximum N.
            StraightPosition=population(:,j)+circshift(distance(:,j).*direction(:,j),-unidrnd(NumCoeffi0+floor(iteration/numShift)));
            Outflag=(StraightPosition>UpperBound(:,j) | StraightPosition<LowerBound(:,j));
            % StraightPosition=StraightPosition-Outflag.*distance(:,j).*direction(:,j);
            StraightPosition=StraightPosition-Outflag.*circshift(distance(:,j).*direction(:,j),-unidrnd(NumCoeffi0+floor(iteration/numShift)));
            [ StraightPosition, outFlagCount ] = getPopInspace( StraightPosition, LowerBound(:,j), UpperBound(:,j), outFlagCount);
            Outflag=(StraightPosition>UpperBound(:,j) | StraightPosition<LowerBound(:,j));
           if Outflag ==0
           else 
               [ StraightPosition, outFlagCount ] = getPopInspace( StraightPosition, LowerBound(:,j), UpperBound(:,j), outFlagCount);
           end
                Straightfunction=strcat(fname,'(StraightPosition)');
                Straightfvalues = eval(Straightfunction);
                [ Straightfvalue1, Straightfvalue2 ] = eval(Straightfunction);
                %Straightfvalues = [ Straightfvalue1, Straightfvalue2 ];
                SamplePosition = [SamplePosition,StraightPosition];
                SampleAngle = [SampleAngle,angle(:,j)];
                %SampleValue = [SampleValue; Straightfvalues];
                SampleValue = [SampleValue; [Straightfvalue1, Straightfvalue2] ];

            %% Look left
            if flagDirec == 1
                i = NDim - direcDul;
                angleRep = repmat(leftangle(:,1), 1, direcDul);
                angleTri = angleRep(end-direcDul+1:end, :);
                angleTril = angleTri.*(tril(ones(direcDul, direcDul), 0));
                direction(i:NDim-1, j) = sin(leftangle(i-1:NDim-2, 1)) .* prod(cos(angleTril))';
                direction(NDim, j) = sin(leftangle(NDim-1, 1));
            else
                direction(1,j)=prod(cos(leftangle(1:NDim-1)));
                for i=2:NDim-1
                    direction(i,j)=sin(leftangle(i-1)).*prod(cos(leftangle(i:NDim-1)));
                end
                direction(NDim,j)=sin(leftangle(NDim-1));
            end

            %LeftPosition=population(:,j)+distance(:,j).*direction(:,j);
            LeftPosition=population(:,j)+circshift(distance(:,j).*direction(:,j),-unidrnd(NumCoeffi0+floor(iteration/numShift)));

            Outflag=(LeftPosition>UpperBound(:,j) | LeftPosition<LowerBound(:,j));
            %LeftPosition=LeftPosition-Outflag.*distance(:,j).*direction(:,j);
            LeftPosition=LeftPosition-Outflag.*circshift(distance(:,j).*direction(:,j),-unidrnd(NumCoeffi0+floor(iteration/numShift)));
                
            [ LeftPosition, outFlagCount ] = getPopInspace( LeftPosition, LowerBound(:,j), UpperBound(:,j), outFlagCount);
Outflag=(LeftPosition>UpperBound(:,j) | LeftPosition<LowerBound(:,j));
if Outflag ==0
else 
    [ LeftPosition, outFlagCount ] = getPopInspace( LeftPosition, LowerBound(:,j), UpperBound(:,j), outFlagCount);
end
                Leftfunction=strcat(fname,'(LeftPosition)');
                %Leftfalues = eval(Leftfunction);
                %实质上是调用了function2(Leftposition)
                [ Leftfvalue1, Leftfvalue2 ] = eval(Leftfunction);
                SamplePosition = [SamplePosition,LeftPosition];
                SampleAngle = [SampleAngle,leftangle(:)];%去掉(:)===========================
                % SampleValue = [SampleValue; Leftfalues'];
                SampleValue = [SampleValue; [Leftfvalue1, Leftfvalue2] ];
                [ LeftPosition, outFlagCount ] = getPopInspace( LeftPosition, LowerBound(:,j), UpperBound(:,j), outFlagCount);

           %% Look right
            if flagDirec == 1
                i = NDim - direcDul;
                angleRep = repmat(rightangle(:,1), 1, direcDul);
                angleTri = angleRep(end-direcDul+1:end, :);
                angleTril = angleTri.*(tril(ones(direcDul, direcDul), 0));
                direction(i:NDim-1, j) = sin(rightangle(i-1:NDim-2, 1)) .* prod(cos(angleTril))';
                direction(NDim, j) = sin(rightangle(NDim-1, 1));
            else
                direction(1,j)=prod(cos(rightangle(1:NDim-1)));
                for i=2:NDim-1
                   direction(i,j)=sin(rightangle(i-1)).*prod(cos(rightangle(i:NDim-1)));
                end
                direction(NDim,j)=sin(rightangle(NDim-1));
            end
            % RightPosition=population(:,j)+distance(:,j).*direction(:,j);
            RightPosition=population(:,j)+circshift(distance(:,j).*direction(:,j),-unidrnd(NumCoeffi0+floor(iteration/numShift)));
            Outflag=(RightPosition>UpperBound(:,j) | RightPosition<LowerBound(:,j));
            % RightPosition=RightPosition-Outflag.*distance(:,j).*direction(:,j);
            RightPosition=RightPosition-Outflag.*circshift(distance(:,j).*direction(:,j),-unidrnd(NumCoeffi0+floor(iteration/numShift)));
            [ RightPosition, outFlagCount ] = getPopInspace( RightPosition, LowerBound(:,j), UpperBound(:,j), outFlagCount);
            Outflag=(RightPosition>UpperBound(:,j) | RightPosition<LowerBound(:,j));
if Outflag ==0
else
    [ RightPosition, outFlagCount ] = getPopInspace( RightPosition, LowerBound(:,j), UpperBound(:,j), outFlagCount);
end
                Rightfunction=strcat(fname,'(RightPosition)');
                %Rightfvalues = eval(Rightfunction);
                [ Rightfvalue1, Rightfvalue2 ] = eval(Rightfunction);
                SamplePosition = [SamplePosition,RightPosition];
                SampleAngle = [SampleAngle,rightangle(:)];
                %SampleValue = [SampleValue; Rightfvalues'];
                SampleValue = [SampleValue; [Rightfvalue1, Rightfvalue2] ];
                % sample 3 points
                %====================================getBestOfSampleValue==========================================
                [fbestdirectionval, bestdirection] = getBestOfSampleValue(SampleValue, SamplePosition,produceri);
                % [fbestdirectionval, bestdirection]=min(SampleValue);
                population(:,j)=SamplePosition(:,bestdirection); 
    % ==========================numObjec====================================
                if bestdirection ~= 1   % if the member find a better place
                    angle(:,j)=SampleAngle(:,bestdirection);
                    oldangle(:,j)=angle(:,j);
                    badcounters(1,j) = 0;
                % badcounter=0;
                else                    % if the member stays
                    badcounters(1,j) = badcounters(1,j)+1;
                % badcounter=badcounter+1;
                    %==============turingAngleCoefficient=======================================================
                    angle(:,j) =  turingAngleCoefficient*pi/(a^2).*R2+angle(:,j); % Turn pi/(a^2)/2 and sample a new direction
                    %===========================================================================================                
                end
            else
                %===================turingAngleCoefficient==================================================
                angle(1:NDim-1,j) = turingAngleCoefficient*pi/(a^2).*R2+angle(1:NDim-1,j);           
                
                %==========rangersPercent===================================================================
                if rand(1)>rangersPercent   
%% ====================== Scroungers perform scrounging %equation(7) ===========================
                    % ==========================numObjec=bestmember? numObjec与目标函数个数相等========================                
                    % archiveNewI = archiveNew(:,getArchiveI(size(archiveNew,2)));
                    archiveNewI = getArchiveNewI(archiveNew, fvaluesNew, tourSize);
                    %R4 = 0;
                    % perform scrounging
                    distance(:,j) = c1 * R3.*(bestmembers(:,geti(numObjec)) - population(:,j)) + ...
                                    c2 * R4.*(archiveNewI - population(:,j));
                    population(:,j) = population(:,j) + distance(:,j);    
                else
%% ========================= Rangers perform ranging ==========================================
                    % distance(:,j)=r*repmat(b*R1,NDim,1);%==============b=====================
                    distance(:,j)=r*b*randn(NDim,1);
                    % direction calculation                    
                if flagDirec == 1
                    i = NDim - direcDul;
                    angleRep = repmat(angle(:,j), 1, direcDul);
                    angleTri = angleRep(end-direcDul+1:end, :);
                    angleTril = angleTri.*(tril(ones(direcDul, direcDul), 0));
                    direction(i:NDim-1, j) = sin(angle(i-1:NDim-2, j)) .* prod(cos(angleTril))';
                    direction(NDim, j) = sin(angle(NDim-1, j));
                else
                    direction(1,j)=prod(cos(angle(1:NDim-1,j)));
                    for i=2:NDim-1
                       direction(i,j)=sin(angle(i-1,j)).*prod(cos(angle(i:NDim-1,j)));
                    end
                    direction(NDim,j)=sin(angle(NDim-1,j));
                end
                %population(:,j) = population(:,j) + distance(:,j).*direction(:,j);
                population(:,j) = population(:,j) + circshift(distance(:,j).*direction(:,j),-unidrnd(NumCoeffi0+floor(iteration/numShift)));
                end
            end
         end
        %fprintf('%3.0f  ',iteration); % for tic,toc
        % 求余数 r = rem(a,b) returns the remainder after division of a by b, where a is the dividend and b is the divisor. 
        if rem(iteration,2)==0    fprintf('%3.0f  ',iteration);    end   % for tic,toc
        if iteration==50          fprintf('\n');                   end
        % Prevent members from flying outside search space
        OutFlag = population<LowerBound | population>UpperBound;
        %population = population - OutFlag.*distance.*direction;
        population = population - OutFlag.*circshift(distance.*direction,-unidrnd(NumCoeffi0+floor(iteration/numShift)));
    	[ population, outFlagCount ] = getPopInspace( population, LowerBound, UpperBound, outFlagCount);
        
        % Evaluate the new swarm
        fvalMulti = eval(exexutefunction);
        [fvalue1, fvalue2] = eval(exexutefunction);     %fvalue = eval(exexutefunction);
        fvalMulti = [fvalue1; fvalue2];
        %=========getArchive=================================================
        % zjp 2011-01-31
        [ archiveNew, fvaluesNew, fvaluesAll, archiveAll, fbestvals(1,1:numObjec) indexs(1,1:numObjec) ] = ...
            getArchiveProducer( fvalMulti', population, fvaluesOld, archiveOld, fvaluesAll, archiveAll, numObjec );
        archiveOld = archiveNew;
        fvaluesOld = fvaluesNew;

%       ==============save for plot convergence=====================================
%         if mod(iteration, 10) == 0
            fvaluesNewIter(fiter).fvaluesNewIter = fvaluesNew;
            fvaluesNewIter(fiter).archiveNewIter = archiveNew;%======================
            fiter = fiter + 1;
%         end
        bestmembers = population(:,indexs);
        Bests(iteration,:) = fbestvals(1,1:numObjec);
        %P=[300,0.8,0.8];
        P=[300,0.8]; % 二维参考点
        [hv(iteration,:),hvcounts]= hvolume2d_hv(fvaluesNew,P,hvcounts); % hv指标
%% =====================hookejeeves模式搜索(单独GSO搜索中不使用)============================================
%         if iteration>6
%             dev1=0;
%             dev2=0;
%             for I=(iteration-4):iteration
%                 dev1=sum(Bests(I,1)+dev1);
%                 dev2=sum(Bests(I,2)+dev2); 
%             end
%             dev1=dev1/5;
%             dev2=dev2/5;
%             dev1=abs(dev1-Bests(iteration,1));
%             dev2=abs(dev2-Bests(iteration,2));
%             if dev1<0.0001 && dev2<0.0001
%             %iteration>6 && Bests(iteration,1)==Bests((iteration-1),1) && Bests(iteration,1)==Bests((iteration-2),1) && Bests(iteration,1)==Bests((iteration-3),1) && Bests(iteration,1)==Bests((iteration-4),1) && Bests(iteration,1)==Bests((iteration-5),1) && Bests(iteration,2)==Bests((iteration-1),2)&& Bests(iteration,2)==Bests((iteration-2),2) && Bests(iteration,2)==Bests((iteration-3),2) && Bests(iteration,2)==Bests((iteration-4),2) && Bests(iteration,2)==Bests((iteration-5),2)
%                N=24;
%                X1= population(:,indexs(1,1)); %Pareto最优解1
%                X2= population(:,indexs(1,2)); %Parero最优解2
%                StepSize=ones(24,1).*0.1;
%                MinStepSize=0.01.*ones(24,1);
%                Eps_Fx=0.01;
%                MaxIter1=200;
%                [X1,X2,BestF1,BestF2,Iters] = hookejeeves(N, X1, X2,StepSize, MinStepSize, Eps_Fx, MaxIter1,indexs);
%                bestmembers=[X1,X2];
%                %Bests(iteration,:)=[BestF1,BestF2];
%                population(:,indexs(1,1))=X1;
%                population(:,indexs(1,2))=X2;
%                [fvalue1, fvalue2] = eval(exexutefunction);     %fvalue = eval(exexutefunction);
%                fvalMulti = [fvalue1; fvalue2];
%                [ archiveNew, fvaluesNew, fvaluesAll, archiveAll, fbestvals(1,1:numObjec),indexs(1,1:numObjec) ] = ...
%                      getArchiveProducer( fvalMulti', population, fvaluesOld, archiveOld, fvaluesAll, archiveAll, numObjec );
%                 Bests(iteration,:) = fbestvals(1,1:numObjec);
%                 [hv(iteration,:),hvcounts]= hvolume2d_hv(fvaluesNew,P,hvcounts);
%             end
%         end
    end
end
%% scrounger following producer
function i = geti(num)
    % 方案一：scrounger随机取一个producer作为跟随对象
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
    % 方案二：scrounger选取欧式距离最近的producer作为跟随对象
    % 待添加
end
%% get members in search space
function [ population, outFlagCount ] = getPopInspace( population, LowerBound, UpperBound, outFlagCount)
%getPopInspace 
% Prevent members from flying outside search space
    OutFlag = population<LowerBound | population>UpperBound;
    aa = 1;
    if ( OutFlag==0 )
    else
        if aa == -1      %不做处理
        elseif aa == 0   %方案一
        [ind1,ind2] = find(population>UpperBound);
        [ind3,ind4] = find(population<LowerBound);
        population = max(population, LowerBound);
        population = min(population, UpperBound);
        outFlagCount = outFlagCount + 1;
        elseif aa == 1  %方案二
            population = population - LowerBound;
            rrem = rem(population,UpperBound-LowerBound); %rem(a,b) returns the remainder after division of a by b, where a is the dividend and b is the divisor
            outFlagL = population < 0;
            % population = population + LowerBound + 2*outFlagL.*(-rrem);
            population = population + LowerBound + outFlagL.*(-rrem) - population.*outFlagL;
            outFlagU = population > UpperBound;
            population = population - outFlagU.*(rrem) - population.*outFlagU + UpperBound.*outFlagU;
        end
    %     population
    %     pause
    end
end