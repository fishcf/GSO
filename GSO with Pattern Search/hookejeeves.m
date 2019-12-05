function [X1,X2,BestF1,BestF2,Iters] = hookejeeves(N, X1, X2,StepSize, MinStepSize, Eps_Fx, MaxIter1,indexs)
% Function HOOKEJEEVS performs multivariate optimization using the
% Hooke-Jeeves search method.
%
% Input
%
% N - number of variables
% X - array of initial guesses
% StepSize - array of search step sizes
% MinStepSize - array of minimum step sizes
% Eps_Fx - tolerance for difference in successive function values
% MaxIter - maximum number of iterations
% myFx - name of the optimized function
%
% Output
%
% X - array of optimized variables
% BestF - function value at optimum
% Iters - number of iterations
%


%% ====================BestF1=====================================
format long;
Xnew = X1;
BestF1 = function3(Xnew,indexs);
LastBestF1 = 100 * BestF1 + 100;
bGoOn = true;
Iters = 0; Flag1=0; Flag2=0;
while bGoOn
    Iters = Iters + 1;
    if Iters > MaxIter1
        break;
    end
    X1= Xnew; %更新X1
    for i=1:N
        bMoved(i) = 0;
        bGoOn2 = true;
        while bGoOn2
            xx = Xnew(i);
            Xnew(i) = xx + StepSize(i); %正向搜索
            F1= function3(Xnew,indexs);
            if sum(F1)<sum(BestF1) %如果计算值小于原有最优值，则更新
                BestF1 = F1;
                bMoved(i) = 1;
                Flag1=1;
            else %否则反向搜索
                Xnew(i) = xx - StepSize(i);
                F1= function3(Xnew,indexs);
                if sum(F1)<sum(BestF1) %如果计算值小于原有最优值，则更新
                    BestF1 = F1;
                    bMoved(i) = 1;
                    Flag1=1;
                else %否则恢复Xnew的值
                    Xnew(i) = xx;
                    bGoOn2 = false;
                end
            end
        end
    end
    bMadeAnyMove = sum(bMoved);
    if bMadeAnyMove > 0
        DeltaX = Xnew - X1;
        lambda = 1.1;
        %lambda = linsearch(X, N, lambda, DeltaX, myFx);
        Xnew = X1+ lambda * DeltaX;
    end
    BestF1 = function3(Xnew,indexs);
    % reduce the step size for the dimensions that had no moves
    for i=1:N
        if bMoved(i) == 0
            StepSize(i) = StepSize(i) / 2;
        end
    end
    if abs(sum(BestF1) - sum(LastBestF1)) < Eps_Fx
        break
    end
    LastBestF1 = BestF1;
    bStop = true;
    for i=1:N
        if StepSize(i) >= MinStepSize(i)
            bStop = false;
        end
    end
    bGoOn = ~bStop;
end
%% =================================BestF2===============================
Xnew = X2;
BestF2 = function4(Xnew,indexs);
LastBestF2 = 100 * BestF2 + 100;
bGoOn = true;
while bGoOn
    Iters = Iters + 1;
    if Iters > MaxIter1
        break;
    end
    X2= Xnew; % 更新X2
    for i=1:N
        bMoved(i) = 0;
        bGoOn2 = true;
        while bGoOn2
            xx = Xnew(i);
            Xnew(i) = xx + StepSize(i);
            F2= function4(Xnew,indexs);
            if sum(F2)<sum(BestF2)
                BestF2 = F2;
                bMoved(i) = 1;
                Flag2=1;
            else
                Xnew(i) = xx - StepSize(i);
                F2= function4(Xnew,indexs);
                
                if sum(F2)<sum(BestF2)
                    BestF2 = F2;
                    bMoved(i) = 1;
                    Flag2=1;
                else
                    Xnew(i) = xx;
                    bGoOn2 = false;
                end
            end
        end
    end
    bMadeAnyMove = sum(bMoved);
    if bMadeAnyMove > 0
        DeltaX = Xnew - X2;
        lambda = 1.1;
        %     lambda = linsearch(X, N, lambda, DeltaX, myFx);
        Xnew = X2 + lambda * DeltaX;
    end
    BestF2 = function4(Xnew,indexs);
    % reduce the step size for the dimensions that had no moves
    for i=1:N
        if bMoved(i) == 0
            StepSize(i) = StepSize(i) / 2;
        end
    end
    if abs(sum(BestF2) - sum(LastBestF2)) < Eps_Fx
        break
    end
    LastBestF2 = BestF2;
    bStop = true;
    for i=1:N
        if StepSize(i) >= MinStepSize(i)
            bStop = false;
        end
    end
    bGoOn = ~bStop;
end
fprintf('Iters== %d, Flag1==%d, Flag2==%d\n',Iters,Flag1,Flag2);
format; %恢复原格式
end
%% function y = myFxEx(N, X, DeltaX, lambda, myFx)
%
%   X = X + lambda * DeltaX;
%   y = feval(myFx, X, N);
%
% % end
%
% function lambda = linsearch(X, N, lambda, D, myFx)
%
%   MaxIt = 100;
%   Toler = 0.000001;
%
%   iter = 0;
%   bGoOn = true;
%   while bGoOn
%     iter = iter + 1;
%     if iter > MaxIt
%       lambda = 0;
%       break
%     end
%
%     h = 0.01 * (1 + abs(lambda));
%     f0 = myFxEx(N, X, D, lambda, myFx);
%     fp = myFxEx(N, X, D, lambda+h, myFx);
%     fm = myFxEx(N, X, D, lambda-h, myFx);
%     deriv1 = (fp - fm) / 2 / h;
%     deriv2 = (fp - 2 * f0 + fm) / h ^ 2;
%     diff = deriv1 / deriv2;
%     lambda = lambda - diff;
%     if abs(diff) < Toler
%       bGoOn = false;
%     end
%   end