此程序旨在测试 "GSO算法" 求解风电调度问题的性能，主函数为 main.m
    1.目的：测试"单独使用GSO算法"的搜索速度（评价指标：hv评价函数的评价次数和指标大小（越大越好））
    2.特点：定义变量 hviters 来保存每次迭代的hv指标
          hviters(1,i)=hvcounts;  ------------------------ hv函数评级次数
          hviters(2,i)=hv(end,1); ------------------------ hv指标
    3.改变风电接入百分比的方法：
          GSOMP_3.m 文件的第97行 windenvironment= Wind_P';
          通过在Wind_P’前乘以一个 风电接入系数  可以改变 windenvironment
