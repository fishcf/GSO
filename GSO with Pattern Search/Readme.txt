此程序旨在探究GSO算法与pattern search算法的配合使用，没有main.m文件，而是从以下两个文件来测试算法性能： 
                 
1.test_GSO_GSOhj.m
  测试“单独使用GSO”和“GSO+hookejeeves”算法的搜索速度（评价指标：hv评价函数的指标大小和评估次数）

2.test_pattern_search
  （1）pattern search 搜索步长的分维计算（不同的维度设置不同的搜索步长向量）
  （2）GSO 与 pattern search的配合（GSO迭代次数为多少时pattern search介入效果最佳）

