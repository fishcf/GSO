1. 此程序是基于分布估计算法（EDA）改进的GSO算法，主函数是 mian.m

2. 设计思路：
   
（1）首先调用GSO算法进行全局搜索，当GSO算法搜索的结果停滞不前时，调用正态分布的EDA算法进行局部搜索，函数是 eda_norm.m
 
（2）若EDA算法产生的种群适应度优于原种群，则更新种群及适应度，若不及原种群，则恢复为原种群；