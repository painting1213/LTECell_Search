# LTE Cell Research V1
* Step 1: data pre process
* Step 2: PSS Sync
* Step 3: SSS Sync
*********
# 程序思路
## 数据预处理
1. 打开文件，文件有两列数据：I路，Q路
2. ss_cache = 对数据1/16重采样
3. ss_cache_f = 对sscache进行滤波(卷积)，给定filter_coffe系数，conv_same
4. ss_cache_d = 对卷积后的结果进行1/2抽取，
**********
## PSS同步
1. ss_cs = sscache_d(1: n_samps_per_subframe*n_subframes_cs)    
2. corr = find_pss(ss_cs, 子帧数量=6，每个子帧的采样点1920/2)，5696*3 matrix
3. 找corr三列每一列最大的绝对值val0 val1 val2，和该最大值的索引值idx0 idx1 idx2
4. Nid_2 = 最大值所在列值i
5. peak_idx = pss_peak_idx * 2
***********
## SSS同步
1. find_sss=(ss_cache, Nid_2, peak_idx)
2. 得到Nid_1, sf
***********


Data: Feb.26 2020
***********