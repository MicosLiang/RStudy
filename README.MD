# Liang的R语言学习仓库

本仓库储存了本人在R语言学习中写的一些脚本，主要是生信方向。尽量不使用第三方库，以保证兼容性。

欢迎各位同辈前辈批评指正，如有问题可通过issue或邮箱联系。

建议clone到RStudio的工作空间下，因为这是我的工作环境。

R version == 4.1.3

联系邮箱：liangjh21@lzu.edu.cn AND liang@iloveyou.uno
个人网站：www.iloveyou.uno

（如果您觉得这个仓库对您有帮助，能不能帮我点个免费的小星星呢/dog）

------

In this repository I push some scripts that I write when I studying R language, mainly focus on bioinformatics. I do my best to avoid using R packages to improve the stability.

I want your suggestions on my code,  you can also put an issue or e-mail me if you have any questions.

I suggest you clone it into your RStudio workspace, because I do so.

R version == 4.1.3

e-mail address :  liangjh21@lzu.edu.cn AND liang@iloveyou.uno
blog : www.iloveyou.uno

(If this repository helped you, I will be very happy if you give me a star : - D) 

I only have poor Chinglish, sorry for this :(

------



###### mcp_ova_iris

​	一个感知器并用OvA技术拓展到三元分类，以R语言内置iris3数据集中的80%作为训练集进行训练，并对剩下20%进行预测。



###### Needleman-Wunsch

  动态规划的双序列比对



###### Decision-Tree

  一个简单的决策树，可处理离散变量。



###### BP

  一个BP神经网络，可以自定义隐藏层神经元数，可用于分类任务。



###### NJ-buildTree

  NJ（邻接法）构建系统进化树。
  已弃用



###### ploy_fitting

  图像边缘的多边形拟合（未完成）


###### China_chess

  纯R语言实现蒙特卡洛树搜索的象棋程序，执行效率不高

###### five_chess_game

  用了Rcpp的蒙特卡洛树的五子棋，做的比较失败

###### TSP

  TSP问题，写了两种方法，蚁群算法和模拟退火，效果比较好
  
###### sequence_alignment

  双序列比对和多序列比对
  
###### ML-buildTree

  各种方法构建系统发育树