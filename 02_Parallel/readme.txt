testStiffODE.mlx has the general information of the problem and it's the main script. 
Change the value of gSize in the range [5-10] locally and for to parfor to show speed up.


Make sure to start pool prior to running the scripts in order to get accurate timings of computations.

batchODE_Script.mlx can be used to show submitting jobs locally or to a cluster.

What is showcased:

  parfor
  parallel.pool.DataQueue
  write code once and use in multiple parallel environments
  scale speed-up with more hardware


This code is for R2022a or later releases.  

