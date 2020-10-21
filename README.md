

# Online Convex Optimization Over Erdős-Rényi Random Networks

[^_^]: # (This repository is the official implementation of [Online Convex Optimization Over Erdős-Rényi Random Networks].)


## Requirements

The code was verified on MatlabR2019b, libsvm-3.24. Not sure on other version.

Libsvm can be installed [here](https://www.csie.ntu.edu.tw/~cjlin/libsvm/).



## Usage

Setups on numerical experiments in this paper:


| Index        | Base Graphs | Losses  | Probability |  Time Horizon *T*
| :------: | :------: | :------: | :--------------: | :--------------: |
| 1   |     complete      |      convex      | 0.2 | 1:200
| 2   |     complete      |      strongly convex      | 0.2 | 1:200
| 3(a)   |     complete      |      convex      | 0.1:0.9 | 200
| 3(b)   |     star      |      convex      | 0.1:0.9 | 200
| 3(c)   |     *k*-regular(*k*=5)      |      convex      | 0.1:0.9 | 200
| 4(a)   |     complete      |      strongly convex      | 0.1:0.9 | 200
| 4(b)   |     star      |      strongly convex      | 0.1:0.9 | 200
| 4(c)   |     *k*-regular(*k*=5)      |      strongly convex      | 0.1:0.9 | 200

For experiment 1 and 2, run "body_fat_main.m". For experiment 3 and 4, run "a_influence_probability.m".

We also provide the code in experiment 5, "a_influence_n.m" and "a_influence_d.m", for researching the impact of network size *N* and decision variable dimension *d* on preformance of algorithms respectively.

Setups on experiment 5:

| Index        | Base Graphs | Losses |  Time Horizon *T* | node number *N*  | vector dimension *d*
| :------: | :------: | :------: | :--------------: | :--------------: | :--------------: |
| 5(a)   |     *k*-regular(*k*=3)      |      convex    | 200 | 5:80 | 14
| 5(b)   |     *k*-regular(*k*=3)      |     strongly convex    | 200 | 5:80 | 14
| 5(c)   |     *k*-regular(*k*=3)      |      convex    | 200 | 20 | 5:100
| 5(d)   |     *k*-regular(*k*=3)      |     strongly convex    | 200 | 20 | 5:100

## Data Preparation

Bodyfat dataset can be downloaded from LIBSVM library ( https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/ ).


## Results


Results are presented:


![avatar](https://raw.githubusercontent.com/TJ2020Lab/Online-Convex-Optimization/main/pic/1%262.png)


![avatar](https://raw.githubusercontent.com/TJ2020Lab/Online-Convex-Optimization/main/pic/3.png)


![avatar](https://raw.githubusercontent.com/TJ2020Lab/Online-Convex-Optimization/main/pic/4.png)


![avatar](https://raw.githubusercontent.com/TJ2020Lab/Online-Convex-Optimization/main/pic/5.png)


## License

This code is published under the BSD [license](http://strategic.mit.edu/docs/matlab_networks/license.txt).


## Acknowledgement

This project is modified from code in  [Distributed Online Optimization with Long-Term Constraints](https://arxiv.org/abs/1912.09705). Thanks to Prof. Yuan for providing the prototype code.

Octave-networks-toolbox (https://zenodo.org/record/10778) is used in this code. Thanks to the developers.

## Citation

```
@inproceedings{lei2020online,
  title={Online Convex Optimization Over Erdős-Rényi Random Networks},
  author={Lei, Jinlong and Yi, Peng and Hong, Yiguang and Chen, Jie and Shi, Guodong},
  booktitle={NeurIPS},
  year={2020}
}
```
