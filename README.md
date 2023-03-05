# Data-Driven-GGM-ET-PHD
Matlab code implementation of a data drieven GGM-ET-PHD filter. A Gaussian Mixtrue Model of car extended measurement model is learned from nuScenes dataset[1]. The original idea is from a series of articles [2-4]. GGM-ET-PHD filter is then used to track multiple extended targets. The simulation scenario is created by MATLAB® Driving Toolbox. The proposed method outperformes the ERHM-PHD filter[5] and is as good as the traditional GGIW-PHD filter[6].
## Reference
[1] H. Caesar et al., “nuScenes: A multimodal dataset for autonomous driving.” arXiv, May 05, 2020. doi: 10.48550/arXiv.1903.11027.

[2] G. Yao, P. Wang, K. Berntorp, H. Mansour, P. Boufounos, and P. V. Orlik, “Extended object tracking with automotive radar using B-spline chained ellipses model,” presented at the ICASSP, IEEE International Conference on Acoustics, Speech and Signal Processing - Proceedings, 2021, vol. 2021-June, pp. 8408–8412. doi: 10.1109/ICASSP39728.2021.9415080.

[3] H. Kaulbersch, J. Honer, and M. Baum, “EM-based Extended Target Tracking with Automotive Radar using Learned Spatial Distribution Models,” in 2019 22th International Conference on Information Fusion (FUSION), Jul. 2019, pp. 1–8. doi: 10.23919/FUSION43075.2019.9011179.

[4] Y. Xia, P. Wang, K. Berntorp, H. Mansour, P. Boufounos, and P. V. Orlik, “Extended object tracking using hierarchical truncation model with partial-view measurements,” presented at the Proceedings of the IEEE Sensor Array and Multichannel Signal Processing Workshop, 2020, vol. 2020-June. doi: 10.1109/SAM48682.2020.9104388.

[5] 李翠芸, 林锦鹏, and 姬红兵, “一种基于椭圆RHM的扩展目标Gamma高斯混合CPHD滤波器,” 控制与决策, vol. 30, no. 9, pp. 1551–1558, 2015, doi: 10.13195/j.kzyjc.2014.0877.

[6] K. Granstrom, C. Lundquist, and O. Orguner, “Extended Target Tracking using a Gaussian-Mixture PHD Filter,” IEEE Trans. Aerosp. Electron. Syst., vol. 48, no. 4, pp. 3268–3286, Oct. 2012, doi: 10.1109/TAES.2012.6324703.
