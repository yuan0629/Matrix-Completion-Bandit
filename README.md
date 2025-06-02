# Matrix-Completion-Bandit
Code to reproduce simulation results, figures, and real data analysis results for the paper "Online Policy Learning and Inference by Matrix Completion" by Congyuan Duan, Jingyang Li, and Dong Xia.


## Organization

### Simulation
Simulation-Inference.R produces the empirical distributions under different settings in Figures 2 and 4-5 in Section 5, and Figures 9-11 in Appendix A.1.

Simulation-Regret.R produces the empirical cumulative regret under $\gamma=\frac{1}{3}$ and $\gamma=\frac{1}{4}$ in Figure 3 in Section 5.

### Simulation-Plot
Inference/Result data contains all the result data from running Simulation-Inference.R under different settings to produce Figures 2 and 4-5 in Section 5, and Figures 9-11 in Appendix A.1. Inference/Simulation-Inference-Plot.R is the code to produce the figures. 

Regret/Result data contains all the result data from running Simulation-Regret.R under under $\gamma=\frac{1}{3}$ and $\gamma=\frac{1}{4}$ to produce Figure 3 in Section 5. Regret/Simulation-Regret-Plot.R is the code to produce the figures. 

### Real Data