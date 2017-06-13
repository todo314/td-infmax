Time-Decaying Influence Maximization
========================

This is a scalable algorithm for influence maximization under the time-varying independent cascade and linear threshold models.

## Usage
Given a graph with edge probability functions and edge length likelihoods, it selects a seed set of a given size.

    $ make
    $ ./imm -graph=GRAPH -k=K -alg=ALG -eps=EPS -ell=ELL -model=MODEL -delay=DELAY -prob=PROB -numMC=NUMMC
* graph: input file (see below)
* k: seed set size
* alg: algorithm (see Algorithms)
* eps: algorithm parameter (eps=0.5 in the paper)
* ell: algorithm parameter (ell=1 in the paper)
* model: diffusion model (see Model parameters)
* delay: edge length likelihoods (see Model parameters)
* prob: edge probability functions (see Model parameters)
* numMC: number of Monte-Carlo simulations for estimating the influence spread of the obtained solution (numMC=10000 in the paper)

### Algorithms
You can choose the proposed algorithm "IMM for TVIC and TVLT" or
three baselines: "IMM-CTIC", "IMM-IC", and "IMM-LT" (see Sec. 6 for further details).
For every algorithm, you have to specify the values of eps ($\eps$ in the paper) and ell ($\ell$ in the paper).

* IMM for TVIC and TVLT (alg=imm): 
This is the proposed algorithm for time-varying diffusion models (see Sec. 5 for further details).

* IMM-CTIC (alg=imm-ctic): 
An adaption of a sketching method for the continuous-time independent cascade model.
Since it takes care of "deadlines" rather than time-decaying edge probabilities,
you are required to specify the "deadline" (deadline=1 in the paper).

* IMM-IC (alg=imm-ic): 
An adaption of a sketching method for the independent cascade model.

* IMM-LT (alg=imm-lt): 
An adaption of a sketching method for the linear threshold model.

### Model parameters
You can choose one of the following four settings of diffusion model (model), edge probability functions (prob), and edge length likelihoods (delay) (See Sec. 6 for further details):

<table>
	<tr>
	<td>model</td>
	<td>prob</td>
	<td>delay</td>
	<td>Description</td>
	<td>Probability function</td>
	<td>Length likelihood</td>
	</tr>
	<tr>
	<td>tvic</td>
	<td>exponential</td>
	<td>weibull</td>
	<td>Weighted exponential IC</td>
	<td><img src="https://latex.codecogs.com/gif.latex?p_{uv}(t) = \frac{1}{|N^-(v)|} \exp(- c_{uv}t)"/></td>
	<td><img src="https://latex.codecogs.com/gif.latex?f_e(\delta) = \frac{\alpha_e}{\beta_e} \cdot \left( \frac{\delta}{\beta_e} \right)^{\alpha_e-1} \cdot \exp\left(-\left( \frac{\delta}{\beta_e} \right)^{\alpha_e}\right)"/></td>
	</tr>
	<tr>
	<td>tvic</td>
	<td>inverse</td>
	<td>weibull</td>
	<td>Weighted reciprocal IC</td>
	<td><img src="https://latex.codecogs.com/gif.latex?p_{uv}(t) = \frac{1}{|N^-(v)| c_{uv} t}"/></td>
	<td><img src="https://latex.codecogs.com/gif.latex?f_e(\delta) = \frac{\alpha_e}{\beta_e} \cdot \left( \frac{\delta}{\beta_e} \right)^{\alpha_e-1} \cdot \exp\left(-\left( \frac{\delta}{\beta_e} \right)^{\alpha_e}\right)"/></td>
	</tr>
	<tr>
	<td>tvlt</td>
	<td>exponential</td>
	<td>weibull</td>
	<td>Weighted exponential LT</td>
	<td><img src="https://latex.codecogs.com/gif.latex?q_{uv}(t) = \frac{1}{|N^-(v)|} \exp(- c_{uv}t)"/></td>
	<td><img src="https://latex.codecogs.com/gif.latex?f_e(\delta) = \frac{\alpha_e}{\beta_e} \cdot \left( \frac{\delta}{\beta_e} \right)^{\alpha_e-1} \cdot \exp\left(-\left( \frac{\delta}{\beta_e} \right)^{\alpha_e}\right)"/></td>
	</tr>
	<tr>
	<td>tvlt</td>
	<td>inverse</td>
	<td>weibull</td>
	<td>Weighted reciprocal LT</td>
	<td><img src="https://latex.codecogs.com/gif.latex?q_{uv}(t) = \frac{1}{|N^-(v)| c_{uv} (t+1)}"/></td>
	<td><img src="https://latex.codecogs.com/gif.latex?f_e(\delta) = \frac{\alpha_e}{\beta_e} \cdot \left( \frac{\delta}{\beta_e} \right)^{\alpha_e-1} \cdot \exp\left(-\left( \frac{\delta}{\beta_e} \right)^{\alpha_e}\right)"/></td>
	</tr>
</table>

### Example
Running IMM for TVIC and TVLT, IMM-CTIC, IMM-IC for the TV-IC model.

    $ ./imm.exe -graph=.test.tsv -k=100 -alg=imm -eps=0.5 -ell=1 -model=tvic -delay=weibull -prob=exponential -numMC=10000
    $ ./imm.exe -graph=.test.tsv -k=100 -alg=imm-ctic -deadline=1 -eps=0.5 -ell=1 -model=tvic -delay=weibull -prob=exponential -numMC=10000
    $ ./imm.exe -graph=.test.tsv -k=100 -alg=imm-ic -eps=0.5 -ell=1 -model=tvic -delay=weibull -prob=exponential -numMC=10000
Running IMM for TVIC and TVLT, IMM-LT for the TV-LT model.

    $ ./imm.exe -graph=.test.tsv -k=100 -alg=imm -eps=0.5 -ell=1 -model=tvlt -delay=weibull -prob=exponential -numMC=10000
    $ ./imm.exe -graph=.test.tsv -k=100 -alg=imm-lt -eps=0.5 -ell=1 -model=tvlt -delay=weibull -prob=exponential -numMC=10000

### Format of the input graph
    u_1 v_1 alpha_1 beta_1 c_1
    ...
    u_i v_i alpha_i beta_i c_i
    ...
    u_m v_m alpha_m beta_m c_m
    
* The i-th line means an directed edge (u_i, v_i) with model parameters alpha_i, beta_i, and c_i.
* Vertices should be described by integers starting from zero.

## References

Naoto Ohsaka, Yutaro Yamaguchi, Naonori Kakimura, and Ken-ichi Kawarabayashi. **[Maximizing Time-Decaying Influence in Social Networks](https://link.springer.com/chapter/10.1007/978-3-319-46128-1_9)**.
European Conference on Machine Learning and Principles and Practice of Knowledge Discovery in Databases (ECML PKDD), pages 132--147, 2016.
