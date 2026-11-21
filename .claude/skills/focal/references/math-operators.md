# Focal math operators

Definitions of the operators Focal implements, matched to the code in
`focal/temporal/_forces.py` and `focal/base/enrichment.py`.

## 1. The force operator

For TF $i$, target $j$, pseudotime $t$:

$$
F(i, j, t) = \text{sign}(\beta_{ij}(t)) \cdot \exp\!\Big(\log_{10}(|\beta_{ij}(t)| + \epsilon) + \log_{10}(E_{\text{TF}_i}(t) + \epsilon)\Big)
$$

$\beta_{ij}(t)$ is the dictys edge weight (signed), $E_{\text{TF}_i}(t)$ the log-CPM expression
of TF $i$, $\epsilon = 10^{-10}$. Equivalently
$F = \text{sign}(\beta)\,\big((|\beta|+\epsilon)(E+\epsilon)\big)^{1/\ln 10}$: a sign-preserving
compression of $\beta \cdot E$.

- **Sign**: $\beta > 0$ activating ($F > 0$), $\beta < 0$ repressive ($F < 0$); the sign rides
  through unchanged, so activating and repressive tails can be thresholded separately.
- **Units**: none. $F$ is a ranking statistic compared within one episode; percentile thresholds
  are unit-invariant.

## 2. The ORA operator

Universe $\Omega$ ($N$ = genes in the episodic GRN), program genes that are active targets
($K$), TF out-degree $n$, overlap $k$:

$$
p = P_{H_0}[X \geq k] = \texttt{hypergeom.sf}(k-1, N, K, n), \qquad ES = \frac{k}{nK/N}
$$

$H_0$: the program is a uniformly random $K$-subset of $\Omega$, independent of the TF's targets.
This is the one-sided Fisher exact test on the 2x2 table with both margins fixed; always use the
exact hypergeometric (overlaps are often 2–3 genes, so $\chi^2$ is invalid). The test is symmetric
in the two sets; Focal frames it as "does the program concentrate on this TF's downstream targets".

## 3. Episodic invariance

An edge $(i, j)$ enters an episodic GRN only if, within the episode:

1. it is present (non-zero row sum) and the regulator is not `ZNF*` / `ZBTB*`;
2. it has >= 3 non-zero points;
3. a one-sample t-test of the non-zero values vs 0 gives p < 0.05, then p < `pval_threshold`
   (default 1e-3);
4. all non-zero values share one sign (direction invariance).

Then the **top-percentile** step keeps the strongest of these coherent edges by |mean force|
(or separate activating / repressive tails). Invariance runs first, magnitude second.

## 4. Multiple testing

Focal applies **no FDR correction**; enrichment tables carry raw p-values and the notebooks use
raw-p display cutoffs. Say so in Methods. If a correction is needed, apply BH (or BY under
arbitrary dependence) per (program, episode) on the output table yourself.
