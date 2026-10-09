# pHigherGreens

**SageMath computations for *p-adic Higher Green's Functions for Stark-Heegner cycles.***

This repository contains code that acommpanies the paper
> Hazem Hassan, p-adic higher Green's functions for Stark--Heegner cycles.

It contains the SageMath implementation of `Algorithm 5.7` of the paper, which compute values of p-adic higher Green's functions. Furthermore, this repository shows the calculations leading to Examples 5.9 and 5.10 of the paper.

For now this repository contains only enough to be able to reproduce the examples from the paper. Certainly one can use this code to explore further examples, but it is not ideal. I am working on improving the speed of the algorithm and streamlining this code so that it is easier to use.

## Requirements
- [SageMath](https://www.sagemath.org/): The computations were computed using **SageMath 10.1**, and relies on several SageMath libraries.

## Files

| File | Purpose |
| --- | --- |
| `LMMC_Main.ipynb` | Object-oriented implementation of the mathematical objects and routines required for the computation of `Algorithm 5.7`. This is the main file.  |
|`Example 5.9.ipynb` | Contains the computation for `Example 5.9`. Relies on stored precomputed data.|
|`Example 5.10.ipynb` | Contains the computation for `Example 5.10`. Relies on stored precomputed data.|
|`Example 5.9 GEN.ipynb` | Generates the cocycle date required for `Example 5.9`. Takes a significant amount of time to run (10 - 12 hours on my device). |
|`Example 5.10 GEN.ipynb` | Generates the cocycle date required for `Example 5.10`. Takes a significant amount of time to run (10 - 12 hours on my device). |
|`*.json` | Precomputed cocycles stored in JSON format. |

## Reproducing the examples
Make sure you have installed Jupyter and SageMath and clone this repository

```bash
git clone https://github.com/HHassan744/pHigherGreens.git
```
Open either of `Example 5.9.ipynb` or `Example 5.10.ipynb` in SageMath, and run the code cells in order.

## Understanding the Example files
The files start by specifying some parameters. For example in `Example 5.10`

```sage
k = 4       
p = 3       
D1 = 5      
D2 = 11   
pprec = 150

%run LMMC_Main.ipynb
```
Here 
- `k` is the weight
- `p` is the prime
- `D1` and `D2` are the discriminants of the divisors
- `pprec` is the p-adic precision

Then it runs 
```sage
LMMC8 = read("J8F1p3pr150.json")
LMMC32 = read("J32F1p3pr150.json")
LMMCcomb = 7 * LMMC8 - 4 * LMMC32
```

The files `J8F1p3pr150.json` and `J32F1p3pr150.json` contain data for the `MittagLeffler` objects associated to the cocycles of weight 4 and RM-points `sqrt(2)` and `2sqrt(2)`, respectively. LMMCcomb is then the linear combination associated to a divisor of strong degree zero.

The next lines are a sanity check 
```sage
LMMCcomb.term2()
LMMCcomb.term3()
LMMCcomb.fix_poly_invariance()
(LMMCcomb + LMMCcomb.mobius_action(Sm)).poly()
(LMMCcomb + LMMCcomb.mobius_action(Um) + LMMCcomb.mobius_action(Um**2)).poly()
```
Both `LMMCcomb.term2()` and `LMMCcomb.term3()` check that the non-polynomial part of the Mittag-Leffler expansion satisfies the two-term and three-term relation, as they should for a divisor of strong degree zero. 

`LMMCcomb.fix_poly_invariance()` applies the step of fixing the polynomial part described in `Algorithm 5.7`.

The next two lines check that the polynomial parts now also satisfy the two-term and three-term relation, as is expected for a principal divisor.

The `LMMC` objects only hold the information of level 1 and higher of the cocycle. The level 0 parts are stored differently as they do not give power series. Instead they are stored as a PolyLog object, which is much quicker to compute. 

```sage
PL8 = Level0PolyLog(sqrt(2))
PL32 = Level0PolyLog(2*sqrt(2))
PL832 = 7 * PL8 - 4 * PL32
PL832.fix_poly_invariance()
```
The lines above compute the `PolyLog` objects associated to our particular divisor and again apply the algorithm fixing the polynomial part.

Now it is time to evaluate these objects (here at `4*sqrt(2) + 5`)

```sage
s = 4 * rd1 + 5
ML_to_eval = MittagLeffler([PS(0) for i in range(p+1)])
PL_to_eval = PolyLog({})

for M in automorph_sequence_pos(s):
    ML_to_eval += LMMCcomb.mobius_action(M**(-1))
    PL_to_eval += PL832.mobius_action(M**(-1))
```
The above lines compute `J(\gamma_{\sigma})` by using the continued fraction expansion. 

```sage
final_answer_832 = 0

for i in range(p):
    final_answer_832 += ML_to_eval.PowSeries[i](1/(hom1(s-i)))/hom1(s-aut1(s))

final_answer_832 += ML_to_eval.PowSeries[p](hom1(s))/hom1(s-aut1(s))
final_answer_832 += (PL_to_eval.evaluate(hom1(s))/(hom1(s-aut1(s))))/4

ML_to_eval_diff = ML_to_eval.differentiate()

for i in range(p):
    final_answer_832 += (-1) * ML_to_eval_diff.PowSeries[i](1/hom1(s-i))/2
final_answer_832 += (-1) * ML_to_eval_diff.PowSeries[p](hom1(s))/2
final_answer_832 += ((-1) * PL_to_eval.derivative().evaluate(hom1(s)))/8
```
Computes the value using the formula described in `Algorithm 5.7`

The next few lines just prepare the field where the values are expected to take place and the expected primes to appear in the factorization 

```sage
QF3.<rd3> = QuadraticField(-1)
hom3 = QF3.hom([K(-1).square_root()])
primes = [5, 13, 19, 29]
pr = []
for i in primes:
    for j in QF3.primes_above(i):
        pr.append(j.gen(1))
pr_log = [hom3(i).log(0) for i in pr]
```

Finally, the `lindep` function uses the `LLL` algorithm to find linear relationships between the logarithms of prime factors in the field `QF3` and the final answer.
The output in this case has in its first line 

```sage
[-624   624   -90   90   0   273   -273   192   0   0]
```
These are coefficients that weight the different values for a linear combination that equals zero. In this line, the final two entries correspond to the error terms of the linear combination (there are two entries because the values take place in `Q_{9}`). The third last entry is the coefficient of `final_answer_832`, which is the value of the higher Green's functions. All other coefficients are the weights of the prime factors. 

## Understanding the cocycle generation files
For example `Example 5.10 GEN.ipynb`
Starts as before by specifying the parameters and running the main file. 

Then computes the Mittag Leffler expansion of the level 1 rational function on the standard affinoid. This is stored in the `MittagLeffler` object `LMMC32`. 
```sage
t = 2*sqrt(2)
J32 = Level1MittagLeffler(t)
```

Then the following lines of code apply the recursion described in `Algorithm 5.7`
```sage
for i in range(pprec):
    previous = (-1) * previous.Next()
    LMMC32 += previous
    print('completed level ' + str(i) + " at " + str(datetime.datetime.now()))
```
This recursion is quite slow and may take several hours to compute. 

Finally, the object is stored in a JSON file.
```sage write(LMMC32, "J32F1p3pr150.json")```
