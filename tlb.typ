#set math.equation(numbering: "(1)")
#set text(font: "New Computer Modern")
#set par(justify: true)

= Bounding the Log Likelihood Difference

Consider the difference in log-likelihood between two waveforms $g$ and $h$ given data $d$ (we will imagine that $g$ and $h$ are two different _models_ for the waveform that is actually contained in $d$, but our argument here is fully general).  We can write 
$ Delta log cal(L) := log cal(L)_h - log cal(L)_g = -1/2 (<d-h|d-h> - <d - g | d-g >). $ <logl-difference>
A sufficient condition that the two models are equivaent in a parameter analysis in some region in parameter space is that $| Delta log cal(L)| << 1$ for parameter values in that region; we may refer to this as "pointwise equivalence" because the log-likelihood difference is evaluated at corresponding points in parameter space.  

Applying the law of cosines, @law-of-cosines, with $a = d-g$, $b = g-h$, and $c = d-h$, the log likelihood becomes 
$ Delta log cal(L) = -1/2 <g-h|g-h> - <d-g|g-h> $ <law-of-cosines-logl>

Let us suppose that the two waveforms of interest, $h$ and $g$, are close, so that we can write $g = h + delta h$ with $|delta h| << 1$.  Then we see that the first term in @law-of-cosines-logl is second order in $|delta h|$; we will neglect it from here onward.  Expanding the second term, and preserving only to first-order in $|delta h|$, we have 
$ |Delta log cal(L)| #sym.tilde.eq |<d-h | delta h>|. $ <logl-difference-data-dependent>
In words: the change in the log likelihood is the projection of the residuals onto the wavefrom difference.

At this point, due to the presence of $d$ in @logl-difference-data-dependent we cannot say much more in general (it is, in principle, possible for the projection of $d-h$ onto $delta h$ to be essentially arbitrarily large).  But if we are willing to assume that the data contain a waveform $H$ that is reasonably close to $h$ (and $g$)
$ d = H + n, $
plus noise $n$ that is consistent with the spectral density used to define our inner product, then 
$ | Delta log cal(L) | #sym.tilde.eq |<n | delta h> + < H | delta h> - <h | delta h> |. $ <logl-noise-term-included>
(Note that it is not the case that $< H - h | H - h> << 1$, since the difference between the true waveform and a fitted waveform can be---in fact, must be---$cal(O)(1)$ for waveforms $h$ with good posterior support in a parameter estimation.)  

The first term in the right hand side of @logl-noise-term-included is a random variable with zero mean and variance $< delta h | delta h> << 1$, and is therefore $cal(O)(delta h)$ in magnitude.  We will see below that the other two terms are $cal(O)(rho delta h)$ (see @rho-definition); presuming that $rho >> 1$, we will ignore the first term.  Applying the triangle inequality to the remainder of @logl-noise-term-included, and using the Cauchy-Schwarz inequality, we have
$ | Delta log cal(L) | #sym.lt.tilde (sqrt(<H|H>) + sqrt(<h|h>)) sqrt(<delta h | delta h>) $
(we use $#sym.lt.tilde$ to remind the reader that we are ignoring sub-leading-order contributions in $rho delta h$).  From here we will assume that $h$ is a sufficiently good model for $H$ that the leading order behavior of $<H|H>$ is the same as $<h|h>$, both scaling as 
$ rho^2 := <h | h> = <H|H> + cal(O)(rho). $ <rho-definition>
Under this assumption, 
$ Delta log cal(L) #sym.lt.tilde 2 rho sqrt(<delta h | delta h>) = 2 rho^2 (| delta h |)/(| h |). $

We can relate $< delta h | delta h > = | delta h |^2$ to the mismatch between $g$ and $h$, which is defined by 
$ cal(M) := 1 - (<g | h>) / (sqrt(<g|g>) sqrt(<h|h>)). $
To linear order in $delta h$ this becomes 
$ cal(M) #sym.tilde.eq 2 (<h| delta h>)/(<h|h>). $
Again applying Cauchy-Schwarz, we have 
$ | cal(M) | #sym.lt.tilde 2 (sqrt(<delta h | delta h>))/(sqrt(<h|h>)). $
This is the *wrong direction* to claim with certainty that
$ Delta log cal(L) #sym.lt.tilde rho^2 cal(M), $
but assuming there is no special orientational issues (i.e. that the projection of $delta h$ onto $h$ is not particularly small compared to their magnitudes), we can say that 
$ Delta log cal(L) ~ rho^2 cal(M); $
but in any case, it is *bounded* by 
$ Delta log cal(L) #sym.lt.tilde 2 rho^2 (|delta h|)/(|h|). $

= Definitions

== Mismatch
The mismatch $cal(M)$ between two waveforms $g$ and $h$ is given by 

$ cal(M) := 1 - (<g|h>)/(sqrt(<g|g>) sqrt(<h|h>)) = 1 - <hat(g) | hat(h) > $ <mismatch>

== Log Likelihood
The log-likelihood for some waveform $h$ given data $d$ is 

$ log cal(L) := -(1)/(2) <d-h|d-h> $

== Triangle Inequality
The triangle inequality applies to any vector space and states that for any vectors $a$ and $b$
$ || a + b || <= ||a|| + ||b|| $  <triangle>
Equivalently, for vectors $c$ and $d$
$ ||c|| - ||d|| <= ||c-d|| $ <triangle-reformulated>
where $||x||$ is the norm of $x$.
(this latter follows from the former by setting $a = c-d$ and $b = d$).

We can also square @triangle to obtain 
$ || a + b ||^2 <= ||a||^2 + 2 ||a|| ||b|| + ||b||^2, $ <triangle-squared>
or
$ ||c||^2 - ||d||^2 <= ||c-d||^2 + 2 ||c-d|| ||d||. $ <triangle-reformulated-squared>

== Cauchy-Schwarz Inequality

The Cauchy-Schwarz inequality for an inner product states that for all vectors $a$ and $b$
$ |<a | b> | <= sqrt(<a|a>) sqrt(<b|b>) $

== The Law of Cosines

For vectors $a$, $b$, and $c$ with $c = a+b$ we have 
$ <c|c> = <a|a> + <b|b> + 2 <a|b> $ <law-of-cosines>

