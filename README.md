# Implementation of the Model Predictive Control for the regulation of the intestinal bacterial overgrowth

### Introduction
The gut microbiome and the influence of its composition on the overall health and well-being has been getting an increasing attention.

Systems biology and methods of the control theory can give us additional insight into this interesting world of small bugs. 
In these scripts, I outline a way in which the optimal control theory could help us to investigate and deal with
the malfunctioning dynamics of the gut ecosystem. 

Specifically, we look at the loss of microbiota diversity (LOMS), which seems to be the most constant finding of the intestinal dysbiosis,
which often seems to be associated with the modern Western lifestyle and diet [1]. A hypothesis suggests that the LOMS is associated with
the decrease of the population of the bacterial predators in the gut. As the health issues caused by the negative perturbations to the
gut microbiome populations can be quite detrimental, it is worth to think about how we could restore the balance of the bacteria in the gut,
so that the population of the predator bacteria increase. This could perhaps be accomplished by the dietary interventions.

While the dynamics of the relationship between the small monsters in our gut is quite intricate, here I consider a simplified situation:
We model the dynamics between two bacterial species in the gut, using the Lotka-Volterra predator-prey model adopted from [2]. The states represent
the bacterial masses of predator (i. e. *Bdellovibrio bacteriovorus*) and prey (i. e. *Prevotella intermedia*) bacteria.
The input to the system is the change in the basal feeding rate of the prey bacterial species, 
while it is assumed that the predator bacteria are only fed by the prey bacteria. The system equations are as follows

\begin{align}
\frac{dB}{dt} &= (r+u)B \left( 1- \frac{B}{k} \right) - \frac{aBP}{c+B}, B \geq 0 \\ 
\frac{dP}{dt} &= b \frac{aBP}{c+B} - dP, P \geq 0
\end{align}

We purposefully start with the very low mass of the predator bacteria (5 g) and high mass of the prey bacteria (100 g), out of the equilibrium
and we simulate the system with the given parameters 
The time evolution of the system is shown in the figure below, where it is apparent that the system oscillates with a period of roughly 12 days.

## References
[1] Mosca, A., Leclerc, M., & Hugot, J. P. (2016). Gut Microbiota Diversity and Human Diseases: Should We Reintroduce Key Predators in Our Ecosystem? Frontiers in Microbiology, 7, 455. https://doi.org/10.3389/fmicb.2016.00455 <br>
[2] Åström, K. J., & Murray, R. M. (2008). Feedback Systems (1st ed.). Princeton: Princeton University Press. Retrieved from http://press.princeton.edu/titles/8701.html.
