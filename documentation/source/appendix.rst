.. _chap:appendix:

Appendix
========

This appendix contains lists of all probability density distributions and all kernels that GRAMPC-S supports.

.. _sec:appendix_pdf:

List of probability density functions
-------------------------------------

The following list contains all probability density functions that are implemented in GRAMPC-S. Note that it is possible to implement your own probability density functions. These must be declared as derived from the ``Distribution`` class.

.. list-table:: Overview of the available probability distributions and associated parameters.
    :name: tab:pdf
    :widths: auto
    :header-rows: 1
    
    * - Name 
      - Probability density function 
      - Parameters
    * - Gaussian distribution 
      - :math:`p(\vm x) = \left(2 \pi\right)^{-\frac{d}{2}} \det(\vm \Sigma)^{-\frac{1}{2}} \exp \left(-\frac{1}{2} (\vm x - \vm \mu)^T \vm \Sigma^{-1} (\vm x - \vm \mu)\right)` 
      - :math:`\vm \mu \in \mathbb{R}^d` :math:`\vm \Sigma \in \mathbb{R}^{d \times d}`  
    * - Beta distribution 
      - :math:`p(x) = \frac{1} {B(p, q)} \, x^{p-1} \,(1 - x)^{q - 1}` 
      - :math:`p>0\,, \; q >0` 
    * - Chi-squared distribution 
      - :math:`p(x) = \frac{x^{\frac{n}{2} - 1} \,\exp\left(-\frac{x}{2}\right)}{\Gamma(\frac{n}{2}) \, 2^{\frac{n}{2}}}` 
      -  :math:`n>0`
    * - Exponential distribution 
      - :math:`p(x) = \lambda \exp(-\lambda \, x)` 
      - :math:`\lambda >0`
    * - Extreme value distribution 
      - :math:`p(x) = \frac{1}{b} \exp\left(\frac{a - x}{b} - \exp\left(\frac{a - x}{b}\right)\right)`
      - :math:`a \in \mathbb{R}\,,\; b>0` 
    * - F-distribution 
      - :math:`p(x) = \frac{\Gamma\left(\frac{m+n}{2}\right)}{\Gamma\left(\frac{m}{2} \right) \Gamma\left( \frac{n}{2}\right)} \left(\frac{m}{n}\right)^\frac{m}{2}  x^{\frac{m}{2} - 1} (1 + \frac{m}{n}  x)^{-\frac{m+n}{2}}` 
      - :math:`m>0\,,\; n >4` 
    * - Gamma distribution 
      - :math:`p(x) = \frac{\exp\left(-\frac{x}{\beta}\right)}{\beta^\alpha  \, \Gamma(\alpha)} \, x^{\alpha - 1}` 
      - :math:`\alpha > 0\,,\; \beta > 0` 
    * - Log-normal distribution 
      - :math:`p(x) = \frac{1}{\sigma x \sqrt{2 \pi}}  \exp\left(- \frac{(\ln(x)-\mu)^2} {2 \sigma^2}\right)` 
      - :math:`\mu \in \mathbb{R} \,,\; \sigma > 0`
    * - Piecewise constant distribution 
      - Histogram of an arbitrary distribution 
      - 
    * - Student's t-distribution 
      - :math:`p(x) = \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sigma \sqrt{\nu \pi} \Gamma\left(\frac{\nu}{2}\right)} \left(1 + \left(\frac{x-\mu}{\sigma}\right)^2 \frac{1}{\nu}\right)^{-\frac{\nu+1}{2}}` 
      - :math:`\mu \in \mathbb{R} \,,\; \sigma > 0` :math:`\nu > 2` 
    * - Uniform distribution 
      - :math:`p(x) = \frac{1}{b - a}` 
      - :math:`a < b`
    * - Weibull distribution 
      - :math:`p(x) = \frac{a}{b} \left(\frac{x}{b}\right)^{a-1} \exp(-\left(\frac{x}{b}\right)^a)` 
      - :math:`a> 0\,,\; b > 0` 
    * - Product of uncorrelated distributions 
      - :math:`\displaystyle p(\vm x) = \prod_{i=1}^d p_i(x_i)`
      - 


.. _sec:appendix_kernels:

List of kernels for Gaussian processes
--------------------------------------

The following list contains all kernels that are implemented in GRAMPC-S.
Note that it is possible to implement your own kernels.
These must be declared as derived from the ``StationaryKernel`` class.

.. list-table:: Overview of the available kernel functions and associated parameters.
    :name: tab:kernels
    :widths: auto
    :header-rows: 1
    
    * - Name 
      - Kernel 
      - Parameters
    * - Squared exponential kernel 
      - :math:`k(\tau) = \sigma^2 \exp\left(-\frac12 \sum\limits_i \left(\frac{\tau_i^2}{l_i^2}\right)\right)` 
      - :math:`\sigma`, :math:`\vm l`
    * - Periodic kernel 
      - :math:`k(\tau) = \sigma^2 \exp\left(-2 \sum\limits_i \left(\frac{\sin^2\left(\pi \frac{\tau_i}{p_i}\right)}{l_i^2}\right)\right)` 
      - :math:`\sigma`, :math:`\vm l`, :math:`\vm p`
    * - Locally periodic kernel 
      - :math:`k(\tau) = \sigma^2 \exp\left(-2 \sum\limits_i\left(\frac{\sin^2(\pi \frac{\tau_i}{p_i})}{l_i^2}\right)\right) \exp\left(-\frac12 \sum\limits_i\left(\frac{\tau_i^2}{l_i^2}\right)\right)` 
      - :math:`\sigma`, :math:`\vm l`, :math:`\vm p`
    * - Matern 3/2 kernel 
      - :math:`k(\tau) = \sigma^2 \left(1 + \sqrt{3 \sum\limits_i\left(\frac{\tau_i^2}{l_i^2}\right)} \right) \exp\left(-  \sqrt{ 3\sum\limits_i\left(\frac{\tau_i^2}{l_i^2}\right)} \right)` 
      - :math:`\sigma`, :math:`\vm l`
    * - Matern 5/2 kernel 
      - :math:`k(\tau) = \sigma^2 \left(1 + \sqrt{5 \sum\limits_i\left(\frac{\tau_i^2}{l_i^2}\right)} + \frac{5}{3} \sum\limits_i\left(\frac{\tau_i^2}{l_i^2}\right)\right) \exp\left(- \sqrt{5 \sum\limits_i\left(\frac{\tau_i^2}{l_i^2}\right)}\right)` 
      - :math:`\sigma`, :math:`\vm l`
    * - Sum of kernels 
      - :math:`k(\tau) = k_1(\tau) + k_2(\tau)` 
      - 
    * - Product of kernels 
      - :math:`k(\tau) = k_1(\tau) k_2(\tau)` 
      - 

