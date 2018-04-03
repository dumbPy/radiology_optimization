---


---

<h1 id="cancer-radiology-optimization-tomotherapy">Cancer Radiology Optimization (Tomotherapy)</h1>
<p>This Project deals with Python Implementation of the Paper <a href="https://doi.org/10.1137/S0036144598342032"><em>Optimizing the Delivery of Radiation Therapy to Cancer Patients</em> by David M Sheperd et. al. 1999</a> and then solving the problem given below using the Optimization Technique described in the above paper.<br>
The Radiotherapy technique described in the above paper is also called as <em><strong>Tomotherapy</strong></em></p>
<p>AMPL API <a href="https://github.com/ampl/amplpy">(amplpy)</a> was used for optimization routine. Two different optimization models were written.<a href="https://github.com/dumbPy/radiology_optimization/blob/master/models/model.mod">[1]</a><a href="https://github.com/dumbPy/radiology_optimization/blob/master/models/model_linear_norm.mod">[2]</a><br>
Suppose we want to design an IMRT setup to obtain a dose corresponding to the following figure.</p>
<p><img src="https://lh3.googleusercontent.com/pVkc8Or2Tj2ninBZCKhSbRMXABF4YZkXpUdC7gaBrtRSJRX-T8m47gbCuDX-CXcxqMaa7P4yxqU" alt="Image"></p>
<p>Assume that we want to achieve close to 10 units of radiation in the red region (denoting a tumor) and be as much below 4 units as possible in the remaining portion. You may assume the radiation intensity diminishes to half every 4 units (which means that this relation is not linear).<br>
Assuming that each beam is 1 unit wide, find a suitable configuration of beam weights and angles that can realize a plan as close as possible to the desired treatment. You may further assume that each beam is composed of 20 small beamlets each of whose intensities can be controlled separately.<br>
Solve one of the LPs or QPs proposed described by Shephard et al and find the optimal configuration. Describe the configuration in a readable manner. Also draw a pictorial representation of the intensity pattern achieved by your proposed solution. Comment on its quality. Clearly state all your assumptions and choices of parameters.</p>
<h2 id="notebook-implementation">Notebook Implementation</h2>
<p>2 Different Jupyter Note Implementations have been tried as below.</p>
<ul>
<li><a href="https://github.com/dumbPy/radiology_optimization/blob/master/tomotherapy_no_hard_cap.ipynb">No Hard Cap</a> - This Implementation uses 1st model above. It puts no cap on how much the radiation can reach in tumour region</li>
<li><a href="https://github.com/dumbPy/radiology_optimization/blob/master/tomotherapy_with_linearized_norm_error.ipynb">Linearized Norm Error</a> - This Implementation uses 2nd model above and is more closer to the optimization equation discussed in paper. It reduces |error| and hence tries to keep radiation in tumor region near the desired dosage level of 10</li>
</ul>
<p>Python <em>.py</em> scripts are also available <a href="https://github.com/dumbPy/radiology_optimization/tree/master/python%20scripts">here</a></p>

