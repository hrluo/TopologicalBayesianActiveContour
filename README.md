


# TopologicalBayesianActiveContour

**Content**
This is the code repository for the research publication "Combining Geometric and Topological Information in Image Segmentation" by [Justin D. Strait](https://jdstrait.weebly.com/) and [Hengrui Luo](https://hrluo.github.io/). 
The manuscript of this paper can be accessed at https://arxiv.org/abs/1910.04778. 

 - The TOP code we store in [TOP folder](https://github.com/hrluo/TopologicalBayesianActiveContour/tree/master/TOP) folder for comparison convenience is associated with the paper [A Topological Approach to Hierarchical Segmentation using Mean Shift](https://people.csail.mit.edu/sparis/)  and Section 1 and 2 of [our paper](https://arxiv.org/abs/1910.04778). The copyright of this code belongs to the original author.
 - The TAC code we store in [TAC folder](https://github.com/hrluo/TopologicalBayesianActiveContour/tree/master/TAC) folder stores the code for both BAC (Bayesian Active Contour) and TAC (Combined method) described in details in [our paper](https://arxiv.org/abs/1910.04778). The copyright of this code belongs to [Justin D. Strait](https://jdstrait.weebly.com/) and [Hengrui Luo](https://hrluo.github.io/). 

**Structure**
Below we describe the files in the [TAC folder](https://github.com/hrluo/TopologicalBayesianActiveContour/tree/master/TAC) folder for your reference.
 - Bayesian Active Contour (BAC).
	 **The full details of BAC method is described both in the paper [Bayesian Active Contours with Affine-Invariant, Elastic Shape Prior](https://ieeexplore.ieee.org/document/6909441) and Section 1 and 2 of [our paper](https://arxiv.org/abs/1910.04778).**
	 - `curve_to_q.m` The function that converts the coordinates of points sampled from a curve into unit-norm SRVF of curve.
	 - `DynamicProgrammingQ.c` The function that finds a minimum-cost path used for elastic shape analysis. The original code is by J. D. Tucker based on equation (5.2) in his thesis, the current version is modified by Justin D. Strait and Hengrui Luo.
	 - `ElasticShooting.m` The function that computes the SRVF of resulting curve after shooting a vector v on tangent space at a direction q1.
	 - `ElasticShootingVector.m` The function that computes the geodesic distance (after optimal re-prameterizatoin and rotation) between two curves represented in the  SRVF form and the optimal re-parameterization and rotation between them.
 		 - `Find_Rotation_and_Seed_unique.m`   The function that computes the optimal reparameterizatoin between two curves (after optimally rotated) represented in the  SRVF form assuming the starting points of these two curves are fixed. 
			 - `Find_Best_Rotation.m`  The function that computes the optimal rotation between two curves represented in the  SRVF form assuming the starting points of these two curves are fixed.
	- `FindElasticMean.m` Compute the the mean shape $\bar{q}$ on shape space (and its SRVF form). 
	- `FindElasticCovariance.m` Compute the covariance matrix on the tangent space (to the shape space) at the mean shape $\bar{q}$.
	- `Form_Basis_Normal_A.m` Find an orthonormal basis of vectors for the normal space of a curve (in SRVF), with respect to inner product
	- `Gram_Schmidt.m` Find an orthonormal basis constructed via Gram-Schmidt on a given space. 
	- `Group_Action_by_Gamma_Coord.m` Apply reparameterization to a curve in its own parameterization.
	- `Group_Action_by_Gamma_Coord_Q.m`  Apply reparameterization to a curve in its SRVF parameterization.
	- `ImageEnergy.m` Compute the $E_{energy}$ term in the energy functional used for BAC method.
		- `PixelDensityEst.m` Estimate the  density of interior  and exterior pixel values of an image.
		- `TrainingPixelDensity.m`   Estimate the  density of interior  and exterior pixel values of training images.
	- `ImageUpdate.m` Compute the needed gradient of $E_{energy}$ along current contour in order to update the energy functional.
	- `InnerProd_Q.m` Compute the inner product of two curves in SRVF form.
	- `invertGamma.m` Compute the inverse of parameterization function.
	- `MainScript.m` The script that we use to run all examples in our paper. The data, results, and TOP initializations called and run in this are found in the `Manuscripts` folder off the main path of the repository.
	- `OutwardUnitNormal.m` Compute the outward normal vector to a parameterized curve.
	- `Parallel_Transport_C.m` Parallel transport a tangent vector along a parameterized curve.
	- `PriorEnergy.m` Compute the $E_{prior}$ term in the energy functional used for BAC method.
	- `PriorUpdate.m` Compute the needed gradient of $E_{prior}$ along current contour in order to update the energy functional.
	- `Project_Tangent.m` Take the projection of a tangent vector into tangent space to the shape space at a point $q$.
	- `ProjectC.m` Take the projection of SRVF $q$ into space of closed curves in SRVF form.
	- `q_to_curve.m` The function that converts unit-norm SRVF of curve into the coordinates of points sampled from a curve.
	- `ReSampleCurve.m` Resample points from a curve with equally spaced points in parameterization with linear interpolation.
	- `ShiftF.m` Shift the curve by a certain amount of unit in parameterization,
	- `SmoothEnergy.m` Compute the $E_{smooth}$ term in the energy functional used for BAC method.
	- `SmoothUpdate.m` Compute the needed gradient of $E_{smooth}$ along current contour in order to update the energy functional.
 - Topological Bayesian Active Contour (TAC).
	 **The full details of TAC method, or combined method is described in [our paper](https://arxiv.org/abs/1910.04778).**
	 - `SegDistTop.m` Take the segmentation results and the truth object boundary and calculate a set of performance evaluation measures.
		 -  `HammingDist.m` Compute the Hamming distance between the truth and segmentation.
		 -  `HausdorffDist.m` Compute the Hausdorff distance between the truth and segmentation.
		 - `JaccardDist.m` Compute the Jaccard distance between the truth and segmentation.
	 - `NeuronScript.m` The script that we use TAC for analyzing the neuron cellular dataset, with additional functionality 
 		 - Gaussian filtering of the original neural cellular image.
		 - Object selection according to the elastic shape distance.
	 - `SimScript.m` The script that we use TAC for analyzing the synthetic MPEG-7 datasets.
	 - `TOPBACSegT.m`  Perform TOP-BAC segmentation with use of training data (e.g., synthetic simulations, skin lesion data). We provide following initialization modes:
	 1. draw contour initializations by hand (default)
	 2. import initialization mask from external file and automatically select n_curve contours with the roughest  estimate of area within
	 3. import initialization mask from external file and cycle through contours until n_curve initializations accepted by user
	 4. input initialization curve from external file
	 - `TOPBACSegNT.m`  Perform TOP-BAC segmentation without use of training data (e.g., synthetic simulations, skin lesion data), with an option of pooling the density estimate of objects.

**Abstract**
A fundamental problem in computer vision is boundary estimation, where the goal is to delineate the boundary of objects in an image.  In this paper, we propose a method which jointly incorporates geometric and topological information within an image to simultaneously estimate boundaries for objects within images with more complex topologies. This method uses a topological method to assist with initializing the active contour algorithm, a method which combines pixel clustering, boundary smoothness, and potential prior shape information to produce an estimated object boundary. Active contour methods are known to be extremely sensitive to algorithm initialization, relying on the user to provide a reasonable starting curve to the algorithm. In the presence of images featuring objects with complex topological structures, such as objects with holes or multiple objects, the user must initialize separate curves for each boundary of interest.  We propose the use of a topological method to provide an automatic initialization in such settings. Our proposed topologically-aware method can produce a smart initialization, freeing up the user from potential pitfalls associated with objects of complex topological structure. We provide a detailed simulation study comparing our initialization to the one obtained from other standard segmentation algorithms. The method is demonstrated on artificially-constructed image datasets from computer vision, as well as applications to skin lesion and neural cellular images, for which multiple topological features can be identified.

**Citation**
We provided MATLAB code for reproducible and experimental purposes under [LICENSE](https://github.com/hrluo/TopologicalBayesianActiveContour).
Please cite our paper using following BibTeX item:

    @article{luo2019combining,
	          title={Combining Geometric and Topological Information in Image Segmentation}, 
	          author={Hengrui Luo and Justin Strait},
	          year={2019},
	          eprint={1910.04778},
	          archivePrefix={arXiv},
	          primaryClass={eess.IV}
    }

Thank you again for the interest and please reach out if you have further questions.
