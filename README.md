


# TopologicalBayesianActiveContour

**Content**
This is the code repository for the research publication "Combining Geometric and Topological Information in Image Segmentation" by [Justin D. Strait](https://jdstrait.weebly.com/) and [Hengrui Luo](https://hrluo.github.io/). 
The manuscript of this paper can be accessed at https://arxiv.org/abs/1910.04778. 

 - The TOP code we store in [TOP folder](https://github.com/hrluo/TopologicalBayesianActiveContour/tree/master/TOP) folder for comparison convenience is associated with the paper [A Topological Approach to Hierarchical Segmentation using Mean Shift](https://people.csail.mit.edu/sparis/). The copyright of this code belongs to the original author.
 - The TAC code we store in [TAC folder](https://github.com/hrluo/TopologicalBayesianActiveContour/tree/master/TAC) folder stores the code for both BAC (Bayesian Active Contour) and TAC (Combined method) described in details in [our paper](https://arxiv.org/abs/1910.04778). The copyright of this code belongs to [Justin D. Strait](https://jdstrait.weebly.com/) and [Hengrui Luo](https://hrluo.github.io/). 

**Structure**

**Abstract**
A fundamental problem in computer vision is image segmentation, where the goal is to delineate the boundary of an object in the image. The focus of this work is on the segmentation of grayscale images and its purpose is two-fold. First, we conduct an in-depth study comparing active contour and topology-based methods in a statistical framework, two popular approaches for boundary detection of 2-dimensional images. Certain properties of the image dataset may favor one method over the other, both from an interpretability perspective as well as through evaluation of performance measures. Second, we propose the use of topological knowledge to assist an active contour method, which can potentially incorporate prior shape information. The latter is known to be extremely sensitive to algorithm initialization, and thus, we use a topological model to provide an automatic initialization. In addition, our proposed model can handle objects in images with more complex topological structures, including objects with holes and multiple objects within one image. We demonstrate this on artificially-constructed image datasets from computer vision, as well as real medical image data.

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
