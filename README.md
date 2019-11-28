# Signal Dissection by Correlation Maximization (SDCM)

M. Grau, G. Lenz and P. Lenz, "Dissection of gene expression datasets into clinically
relevant interaction signatures via high-dimensional correlation maximization", 
Nature Communications, 10, 5417 (2019) doi:10.1038/s41467-019-12713-5

## Documentation

**tl;dr** Download and extract the release ZIP or clone this repository. 
Open `selfTest_lowDim.m` or `selfTest_highDim.m` in your Matlab editor. 
Press F5 to run. Then read the in-code docu/comments/config explanations.

### About the method
SDCM is a novel first-principles method that overcomes limitations of principal components 
analysis (PCA), namely linearity and orthogonality. Via its signature functional, it searches
for and maximizes correlation and consistency in the data. It models, 
discovers and dissects overlapping signatures from real interactions via its nonlinear signal model
based on bimonotonic functions. For data focusing, it utilizes weight hyper cones. 
SDCM has already discovered novel gene expression signatures with strong impact on
patient survival. Like PCA, it is mathematically generic and ready to use for new datasets.

### Research Article
#### Abstract
Gene expression is controlled by many simultaneous interactions, frequently measured 
collectively in biology and medicine by high-throughput technologies. It is a highly 
challenging task to infer from these data the generating effects and cooperating genes. 
Here, we present an unsupervised hypothesis-generating learning concept termed signal 
dissection by correlation maximization (SDCM) that dissects large high-dimensional datasets 
into signatures. Each signature captures a particular signal pattern that was consistently 
observed for multiple genes and samples, likely caused by the same underlying interaction. 
A key difference to other methods is our flexible nonlinear signal superposition model, 
combined with a precise regression technique. Analyzing gene expression of diffuse large 
B-cell lymphoma, our method discovers previously unidentified signatures that reveal 
significant differences in patient survival. These signatures are more predictive than 
those from various methods used for comparison and robustly validate across technological 
platforms. This implies highly specific extraction of clinically relevant gene interactions.

#### Paper Link and Methods Outline
For a detailed introduction for introduction, definition and discussion of results by SDCM
for Diffuse Large B-Cell Lymphoma, read our paper at https://www.nature.com/articles/s41467-019-12713-5.
In particular, the methods section provides a comprehensive mathematical description 
of SDCM, organized in the following subsections:
  - Conceptual summary
  - Mathematical framework
  - Model
  - Measure for interactions
  - Algorithm Overview
    - Step 1: Initial representative or termination 
    - Step 2: Signature axes via maximization principle
    - Step 3: Bi-monotonic regression and smoothening
    - Step 4: Signature signal E_k and its dissection
  - Signature focus
  - Signature functional
  - Signature strengths
  - Smoothening operator
  - Bi-monotonic regression
  - Algorithmic complexity
  - and Further Methods.

## Dissecting Your Own Data

### System Requirements
- This SDCM software package is a Matlab toolbox developed with Matlab 
  versions 2013b-R2016b (The MathWorks® Inc., Natick, Massachusetts, United States).
- It should also work with any newer Matlab version (self-tests run fine 
  on R2017a and R2018a on Windows, Linux and MAC).
- SDCM internally requires Matlab's Statistics Toolbox, Bioinformatics Toolbox and 
  Matlab's Parallel Computing Toolbox for parallelization. 
- Additional required tools and helper functions are provided in the `\SDCM_Library` 
  folder for standalone deployment and execution.

### Installation and Self-Test Examples
- As SDCM is a Matlab toolbox, just extract this package to a dedicated folder on 
  your disk, start MATLAB and open either the low-dimensional example and self-test
  `selfTest_lowDim.m` or the high-dimensional example `selfTest_highDim.m` in the
  Matlab editor, then press F5 to execute the script. Required library functions
  are automatically appended to your Matlab library path. (If you want to make this 
  persistent, manually save your Matlab path.)
- Self-tests simulate either the 3D data pattern presented in Fig.1 in the paper or 
  a high-dimensional pattern and then apply SDCM to these data. For each detected
  signature, a definition plot is generated. Finally, signature curves for the 3D
  example or, respectively, correlation matrices comparing simulated with detected
  signature axes are shown for method validation against the high-dimensional example. 
  For method comparison, likewise matrices are plotted for PCA results for the 
  same input data. 
  These self-tests for randomly simulated data can be considered successful, 
  if plotted results resemble those provided in `\SDCM_selfTest_expectedResults` 
  that were obtained in a tested computing environment for identical simulation models.
- Installation and running the low dimensional self-test should be completed within 
  minutes. Depending on your computing equipment, running the high-dimensional 
  self-test might take up to half an hour. This self-test can be configured with
  a lower amount of simulated genes and samples, if needed.
  
### Instructions for New Datasets
- Quick start: To dissect your own data with SDCM, we suggest simply using 
  `selfTest_highDim.m` as a template. Just replace the simulated input data with your data.
- For more detailed control over SDCM, execute 
  `configuration = SDCM_defaultConfig(numberOfGenes, numberOfPatients);`
  The resulting structure contains all configurable SDCM parameters, including
  visualization and export options. Parameters are organized hierarchically as a tree.
- Each parameter is explained in `configuration.explanations.*` at the same tree 
  position (or, alternatively, read the source file `SDCM_defaultConfig.m`).

### Support and Hints
- For support, please open a new issue here at GitHub. If applicable, include or upload a 
  ZIP with your input log2(intensities) matrix and input configuration structure saved 
  as .mat file and provide a Matlab script that reproduces or demonstrates the problem. 
  I will try to debug any reproducible problems on the basis of available time resources.
- Questions that are already answered in our paper, its Supplementery Information, by source 
  code comments or by MATLAB's own documentation have to be ignored; sorry.
- SDCM internally makes use of highly optimized Matlab functions like `ifft` for
  inverse Fourier transforms. Underlying numeric libraries might be upgraded 
  between Matlab versions and might work slightly different on different workstations,
  for example depending on availability of CPU features like AVX. As a consequence, 
  numeric results for the same inputs might not be exactly reproducible between 
  different Matlab versions. To mitigate this problem, we recommend using double 
  precision for all new projects (i.e. do not use the memory saver option 
  `config.preprocessing.numericTargetPrecision = 'single'`) and to enable the 
  new option `config.internal.bEnablePostProductionCodeOptimizations = true` 
  which implements stabilization against numeric roundoff errors in 
  `smoothenInSignatureStrengthsSpace.m` by symmetrization. This will probably 
  become the default in future releases, although it is slightly slower.
- The source code is structured like a document where comments are headers for the
  following code lines. The document's hierarchy is realized by indentation and code 
  folding. When working with the source code, it is recommended to configure Matlab's 
  code folding to also include `if`-structures etc., as I often use them to semantically 
  group lines of code. Then your chosen fold-all hotkey and hotkeys for jumping to 
  the previous/next Matlab code cell should give you a high-level
  code overview and help navigating the code quickly.
