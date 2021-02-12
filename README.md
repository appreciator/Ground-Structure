# Ground-Structure
Generating a ground structure with geometrically orthotropic symmetry for topology optimization  
  
Please cite the following paper when using these Matlab scripts:   
Li, Z., Luo, Z., Zhang, L.-C., & Wang, C.-H. (2021). Topological design of pentamode lattice metamaterials using a ground structure method. Materials & Design, 202, 109523. https://doi.org/10.1016/j.matdes.2021.109523

Example:  
[barNode,nodeCoor,nDesVar] = GenerateGS(5,5,5);  
designVar=zeros(1,405); designVar(300)=1;
isActive = ActiveBar(barNode,nodeCoor,nDesVar,designVar);
