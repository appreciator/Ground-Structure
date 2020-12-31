# Ground-Structure
Generating a ground structure with geometrically orthotropic symmetry for topology optimization  
  
Matlab scripts used for the paper submmited to JMADE:   
Topological Design of Pentamode Lattice Metamaterials Using A Ground Structure Method 

Example:  
[barNode,nodeCoor,nDesVar] = GenerateGS(5,5,5);  
designVar=zeros(1,405); designVar(300)=1;  
isActive = ActiveBar(barNode,nodeCoor,nDesVar,designVar);
