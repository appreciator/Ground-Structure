function [barNode,nodeCoor] = GenerateGS(nelx,nely,nelz)
%Generating a ground structure with geometrically orthotropic symmetry
%Developed by Zuyu Li and Zhen Luo
%Email:zuyu.li@student.uts.edu.au; zhen.luo@uts.edu.au
%December 27, 2020
%% GENERATION OF BASE MESH
nelx = ceil(nelx/2)*2; nely = ceil(nely/2)*2; nelz = ceil(nelz/2)*2;
[X,Y,Z] = meshgrid(-nelx/2:1:0,-nely/2:1:0,-nelz/2:1:0);
nodeCoor = [X(:),Y(:),Z(:);X(:),-Y(:),Z(:);-X(:),Y(:),Z(:);-X(:),-Y(:),Z(:)];
nodeCoor = [nodeCoor;[nodeCoor(:,1:2),-nodeCoor(:,3)]];
nNode = size(nodeCoor,1);
%% GENERATION OF GROUND STRUCTURE
nodeXL = (nely/2+1)*(nelx/2+1)*nelz/2+(nely/2+1:nely/2+1:(nely/2+1)*(nelx/2+1))';
nNodeXL = numel(nodeXL);
[barXL_j,barXL_i] = find(tril(ones(nNodeXL,nNodeXL))-eye(nNodeXL));

nodeYD = (nely/2+1)*(nelx/2+1)*nelz/2+(nely/2+1)*nelx/2+(1:nely/2+1)';
nNodeYD = numel(nodeYD);
[barYD_j,barYD_i] = find(tril(ones(nNodeYD,nNodeYD))-eye(nNodeYD));

nodeZN = ((nely/2+1)*(nelx/2+1):(nely/2+1)*(nelx/2+1):(nely/2+1)*(nelx/2+1)*(nelz/2+1))';
nNodeZN = numel(nodeZN);
[barZN_j,barZN_i] = find(tril(ones(nNodeZN,nNodeZN))-eye(nNodeZN));

nodeLDXY = (nely/2+1)*(nelx/2+1)*nelz/2+(1:(nely/2+1)*nelx/2)'; 
nodeLDXY(nely/2+1:nely/2+1:end)=[];

nodeLDYZ = bsxfun(@plus,(nely/2+1)*nelx/2+(1:nely/2)',(nely/2+1)*(nelx/2+1)*(0:1:(nelz/2-1)));
nodeLDYZ = nodeLDYZ(:);

nodeLDXZ = bsxfun(@plus,((nely/2+1):(nely/2+1):(nely/2+1)*nelx/2)',(nely/2+1)*(nelx/2+1)*(0:1:(nelz/2-1)));
nodeLDXZ = nodeLDXZ(:);


barNodeA = [nodeXL(1:end-1),nodeXL(1:end-1)+nNode/4;
            nodeYD(1:end-1),nodeYD(1:end-1)+nNode/8;
            nodeZN(1:end-1),nodeZN(1:end-1)+nNode/2];%Œﬁ÷ÿ∏¥
        
barNodeB1 = [nodeXL(barXL_i),nodeXL(barXL_j);nodeYD(barYD_i),nodeYD(barYD_j);
             nodeZN(barZN_i),nodeZN(barZN_j);];%X2         
barNodeB2 = nNode*[repmat([1/4,1/4],nelx^2/8+nelx/4,1);repmat([1/8,1/8],nely^2/8+nely/4,1);
            repmat([1/2,1/2],nelz^2/8+nelz/4,1)]+barNodeB1;%X2         
         
barNodeB1XY = [nodeLDXY,nodeLDXY+nNode/4;nodeLDXY,nodeLDXY+nNode/8];%X2
barNodeB2XY = nNode*[repmat([1/8,1/8],nelx*nely/4,1);repmat([1/4,1/4],nelx*nely/4,1)]+barNodeB1XY;%X2 

barNodeB1YZ = [nodeLDYZ,nodeLDYZ+nNode/8;nodeLDYZ,nodeLDYZ+nNode/2];%X2
barNodeB2YZ = nNode*[repmat([1/2,1/2],nely*nelz/4,1);repmat([1/8,1/8],nely*nelz/4,1)]+barNodeB1YZ;%X2
         
barNodeB1XZ = [nodeLDXZ,nodeLDXZ+nNode/4;nodeLDXZ,nodeLDXZ+nNode/2];%X2
barNodeB2XZ = nNode*[repmat([1/2,1/2],nelx*nelz/4,1);repmat([1/4,1/4],nelx*nelz/4,1)]+barNodeB1XZ;%X2                          

nodeXY = (nely/2+1)*(nelx/2+1)*nelz/2+(1:(nely/2+1)*(nelx/2+1))'; 

nodeYZ = bsxfun(@plus,(nely/2+1)*nelx/2+(1:(nely/2+1))',(nely/2+1)*(nelx/2+1)*(0:1:nelz/2));
nodeYZ = nodeYZ(:);

nodeXZ = bsxfun(@plus,((nely/2+1):(nely/2+1):(nely/2+1)*(nelx/2+1))',(nely/2+1)*(nelx/2+1)*(0:1:nelz/2));
nodeXZ = nodeXZ(:);

[bar_j,bar_i] = find(tril(ones((nely/2+1)*(nelx/2+1),(nely/2+1)*(nelx/2+1)))-eye((nely/2+1)*(nelx/2+1)));
barNodeC1XY = setdiff([nodeXY(bar_i),nodeXY(bar_j)],barNodeB1,'rows');%X4
barNodeC2XY = nNode/8+barNodeC1XY;%X4
barNodeC3XY = nNode*2/8+barNodeC1XY;%X4
barNodeC4XY = nNode*3/8+barNodeC1XY;%X4

[bar_j,bar_i] = find(tril(ones((nely/2+1)*(nelz/2+1),(nely/2+1)*(nelz/2+1)))-eye((nely/2+1)*(nelz/2+1)));
barNodeC1YZ = setdiff([nodeYZ(bar_i),nodeYZ(bar_j)],barNodeB1,'rows');%X4
barNodeC2YZ = nNode/8+barNodeC1YZ;%X4
barNodeC3YZ = nNode*4/8+barNodeC1YZ;%X4
barNodeC4YZ = nNode*5/8+barNodeC1YZ;%X4

[bar_j,bar_i] = find(tril(ones((nelz/2+1)*(nelx/2+1),(nelz/2+1)*(nelx/2+1)))-eye((nelz/2+1)*(nelx/2+1)));
barNodeC1XZ = setdiff([nodeXZ(bar_i),nodeXZ(bar_j)],barNodeB1,'rows');%X4
barNodeC2XZ = nNode*2/8+barNodeC1XZ;%X4
barNodeC3XZ = nNode*4/8+barNodeC1XZ;%X4
barNodeC4XZ = nNode*6/8+barNodeC1XZ;%X4

nodeLDN = 1:(nely/2+1)*nelx/2;
nodeLDN(nely/2+1:nely/2+1:end)=[];
nodeLDN = bsxfun(@plus,nodeLDN',(nely/2+1)*(nelx/2+1)*(0:1:(nelz/2-1)));
nodeLDN = nodeLDN(:);
barNodeD1 = [nodeLDN,nodeLDN+nNode*2/8;nodeLDN,nodeLDN+nNode/8;nodeLDN,nodeLDN+nNode*4/8];%X4
barNodeD2 = [nodeLDN+nNode/8,nodeLDN+nNode*2/8+nNode/8;nodeLDN+nNode*2/8,nodeLDN+nNode/8+nNode*2/8;nodeLDN+nNode/8,nodeLDN+nNode*4/8+nNode/8];%X4
barNodeD3 = [nodeLDN+nNode*4/8,nodeLDN+nNode*2/8+nNode*4/8;nodeLDN+nNode*4/8,nodeLDN+nNode/8+nNode*4/8;nodeLDN+nNode*2/8,nodeLDN+nNode*4/8+nNode*2/8];%X4
barNodeD4 = [nodeLDN+nNode*5/8,nodeLDN+nNode*2/8+nNode*5/8;nodeLDN+nNode*6/8,nodeLDN+nNode/8+nNode*6/8;nodeLDN+nNode*3/8,nodeLDN+nNode*4/8+nNode*3/8];%X4

[bar_j,bar_i] = find(tril(ones(nNode/8,nNode/8))-eye(nNode/8));
barNode = setdiff([bar_i,bar_j],[barNodeB1;barNodeC1XY;barNodeC1YZ;barNodeC1XZ],'rows');
barNode = repmat(barNode,8,1)+kron(nNode/8*(0:7)',ones(size(barNode)));

barNodeB1 = [barNodeB1;barNodeB1XY;barNodeB1YZ;barNodeB1XZ];
barNodeB2 = [barNodeB2;barNodeB2XY;barNodeB2YZ;barNodeB2XZ];

barNodeC1 = [barNodeC1XY;barNodeC1YZ;barNodeC1XZ];
barNodeC2 = [barNodeC2XY;barNodeC2YZ;barNodeC2XZ];
barNodeC3 = [barNodeC3XY;barNodeC3YZ;barNodeC3XZ];
barNodeC4 = [barNodeC4XY;barNodeC4YZ;barNodeC4XZ];
%% Merge and Renumber Nodes
[nodeCoor,~,indexn] = unique(nodeCoor(:,[3,1,2]),'rows');
nodeCoor = nodeCoor(:,[2,3,1]);
barNodeA = indexn(barNodeA);
barNodeB1 = indexn(barNodeB1);
barNodeB2 = indexn(barNodeB2);
barNodeC1 = indexn(barNodeC1);
barNodeC2 = indexn(barNodeC2);
barNodeC3 = indexn(barNodeC3);
barNodeC4 = indexn(barNodeC4);
barNodeD1 = indexn(barNodeD1);
barNodeD2 = indexn(barNodeD2);
barNodeD3 = indexn(barNodeD3);
barNodeD4 = indexn(barNodeD4);
barNode = sort([barNodeA;barNodeB1;barNodeB2;barNodeC1;barNodeC2;barNodeC3;barNodeC4;barNodeD1;barNodeD2;barNodeD3;barNodeD4;indexn(barNode)],2);
%% Information of Ground Structure
nBar = size(barNode,1);
nBarA = size(barNodeA,1);
nBarB = size(barNodeB1,1);
nBarC = size(barNodeC1,1);
nBarD = size(barNodeD1,1);
end

function isActive = ActiveBar(designVar,nBar,nBarA,nBarB,nBarC,nBarD)
%% INDEX OF ACTIVE BARS
isActive = zeros(nBar,1,'logical');
isActive(1:nBarA) = designVar(1:nBarA);
isActive(nBarA+1:nBarA+nBarB*2) = repmat(designVar(nBarA+1:nBarA+nBarB)',2,1);
isActive(nBarA+nBarB*2+1:nBarA+nBarB*2+nBarC*4) = repmat(designVar(nBarA+nBarB+1:nBarA+nBarB+nBarC)',4,1);
isActive(nBarA+nBarB*2+nBarC*4+1:nBarA+nBarB*2+nBarC*4+nBarD*4) = repmat(designVar(nBarA+nBarB+nBarC+1:nBarA+nBarB+nBarC+nBarD)',4,1);
isActive(nBarA+nBarB*2+nBarC*4+nBarD*4+1:end) = repmat(designVar(nBarA+nBarB+nBarC+nBarD+1:end)',8,1);
end