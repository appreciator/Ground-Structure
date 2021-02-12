function [barNode,nodeCoor,nDesVar] = GenerateGS(nNodeX,nNodeY,nNodeZ)
%Generating a ground structure with geometrically orthotropic symmetry
%
%Inputs:
%  nNodeX - Number of node columns along the X-axis. It must be an odd number >1 .
%  nNodeY - Number of node columns along the Y-axis. It must be an odd number >1 .
%  nNodeZ - Number of node columns along the Z-axis. It must be an odd number >1 .
%Outputs:
%  barNode - Node numbers of each bar element. It is an M-by-2 matrix (M is the number of bars).
%  nodeCoor - Coordinates of each node. It is an N-by-3 matrix (N is the number of nodes).
%  nDesVar - Number of design variables corresponding to each bar group. It is an 1-by-5 vector.
%
%Please cite the following paper when using this Matlab script:
%https://doi.org/10.1016/j.matdes.2021.109523
%
%Copyright (c) 2020 Zuyu Li, Zhen Luo
%Email: lizuyu0091@163.com; zhen.luo@uts.edu.au
%ALL USAGE OF GENERATEGS IS SUBJECT TO LICENSE

%% Generation of Base Mesh
nEX = nNodeX-1; nEY = nNodeY-1; nEZ = nNodeZ-1; 
nEX = ceil(nEX/2)*2; nEY = ceil(nEY/2)*2; nEZ = ceil(nEZ/2)*2;
[X,Y,Z] = meshgrid(-nEX/2:1:0,-nEY/2:1:0,-nEZ/2:1:0);
nodeCoor = [X(:),Y(:),Z(:);X(:),-Y(:),Z(:);-X(:),Y(:),Z(:);-X(:),-Y(:),Z(:)];
nodeCoor = [nodeCoor;[nodeCoor(:,1:2),-nodeCoor(:,3)]];
nNode = size(nodeCoor,1);
%% Generation of Ground Structure
nodeXL = (nEY/2+1)*(nEX/2+1)*nEZ/2+(nEY/2+1:nEY/2+1:(nEY/2+1)*(nEX/2+1))';
nNodeXL = numel(nodeXL);
[barXL_j,barXL_i] = find(tril(ones(nNodeXL,nNodeXL))-eye(nNodeXL));

nodeYD = (nEY/2+1)*(nEX/2+1)*nEZ/2+(nEY/2+1)*nEX/2+(1:nEY/2+1)';
nNodeYD = numel(nodeYD);
[barYD_j,barYD_i] = find(tril(ones(nNodeYD,nNodeYD))-eye(nNodeYD));

nodeZN = ((nEY/2+1)*(nEX/2+1):(nEY/2+1)*(nEX/2+1):(nEY/2+1)*(nEX/2+1)*(nEZ/2+1))';
nNodeZN = numel(nodeZN);
[barZN_j,barZN_i] = find(tril(ones(nNodeZN,nNodeZN))-eye(nNodeZN));

nodeLDXY = (nEY/2+1)*(nEX/2+1)*nEZ/2+(1:(nEY/2+1)*nEX/2)'; 
nodeLDXY(nEY/2+1:nEY/2+1:end)=[];

nodeLDYZ = bsxfun(@plus,(nEY/2+1)*nEX/2+(1:nEY/2)',(nEY/2+1)*(nEX/2+1)*(0:1:(nEZ/2-1)));
nodeLDYZ = nodeLDYZ(:);

nodeLDXZ = bsxfun(@plus,((nEY/2+1):(nEY/2+1):(nEY/2+1)*nEX/2)',(nEY/2+1)*(nEX/2+1)*(0:1:(nEZ/2-1)));
nodeLDXZ = nodeLDXZ(:);

%Group A
barNodeA = [nodeXL(1:end-1),nodeXL(1:end-1)+nNode/4;
            nodeYD(1:end-1),nodeYD(1:end-1)+nNode/8;
            nodeZN(1:end-1),nodeZN(1:end-1)+nNode/2];

%Group B        
barNodeB1 = [nodeXL(barXL_i),nodeXL(barXL_j);nodeYD(barYD_i),nodeYD(barYD_j);
             nodeZN(barZN_i),nodeZN(barZN_j);];      
barNodeB2 = nNode*[repmat([1/4,1/4],nEX^2/8+nEX/4,1);repmat([1/8,1/8],nEY^2/8+nEY/4,1);
            repmat([1/2,1/2],nEZ^2/8+nEZ/4,1)]+barNodeB1;               
barNodeB1XY = [nodeLDXY,nodeLDXY+nNode/4;nodeLDXY,nodeLDXY+nNode/8];
barNodeB2XY = nNode*[repmat([1/8,1/8],nEX*nEY/4,1);repmat([1/4,1/4],nEX*nEY/4,1)]+barNodeB1XY;
barNodeB1YZ = [nodeLDYZ,nodeLDYZ+nNode/8;nodeLDYZ,nodeLDYZ+nNode/2];
barNodeB2YZ = nNode*[repmat([1/2,1/2],nEY*nEZ/4,1);repmat([1/8,1/8],nEY*nEZ/4,1)]+barNodeB1YZ;         
barNodeB1XZ = [nodeLDXZ,nodeLDXZ+nNode/4;nodeLDXZ,nodeLDXZ+nNode/2];
barNodeB2XZ = nNode*[repmat([1/2,1/2],nEX*nEZ/4,1);repmat([1/4,1/4],nEX*nEZ/4,1)]+barNodeB1XZ;   
barNodeB1 = [barNodeB1;barNodeB1XY;barNodeB1YZ;barNodeB1XZ];
barNodeB2 = [barNodeB2;barNodeB2XY;barNodeB2YZ;barNodeB2XZ];

%Group C
nodeXY = (nEY/2+1)*(nEX/2+1)*nEZ/2+(1:(nEY/2+1)*(nEX/2+1))'; 
nodeYZ = bsxfun(@plus,(nEY/2+1)*nEX/2+(1:(nEY/2+1))',(nEY/2+1)*(nEX/2+1)*(0:1:nEZ/2));
nodeYZ = nodeYZ(:);
nodeXZ = bsxfun(@plus,((nEY/2+1):(nEY/2+1):(nEY/2+1)*(nEX/2+1))',(nEY/2+1)*(nEX/2+1)*(0:1:nEZ/2));
nodeXZ = nodeXZ(:);
[bar_j,bar_i] = find(tril(ones((nEY/2+1)*(nEX/2+1),(nEY/2+1)*(nEX/2+1)))-eye((nEY/2+1)*(nEX/2+1)));
barNodeC1XY = setdiff([nodeXY(bar_i),nodeXY(bar_j)],barNodeB1,'rows');
barNodeC2XY = nNode/8+barNodeC1XY;
barNodeC3XY = nNode*2/8+barNodeC1XY;
barNodeC4XY = nNode*3/8+barNodeC1XY;
[bar_j,bar_i] = find(tril(ones((nEY/2+1)*(nEZ/2+1),(nEY/2+1)*(nEZ/2+1)))-eye((nEY/2+1)*(nEZ/2+1)));
barNodeC1YZ = setdiff([nodeYZ(bar_i),nodeYZ(bar_j)],barNodeB1,'rows');
barNodeC2YZ = nNode/8+barNodeC1YZ;
barNodeC3YZ = nNode*4/8+barNodeC1YZ;
barNodeC4YZ = nNode*5/8+barNodeC1YZ;
[bar_j,bar_i] = find(tril(ones((nEZ/2+1)*(nEX/2+1),(nEZ/2+1)*(nEX/2+1)))-eye((nEZ/2+1)*(nEX/2+1)));
barNodeC1XZ = setdiff([nodeXZ(bar_i),nodeXZ(bar_j)],barNodeB1,'rows');
barNodeC2XZ = nNode*2/8+barNodeC1XZ;
barNodeC3XZ = nNode*4/8+barNodeC1XZ;
barNodeC4XZ = nNode*6/8+barNodeC1XZ;
barNodeC1 = [barNodeC1XY;barNodeC1YZ;barNodeC1XZ];
barNodeC2 = [barNodeC2XY;barNodeC2YZ;barNodeC2XZ];
barNodeC3 = [barNodeC3XY;barNodeC3YZ;barNodeC3XZ];
barNodeC4 = [barNodeC4XY;barNodeC4YZ;barNodeC4XZ];

%Group D
nodeLDN = 1:(nEY/2+1)*nEX/2;
nodeLDN(nEY/2+1:nEY/2+1:end)=[];
nodeLDN = bsxfun(@plus,nodeLDN',(nEY/2+1)*(nEX/2+1)*(0:1:(nEZ/2-1)));
nodeLDN = nodeLDN(:);
barNodeD1 = [nodeLDN,nodeLDN+nNode*2/8;nodeLDN,nodeLDN+nNode/8;nodeLDN,nodeLDN+nNode*4/8];
barNodeD2 = [nodeLDN+nNode/8,nodeLDN+nNode*2/8+nNode/8;
             nodeLDN+nNode*2/8,nodeLDN+nNode/8+nNode*2/8;
             nodeLDN+nNode/8,nodeLDN+nNode*4/8+nNode/8];
barNodeD3 = [nodeLDN+nNode*4/8,nodeLDN+nNode*2/8+nNode*4/8;
             nodeLDN+nNode*4/8,nodeLDN+nNode/8+nNode*4/8;
             nodeLDN+nNode*2/8,nodeLDN+nNode*4/8+nNode*2/8];
barNodeD4 = [nodeLDN+nNode*5/8,nodeLDN+nNode*2/8+nNode*5/8;
             nodeLDN+nNode*6/8,nodeLDN+nNode/8+nNode*6/8;
             nodeLDN+nNode*3/8,nodeLDN+nNode*4/8+nNode*3/8];

%Group E
[bar_j,bar_i] = find(tril(ones(nNode/8,nNode/8))-eye(nNode/8));
barNodeE = setdiff([bar_i,bar_j],[barNodeB1;barNodeC1XY;barNodeC1YZ;barNodeC1XZ],'rows');
barNodeE = repmat(barNodeE,8,1)+kron(nNode/8*(0:7)',ones(size(barNodeE)));
%% Merge and Renumber Nodes
[nodeCoor,~,indexn] = unique(nodeCoor(:,[3,1,2]),'rows');
nodeCoor = nodeCoor(:,[2,3,1]);
barNodeA = indexn(barNodeA); 
barNodeB = indexn([barNodeB1;barNodeB2]);
barNodeC = indexn([barNodeC1;barNodeC2;barNodeC3;barNodeC4]);
barNodeD = indexn([barNodeD1;barNodeD2;barNodeD3;barNodeD4]);
barNodeE = indexn(barNodeE);
barNode = sort([barNodeA;barNodeB;barNodeC;barNodeD;barNodeE],2);
%% Information of Design Variables
nDesVar = [size(barNodeA,1),size(barNodeB1,1),size(barNodeC1,1),size(barNodeD1,1),size(barNodeE,1)/8];
%% Visualization of Ground Structure
figure; hold on; view(3); axis equal; axis off; axis vis3d;
set(gcf,'Color','w','Position',[100,100,1200,600]);
nodeX = nodeCoor(:,1); nodeY = nodeCoor(:,2); nodeZ = nodeCoor(:,3);
line(nodeX(barNode)',nodeY(barNode)',nodeZ(barNode)','Color','r','LineWidth',1,...
    'Marker','o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',4);
end
