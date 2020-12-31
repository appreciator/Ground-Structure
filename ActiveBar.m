function isActive = ActiveBar(barNode,nodeCoor,nDesVar,designVar)
%Extracting active bars in a ground structure
%Inputs:
%  barNode - Node numbers of each bar element. It is an M-by-2 matrix (M is the number of bars).
%  nodeCoor - Coordinates of each node. It is an N-by-3 matrix (N is the number of nodes).
%  nDesVar - Number of design variables of each bar group. It is an 1-by-5 vector.
%  designVar - A vector of design variables. numel(designVar)==sum(nDesVar)
%Outputs:
%  isActive - Active state of each bar element. It is an M-by-1 boolean vector.
%Copyright (c) 2020-2021 Zuyu Li, Zhen Luo
%Email: lizuyu0091@163.com; zhen.luo@uts.edu.au
%ALL USAGE OF GENERATEGS IS SUBJECT TO LICENSE
%% Active states of bars
if numel(designVar)~=sum(nDesVar)
    error('Wrong number of design variables');
end
designVar = designVar(:)';%Row vector
nBar = size(barNode,1);
nDesA = nDesVar(1); nDesB = nDesVar(2); nDesC = nDesVar(3); nDesD = nDesVar(4);
isActive = zeros(nBar,1,'logical');
isActive(1:nDesA) = designVar(1:nDesA);
isActive(nDesA+1:nDesA+nDesB*2) = repmat(designVar(nDesA+1:nDesA+nDesB)',2,1);
isActive(nDesA+nDesB*2+1:nDesA+nDesB*2+nDesC*4) = repmat(designVar(nDesA+nDesB+1:nDesA+nDesB+nDesC)',4,1);
isActive(nDesA+nDesB*2+nDesC*4+1:nDesA+nDesB*2+nDesC*4+nDesD*4) = repmat(designVar(nDesA+nDesB+nDesC+1:nDesA+nDesB+nDesC+nDesD)',4,1);
isActive(nDesA+nDesB*2+nDesC*4+nDesD*4+1:end) = repmat(designVar(nDesA+nDesB+nDesC+nDesD+1:end)',8,1);
%% Visualization of Ground Structure
figure; hold on; view(3); axis equal; axis off; axis vis3d;
set(gcf,'Color','w','Position',[100,100,1200,600]);
nodeX = nodeCoor(:,1); nodeY = nodeCoor(:,2); nodeZ = nodeCoor(:,3);
line(nodeX(barNode(isActive,:))',nodeY(barNode(isActive,:))',nodeZ(barNode(isActive,:))','Color','r','LineWidth',1,...
    'Marker','o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',4);
end