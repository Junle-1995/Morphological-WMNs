function plotCouplingSeq(curve1,curve2, couplingSeq, pairDist, curveLabels, figOrAxesHandle)
% Make a nice plot illustrating the Frechet couplings and highlighting
% where the maximal Frechet distance is

if nargin<5 || isempty(curveLabels)
    curveLabels = {'curve1','curve2'};
end
if nargin<6 || isempty(figOrAxesHandle)
    figHandle = figure;
    axHandle = axes;
else
    switch figOrAxesHandle.Type
        case 'axes'
            figHandle = ancestor(figOrAxesHandle, 'figure');
            axHandle = figOrAxesHandle;
        case 'figure'
            figHandle = figOrAxesHandle;
            axHandle = axes;
        otherwise
            error('Input "figOrAxesHandle" does not appear to be either!')
    end
end


figure(figHandle);
axes(axHandle);
cla;


x1 = curve1(couplingSeq(:,1),1);
x2 = curve2(couplingSeq(:,2),1); 
y1 = curve1(couplingSeq(:,1),2);
y2 = curve2(couplingSeq(:,2),2);

ind = sub2ind(size(pairDist), couplingSeq(:,1), couplingSeq(:,2));
dists = pairDist(ind);

nans = nan(size(x1));

link_xs = a2v([x1,x2,nans]');
link_ys = a2v([y1,y2,nans]');

links_lh = line(link_xs, link_ys, 'LineWidth',0.5, 'Color',[0.75,0.75,0.75]);

curve1_lh = line(curve1(:,1), curve1(:,2), 'LineWidth',1, 'Color',[0.7,0,0]);
curve2_lh = line(curve2(:,1), curve2(:,2), 'LineWidth',1, 'Color',[0,0,0.7]);

maxMask = dists==max(dists);
maxLink_xs =  a2v([x1(maxMask) x2(maxMask) nans(maxMask)]');
maxLink_ys = a2v([y1(maxMask) y2(maxMask) nans(maxMask)]');
frechet_lh = line(maxLink_xs, maxLink_ys, 'LineWidth',1.5,'color','g');

legend('coupling', curveLabels{:}, 'Frechet distance')
axis equal
grid on
