% <strong>mplot2(t,x,axislabels*,titles*)</strong>
% 
% This function plots several line plots in one figure (subplots). Args with * are
% optional.
% 
% <strong>t</strong> is a 1xN line vector.
%
% <strong>X</strong> is a MxN matrix where line vectors 1-M give the 
% multiple lines to plot.
%
% <strong>axislabels</strong> is a Mx2 string line vector, with the desired 
% labels for each of the axis [x,y]. Defaults to ["t","x1";...;"t","xM"], but 
% requires [] to be passed to the function. Accepts LaTeX.
%
% <strong>titles</strong> works similarly to <strong>axislabels</strong>
% but is only Mx1. Contains the titles of each subplot.

% mplot2 by R.Campos (2025)
% ricardo.jpcampos@ua.pt

function mplot2(t,X,axislabels,titles)
    rows=size(X,1);
    if ~exist('axislabels','var') || isempty(axislabels)
      axislabels=strings(rows,2);
      axislabels(:,1)="$t$";
      for i=1:rows
          axislabels(i,2)=sprintf("$x%d$",i);
      end
    end
    if ~exist('titles','var') || isempty(titles)
        titles=strings(rows,1);
        for i=1:rows
            titles(i,1)=sprintf('Plot %d',i);
        end
    end
    for i=1:rows
        subplot(rows,1,i)
        plot(t,X(i,:))
        xlabel(axislabels(i,1),'interpreter','latex')
        ylabel(axislabels(i,2),'interpreter','latex')
        title(titles(i,1),'interpreter','latex')
    end
end