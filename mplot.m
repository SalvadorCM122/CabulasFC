% <strong>mplot(t,x,labels*,Title*)</strong>
%
% This function plots several line plots in one figure. Args with * are
% optional.
% 
% <strong>t</strong> is a 1xN line vector.
%
% <strong>x</strong> is a MxN matrix where line vectors 1-M give the 
% multiple lines to plot.
% 
% <strong>labels</strong> is a 1xM string line vector, with the desired 
% labels for each of the lines. Defaults to ["x1",...,"xm"], but 
% requires [] to be passed to the function. Accepts LaTeX.
% 
% <strong>Title</strong> is the title of the plot. Defaults to x(t). It is
% rendered with LaTeX.
%
% Examples: 
% mplot(t,x)                  uses default labels and title.
% mplot(t,x, [], Title)       uses default labels and custom title
% mplot(t,x,labels,Title)     all user-defined

% mplot by R.Campos (2025)
% ricardo.jpcampos@ua.pt

function mplot(t,x,labels,Title)
    if ~exist('labels','var') || isempty(labels)
       labels=strings(1,size(x,1));
       for i=1:size(x,1)
            label=sprintf('x%d',i);
            labels(i)=label;
       end
    end
    if ~exist('Title','var') || isempty(Title)
        Title='$\textbf{x(t)}$';
    end
    figure 'Name' 'posplot Figure'
    for i = 1:size(x,1)
        plot(t,x(i,:),'DisplayName',labels(i))
        hold on
    end
    hold off
    legend()
    drawnow
    title(Title,'interpreter','latex')
    xlabel('$\textbf{t}$','interpreter','latex')
end