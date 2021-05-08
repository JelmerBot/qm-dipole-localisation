function add_subplot_grid_labels(rows, cols, row_labels, col_labels)

% Resize the figure before using this!
% Add a suptitle after using this

% Setup some plots
m=rows;n=cols;  % Number of plots to make
padding = [0.3 0.1]; % Determines space between labels and plots
handle = gcf;

tmp_fig = figure;
set(tmp_fig, 'Units',  get(handle, 'Units'))
set(tmp_fig, 'Position', get(handle, 'Position'))

for j=1:m
    for k=1:n
        figure(tmp_fig.Number)
%         new_ax = subplot(m,n,sub2ind([n,m],k,j));
               
        % Add labels
        if j==1 && ~isempty(col_labels{k})% Top labels
            figure(handle.Number);
            subplot(m,n,sub2ind([n, m], k, j));
            target_ax = gca;
            copyobj(target_ax, tmp_fig)
            figure(tmp_fig.Number)
            view(2)
            
            htmp=xlabel(gca, col_labels{k});
            htmp.Units='normalized';
            htmp.Position(2)= 1+padding(2);
            
            copyobj(htmp,target_ax);
        end
        if k==1 && ~isempty(row_labels{j}) % Left Labels
            figure(handle.Number);
            subplot(m,n,sub2ind([n, m], k, j));
            target_ax = gca;
            copyobj(target_ax, tmp_fig)
            figure(tmp_fig.Number)
            view(2)
            
            htmp=ylabel(gca, row_labels{j});
            htmp.Units='normalized';
            htmp.Position(1)= -padding(1);
            
            copyobj(htmp,target_ax);
        end
    end
end
close(tmp_fig); % Close the temporary axes
end