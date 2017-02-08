function [serr, fig, h] = sub_conv(AS,k,sstype,Nbs,fig)
W1       = AS.sub.W1;
serr.err = zeros(length(k),3);
serr.eig = serr.err;
err      = zeros(Nbs,1);
eigr     = err;
N        = size(AS.X,1);

for i=1:length(k)
    for ii=1:Nbs
        ind   = randi(N,k(i),1);
        sub   = compute(AS.X(ind,:),AS.F(ind),[],[],sstype,0,0);
       
        if size(W1,2) ~= sub.W1
            W = [sub.W1, sub.W2];
            sub.W1 = W(:,1:size(W1,2));
            sub.W2 = W(:,size(W1,2)+1:end);
        end
        
        eigr(ii)  = norm(AS.sub.eigenvalues - sub.eigenvalues);
        err(ii)   = norm(W1'*sub.W2);
    end
    serr.eig(i,:)   = [min(eigr),mean(eigr),max(eigr)];
    serr.err(i,:) = [min(err),mean(err),max(err)];
end

%% Plot
% Get plotting options.
opts = plot_opts([]);
opts.fontsize = 16;

if isempty(fig)
    fig = figure();
else
    fig = figure(fig.Number);
end

h = loglog(k,serr.err(:,2), ...
         'marker', opts.marker, ...
         'markersize', opts.markersize, ...
         'linewidth', opts.linewidth);
hold on

% Plot bootstrap errors.
if max(abs(serr.err(:,3) - serr.err(:,1))) > 1e-7
    fillh = fill([k, k(end:-1:1)], [serr.err(:,1)', fliplr(serr.err(:,3)')], h.Color);
    hAnnotation = get(fillh,'Annotation');
    hLegendEntry = get(hAnnotation','LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off')
    set(fillh,'facealpha',.25);
end
h = loglog(k,serr.err(:,2), ...
         'color',h.Color,...
         'MarkerEdgeColor',h.Color,...
         'MarkerFaceColor',h.Color,...
         'marker', opts.marker, ...
         'markersize', opts.markersize, ...
         'linewidth', opts.linewidth);
title('Convergence','fontsize',14)
ylabel('Subspace Error','fontsize',14)
xlabel('Sample Size','fontsize',14)
grid on;
hold on;

set(gca, ...
    'XLim', [k(1), k(end)], ...
    'XTick', k, ...
    'XScale', 'Log', ...
    'YScale', 'Log', ...
    'fontsize', opts.fontsize)

