% choose network name and threshold for showing connections

% name
nm = 'chondrogenic';
full_nm = 'chondrogenic';

% threshold
tr = 0.1;

% amount of network measures to show
N = 5;

% load files
flnm = [nm,'_GENIEnet.mat'];
gnm = [nm,'_GENIEgenes.mat'];
load(flnm)
load(gnm)

% construct and plot graph
x(isnan(x)) = 0;

A = x > tr;

G = digraph(x);
Gp = digraph(A);
plot(Gp,'-.sr','NodeLabel',g,'LineWidth',1)
title(['Most important edges in network generated for ', full_nm, ' cells'])
set(gcf,'Position',[100 100 900 900])
print(gcf,[nm,'_net.png'],'-r700','-dpng');

% network measures
pr = centrality(G,'pagerank');
show_rank(pr,g,'pagerank',N)
bc = centrality(G,'betweenness');
show_rank(bc,g,'out closeness',N)

function show_rank(x,g,nm,N)
    [~,I] = maxk(x,N);
    pr_gene_list = g(I);
    disp(['Genes ranked by ', nm,' for GENIE network']);
    for k=1:length(pr_gene_list)
        nm = pr_gene_list(k);
        disp(nm{1});  
    end
end


