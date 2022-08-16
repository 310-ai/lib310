from typing import Optional

import scanpy as sc


def _construct_adata(model_results, color: Optional[str] = None, **kwargs):
    adata = sc.AnnData(X=model_results.get_hidden_states(layer=kwargs.get('layer', -1)))
    adata.obs = model_results.obs.copy()

    if color:
        adata.obs['color'] = model_results.obs[color].values
    else:
        pass

    return adata


def convert_column_to_list(adata, column):
    adata.obs[column] = adata.obs[column].apply(lambda cell: ''.join(c for c in cell if c not in "\'\"[]").split(', '))


def umap(model_results, color: Optional[str] = None, plot_type: Optional[str] = None, **kwargs):
    sc.settings.set_figure_params(dpi=kwargs.pop('dpi', 100))
    if plot_type == '1_vs_all':
        adata = _construct_adata(model_results, color=None, layer=kwargs.pop('layer', -1))
    else:
        adata = _construct_adata(model_results, color, layer=kwargs.pop('layer', -1))

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    if plot_type == '1_vs_all':
        convert_column_to_list(adata, color)

        unique_values = set()
        for group in adata.obs[color].to_list():
            unique_values |= set(group)
        unique_values = sorted(list(unique_values))

        for term in unique_values:
            adata.obs[term] = 'Others'
            adata.obs.loc[adata.obs[color].apply(lambda x: term in x), term] = term
            sc.pl.umap(adata, color=term, frameon=False, **kwargs)
    else:
        sc.pl.umap(adata, color='color' if color else None, **kwargs)




def draw_term_family(id):
    import lib310
    import graphviz

    results = lib310.db.fetch(query="""SELECT
      pathList
    FROM
      `pfsdb3.1_go.term_term_path`
    WHERE
      child_id = {id}
    """.format(id=id))
    import graphviz
    from collections import defaultdict

    relationtype = {
        'is_a': {'color': 'black', 'style': 'solid'},
        'part_of': {'color': 'blue', 'style': 'dashed'},
        'regulates': {'color': 'darkorange', 'style': 'bold'},
        'positively_regulates': {'color': 'chartreuse1', 'style': 'bold'},
        'negatively_regulates': {'color': 'red', 'style': 'bold'},
        'occurs_in': {'color': 'darkslategray4', 'style': 'bold'},
    }

    defultrel = {'color': 'gray', 'style': 'dotted'}

    u = graphviz.Digraph('unix', filename='unix.gv')
    u.attr(size='100,100')
    u.attr('node', shape='record')

    edgeList = defaultdict(set)
    rel = defaultdict(dict)
    data = results[['pathList']].values.tolist()
    for row in data:
        paths = row
        pathList = paths[0].strip('][').split(', ')
        x = 0
        while (x + 2 < len(pathList)):
            edgeList[pathList[x]].add(pathList[x + 2])
            rel[pathList[x]][pathList[x + 2]] = pathList[x + 1]
            x = x + 2
    for x in edgeList:
        for y in edgeList[x]:
            z = rel[x][y].replace(' " ', " ")
            z = z.strip(' " " ')
            z = z.replace("'", '')
            if z in relationtype.keys():
                u.edge(x[4:11], y[4:11], color=relationtype[z]['color'], style=relationtype[z]['style'])
            else:
                z = defultrel
                u.edge(x[4:11], y[4:11], color=z['color'], style=z['style'])
    u.view()
    u


def show_umap(features):
    from umap import UMAP
    import plotly.express as px

    umap_2d = UMAP(n_components=2, init='random', random_state=0)
    umap_3d = UMAP(n_components=3, init='random', random_state=0)

    proj_2d = umap_2d.fit_transform(features)
    proj_3d = umap_3d.fit_transform(features)

    fig_2d = px.scatter(
        proj_2d, x=0, y=1,
        color=0,
    )
    fig_3d = px.scatter_3d(
        proj_3d, x=0, y=1, z=2,
        color=0,
    )
    fig_3d.update_traces(marker_size=5)

    fig_2d.show()
    fig_3d.show()


def pca(model_results, color: Optional[str] = None, plot_type: Optional[str] = None, **kwargs):
    sc.settings.set_figure_params(dpi=kwargs.pop('dpi', 100))
    if plot_type == '1_vs_all':
        adata = _construct_adata(model_results, color=None, layer=kwargs.pop('layer', -1))
    else:
        adata = _construct_adata(model_results, color=color, layer=kwargs.pop('layer', -1))

    sc.pp.pca(adata, n_comps=50)

    if plot_type == '1_vs_all':
        convert_column_to_list(adata, color)

        unique_values = set()
        for group in adata.obs[color].to_list():
            unique_values |= set(group)
        unique_values = sorted(list(unique_values))

        for term in unique_values:
            adata.obs[term] = 'Others'
            adata.obs.loc[adata.obs[color].apply(lambda x: term in x), term] = term
            sc.pl.pca(adata, color=term, **kwargs)
    else:
        sc.pl.pca(adata, color='color' if color else None, **kwargs)


def tsne(model_results, color: Optional[str] = None, plot_type: Optional[str] = None, **kwargs):
    sc.settings.set_figure_params(dpi=kwargs.get('dpi', 100))
    if plot_type == '1_vs_all':
        adata = _construct_adata(model_results, color=None, layer=kwargs.pop('layer', -1))
    else:
        adata = _construct_adata(model_results, color, layer=kwargs.pop('layer', -1))

    sc.pp.neighbors(adata)
    sc.tl.tsne(adata)

    if plot_type == '1_vs_all':
        convert_column_to_list(adata, color)

        unique_values = set()
        for group in adata.obs[color].to_list():
            unique_values |= set(group)
        unique_values = sorted(list(unique_values))

        for term in unique_values:
            adata.obs[term] = 'Others'
            adata.obs.loc[adata.obs[color].apply(lambda x: term in x), term] = term
            sc.pl.tsne(adata, color=term, **kwargs)
    else:
        sc.pl.tsne(adata, color='color' if color else None, frameon=False, **kwargs)
