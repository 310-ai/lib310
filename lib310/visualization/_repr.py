from collections import defaultdict
from typing import Optional

import scanpy as sc



class Graph:
    def __init__(self, vertices):
        self.graph = defaultdict(list)
        self.ids = defaultdict(list)
        self.V = vertices
        self.visit = {}
        self.label = {}
        self.paths = defaultdict(list)

    def add_edge(self, u, v):
        self.graph[u].append(v)
        self.ids[u] = u
        self.ids[v] = v
        self.visit[u] = False
        self.visit[v] = False

    def dfs_util(self, v, label):
        if self.label.get(v) is None:
            self.label[v] = label
        else:
            self.label[v] = max(self.label[v],label)

        for x in self.graph[v]:
            self.dfs_util(x,self.label[v]+1)



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
    query = """SELECT
          pathList
        FROM
          `pfsdb3.1_go.term_term_path`
        WHERE
          child_id = "{id}";""".format(id=id)
    results = lib310.db.fetch(query=query)



    results = lib310.db.fetch(query="""SELECT
      pathList
    FROM
      `pfsdb3.1_go.term_term_path`
    WHERE
      child_id = "{id}"
    """.format(id=id))
    # import graphviz
    from collections import defaultdict

    g = Graph(len(results))
    relationtype = {
        'is_a': {'color': 'black', 'style': 'bold'},
        'part_of': {'color': 'blue', 'style': 'dashed'},
        'regulates': {'color': 'darkorange', 'style': 'bold'},
        'positively_regulates': {'color': 'darkgreen', 'style': 'bold'},
        'negatively_regulates': {'color': 'red', 'style': 'bold'},
        'occurs_in': {'color': 'darkslategray4', 'style': 'bold'},
    }

    defultrel = {'color': 'gray', 'style': 'dotted'}

    nodes = {}
    u = graphviz.Digraph('unix', filename='unix.gv')
    u.attr(size='1000,1000')
    u.attr(ranksep="2")
    u.attr('node', shape='plaintext')
    u.attr('node', font='helvetica')

    res = lib310.db.fetch(query="""SELECT go_id,name FROM `pfsdb3.1_go.terms`  """)
    node_names = {}
    for row in res.values.tolist():
        node_names[row[0]] = row[1]
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
            X = "" + x[4:11]
            Y = "" + y[4:11]
            if X not in nodes:
                nodes[X] = node_names["GO:" + X]
            if Y not in nodes:
                nodes[Y] = node_names["GO:" + Y]
            if z in relationtype.keys():
                u.edge(x[4:11], y[4:11], color=relationtype[z]['color'], style=relationtype[z]['style'])
                g.add_edge(y[4:11], x[4:11])
            else:
                z = defultrel
                u.edge(x[4:11], y[4:11], color=z['color'], style=z['style'])
                # g.add_edge( y[4:11],x[4:11])

    g.dfs_util(id[3:], 1)

    i = 1
    while i < 1000:
        with u.subgraph() as s:
            s.attr(rank='same')

            for x in g.label.keys():
                if g.label[x] == i:
                    # s.node(x)
                    name = "GO:" + str(x)
                    s.node(x, label='''<
                    <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="2" CELLPADDING="3">
                      <TR>

                        <TD BGCOLOR="deepskyblue3" COLSPAN="1">{id}</TD>
                      </TR>

                      <TR>
                        <TD  COLSPAN="1">{name}</TD>
                      </TR>
                    </TABLE>>'''.format(id=x, name=node_names[name]))

                    # print(str(x)+': '+str(i))

        i = i + 1

    for node in nodes:
        if node not in g.label.keys():
            u.node(node, label='''<
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="2" CELLPADDING="3">
          <TR>

            <TD BGCOLOR="deepskyblue3" COLSPAN="1">{id}</TD>
          </TR>

          <TR>
            <TD  COLSPAN="1">{name}</TD>
          </TR>
        </TABLE>>'''.format(id=node, name=node_names["GO:" + node]))
    return u


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


def visualize_3model(fpath, chain_id):
        import py3Dmol

        with open(fpath) as ifile:
            system = "".join([x for x in ifile])

        view = py3Dmol.view(width=600, height=400)
        view.addModelsAsFrames(system)
        view.setStyle({'model': -1, 'chain': chain_id}, {"cartoon": {'color': 'spectrum'}})
        view.setBackgroundColor(0xe6ecf5)
        view.zoomTo()
        view.show()