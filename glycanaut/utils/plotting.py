import pandas as pd
import networkx as nx

import plotly.graph_objects as go
from plotly.graph_objs import Figure

from typing import Tuple, List


def plot_mass_spectrum(df: pd.DataFrame) -> Figure:
    """
    Plot a mass spectrum from a DataFrame containing "m/z" and "Intensity" columns.
    """
    fig = go.Figure()
    if df.empty:
        return fig

    for x, y in zip(df["m/z"], df["Intensity"]):
        fig.add_trace(
            go.Scatter(
                x=[x, x],
                y=[0, y],
                mode="lines+text",
                line=dict(color="grey", width=2),
                text=[None, f"{x:.2f}"],
                textposition="top center",
                showlegend=False,
            )
        )

    fig.update_layout(
        xaxis_title="m/z",
        yaxis_title="Intensity",
        template="simple_white",
        height=300,
        margin=dict(l=40, r=40, t=20, b=40),
    )

    fig.update_traces(hoverinfo="skip")

    return fig


def plot_peak_diff_histogram(
    df_diffs: pd.DataFrame, df_diffs_assigned: pd.DataFrame
) -> Figure:
    """
    Plot histogram of peak differences and annotate wherever assigned.
    """
    df_diffs["Peak Difference (rounded)"] = df_diffs["Peak Difference"].round(1)
    counts = (
        df_diffs["Peak Difference (rounded)"].value_counts().sort_index().reset_index()
    )
    counts.columns = ["Peak Difference", "Count"]

    assigned_info = (
        df_diffs.groupby("Peak Difference (rounded)")
        .agg(Assigned=("Assigned Symbol", "first"))
        .reset_index()
    )
    assigned_info["Peak Difference"] = assigned_info["Peak Difference (rounded)"]
    counts = counts.merge(assigned_info, on="Peak Difference", how="left")

    assigned = counts[
        counts["Peak Difference"].isin(df_diffs_assigned["Peak Difference"].round(1))
    ]
    unassigned = counts[
        ~counts["Peak Difference"].isin(df_diffs_assigned["Peak Difference"].round(1))
    ]

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=assigned["Peak Difference"],
            y=assigned["Count"],
            name="Assigned",
            marker_color="#669673",
            text=assigned["Assigned"],
            hovertemplate="Peak Difference: %{x}<br>"
            + "Count: %{y}<br>"
            + "Assigned: %{text}<extra></extra>",
        )
    )

    fig.add_trace(
        go.Bar(
            x=unassigned["Peak Difference"],
            y=unassigned["Count"],
            name="Not assigned",
            marker_color="#d3e3d7",
            hovertemplate="Peak Difference: %{x}<br>"
            + "Count: %{y}<br>"
            + "Assigned: None<extra></extra>",
        )
    )

    fig.update_layout(
        xaxis_title="Peak Differences",
        yaxis_title="Count",
        margin=dict(l=40, r=40, t=20, b=40),
        showlegend=True,
    )

    return fig


def plot_peak_diff_graph(df_assigned: pd.DataFrame) -> Figure:
    """
    Plot a network graph of assigned peak differences.
    """

    if df_assigned.empty:
        fig = go.Figure()
        fig.update_layout(
            showlegend=False,
            xaxis_visible=False,
            yaxis_visible=False,
            margin=dict(l=20, r=20, t=20, b=20),
        )
        return fig

    G = nx.from_pandas_edgelist(
        df_assigned[df_assigned["Type"] != "Modification"],
        "Peak 1",
        "Peak 2",
        edge_attr=["Peak Difference", "Assigned Symbol"],
    )
    pos = nx.spring_layout(G)

    smallest_peak = min(G.nodes)
    largest_peak = max(G.nodes)
    backbone = nx.shortest_path(G, source=largest_peak, target=smallest_peak)
    backbone_edges = list(zip(backbone[:-1], backbone[1:]))

    edge_x, edge_y, edge_colors, opacities = [], [], [], []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
        edge_colors.append(
            "#669673"
            if (u, v) in backbone_edges or (v, u) in backbone_edges
            else "grey"
        )
        opacities.append(
            0.8 if (u, v) in backbone_edges or (v, u) in backbone_edges else 0.2
        )

    edge_traces = []
    for i in range(len(G.edges())):
        edge_traces.append(
            go.Scatter(
                x=edge_x[i * 3 : i * 3 + 2],
                y=edge_y[i * 3 : i * 3 + 2],
                line=dict(width=1, color=edge_colors[i]),
                hoverinfo="none",
                opacity=opacities[i],
                mode="lines",
            )
        )

    node_x, node_y, node_text = zip(
        *[(pos[n][0], pos[n][1], str(n)) for n in G.nodes()]
    )
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        marker=dict(size=3, color="#669673"),
        text=node_text,
        textposition="top center",
    )

    fig = go.Figure(data=edge_traces + [node_trace])
    for u, v, attr in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        label = attr.get("Assigned Symbol", "")
        fig.add_annotation(
            x=(x0 + x1) / 2,
            y=(y0 + y1) / 2,
            text=label,
            showarrow=False,
            font=dict(size=30, color="grey"),
        )
    fig.update_layout(
        showlegend=False,
        xaxis_visible=False,
        yaxis_visible=False,
        margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig


def generate_glycan_graph(G: nx.Graph) -> Tuple[List, List]:

    # Example node / edge list:
    # nodes = [
    #     {"x": 0, "y": 0, "label": "GlcNAc", "color": "orange"},
    #     {"x": 1, "y": 0, "label": "Man", "color": "green"},
    #     {"x": 2, "y": 0, "label": "Man", "color": "green"},
    # ]
    # edges = [(0, 1), (1, 2)]

    nodes = []
    edges = []
    smallest_peak = min(G.nodes)
    largest_peak = max(G.nodes)
    backbone = nx.shortest_path(G, source=largest_peak, target=smallest_peak)
    branches = {}
    for i in range(len(backbone) - 1):
        node = backbone[i]
        nodes.append(
            {
                "x": i,
                "y": 0,
                "label": G.edges[backbone[i], backbone[i + 1]]["Assigned Symbol"],
            }
        )
        deg = G.degree(node)
        if deg > 2:
            neighbors = list(G.neighbors(node))
            for n in neighbors:
                if n not in backbone:
                    for path in nx.all_simple_paths(G, source=node, target=n):
                        for j in range(len(path) - 1):
                            nodes.append(
                                {
                                    "x": i,
                                    "y": j + 1,
                                    "label": G.edges[path[j], backbone[j + 1]][
                                        "Assigned Symbol"
                                    ],
                                }
                            )

    return nodes, edges


def plot_glycan(nodes: List, edges: List) -> Figure:
    fig = go.Figure()

    for i, j in edges:
        fig.add_shape(
            type="line",
            x0=nodes[i]["x"],
            y0=nodes[i]["y"],
            x1=nodes[j]["x"],
            y1=nodes[j]["y"],
            line=dict(color="black", width=2),
            layer="below",
        )

    fig.add_trace(
        go.Scatter(
            x=[n["x"] for n in nodes],
            y=[n["y"] for n in nodes],
            mode="markers+text",
            marker=dict(
                size=40,
                color=[n["color"] for n in nodes],
                line=dict(color="black", width=1),
            ),
            text=[n["label"] for n in nodes],
            textposition="bottom center",
            hovertext=[n["label"] for n in nodes],
            hoverinfo="text",
        )
    )

    fig.update_layout(
        showlegend=False,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        width=600,
        height=400,
    )

    return fig