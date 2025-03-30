import dash
from dash import dcc, html, Input, Output, ctx
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import geopandas as gpd
import json
from plotly.subplots import make_subplots
import dash_bootstrap_components as dbc

# Load preprocessed data (you'll need to adjust paths)
df = pd.read_csv(
    "data/processed_philippine_cities_monthly.csv",
)
gdf = gpd.read_file("data/phl_adm_simple_maps.gpkg")
geojson_data = json.loads(gdf.to_json())

# identify quantitative columns for aggregation
quantitative_columns = [
    "temperature_2m_mean",
    "apparent_temperature_mean",
    "wind_speed_10m_max",
    "shortwave_radiation_sum",
    "HLI",
]

def process_spatial_aggregation(df, group_cols, gdf, crs, adm_col="adm1"):
    # Aggregate by specified grouping columns
    df_agg = df.groupby(group_cols).mean().reset_index()

    # Convert to GeoDataFrame and reproject
    df_agg["geometry"] = gpd.points_from_xy(df_agg["longitude"], df_agg["latitude"])
    gdf_agg = gpd.GeoDataFrame(df_agg, geometry="geometry", crs=crs).to_crs(crs)

    # Perform spatial join with administrative boundaries
    merged_gdf = gpd.sjoin(gdf_agg, gdf[["name", "geometry"]], how="left", predicate="intersects")
    merged_gdf = merged_gdf.rename(columns={"name": adm_col})

    df_final =  merged_gdf.groupby([adm_col] + [col for col in group_cols if col != "city_name"])[quantitative_columns].mean().reset_index()

    df_final["excl_CAR"] = np.where(df_final[adm_col] == "Cordillera Administrative Region", 1, 0)

    return df_final

# Process monthly decadal aggregation
gdf_month_decadal_adm1 = process_spatial_aggregation(df, ["decade", "month", "city_name"], gdf, gdf.crs)
gdf_month_decadal_adm1["decade"] = gdf_month_decadal_adm1["decade"].astype(int, errors="ignore")

# Process decadal aggregation
df["decade"] = df["year"] // 10 * 10  # Ensure decade column is created
gdf_decadal_adm1 = process_spatial_aggregation(df, ["decade", "city_name"], gdf, gdf.crs)
gdf_decadal_adm1["decade"] = gdf_decadal_adm1["decade"].astype(int, errors="ignore")

# create dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server  # Required for deployment with Gunicorn

app.layout = dbc.Container(
    [
        # Title
        dbc.Row(
            dbc.Col(html.H1("Philippine Climate Data Visualization - WIP", className="text-center mb-4"))
        ),
        # Main content (Map + Charts)
        dbc.Row(
            [
                # Left Column: Map
                dbc.Col(
                    [
                        html.Div(
                            dcc.Slider(
                                id="year-slider",
                                min=int(gdf_decadal_adm1["decade"].min()),
                                max=int(gdf_decadal_adm1["decade"].max()),
                                value=int(gdf_decadal_adm1["decade"].min()),
                                marks={str(int(decade)): str(int(decade)) for decade in gdf_decadal_adm1["decade"].unique()},
                                step=None,
                            ),
                            className="slider-container",  # Assign a CSS class
                        ),
                        html.Div(
                            [
                                html.Div(
                                    dcc.Graph(id="choropleth-map", responsive=True),
                                    className="map-frame",
                                ),
                                html.Div(
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                html.Button("Reset Selection", id="reset-button", n_clicks=0, className="btn btn-secondary"),
                                                width="auto",
                                            ),
                                            dbc.Col(
                                                dcc.Checklist(
                                                    id="exclude-car-toggle",
                                                    options=[{"label": " Exclude CAR", "value": "exclude"}],
                                                    value=[],
                                                    inline=True,
                                                    inputStyle={"margin-right": "5px"},
                                                ),
                                                width="auto",
                                                className="d-flex align-items-center",
                                            ),
                                        ],
                                        className="d-flex align-items-center",
                                    ),
                                    className="map-controls",  # <--- Add this class
                                ),
                            ],
                            className="map-container",
                        ),
                    ],
                    width=6,
                    xs=12,
                    sm=12,
                    md=6,
                    lg=6,
                    
                ),

                # Right Column: Charts
                dbc.Col(
                    [
                        html.Div(dcc.Graph(id="line-chart", responsive=True), className="chart-container"),
                        html.Div(dcc.Graph(id="bar-chart", responsive=True), className="chart-container"),
                        html.Div(dcc.Graph(id="line-chart-hli-monthly", responsive=True), className="chart-container"),
                    ],
                    width=6,
                    xs=12,
                    sm=12,
                    md=6,
                    lg=6,

                )
            ]
        ),
    ],
    fluid=True,  # Makes the container full-width
    className="pb-5 bg-light",  # Adds padding and background color
)

def get_selected_region(clickData):
    if not clickData or "points" not in clickData or len(clickData["points"]) == 0:
        return "All Regions"
    return clickData["points"][0]["location"]

# Helper function for layout
def apply_chart_layout(fig, title, x_label="Year", y_label="Value", df=None, x_col=None, y_col=None):
    fig.update_layout(
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        autosize=True,  # Ensures the figure resizes dynamically
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.5,
            xanchor="center",
            x=0.5,
        ),
        margin=dict(b=80),  # Adjust margin for legend placement
        xaxis=dict(
            showgrid=False,
            zeroline=False,
        ),
        yaxis=dict(
            showgrid=True,  # Keeps grid lines visible
            zeroline=True,  # Keeps the zero line visible
        ),
    )

    if df is not None and x_col is not None:
        if x_col in df.columns:
            fig.update_layout(
                xaxis=dict(
                    type="category" if df[x_col].dtype == 'O' else "linear",  # Categorical vs numeric handling
                    tickmode="array",
                    tickvals=df[x_col].unique().tolist(),  # Ensure all x-axis values are shown
                )
            )

    if df is not None and y_col is not None:
        if y_col in df.columns:
            fig.update_layout(
                yaxis=dict(
                    autorange=True,
                    tickmode="array",
                    tickvals=df[y_col].unique().tolist() if df[y_col].dtype == 'O' else None,  # Only force ticks for categorical values
                )
            )

# calculate the aggregate on all regions
all_regions_decadal_avg = gdf_decadal_adm1.groupby(["decade"])[quantitative_columns].mean().reset_index()
all_regions_monthly_avg = gdf_month_decadal_adm1.groupby(["month", "decade"])[quantitative_columns].mean().reset_index()

# Callbacks
@app.callback(
    [
        Output("choropleth-map", "figure"),
        Output("exclude-car-toggle", "value"),  # Ensure checkbox resets visually
    ],
    [
        Input("year-slider", "value"),
        Input("choropleth-map", "clickData"),
        Input("reset-button", "n_clicks"),  # Reset trigger
        Input("exclude-car-toggle", "value"),
    ],
)
def update_choropleth(selected_year, clickData, n_clicks, exclude_car):
    triggered_id = ctx.triggered_id

    # Reset checkbox when reset button is clicked
    if triggered_id == "reset-button":
        exclude_car = []

    global_hli_min = gdf_decadal_adm1["HLI"].min()
    global_hli_max = gdf_decadal_adm1["HLI"].max()

    # Filter data for the selected year
    filtered_df = gdf_decadal_adm1[gdf_decadal_adm1["decade"] == selected_year]

    # Apply exclusion filter if checkbox is checked
    if "exclude" in exclude_car:
        filtered_df = filtered_df[filtered_df["excl_CAR"] == 0]

    # Determine the selected region, reset if the reset button is clicked
    if triggered_id == "reset-button":
        selected_region = "All Regions"  # Force reset to ALL REGIONS
    else:
        selected_region = get_selected_region(clickData)  # Only run if reset wasn't clicked

    # Create Choropleth Map
    fig = px.choropleth(
        filtered_df,
        geojson=geojson_data,
        locations="adm1",
        featureidkey="properties.name",
        color="HLI",
        color_continuous_scale="thermal",
        range_color=[global_hli_min, global_hli_max],
    )

    # Reset all borders to white if reset button is clicked
    if selected_region == "All Regions":
        print("Resetting map colors to white.")
        fig.update_traces(
            marker_line_color="white"
        )
    else:
        print(f"Highlighting region: {selected_region}")
        border_colors = [
            "black" if region == selected_region else "white"
            for region in filtered_df["adm1"]
        ]
        fig.update_traces(
            marker_line_color=border_colors
        )

    # Remove background map and maximize plot area
    fig.update_geos(
        fitbounds="locations",
        visible=False,
        projection_type="mercator",
    )

    # Maximize the size of the map in the canvas
    fig.update_layout(
        uirevision=str(n_clicks),
        margin={"r": 0, "t": 0, "l": 0, "b": 0},
        dragmode=False,
    )

    return (fig, exclude_car)

@app.callback(
    [
        Output("line-chart", "figure"),
        Output("bar-chart", "figure"),
        Output("line-chart-hli-monthly", "figure"),
    ],
    [
        Input("year-slider", "value"),
        Input("choropleth-map", "clickData"),
        Input("reset-button", "n_clicks"),
        Input("exclude-car-toggle", "value"),  # Add checkbox input
    ],  
)
def update_charts(selected_year, clickData, _, exclude_car):
    triggered_id = ctx.triggered_id

    # Force reset if reset-button is clicked
    selected_region = "All Regions" if triggered_id == "reset-button" else get_selected_region(clickData)

    # Apply exclusion filter if checkbox is checked
    if "exclude" in exclude_car:
        filtered_gdf_decadal_adm1 = gdf_decadal_adm1[gdf_decadal_adm1["excl_CAR"] != 1]
        filtered_gdf_month_decadal_adm1 = gdf_month_decadal_adm1[gdf_month_decadal_adm1["excl_CAR"] != 1]
    else:
        filtered_gdf_decadal_adm1 = gdf_decadal_adm1
        filtered_gdf_month_decadal_adm1 = gdf_month_decadal_adm1

    # Use precomputed averages instead of redundant groupby calculations
    adm1_df = all_regions_decadal_avg if selected_region == "All Regions" else filtered_gdf_decadal_adm1[filtered_gdf_decadal_adm1["adm1"] == selected_region].copy()
    adm1_month_df = all_regions_monthly_avg if selected_region == "All Regions" else filtered_gdf_month_decadal_adm1[filtered_gdf_month_decadal_adm1["adm1"] == selected_region].copy()

    adm1_df["HLI_5yr_MA"] = adm1_df["HLI"].rolling(window=5, min_periods=1).mean()
    adm1_df["HLI_10yr_MA"] = adm1_df["HLI"].rolling(window=10, min_periods=1).mean()
    adm1_df["HLI_20yr_MA"] = adm1_df["HLI"].rolling(window=20, min_periods=1).mean()

    # LINE CHART: HLI Trends + Moving Averages
    fig_line = go.Figure()

    # Original HLI Trend
    fig_line.add_trace(
        go.Scatter(
            x=adm1_df["decade"],
            y=adm1_df["HLI"],
            mode="lines+markers",
            name="HLI",
            line=dict(width=2),
        )
    )

    # 5-Year Moving Average
    fig_line.add_trace(
        go.Scatter(
            x=adm1_df["decade"],
            y=adm1_df["HLI_5yr_MA"],
            mode="lines",
            name="5-Year MA",
            line=dict(dash="dash", width=2),
        )
    )

    # 10-Year Moving Average
    fig_line.add_trace(
        go.Scatter(
            x=adm1_df["decade"],
            y=adm1_df["HLI_10yr_MA"],
            mode="lines",
            name="10-Year MA",
            line=dict(dash="dash", width=2),
        )
    )

    # 20-Year Moving Average
    fig_line.add_trace(
        go.Scatter(
            x=adm1_df["decade"],
            y=adm1_df["HLI_20yr_MA"],
            mode="lines",
            name="20-Year MA",
            line=dict(dash="dash", width=2),
        )
    )

    # Add Vertical Line at Selected Year
    fig_line.add_shape(
        dict(
            type="line",
            x0=selected_year,  # Position on x-axis
            x1=selected_year,  # Same x position to form a vertical line
            y0=adm1_df["HLI"].min(),  # Start from the lowest value in HLI
            y1=adm1_df["HLI"].max(),  # End at the highest value in HLI
            line=dict(width=0.5, dash="dot"),  # Customize color and style
        )
    )
    
    # Layout settings
    apply_chart_layout(fig_line, f"HLI Trends - {selected_region}", "Year", "Heat Load Index (HLI)", x_col="decade", df=adm1_df)

    # BAR CHART: Temperature and Wind Speed for Selected Year & Region
    # Create figure with secondary y-axis
    fig_bar = make_subplots(specs=[[{"secondary_y": True}]])

    # Add Bar Chart (Primary Y-Axis: Temperature)
    fig_bar.add_trace(
        go.Bar(
            x=adm1_df["decade"],
            y=adm1_df["temperature_2m_mean"],
            name="Mean Temp",
            marker=dict(color="#A7C7E7"),
        ),
        secondary_y=False,
    )

    fig_bar.add_trace(
        go.Bar(
            x=adm1_df["decade"],
            y=adm1_df["apparent_temperature_mean"],
            name="Apparent Temp",
            marker=dict(color="#D1B3E0"),
        ),
        secondary_y=False,
    )

    # Add Scatter Line (Secondary Y-Axis: Wind Speed)
    fig_bar.add_trace(
        go.Scatter(
            x=adm1_df["decade"],
            y=adm1_df["wind_speed_10m_max"],
            mode="markers+lines",
            name="Wind Speed",
            line=dict(dash="dot", width=2, color="#5E548E"),
        ),
        secondary_y=True,  # Assign this trace to the secondary y-axis
    )

    # Add Scatter Line (Secondary Y-Axis: Shortwave Radiation)
    fig_bar.add_trace(
        go.Scatter(
            x=adm1_df["decade"],
            y=adm1_df["shortwave_radiation_sum"],
            mode="markers+lines",
            name="Shortwave Radiation",
            line=dict(dash="dot", width=2,  color="#AC6C57"),
        ),
        secondary_y=True,  # Assign this trace to the secondary y-axis
    )

    # Add Vertical Line at Selected Year in Bar Chart
    fig_bar.add_shape(
        dict(
            type="rect",
            x0=selected_year - 5,  # Adjusting to center around the decade
            x1=selected_year + 5,  # Ensuring the highlight covers the decade properly
            y0=0,
            y1=max(adm1_df["temperature_2m_mean"].max(), adm1_df["apparent_temperature_mean"].max()),  
            fillcolor="rgba(200, 200, 200, 1)",  # Light gray with 30% opacity
            layer="below",  # Ensure it is behind all other elements
            line=dict(width=0),  # No border
        )
    )

    # Remove Y-axis grid lines
    fig_bar.update_layout(
        yaxis=dict(
            showgrid=False,  # Removes primary Y-axis grid
            zeroline=False,  # Removes zero line
        ),
        yaxis2=dict(
            showgrid=False,  # Removes secondary Y-axis grid
            zeroline=False,
        ),
    )

    # Update layout for dual Y-Axis
    apply_chart_layout(fig_bar, f"Temperature & Wind Speed - {selected_region}", "Decade", "Temperature (°C)", x_col="decade", df=adm1_df)

    # LINE CHART: HLI by Month and Decade
    fig_monthly_hli = px.line(
        adm1_month_df,
        x="month",
        y="HLI",
        color="decade",
        labels={"month": "Month", "HLI": "Heat Load Index (HLI)", "decade": "Decade"},
        markers=True,  # Adds markers
    )
    # Set opacity for all lines to 50% by default
    for trace in fig_monthly_hli.data:
        trace.opacity = 0.2  # Default opacity

    # Increase opacity to 100% if the decade matches the selected year
    for trace in fig_monthly_hli.data:
        if str(trace.name) == str(selected_year):  # Match decade with slider
            trace.opacity = 1.0  # Full opacity for matching line

    # Add a title
    apply_chart_layout(fig_monthly_hli, f"Monthly HLI Trends by Decade - {selected_region}", "Month", "HLI", x_col="month", df=adm1_month_df)

    return fig_line, fig_bar, fig_monthly_hli


# Run app
if __name__ == "__main__":
    # app.run_server(debug=True)
    app.run(debug=True)