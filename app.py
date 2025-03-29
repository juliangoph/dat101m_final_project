import dash
from dash import dcc, html, Input, Output, ctx
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import geopandas as gpd
import json
from plotly.subplots import make_subplots
import dash_bootstrap_components as dbc

# Load preprocessed data (you'll need to adjust paths)
df = pd.read_csv(
    "data/processed_philippine_cities_monthly.csv",
)
gdf = gpd.read_file("data/phl_adm_simple_maps.gpkg")

# identify quantitative columns for aggregation
quantitative_columns = [
    "temperature_2m_mean",
    "apparent_temperature_mean",
    "wind_speed_10m_max",
    "shortwave_radiation_sum",
    "HLI",
]

## aggregate by df_month_decadal and convert to geospatial dataset
# aggregate data by month-decade
df_month_decadal = (
    df.groupby(["decade", "month", "city_name"]).mean().reset_index()
)

# convert GeoDataFrame to JSON for choropleth mapping
geojson_data = json.loads(gdf.to_json())

# convert pandas to a gdf
df_month_decadal["geometry"] = gpd.points_from_xy(
    df_month_decadal["longitude"], df_month_decadal["latitude"]
)

# reproject "gdf_month_decadal" to match crs on "gdf"
gdf_month_decadal = gpd.GeoDataFrame(df_month_decadal, geometry="geometry", crs=gdf.crs)
gdf_month_decadal = gdf_month_decadal.to_crs(gdf.crs)

# perform spatial join
merged_gdf_month_decadal = gpd.sjoin(
    gdf_month_decadal, gdf[["name", "geometry"]], how="left", predicate="intersects"
)
merged_gdf_month_decadal = merged_gdf_month_decadal.rename(
    columns={"name": "adm1"}
)  # rename col to "adm1"

# perform aggregation by decade and month
gdf_month_decadal_adm1 = (
    merged_gdf_month_decadal.groupby(["adm1", "decade", "month"])[quantitative_columns]
    .mean()
    .reset_index()
)

gdf_month_decadal_adm1["decade"] = gdf_month_decadal_adm1["decade"].astype(
    int, errors="ignore"
)

## aggregate by df_decadal and convert to geospatial dataset
# aggregate data by year and then by decade
df_yearly = df.groupby(["year", "city_name"]).mean().reset_index()
df_decadal = df_yearly.groupby(["decade", "city_name"]).mean().reset_index()

# convert c to a GeoDataFrame
df_decadal["geometry"] = gpd.points_from_xy(
    df_decadal["longitude"], df_decadal["latitude"]
)

# reproject "gdf_decadal" to match crs on "gdf"
gdf_decadal = gpd.GeoDataFrame(df_decadal, geometry="geometry", crs=gdf.crs)
gdf_decadal = gdf_decadal.to_crs(gdf.crs)

# perform spatial join
merged_gdf = gpd.sjoin(
    gdf_decadal, gdf[["name", "geometry"]], how="left", predicate="intersects"
)
merged_gdf = merged_gdf.rename(columns={"name": "adm1"})  # rename col to "adm1"

# perform aggregation by decade and month
gdf_decadal_adm1 = (
    merged_gdf.groupby(["adm1", "decade"])[quantitative_columns].mean().reset_index()
)

gdf_decadal_adm1["decade"] = gdf_decadal_adm1["decade"].astype(int, errors="ignore")

# create dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server  # Required for deployment with Gunicorn

app.layout = dbc.Container(
    [
        # Title
        dbc.Row(
            dbc.Col(html.H1("Philippine Climate Data Visualization", className="text-center mb-4"))
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
                                    dcc.Graph(id="choropleth-map", style={"height": "87vh"}),
                                    className="map-frame",
                                ),
                                html.Button("Reset Selection", id="reset-button", n_clicks=0, className="btn btn-secondary"),
                            ],
                            className="map-container",
                        ),
                    ],
                    width=6,  # Takes 6/12 of the screen (half)
                ),

                # Right Column: Charts
                dbc.Col(
                    [
                        html.Div(dcc.Graph(id="line-chart", style={"height": "30vh"}), style={"padding-bottom": "15px"}),
                        html.Div(dcc.Graph(id="bar-chart", style={"height": "30vh"}), style={"padding-bottom": "15px"}),
                        html.Div(dcc.Graph(id="line-chart-hli-monthly", style={"height": "30vh"}), style={"padding-bottom": "15px"}),
                    ],
                    width=6,  # Takes 6/12 of the screen (half)
                )
            ]
        ),
    ],
    fluid=True,  # Makes the container full-width
    className="p-4 bg-light",  # Adds padding and background color
)

def get_selected_region(clickData):
    if not clickData or "points" not in clickData or len(clickData["points"]) == 0:
        return "All Regions"
    return clickData["points"][0]["location"]

# Helper function for layout
def apply_chart_layout(fig, title, x_label="Year", y_label="Value"):
    fig.update_layout(
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        legend=dict(
            orientation="h",  # Horizontal legend
            yanchor="top",
            y=-0.5,  # Moves the legend below the x-axis label
            xanchor="center",
            x=0.5,
        ),
        margin=dict(b=80)  # Increases bottom margin to prevent overlap
    )

# calculate the aggregate on all regions
all_regions_decadal_avg = gdf_decadal_adm1.groupby(["decade"])[quantitative_columns].mean().reset_index()
all_regions_monthly_avg = gdf_month_decadal_adm1.groupby(["month", "decade"])[quantitative_columns].mean().reset_index()

# Callbacks
@app.callback(
    Output("choropleth-map", "figure"),
    Input("year-slider", "value"),
    Input("choropleth-map", "clickData")  # Capture clicked region
)
def update_choropleth(selected_year, clickData):
    global_hli_min = gdf_decadal_adm1["HLI"].min()
    global_hli_max = gdf_decadal_adm1["HLI"].max()

    # Filter data for the selected year
    filtered_df = gdf_decadal_adm1[gdf_decadal_adm1["decade"] == selected_year]

    # Get the clicked region (default to None)
    selected_region = get_selected_region(clickData)

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

    # Update Borders: Highlight Selected Region
    fig.update_traces(
        marker_line_color=[
            "black" if region == selected_region else "white"
            for region in filtered_df["adm1"]
        ],
    )

    # Remove background map and maximize plot area
    fig.update_geos(
        fitbounds="locations",
        visible=False,
        projection_type="mercator",
    )

    # Maximize the size of the map in the canvas
    fig.update_layout(
        margin={"r": 0, "t": 0, "l": 0, "b": 0},
        dragmode=False,
    )

    return fig

@app.callback(
    [
        Output("line-chart", "figure"),
        Output("bar-chart", "figure"),
        Output("line-chart-hli-monthly", "figure"),
    ],
    [Input("choropleth-map", "clickData"),
     Input("reset-button", "n_clicks")],  
)
def update_charts(clickData, _):
    triggered_id = ctx.triggered_id

    # Force reset if reset-button is clicked
    selected_region = "All Regions" if triggered_id == "reset-button" else get_selected_region(clickData)

    # Use precomputed averages instead of redundant groupby calculations
    adm1_df = all_regions_decadal_avg if selected_region == "All Regions" else gdf_decadal_adm1[gdf_decadal_adm1["adm1"] == selected_region]
    adm1_month_df = all_regions_monthly_avg if selected_region == "All Regions" else gdf_month_decadal_adm1[gdf_month_decadal_adm1["adm1"] == selected_region]

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

    # Layout settings
    apply_chart_layout(fig_line, f"HLI Trends - {selected_region}", "Year", "Heat Load Index (HLI)")

    # BAR CHART: Temperature and Wind Speed for Selected Year & Region
    # Create figure with secondary y-axis
    fig_bar = make_subplots(specs=[[{"secondary_y": True}]])

    # Add Bar Chart (Primary Y-Axis: Temperature)
    fig_bar.add_trace(
        go.Bar(
            x=adm1_df["decade"],
            y=adm1_df["temperature_2m_mean"],
            name="Mean Temp",
        ),
        secondary_y=False,
    )

    fig_bar.add_trace(
        go.Bar(
            x=adm1_df["decade"],
            y=adm1_df["apparent_temperature_mean"],
            name="Apparent Temp",
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
            line=dict(dash="dot", width=2),
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
            line=dict(dash="dot", width=2),
        ),
        secondary_y=True,  # Assign this trace to the secondary y-axis
    )

    # Update layout for dual Y-Axis
    apply_chart_layout(fig_bar, f"Temperature & Wind Speed - {selected_region}", "Decade", "Temperature (°C)")

    # LINE CHART: HLI by Month and Decade
    fig_monthly_hli = px.line(
        adm1_month_df,
        x="month",
        y="HLI",
        color="decade",
        labels={"month": "Month", "HLI": "Heat Load Index (HLI)", "decade": "Decade"},
        markers=True,  # Adds markers
    )

    # Add a title
    apply_chart_layout(fig_monthly_hli, f"Monthly HLI Trends by Decade - {selected_region}", "Month", "HLI")

    return fig_line, fig_bar, fig_monthly_hli


# Run app
if __name__ == "__main__":
    # app.run_server(debug=True)
    app.run(debug=True)