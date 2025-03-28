import pandas as pd

# Load preprocessed data (you'll need to adjust paths)
df = pd.read_csv(
    "data/philippine_cities_daily_data_1950_to_2025.csv",
    parse_dates=["date"],
)


# add heat load index
df["HLI"] = (df["temperature_2m_mean"] + (0.1 * df["shortwave_radiation_sum"])) - (
    0.05 * df["wind_speed_10m_max"]
)

# Extract month, year and decade from the date column
df["month"] = df["date"].dt.month
df["year"] = df["date"].dt.year
df["decade"] = (df["date"].dt.year // 10) * 10

# Aggregate data monthly
df_monthly = (
    df.groupby(
        [
            "decade",
            "year",
            "month",
            "city_name",
            "latitude",
            "longitude",
        ]
    )
    .agg(
        {
            "temperature_2m_mean": "mean",
            "apparent_temperature_mean": "mean",
            "wind_speed_10m_max": "mean",
            "shortwave_radiation_sum": "mean",
            "HLI": "mean",
        }
    )
    .reset_index()
)

df_monthly.to_csv("data/processed_philippine_cities_monthly.csv", index=False)
