import geopandas as gpd
import pandas as pd
import folium
import branca.colormap as cm
from pathlib import Path
import os

# ================= 0. Configuration & Paths =================
bound_path = Path("CAUTH_MAY_2025_EN_BSC_1271049543973882786.geojson")
brownfield_path = Path("brownfield-land.geojson")
lsoa_path = Path("Lower_layer_Super_Output_Areas_December_2021_Boundaries_EW_BSC_V4_-4299016806856585929.geojson")
fuel_excel_path = Path("Sub-regional_fuel_poverty_statistics_2023.xlsx")
imd_excel_path = Path("Deprivation_2025.xlsx")
police_dir = Path("police_data")

# ================= 1. Process Administrative Boundaries (LCR) =================
print("1/8 Reading administrative boundary data...")
bound_gdf = gpd.read_file(bound_path)
if bound_gdf.crs is None or bound_gdf.crs.to_epsg() != 4326:
    bound_gdf = bound_gdf.to_crs(epsg=4326)

lcr = bound_gdf[bound_gdf["CAUTH25NM"] == "Liverpool City Region"].copy()
if lcr.empty:
    raise ValueError("LCR boundary not found.")

centre = lcr.geometry.unary_union.centroid
centre_lat, centre_lon = centre.y, centre.x

# ================= 2. Process LSOA Spatial Data =================
print("2/8 Reading and filtering LSOA spatial data...")
lsoa_gdf = gpd.read_file(lsoa_path)
if lsoa_gdf.crs is None or lsoa_gdf.crs.to_epsg() != 4326:
    lsoa_gdf = lsoa_gdf.to_crs(epsg=4326)

lsoa_lcr = gpd.sjoin(
    lsoa_gdf,
    lcr[["geometry"]],
    how="inner",
    predicate="intersects"
).drop(columns=['index_right'])

valid_lsoa_codes = set(lsoa_lcr['LSOA21CD'].unique())

# ================= 3. Load Datasets =================
print("3/8 Processing Fuel Poverty Data...")
try:
    fuel_data = pd.read_excel(fuel_excel_path, sheet_name="Table 4", header=2)
    fuel_data = fuel_data.rename(columns={"LSOA Code": "LSOA21CD", "Proportion of households fuel poor (%)": "Fuel_Poverty_Rate"})
    fuel_data["Fuel_Poverty_Rate"] = pd.to_numeric(fuel_data["Fuel_Poverty_Rate"], errors="coerce")
    fuel_data = fuel_data[["LSOA21CD", "Fuel_Poverty_Rate"]]
except:
    fuel_data = pd.DataFrame(columns=["LSOA21CD", "Fuel_Poverty_Rate"])

print("4/8 Processing IMD Deprivation Data...")
try:
    imd_data = pd.read_excel(imd_excel_path, sheet_name="IMD25", header=0)
    imd_col = next((c for c in imd_data.columns if "Decile" in str(c) and "IMD" in str(c)), None)
    if imd_col:
        imd_data = imd_data.rename(columns={"LSOA code (2021)": "LSOA21CD", imd_col: "IMD_Decile"})
        imd_data["IMD_Decile"] = pd.to_numeric(imd_data["IMD_Decile"], errors="coerce")
        imd_data = imd_data[["LSOA21CD", "IMD_Decile"]]
    else:
        imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])
except:
    imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])

print("5/8 Processing Police Data...")
street_counts = {}
stop_search_points = []
all_csvs = list(police_dir.rglob("*.csv"))

for csv_file in all_csvs:
    fname = csv_file.name.lower()
    try:
        if "street" in fname:
            df = pd.read_csv(csv_file, usecols=["LSOA code"])
            df_filtered = df[df["LSOA code"].isin(valid_lsoa_codes)]
            counts = df_filtered["LSOA code"].value_counts()
            for code, count in counts.items():
                street_counts[code] = street_counts.get(code, 0) + count
        elif "stop-and-search" in fname:
            df = pd.read_csv(csv_file, usecols=["Latitude", "Longitude"])
            df = df.dropna(subset=["Latitude", "Longitude"])
            # Rough filter for speed
            mask = (df['Latitude'].between(53.1, 53.8)) & (df['Longitude'].between(-3.4, -2.4))
            df_trimmed = df[mask]
            if not df_trimmed.empty:
                stop_search_points.append(df_trimmed)
    except:
        continue

df_street = pd.DataFrame(list(street_counts.items()), columns=["LSOA21CD", "Crime_Count"])

df_stop_agg = pd.DataFrame(columns=["LSOA21CD", "StopSearch_Count"])
if stop_search_points:
    all_stops = pd.concat(stop_search_points, ignore_index=True)
    gdf_stops = gpd.GeoDataFrame(all_stops, geometry=gpd.points_from_xy(all_stops.Longitude, all_stops.Latitude), crs="EPSG:4326")
    stops_with_lsoa = gpd.sjoin(gdf_stops, lsoa_lcr[['LSOA21CD', 'geometry']], predicate='within')
    stop_counts = stops_with_lsoa['LSOA21CD'].value_counts().reset_index()
    stop_counts.columns = ['LSOA21CD', 'StopSearch_Count']
    df_stop_agg = stop_counts

# ================= 6. Merge & Normalize Data =================
print("6/8 Merging and Normalizing datasets...")

master_gdf = lsoa_lcr.merge(fuel_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(imd_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_street, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_stop_agg, on="LSOA21CD", how="left")

# Fill NaNs
master_gdf["Fuel_Poverty_Rate"] = master_gdf["Fuel_Poverty_Rate"].fillna(0)
master_gdf["IMD_Decile"] = master_gdf["IMD_Decile"].fillna(10) # 10 is least deprived (good)
master_gdf["Crime_Count"] = master_gdf["Crime_Count"].fillna(0)
master_gdf["StopSearch_Count"] = master_gdf["StopSearch_Count"].fillna(0)

# --- Normalization (0.0 = Good/Low Intensity, 1.0 = Bad/High Intensity) ---
# This ensures all layers share the same "Red" scale meaning

# 1. Fuel: Simply min-max scale
f_min, f_max = master_gdf["Fuel_Poverty_Rate"].min(), master_gdf["Fuel_Poverty_Rate"].max()
master_gdf["norm_fuel"] = (master_gdf["Fuel_Poverty_Rate"] - f_min) / (f_max - f_min + 1e-9)

# 2. IMD: Invert (1 is Bad, 10 is Good) -> (11 - Decile) / 10 -> 1 becomes 1.0 (Bad)
master_gdf["norm_imd"] = (11 - master_gdf["IMD_Decile"]) / 10.0

# 3. Crime: Scale to 95th percentile (to handle extreme hotspots)
c_max = master_gdf["Crime_Count"].quantile(0.95)
master_gdf["norm_crime"] = master_gdf["Crime_Count"] / c_max
master_gdf["norm_crime"] = master_gdf["norm_crime"].clip(upper=1.0) # Cap at 1.0

# 4. StopSearch: Scale to 95th percentile
s_max = master_gdf["StopSearch_Count"].quantile(0.95)
master_gdf["norm_stop"] = master_gdf["StopSearch_Count"] / s_max
master_gdf["norm_stop"] = master_gdf["norm_stop"].clip(upper=1.0)

# 5. Composite Score (Average of the available normalized metrics)
# We use this to color the Brownfield sites to reflect their area
master_gdf["Composite_Score"] = master_gdf[["norm_fuel", "norm_imd", "norm_crime"]].mean(axis=1)

# ================= 7. Prepare Brownfield Data (Spatial Join for Color) =================
print("7/8 Processing Brownfield Sites...")
brownfield = gpd.read_file(brownfield_path).to_crs(epsg=4326)

# Join Brownfield points to Master GDF to get the Composite Score of the area
bf_joined = gpd.sjoin(brownfield, master_gdf[["geometry", "Composite_Score", "LSOA21NM"]], predicate="within")
bf_joined["hectares"] = pd.to_numeric(bf_joined["hectares"], errors="coerce")

# ================= 8. Mapping =================
print("8/8 Generating Map...")

m = folium.Map(location=[centre_lat, centre_lon], zoom_start=11, tiles="CartoDB positron")

# --- UNIFIED COLORMAP ---
# Transparent White -> Red. 
# When layers overlap, the Red opacity accumulates, creating a "blended" darker red.
# We use 'Reds' from 0.0 to 1.0.
cmap_unified = cm.LinearColormap(['#ffffff', '#fee5d9', '#fcae91', '#fb6a4a', '#de2d26', '#a50f15'], 
                                 vmin=0.0, vmax=1.0, 
                                 caption="Vulnerability / Intensity Score (0-1)")
m.add_child(cmap_unified)

# --- Custom Panes for Z-Index Control ---
# Default polygon pane z-index is usually around 400.
# We place markers on a custom pane with z-index 600 to ensure they are ALWAYS on top.
folium.map.CustomPane("sites_pane", z_index=600).add_to(m)

# --- Helper Function for Polygon Layers ---
def add_layer(gdf, name, norm_column, real_column, show=False):
    fg = folium.FeatureGroup(name=name, show=show)
    folium.GeoJson(
        gdf,
        style_function=lambda x: {
            # Use the normalized column (0-1) for color
            "fillColor": cmap_unified(x["properties"][norm_column]),
            "color": "transparent",
            "weight": 0,
            # Lower opacity allows blending when multiple layers are on
            "fillOpacity": 0.5, 
        },
        highlight_function=lambda x: {"weight": 1, "color": "#333", "fillOpacity": 0.8},
        tooltip=folium.GeoJsonTooltip(
            fields=["LSOA21NM", real_column],
            aliases=["Area:", f"{name} Value:"],
            sticky=False,
            # Format numbers nicely
            style="font-family: sans-serif;" 
        )
    ).add_to(fg)
    fg.add_to(m)

# Add Layers (All use the same Unified Colormap)
add_layer(master_gdf, "Fuel Poverty", "norm_fuel", "Fuel_Poverty_Rate", show=True)
add_layer(master_gdf, "IMD Deprivation", "norm_imd", "IMD_Decile", show=False)
add_layer(master_gdf, "Crime Density", "norm_crime", "Crime_Count", show=False)
add_layer(master_gdf, "Stop & Search", "norm_stop", "StopSearch_Count", show=False)

# --- Brownfield Sites (On Top Pane) ---
fg_sites = folium.FeatureGroup(name="Brownfield Sites (Composite Color)", show=True)

for _, row in bf_joined.iterrows():
    if row.geometry.is_empty or pd.isna(row["hectares"]): continue
    
    # Get the score (default to 0 if missing)
    score = row["Composite_Score"] if not pd.isna(row["Composite_Score"]) else 0
    
    folium.CircleMarker(
        location=[row.geometry.centroid.y, row.geometry.centroid.x],
        radius=3 + (row["hectares"] * 0.4),
        color="#333",        # Dark grey border
        weight=1,
        fillColor=cmap_unified(score), # Fill with the area's composite color
        fillOpacity=1.0,     # Solid opacity to stand out against background
        popup=f"<b>Site:</b> {row.get('SiteName', 'N/A')}<br>"
              f"<b>Size:</b> {row['hectares']:.2f} ha<br>"
              f"<b>LSOA:</b> {row['LSOA21NM']}<br>"
              f"<b>Area Score:</b> {score:.2f} / 1.0",
    ).add_to(fg_sites)

# IMPORTANT: Add the FeatureGroup to the map, but we don't attach the pane to the Group itself 
# in the standard way for Markers easily in pure Python Folium without loop hacks.
# Wait, Folium CircleMarker accepts a 'pane' argument? Yes!
# Let's rebuild the loop to use the pane correctly.

fg_sites = folium.FeatureGroup(name="Brownfield Sites", show=True)
for _, row in bf_joined.iterrows():
    if row.geometry.is_empty or pd.isna(row["hectares"]): continue
    score = row["Composite_Score"] if not pd.isna(row["Composite_Score"]) else 0
    
    folium.CircleMarker(
        location=[row.geometry.centroid.y, row.geometry.centroid.x],
        radius=4 + (row["hectares"] * 0.4),
        color="black",
        weight=1,
        fillColor=cmap_unified(score),
        fillOpacity=0.9,
        popup=f"Site: {row.get('SiteName', 'N/A')}",
        pane="sites_pane"  # <--- MAGIC: Assign to the high Z-Index pane
    ).add_to(fg_sites)

fg_sites.add_to(m)

# --- Finalize ---
folium.LayerControl(collapsed=False).add_to(m)
outfile = "LCR_Integrated_Map_Unified.html"
m.save(outfile)
print(f"Done! Map saved to {outfile}")