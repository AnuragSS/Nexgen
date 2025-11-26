import geopandas as gpd
import pandas as pd
import folium
import branca.colormap as cm
from pathlib import Path

# ================= 0. Configuration Paths =================
bound_path = Path("CAUTH_MAY_2025_EN_BSC_1271049543973882786.geojson")
brownfield_path = Path("brownfield-land.geojson")
lsoa_path = Path("Lower_layer_Super_Output_Areas_December_2021_Boundaries_EW_BSC_V4_-4299016806856585929.geojson")
excel_path = Path("Sub-regional_fuel_poverty_statistics_2023.xlsx")

# ================= 1. Process Administrative Boundaries (LCR) =================
print("1/5 Reading administrative boundary data...")
bound_gdf = gpd.read_file(bound_path)
if bound_gdf.crs is None or bound_gdf.crs.to_epsg() != 4326:
    bound_gdf = bound_gdf.to_crs(epsg=4326)

# Filter for Liverpool City Region
lcr = bound_gdf[bound_gdf["CAUTH25NM"] == "Liverpool City Region"].copy()
if lcr.empty:
    raise ValueError("LCR boundary not found, please check the data.")

centre = lcr.geometry.unary_union.centroid
centre_lat, centre_lon = centre.y, centre.x

# ================= 2. Process LSOA Spatial Data =================
print("2/5 Reading and filtering LSOA spatial data...")
lsoa_gdf = gpd.read_file(lsoa_path)
if lsoa_gdf.crs is None or lsoa_gdf.crs.to_epsg() != 4326:
    lsoa_gdf = lsoa_gdf.to_crs(epsg=4326)

# Spatial filter: Keep only LSOAs within LCR
lsoa_lcr = gpd.sjoin(
    lsoa_gdf,
    lcr[["geometry"]],
    how="inner",
    predicate="intersects"
)

# ================= 3. Parse and Merge Excel Poverty Data (New) =================
print("3/5 Reading Excel poverty data and merging...")

# Read Excel. Based on screenshot, header is at row 3 (index 2), Sheet name is "Table 4"
try:
    fuel_data = pd.read_excel(excel_path, sheet_name="Table 4", header=2)
except Exception as e:
    print(f"Failed to read Excel: {e}")
    # If Table 4 not found, try default read
    fuel_data = pd.read_excel(excel_path, header=2)

# Key column name mapping (based on screenshot)
# Assuming column names in Excel are as follows (tweak if actual file differs slightly)
lsoa_col_excel = "LSOA Code"
rate_col_excel = "Proportion of households fuel poor (%)"

# Data cleaning: Keep only necessary columns
fuel_data = fuel_data[[lsoa_col_excel, "LSOA Name", rate_col_excel]].copy()
fuel_data.rename(columns={rate_col_excel: "Fuel_Poverty_Rate"}, inplace=True)

# Ensure numeric type
fuel_data["Fuel_Poverty_Rate"] = pd.to_numeric(fuel_data["Fuel_Poverty_Rate"], errors="coerce")

# == Data Merge ==
# Merge GeoDataFrame (lsoa_lcr) with DataFrame (fuel_data)
# Left join: Ensures map shapes remain; areas without data will show NaN
lsoa_merged = lsoa_lcr.merge(
    fuel_data,
    left_on="LSOA21CD", # Code column in GeoJSON
    right_on=lsoa_col_excel, # Code column in Excel
    how="left"
)

# Fill missing values (prevent plotting errors), or handle in style_function later
lsoa_merged["Fuel_Poverty_Rate"] = lsoa_merged["Fuel_Poverty_Rate"].fillna(0)

print(f"   - Merge complete. {len(lsoa_merged)} areas within LCR contain poverty data.")

# ================= 4. Prepare Brownfield Data =================
print("4/5 Processing Brownfield data...")
brownfield = gpd.read_file(brownfield_path)
if brownfield.crs is None or brownfield.crs.to_epsg() != 4326:
    brownfield = brownfield.to_crs(epsg=4326)

brownfield_lcr = gpd.sjoin(
    brownfield,
    lcr[["geometry"]],
    how="inner",
    predicate="within"
)
brownfield_lcr["hectares"] = pd.to_numeric(brownfield_lcr["hectares"], errors="coerce")

# ================= 5. Plotting (Styling) =================
print("5/5 Generating final map...")

m = folium.Map(
    location=[centre_lat, centre_lon],
    zoom_start=11,
    tiles="CartoDB positron" # This basemap makes colored areas stand out more
)

# --- A. Set LSOA Color Map (Poverty Rate) ---
fp_min = lsoa_merged["Fuel_Poverty_Rate"].min()
fp_max = lsoa_merged["Fuel_Poverty_Rate"].max()

# Use YlOrRd (Yellow-Orange-Red) scale, Red indicates high poverty rate
fp_cmap = cm.LinearColormap(
    colors=['#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026'],
    vmin=fp_min,
    vmax=fp_max,
    caption="Proportion of households fuel poor (%)"
)

# --- B. Add LSOA Choropleth Layer ---
fg_lsoa = folium.FeatureGroup(name="Fuel Poverty by LSOA (2023)")

folium.GeoJson(
    lsoa_merged,
    name="LSOA Fuel Poverty",
    style_function=lambda x: {
        "fillColor": fp_cmap(x["properties"]["Fuel_Poverty_Rate"]) if x["properties"]["Fuel_Poverty_Rate"] > 0 else "#cccccc",
        "color": "white",       # Polygon outline color
        "weight": 0.5,          # Very thin outline to avoid clutter
        "fillOpacity": 0.5,     # 50% Opacity (Requirement)
    },
    tooltip=folium.GeoJsonTooltip(
        fields=["LSOA21NM", "Fuel_Poverty_Rate"],
        aliases=["Area Name:", "Fuel Poverty (%):"],
        sticky=False,
        localize=True
    )
).add_to(fg_lsoa)

fg_lsoa.add_to(m)
fp_cmap.add_to(m) # Add color bar legend

# --- C. Add LCR Boundary ---
folium.GeoJson(
    lcr,
    name="LCR Boundary",
    style_function=lambda x: {
        "fillColor": "transparent",
        "color": "black",
        "weight": 2,
        "fillOpacity": 0,
    }
).add_to(m)

# --- D. Add Brownfield Sites ---
fg_sites = folium.FeatureGroup(name="Brownfield Sites")

# Brownfield color map (based on size)
h_min, h_max = brownfield_lcr["hectares"].min(), brownfield_lcr["hectares"].max()
site_cmap = cm.LinearColormap(["blue", "cyan"], vmin=h_min, vmax=h_max) # Use cool tones to differentiate

for _, row in brownfield_lcr.iterrows():
    geom = row.geometry
    if geom is None or geom.is_empty or pd.isna(row["hectares"]):
        continue
    
    lat, lon = geom.centroid.y, geom.centroid.x
    h = row["hectares"]
    radius = 3 + 9 * (h - h_min) / (h_max - h_min + 1e-9)

    popup_text = f"<b>Site:</b> {row.get('SiteName', 'N/A')}<br><b>Hectares:</b> {h}"
    
    folium.CircleMarker(
        location=[lat, lon],
        radius=radius,
        color="#333333",      # Dark border
        weight=1,
        fillColor=site_cmap(h), # Fill color
        fillOpacity=0.9,
        popup=popup_text
    ).add_to(fg_sites)

fg_sites.add_to(m)

folium.LayerControl(collapsed=False).add_to(m)

output_file = "lcr_fuel_poverty_map.html"
m.save(output_file)
print(f"Map saved to: {output_file}")