import geopandas as gpd
import pandas as pd
import folium
import branca.colormap as cm
from pathlib import Path

# ================= 0. Configuration Paths =================
bound_path = Path("CAUTH_MAY_2025_EN_BSC_1271049543973882786.geojson")
brownfield_path = Path("brownfield-land.geojson")
lsoa_path = Path("Lower_layer_Super_Output_Areas_December_2021_Boundaries_EW_BSC_V4_-4299016806856585929.geojson")

# Excel Paths
fuel_excel_path = Path("Sub-regional_fuel_poverty_statistics_2023.xlsx")
imd_excel_path = Path("Deprivation_2025.xlsx")

# ================= 1. Process Administrative Boundaries (LCR) =================
print("1/6 Reading administrative boundary data...")
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
print("2/6 Reading and filtering LSOA spatial data...")
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

# ================= 3. Load & Clean Dataset A: Fuel Poverty =================
print("3/6 Processing Fuel Poverty Data...")
try:
    # Sheet "Table 4", Header at index 2 (Row 3)
    fuel_data = pd.read_excel(fuel_excel_path, sheet_name="Table 4", header=2)
    
    # Select and Rename Columns
    fuel_cols = {
        "LSOA Code": "LSOA21CD", 
        "Proportion of households fuel poor (%)": "Fuel_Poverty_Rate"
    }
    # Keep only columns that exist
    valid_cols = [c for c in fuel_cols.keys() if c in fuel_data.columns]
    fuel_data = fuel_data[valid_cols].rename(columns=fuel_cols)
    
    # Convert to numeric
    fuel_data["Fuel_Poverty_Rate"] = pd.to_numeric(fuel_data["Fuel_Poverty_Rate"], errors="coerce")
    
except Exception as e:
    print(f"Error loading Fuel Poverty data: {e}")
    fuel_data = pd.DataFrame(columns=["LSOA21CD", "Fuel_Poverty_Rate"]) # Empty fallback

# ================= 4. Load & Clean Dataset B: IMD Deprivation =================
print("4/6 Processing IMD Deprivation Data...")
try:
    # Sheet "IMD25", Header at index 0 (Row 1)
    imd_data = pd.read_excel(imd_excel_path, sheet_name="IMD25", header=0)
    
    # Find the specific Decile column
    imd_col_name = "Index of Multiple Deprivation (IMD) Decile (where 1 is most deprived 10% of LSOA)"
    if imd_col_name not in imd_data.columns:
        # Fallback search
        cols = [c for c in imd_data.columns if "Decile" in str(c) and "IMD" in str(c)]
        imd_col_name = cols[0] if cols else None

    if imd_col_name:
        imd_cols = {
            "LSOA code (2021)": "LSOA21CD", 
            imd_col_name: "IMD_Decile"
        }
        imd_data = imd_data[list(imd_cols.keys())].rename(columns=imd_cols)
        imd_data["IMD_Decile"] = pd.to_numeric(imd_data["IMD_Decile"], errors="coerce")
    else:
        raise ValueError("IMD Decile column not found.")
        
except Exception as e:
    print(f"Error loading IMD data: {e}")
    imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])

# ================= 5. Merge Data into Spatial GeoDataFrame =================
print("5/6 Merging datasets...")

# Merge Fuel Poverty
lsoa_merged = lsoa_lcr.merge(fuel_data, on="LSOA21CD", how="left")

# Merge IMD
lsoa_merged = lsoa_merged.merge(imd_data, on="LSOA21CD", how="left")

# Handle NaNs for styling (fill with -1 or similar to color them grey)
lsoa_merged["Fuel_Poverty_Rate"] = lsoa_merged["Fuel_Poverty_Rate"].fillna(0)
lsoa_merged["IMD_Decile"] = lsoa_merged["IMD_Decile"].fillna(-1)

# ================= 6. Mapping =================
print("6/6 Generating Combined Map...")

m = folium.Map(
    location=[centre_lat, centre_lon],
    zoom_start=11,
    tiles="CartoDB positron"
)

# --- Define Color Maps ---

# 1. Fuel Poverty Colormap (Yellow -> Red)
fp_min, fp_max = lsoa_merged["Fuel_Poverty_Rate"].min(), lsoa_merged["Fuel_Poverty_Rate"].max()
cmap_fp = cm.LinearColormap(
    ['#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026'],
    vmin=fp_min, vmax=fp_max,
    caption="Fuel Poverty Rate (%)"
)

# 2. IMD Colormap (Red -> Green, 1 is Bad)
cmap_imd = cm.LinearColormap(
    ['#d73027', '#fc8d59', '#fee08b', '#d9ef8b', '#91cf60', '#1a9850'],
    vmin=1, vmax=10,
    caption="IMD Decile (1=Most Deprived, 10=Least)"
)

# --- Layer 1: Fuel Poverty ---
# We use FeatureGroup with 'overlay=False' logic (handled by LayerControl)
# But standard FeatureGroups overlap. We will let the user toggle.
fg_fp = folium.FeatureGroup(name="Fuel Poverty (2023)", show=True)

folium.GeoJson(
    lsoa_merged,
    name="Fuel Poverty Polygons",
    style_function=lambda x: {
        "fillColor": cmap_fp(x["properties"]["Fuel_Poverty_Rate"]) if x["properties"]["Fuel_Poverty_Rate"] > 0 else "#cccccc",
        "color": "white",
        "weight": 0.5,
        "fillOpacity": 0.5,
    },
    tooltip=folium.GeoJsonTooltip(
        fields=["LSOA21NM", "Fuel_Poverty_Rate"],
        aliases=["Area:", "Fuel Poverty %:"],
        sticky=False
    )
).add_to(fg_fp)
fg_fp.add_to(m)

# --- Layer 2: IMD Deprivation ---
fg_imd = folium.FeatureGroup(name="IMD Deprivation (2025)", show=False) # Hidden by default

folium.GeoJson(
    lsoa_merged,
    name="IMD Polygons",
    style_function=lambda x: {
        "fillColor": cmap_imd(x["properties"]["IMD_Decile"]) if x["properties"]["IMD_Decile"] > 0 else "#cccccc",
        "color": "white",
        "weight": 0.5,
        "fillOpacity": 0.5,
    },
    tooltip=folium.GeoJsonTooltip(
        fields=["LSOA21NM", "IMD_Decile"],
        aliases=["Area:", "IMD Decile:"],
        sticky=False
    )
).add_to(fg_imd)
fg_imd.add_to(m)

# --- Layer 3: Brownfield Sites (Always Top) ---
brownfield = gpd.read_file(brownfield_path).to_crs(epsg=4326)
bf_lcr = gpd.sjoin(brownfield, lcr[["geometry"]], predicate="within")
bf_lcr["hectares"] = pd.to_numeric(bf_lcr["hectares"], errors="coerce")

fg_sites = folium.FeatureGroup(name="Brownfield Sites", show=True)
h_min, h_max = bf_lcr["hectares"].min(), bf_lcr["hectares"].max()
# Blue/Cyan scale for contrast
cmap_sites = cm.LinearColormap(["#2171b5", "#08306b"], vmin=h_min, vmax=h_max) 

for _, row in bf_lcr.iterrows():
    if row.geometry.is_empty or pd.isna(row["hectares"]): continue
    
    folium.CircleMarker(
        location=[row.geometry.centroid.y, row.geometry.centroid.x],
        radius=3 + 9 * (row["hectares"] - h_min) / (h_max - h_min + 1e-9),
        color="#ffffff",
        weight=1,
        fillColor=cmap_sites(row["hectares"]),
        fillOpacity=0.9,
        popup=f"Site: {row.get('SiteName', 'N/A')}<br>Ha: {row['hectares']}"
    ).add_to(fg_sites)

fg_sites.add_to(m)

# --- Add Legends ---
# Note: In standard Folium, legends are static. Both will appear on the map.
cmap_fp.add_to(m)
cmap_imd.add_to(m)

# --- Add Controls ---
# LayerControl allows toggling layers on/off
folium.LayerControl(collapsed=False).add_to(m)

# --- Save ---
outfile = "LCR_Combined_Map.html"
m.save(outfile)
print(f"Map saved successfully to: {outfile}")
print("Open the HTML file and use the Layer Control (top right) to switch between Fuel Poverty and IMD.")