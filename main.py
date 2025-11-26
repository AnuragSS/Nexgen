from pathlib import Path

import branca.colormap as cm
import folium
import geopandas as gpd
import pandas as pd

# ================= 0. Configuration & Paths =================
bound_path = Path("CAUTH_MAY_2025_EN_BSC_1271049543973882786.geojson")
brownfield_path = Path("brownfield-land.geojson")
lsoa_path = Path("Lower_layer_Super_Output_Areas_December_2021_Boundaries_EW_BSC_V4_-4299016806856585929.geojson")
fuel_excel_path = Path("Sub-regional_fuel_poverty_statistics_2023.xlsx")
imd_excel_path = Path("Deprivation_2025.xlsx")
police_dir = Path("police_data")

# ================= 1. Process Administrative Boundaries (LCR) =================
print("1/7 Reading administrative boundary data...")
bound_gdf = gpd.read_file(bound_path)
if bound_gdf.crs is None or bound_gdf.crs.to_epsg() != 4326:
	bound_gdf = bound_gdf.to_crs(epsg=4326)

# Filter for Liverpool City Region
lcr = bound_gdf[bound_gdf["CAUTH25NM"] == "Liverpool City Region"].copy()
if lcr.empty:
	raise ValueError("LCR boundary not found. Please check the Combined Authority GeoJSON.")

centre = lcr.geometry.unary_union.centroid
centre_lat, centre_lon = centre.y, centre.x

# ================= 2. Process LSOA Spatial Data =================
print("2/7 Reading and filtering LSOA spatial data...")
lsoa_gdf = gpd.read_file(lsoa_path)
if lsoa_gdf.crs is None or lsoa_gdf.crs.to_epsg() != 4326:
	lsoa_gdf = lsoa_gdf.to_crs(epsg=4326)

# Spatial filter: Keep only LSOAs within LCR
# We keep LSOA21CD (Code) and LSOA21NM (Name)
lsoa_lcr = gpd.sjoin(
	lsoa_gdf,
	lcr[["geometry"]],
	how="inner",
	predicate="intersects"
).drop(columns=['index_right'])

# Create a set of valid LCR LSOA codes for filtering the massive police data later
valid_lsoa_codes = set(lsoa_lcr['LSOA21CD'].unique())

# ================= 3. Load Dataset A: Fuel Poverty =================
print("3/7 Processing Fuel Poverty Data...")
try:
	fuel_data = pd.read_excel(fuel_excel_path, sheet_name="Table 4", header=2)
	fuel_data = fuel_data.rename(columns={
		"LSOA Code"                             : "LSOA21CD",
		"Proportion of households fuel poor (%)": "Fuel_Poverty_Rate"
	})
	fuel_data["Fuel_Poverty_Rate"] = pd.to_numeric(fuel_data["Fuel_Poverty_Rate"], errors="coerce")
	fuel_data = fuel_data[["LSOA21CD", "Fuel_Poverty_Rate"]]
except Exception as e:
	print(f"Warning: Fuel Poverty data error: {e}")
	fuel_data = pd.DataFrame(columns=["LSOA21CD", "Fuel_Poverty_Rate"])

# ================= 4. Load Dataset B: IMD Deprivation =================
print("4/7 Processing IMD Deprivation Data...")
try:
	imd_data = pd.read_excel(imd_excel_path, sheet_name="IMD25", header=0)
	# Search for Decile column
	imd_col = next((c for c in imd_data.columns if "Decile" in str(c) and "IMD" in str(c)), None)

	if imd_col:
		imd_data = imd_data.rename(columns={"LSOA code (2021)": "LSOA21CD", imd_col: "IMD_Decile"})
		imd_data["IMD_Decile"] = pd.to_numeric(imd_data["IMD_Decile"], errors="coerce")
		imd_data = imd_data[["LSOA21CD", "IMD_Decile"]]
	else:
		imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])
except Exception as e:
	print(f"Warning: IMD data error: {e}")
	imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])

# ================= 5. Load Dataset C: POLICE DATA (Complex) =================
print("5/7 Processing Police Data (This may take a moment)...")

street_counts = {}  # Dictionary to aggregate crimes: {LSOA_CODE: count}
stop_search_points = []  # List to hold lat/lon for spatial join later

# Recursively find all CSV files
all_csvs = list(police_dir.rglob("*.csv"))
print(f"   Found {len(all_csvs)} police data files. parsing...")

for csv_file in all_csvs:
	fname = csv_file.name.lower()

	try:
		# --- CASE 1: Street Crime (Has LSOA Code column) ---
		if "street" in fname:
			# Only read necessary columns to save RAM
			# Note: Column names in your sample are "LSOA code"
			df = pd.read_csv(csv_file, usecols=["LSOA code"])

			# Filter rows where LSOA is in our LCR list
			df_filtered = df[df["LSOA code"].isin(valid_lsoa_codes)]

			# Aggregate
			counts = df_filtered["LSOA code"].value_counts()
			for code, count in counts.items():
				street_counts[code] = street_counts.get(code, 0) + count

		# --- CASE 2: Stop and Search (Has Lat/Lon, NO LSOA Code) ---
		elif "stop-and-search" in fname:
			# Need Latitude and Longitude to spatially join to LSOA
			df = pd.read_csv(csv_file, usecols=["Latitude", "Longitude"])

			# Drop invalid coordinates
			df = df.dropna(subset=["Latitude", "Longitude"])

			# Pre-filter rough bounding box of Liverpool to speed up (optional but good)
			# Rough LCR box: Lat 53.2 to 53.7, Lon -3.3 to -2.5
			mask = (df['Latitude'].between(53.1, 53.8)) & (df['Longitude'].between(-3.4, -2.4))
			df_trimmed = df[mask]

			if not df_trimmed.empty:
				stop_search_points.append(df_trimmed)

	except ValueError as ve:
		# Happens if columns like 'LSOA code' are missing in a specific file
		continue
	except Exception as e:
		print(f"   Skipped {fname}: {e}")

# --- Aggregate Street Crimes to DataFrame ---
df_street = pd.DataFrame(list(street_counts.items()), columns=["LSOA21CD", "Crime_Count"])

# --- Process Stop & Search (Spatial Join) ---
print("   Performing spatial join for Stop & Search data...")
df_stop_agg = pd.DataFrame(columns=["LSOA21CD", "StopSearch_Count"])

if stop_search_points:
	# Combine all individual stop/search dataframes
	all_stops = pd.concat(stop_search_points, ignore_index=True)

	# Convert to GeoDataFrame
	gdf_stops = gpd.GeoDataFrame(
		all_stops,
		geometry=gpd.points_from_xy(all_stops.Longitude, all_stops.Latitude),
		crs="EPSG:4326"
	)

	# Spatial Join: Points WITHIN LSOA Polygons
	# We join stops to the previously filtered 'lsoa_lcr'
	stops_with_lsoa = gpd.sjoin(gdf_stops, lsoa_lcr[['LSOA21CD', 'geometry']], predicate='within')

	# Count by LSOA
	stop_counts = stops_with_lsoa['LSOA21CD'].value_counts().reset_index()
	stop_counts.columns = ['LSOA21CD', 'StopSearch_Count']
	df_stop_agg = stop_counts

# ================= 6. Merge All Data =================
print("6/7 Merging all datasets into Master GeoDataFrame...")

# Start with Spatial LSOAs
master_gdf = lsoa_lcr.merge(fuel_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(imd_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_street, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_stop_agg, on="LSOA21CD", how="left")

# Fill NaNs
master_gdf["Fuel_Poverty_Rate"] = master_gdf["Fuel_Poverty_Rate"].fillna(0)
master_gdf["IMD_Decile"] = master_gdf["IMD_Decile"].fillna(-1)  # -1 indicates no data
master_gdf["Crime_Count"] = master_gdf["Crime_Count"].fillna(0)
master_gdf["StopSearch_Count"] = master_gdf["StopSearch_Count"].fillna(0)

# ================= 7. Mapping =================
print("7/7 Generating Map...")

m = folium.Map(location=[centre_lat, centre_lon], zoom_start=11, tiles="CartoDB positron")

# --- Define Colormaps ---

# 1. Fuel Poverty (Yellow -> Red)
cmap_fp = cm.LinearColormap(['#ffffb2', '#fd8d3c', '#bd0026'],
                            vmin=master_gdf["Fuel_Poverty_Rate"].min(),
                            vmax=master_gdf["Fuel_Poverty_Rate"].max(),
                            caption="Fuel Poverty (%)")

# 2. IMD (Red -> Green)
cmap_imd = cm.LinearColormap(['#d73027', '#fee08b', '#1a9850'],
                             vmin=1, vmax=10,
                             caption="IMD Decile (1=Deprived)")

# 3. Crime Count (White -> Purple)
# Use 95th percentile for max to avoid skewing by city centre hotspots
crime_max = master_gdf["Crime_Count"].quantile(0.95)
cmap_crime = cm.LinearColormap(['#fcfbfd', '#9e9ac8', '#3f007d'],
                               vmin=0, vmax=crime_max,
                               caption="Reported Crimes (2022-2025)")

# 4. Stop & Search (White -> Orange)
ss_max = master_gdf["StopSearch_Count"].quantile(0.95)
cmap_ss = cm.LinearColormap(['#fff5eb', '#fd8d3c', '#7f2704'],
                            vmin=0, vmax=ss_max,
                            caption="Stop & Search Events")


# --- Helper Function for Layers ---
def add_layer(gdf, name, column, cmap, show=False):
	fg = folium.FeatureGroup(name=name, show=show)
	folium.GeoJson(
		gdf,
		style_function=lambda x: {
			"fillColor"  : cmap(x["properties"][column]) if x["properties"][column] > 0 else "#f0f0f0",
			"color"      : "transparent",
			"weight"     : 0,
			"fillOpacity": 0.6,
		},
		highlight_function=lambda x: {"weight": 2, "color": "#666", "fillOpacity": 0.8},
		tooltip=folium.GeoJsonTooltip(
			fields=["LSOA21NM", column],
			aliases=["Area:", f"{name}:"],
			sticky=False
		)
	).add_to(fg)
	fg.add_to(m)
	return cmap


# Add Polygon Layers
add_layer(master_gdf, "Fuel Poverty", "Fuel_Poverty_Rate", cmap_fp, show=True).add_to(m)
add_layer(master_gdf, "IMD Deprivation", "IMD_Decile", cmap_imd, show=False).add_to(m)
add_layer(master_gdf, "Total Crimes", "Crime_Count", cmap_crime, show=False).add_to(m)
add_layer(master_gdf, "Stop & Search", "StopSearch_Count", cmap_ss, show=False).add_to(m)

# --- Brownfield Sites (Points) ---
brownfield = gpd.read_file(brownfield_path).to_crs(epsg=4326)
bf_lcr = gpd.sjoin(brownfield, lcr[["geometry"]], predicate="within")
bf_lcr["hectares"] = pd.to_numeric(bf_lcr["hectares"], errors="coerce")

fg_sites = folium.FeatureGroup(name="Brownfield Sites", show=True)
for _, row in bf_lcr.iterrows():
	if row.geometry.is_empty or pd.isna(row["hectares"]): continue
	folium.CircleMarker(
		location=[row.geometry.centroid.y, row.geometry.centroid.x],
		radius=3 + (row["hectares"] * 0.5),  # Scale radius
		color="black", weight=1,
		fillColor="#00bcd4", fillOpacity=0.8,
		popup=f"Site: {row.get('SiteName', 'N/A')}<br>Ha: {row['hectares']}"
	).add_to(fg_sites)
fg_sites.add_to(m)

# --- Finalize ---
folium.LayerControl(collapsed=False).add_to(m)
outfile = "LCR_Integrated_Map.html"
m.save(outfile)
print(f"Done! Map saved to {outfile}")
