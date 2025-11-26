import json
from pathlib import Path

import folium
import geopandas as gpd
import numpy as np
import pandas as pd
import pgeocode
from branca.element import Element
from folium.plugins import HeatMap

# ================= 0. Configuration =================
bound_path = Path("CAUTH_MAY_2025_EN_BSC_1271049543973882786.geojson")
brownfield_path = Path("brownfield-land.geojson")
lsoa_path = Path("Lower_layer_Super_Output_Areas_December_2021_Boundaries_EW_BSC_V4_-4299016806856585929.geojson")
fuel_excel_path = Path("Sub-regional_fuel_poverty_statistics_2023.xlsx")
imd_excel_path = Path("Deprivation_2025.xlsx")
housing_excel_path = Path("Copy of Housing_Pipeline_Long_List_External_Hackathon.xlsx")  # 确保文件名完全一致
police_dir = Path("police_data")

# ================= 1. Data Processing (Backend) =================
print(">>> Step 1/6: Loading Administrative & Metric Data...")

# 1.1 Boundaries & LSOA
bound_gdf = gpd.read_file(bound_path)
if bound_gdf.crs is not None and bound_gdf.crs.to_epsg() != 4326:
	bound_gdf = bound_gdf.to_crs(epsg=4326)

lcr = bound_gdf[bound_gdf["CAUTH25NM"] == "Liverpool City Region"].copy()
if lcr.empty:
	raise ValueError("ERROR: Could not find 'Liverpool City Region' in boundary file!")

try:
	centre = lcr.geometry.union_all().centroid
except AttributeError:
	centre = lcr.geometry.unary_union.centroid
centre_lat, centre_lon = centre.y, centre.x

lsoa_gdf = gpd.read_file(lsoa_path)
if lsoa_gdf.crs is not None and lsoa_gdf.crs.to_epsg() != 4326:
	lsoa_gdf = lsoa_gdf.to_crs(epsg=4326)

lsoa_lcr = gpd.sjoin(lsoa_gdf, lcr[["geometry"]], how="inner", predicate="intersects").drop(columns=['index_right'])
valid_lsoa_codes = set(lsoa_lcr['LSOA21CD'].unique())

# 1.2 Metrics (Fuel, IMD, Police)
try:
	fuel_data = pd.read_excel(fuel_excel_path, sheet_name="Table 4", header=2)
	fuel_data = fuel_data.rename(columns={
		"LSOA Code"                             : "LSOA21CD",
		"Proportion of households fuel poor (%)": "Fuel_Poverty_Rate"})
	fuel_data["Fuel_Poverty_Rate"] = pd.to_numeric(fuel_data["Fuel_Poverty_Rate"], errors="coerce")
except:
	fuel_data = pd.DataFrame(columns=["LSOA21CD", "Fuel_Poverty_Rate"])

try:
	imd_data = pd.read_excel(imd_excel_path, sheet_name="IMD25", header=0)
	imd_col = next((c for c in imd_data.columns if "Decile" in str(c) and "IMD" in str(c)), None)
	if imd_col:
		imd_data = imd_data.rename(columns={"LSOA code (2021)": "LSOA21CD", imd_col: "IMD_Decile"})
		imd_data["IMD_Decile"] = pd.to_numeric(imd_data["IMD_Decile"], errors="coerce")
	else:
		imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])
except:
	imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])

street_counts = {}
stop_search_points = []
all_csvs = list(police_dir.rglob("*.csv"))

for csv_file in all_csvs:
	fname = csv_file.name.lower()
	try:
		if "street" in fname:
			df = pd.read_csv(csv_file, usecols=["LSOA code"])
			df = df[df["LSOA code"].isin(valid_lsoa_codes)]
			counts = df["LSOA code"].value_counts()
			for code, count in counts.items():
				street_counts[code] = street_counts.get(code, 0) + count
		elif "stop-and-search" in fname:
			df = pd.read_csv(csv_file, usecols=["Latitude", "Longitude"])
			df = df.dropna(subset=["Latitude", "Longitude"])
			mask = (df['Latitude'].between(53.0, 54.0)) & (df['Longitude'].between(-3.5, -2.0))
			df_trimmed = df[mask]
			if not df_trimmed.empty:
				stop_search_points.append(df_trimmed)
	except:
		continue

df_street = pd.DataFrame(list(street_counts.items()), columns=["LSOA21CD", "Crime_Count"])
df_stop_agg = pd.DataFrame(columns=["LSOA21CD", "StopSearch_Count"])
if stop_search_points:
	all_stops = pd.concat(stop_search_points, ignore_index=True)
	gdf_stops = gpd.GeoDataFrame(all_stops,
	                             geometry=gpd.points_from_xy(all_stops.Longitude, all_stops.Latitude),
	                             crs="EPSG:4326")
	stops_with_lsoa = gpd.sjoin(gdf_stops, lsoa_lcr[['LSOA21CD', 'geometry']], predicate='within')
	stop_counts = stops_with_lsoa['LSOA21CD'].value_counts().reset_index()
	stop_counts.columns = ['LSOA21CD', 'StopSearch_Count']
	df_stop_agg = stop_counts

# 1.3 Merge & Normalize
master_gdf = lsoa_lcr.merge(fuel_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(imd_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_street, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_stop_agg, on="LSOA21CD", how="left")

master_gdf["Fuel_Poverty_Rate"] = master_gdf["Fuel_Poverty_Rate"].fillna(0)
master_gdf["IMD_Decile"] = master_gdf["IMD_Decile"].fillna(10)
master_gdf["Crime_Count"] = master_gdf["Crime_Count"].fillna(0)
master_gdf["StopSearch_Count"] = master_gdf["StopSearch_Count"].fillna(0)

f_min, f_max = master_gdf["Fuel_Poverty_Rate"].min(), master_gdf["Fuel_Poverty_Rate"].max()
master_gdf["norm_fuel"] = (master_gdf["Fuel_Poverty_Rate"] - f_min) / (f_max - f_min + 1e-9)
master_gdf["norm_imd"] = (11 - master_gdf["IMD_Decile"]) / 10.0
c_max = master_gdf["Crime_Count"].quantile(0.95)
master_gdf["norm_crime"] = (master_gdf["Crime_Count"] / c_max).clip(upper=1.0)
s_max = master_gdf["StopSearch_Count"].quantile(0.95)
master_gdf["norm_stop"] = (master_gdf["StopSearch_Count"] / s_max).clip(upper=1.0)

cols_to_round = ["norm_fuel", "norm_imd", "norm_crime", "norm_stop", "Fuel_Poverty_Rate"]
master_gdf[cols_to_round] = master_gdf[cols_to_round].round(3)

print("   -> Simplifying LSOA Geometries...")
master_gdf.geometry = master_gdf.geometry.simplify(tolerance=0.0001, preserve_topology=True)

# ================= 2. Process Housing Pipeline (With Debugging) =================
print(">>> Step 2/6: Processing Housing Pipeline Data...")
housing_geo = pd.DataFrame()

try:
	if not housing_excel_path.exists():
		print(f"   [ERROR] Housing Excel file not found at: {housing_excel_path}")
	else:
		housing_df = pd.read_excel(housing_excel_path)
		print(f"   Loaded Housing Excel. Columns found: {housing_df.columns.tolist()}")

		# Column matching
		cols = housing_df.columns.tolist()


		def pick_col(substring):
			for c in cols:
				if substring.lower() in c.lower(): return c
			return None


		# Map columns
		col_map = {
			pick_col("Project Name")      : "project_name",
			pick_col("LA")                : "la",
			pick_col("Pipeline First Cut"): "pipeline_stage",
			pick_col("Postcode")          : "postcode",
			pick_col("Developer")         : "developer",
			pick_col("One Line Summary")  : "summary",
			pick_col("Number Of Units")   : "units",
			pick_col("Site Size")         : "site_size_ha"
		}

		# Check mapped columns
		print(f"   Mapped columns: {col_map}")

		# Filter out missing columns
		col_map = {k: v for k, v in col_map.items() if k is not None}

		housing_clean = housing_df[col_map.keys()].rename(columns=col_map)

		# Fill numbers
		if "units" in housing_clean.columns:
			housing_clean["units"] = pd.to_numeric(housing_clean["units"], errors="coerce").fillna(0)
		else:
			housing_clean["units"] = 0

		if "site_size_ha" in housing_clean.columns:
			housing_clean["site_size_ha"] = pd.to_numeric(housing_clean["site_size_ha"], errors="coerce").fillna(0)
		else:
			housing_clean["site_size_ha"] = 0

		# Geocode
		if "postcode" in housing_clean.columns:
			print("   -> Geocoding Housing Postcodes (pgeocode)...")
			nomi = pgeocode.Nominatim("GB")
			housing_clean["postcode_clean"] = housing_clean["postcode"].astype(str).str.strip()

			# Batch geocode
			unique_pcs = housing_clean["postcode_clean"].unique()
			print(f"      Attempting to geocode {len(unique_pcs)} unique postcodes...")

			geo_results = nomi.query_postal_code(unique_pcs)

			pc_to_lat = dict(zip(unique_pcs, geo_results.latitude))
			pc_to_lon = dict(zip(unique_pcs, geo_results.longitude))

			housing_clean["lat"] = housing_clean["postcode_clean"].map(pc_to_lat)
			housing_clean["lon"] = housing_clean["postcode_clean"].map(pc_to_lon)

			housing_geo = housing_clean.dropna(subset=["lat", "lon"]).copy()
			print(f"      [SUCCESS] Geocoded {len(housing_geo)} / {len(housing_clean)} sites.")

			if len(housing_geo) > 0:
				# Create Bands
				housing_geo["units_band"] = pd.cut(
					housing_geo["units"],
					bins=[0, 20, 50, 100, 200, np.inf],
					labels=["0–20 homes", "21–50 homes", "51–100 homes", "101–200 homes", "200+ homes"],
					include_lowest=True
				)

				housing_geo["size_band"] = pd.cut(
					housing_geo["site_size_ha"],
					bins=[0, 0.5, 2, 5, 10, np.inf],
					labels=["<0.5 ha", "0.5–2 ha", "2–5 ha", "5–10 ha", "10+ ha"],
					include_lowest=True
				)
			else:
				print(
					"      [WARNING] No sites were successfully geocoded. Check postcode format or pgeocode installation.")
		else:
			print("   [ERROR] 'Postcode' column not found in Excel.")

except Exception as e:
	print(f"   [CRITICAL ERROR] processing Housing Pipeline: {e}")
	housing_geo = pd.DataFrame()

# ================= 3. Process Brownfield Data =================
print(">>> Step 3/6: Preparing Brownfield Data...")
brownfield = gpd.read_file(brownfield_path).to_crs(epsg=4326)
bf_joined = gpd.sjoin(brownfield, master_gdf[["geometry", "LSOA21CD"]], predicate="within")
bf_joined["hectares"] = pd.to_numeric(bf_joined["hectares"], errors="coerce").fillna(0)

bf_data_list = []
for _, row in bf_joined.iterrows():
	if row.geometry.is_empty: continue
	bf_data_list.append({
		"lat"      : row.geometry.centroid.y,
		"lng"      : row.geometry.centroid.x,
		"ha"       : round(row["hectares"], 2),
		"name": row.get("name", row.get("SiteName", "Unknown Site")),
		"address"  : row.get("site-address", "No Address Provided"),
		"status"   : row.get("planning-permission-status", "N/A"),
		"dwellings": row.get("minimum-net-dwellings", "N/A"),
		"lsoa"     : row["LSOA21CD"]
	})

# ================= 4. Generating Map =================
print(">>> Step 4/6: Initializing Map...")

m = folium.Map(location=[centre_lat, centre_lon], zoom_start=11, tiles="CartoDB positron")

# 4.1 Define Z-Index Panes
folium.map.CustomPane("poly_pane", z_index=400).add_to(m)
folium.map.CustomPane("housing_pane", z_index=620).add_to(m)
folium.map.CustomPane("brownfield_pane", z_index=650).add_to(m)

# 4.2 Add Housing Pipeline Layers (Standard Folium)
if not housing_geo.empty:
	print(f">>> Step 5/6: Adding {len(housing_geo)} Housing Markers to Map...")

	stage_colors = {"Short": "green", "Medium": "orange", "Long": "red"}
	size_colors = {
		"<0.5 ha": "#999999", "0.5–2 ha": "#377eb8", "2–5 ha": "#4daf4a",
		"5–10 ha": "#984ea3", "10+ ha": "#e41a1c"
	}

	# A. Marker Layers (Grouped by Unit Band)
	unique_bands = housing_geo["units_band"].unique()
	# Filter out NaNs explicitly for the loop
	unique_bands = [b for b in unique_bands if not pd.isna(b)]

	fg_by_units = {lab: folium.FeatureGroup(name=f"Housing: {lab}", show=True) for lab in unique_bands}

	for _, row in housing_geo.iterrows():
		ub = row["units_band"]
		if pd.isna(ub) or ub not in fg_by_units: continue

		radius = 4 + (np.sqrt(row["units"]) if row["units"] > 0 else 0)
		fill_color = stage_colors.get(str(row["pipeline_stage"]).strip(), "blue")
		border_color = size_colors.get(str(row["size_band"]), "#333")

		popup_html = f"""
        <div style="font-family:sans-serif; font-size:12px; width:200px;">
            <h4 style="margin:0 0 5px; border-bottom:2px solid {fill_color};">{row['project_name']}</h4>
            <b>Units:</b> {int(row['units'])}<br>
            <b>Stage:</b> {row['pipeline_stage']}<br>
            <b>Size:</b> {row['site_size_ha']} ha<br>
            <b>Developer:</b> {row['developer']}<br>
            <hr style="margin:5px 0;">
            <i>{row['summary']}</i>
        </div>
        """

		folium.CircleMarker(
			location=[row["lat"], row["lon"]],
			radius=radius,
			color=border_color, weight=2,
			fill=True, fill_color=fill_color, fill_opacity=0.8,
			popup=folium.Popup(popup_html, max_width=300),
			tooltip=f"{row['project_name']} ({int(row['units'])} units)",
			pane="housing_pane"
		).add_to(fg_by_units[ub])

	for fg in fg_by_units.values():
		fg.add_to(m)

	# B. Heatmaps
	h_units = housing_geo[housing_geo["units"] > 0].copy()
	if not h_units.empty:
		h_units["w"] = h_units["units"] / h_units["units"].max()
		hm_layer = folium.FeatureGroup(name="Heatmap: Housing Units", show=False)
		HeatMap(
			data=list(zip(h_units["lat"], h_units["lon"], h_units["w"])),
			radius=20, blur=15, min_opacity=0.3
		).add_to(hm_layer)
		hm_layer.add_to(m)
else:
	print(">>> Step 5/6: Skipping Housing layers (No valid data found).")

# 4.3 Prepare Data for Custom JS
geojson_str = master_gdf.to_json()

# 4.4 The HTML/CSS/JS Interface
custom_ui = f"""
<style>
    #controls {{
        position: absolute; bottom: 30px; right: 10px; width: 300px;
        background: rgba(255, 255, 255, 0.95); backdrop-filter: blur(5px);
        padding: 15px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.3);
        z-index: 1000; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }}
    #controls h4 {{ margin: 0 0 10px 0; color: #333; }}
    .slider-container {{ margin-bottom: 12px; }}
    .slider-label {{ display: flex; justify-content: space-between; font-size: 12px; font-weight: 600; color: #555; }}
    input[type=range] {{ width: 100%; margin: 5px 0; cursor: pointer; }}
    .legend-gradient {{
        height: 10px; background: linear-gradient(to right, #ffffff, #fee5d9, #de2d26, #a50f15);
        border-radius: 5px; margin-top: 10px; border: 1px solid #ccc;
    }}
    .legend-labels {{ display: flex; justify-content: space-between; font-size: 10px; color: #666; margin-top: 2px; }}
    
    /* Popup Styling */
    .popup-table td {{ padding: 2px 5px; border-bottom: 1px solid #eee; }}
    .popup-header {{ margin: 0 0 5px 0; border-bottom: 2px solid #a50f15; padding-bottom: 3px; font-size: 14px; color: #333; }}
</style>

<div id="controls">
    <h4>Layer Weights (LSOA Scoring)</h4>
    <div class="slider-container"><div class="slider-label"><span>Fuel Poverty</span><span id="val_fuel">1.0</span></div><input type="range" id="w_fuel" min="0" max="5" step="0.1" value="1.0" oninput="updateWeights()"></div>
    <div class="slider-container"><div class="slider-label"><span>IMD Deprivation</span><span id="val_imd">1.0</span></div><input type="range" id="w_imd" min="0" max="5" step="0.1" value="1.0" oninput="updateWeights()"></div>
    <div class="slider-container"><div class="slider-label"><span>Crime Density</span><span id="val_crime">1.0</span></div><input type="range" id="w_crime" min="0" max="5" step="0.1" value="1.0" oninput="updateWeights()"></div>
    <div class="slider-container"><div class="slider-label"><span>Stop & Search</span><span id="val_stop">1.0</span></div><input type="range" id="w_stop" min="0" max="5" step="0.1" value="1.0" oninput="updateWeights()"></div>
    <div class="legend-gradient"></div>
    <div class="legend-labels"><span>Safe (Low)</span><span>Vulnerable (High)</span></div>
</div>

<script>
    console.log("Initializing Custom Map Script...");
    
    // 1. DATA INJECTION
    try {{
        var lsoaData = {geojson_str};
        var brownfieldData = {json.dumps(bf_data_list)};
    }} catch(e) {{
        console.error("Data Injection Failed:", e);
    }}
    
    // 2. STATE
    var mapInstance = null;
    var weights = {{ fuel: 1.0, imd: 1.0, crime: 1.0, stop: 1.0 }};
    var geojsonLayer = null;
    var bfLayerGroup = null;
    var lsoaScores = {{}}; 

    // 3. COLOR SCALE
    function getColor(score) {{
        if(isNaN(score)) return "#cccccc";
        if(score > 1) score = 1;
        if(score < 0) score = 0;
        
        var r = Math.round(255 + (165 - 255) * score);
        var g = Math.round(255 + (15 - 255) * score);
        var b = Math.round(255 + (21 - 255) * score);
        return "rgb(" + r + "," + g + "," + b + ")";
    }}

    function findFoliumMap() {{
        for (var key in window) {{
            if (window[key] instanceof L.Map) return window[key];
        }}
        return null;
    }}

    // 4. INIT LAYERS
    function initMap() {{
        mapInstance = findFoliumMap();
        if (!mapInstance) {{
            setTimeout(initMap, 200);
            return;
        }}

        // Create Layer Group for Brownfields
        bfLayerGroup = L.layerGroup().addTo(mapInstance);

        // Add Polygons
        geojsonLayer = L.geoJson(lsoaData, {{
            pane: 'poly_pane',
            style: function(feature) {{
                return {{
                    fillColor: '#999',
                    weight: 0.5,
                    opacity: 0.5,
                    color: '#ccc',
                    fillOpacity: 0.6 
                }};
            }},
            onEachFeature: function(feature, layer) {{
                var p = feature.properties;
                var content = `
                    <div style="font-family: sans-serif; font-size: 12px; min-width: 220px;">
                        <h4 class="popup-header">${{p.LSOA21NM}}</h4>
                        <table class="popup-table" style="width:100%; border-collapse:collapse;">
                            <tr><td style="color:#666;">Code</td><td>${{p.LSOA21CD}}</td></tr>
                            <tr><td>Fuel Poverty</td><td><b>${{p.Fuel_Poverty_Rate}}%</b> <small>(${{p.norm_fuel}})</small></td></tr>
                            <tr><td>IMD Decile</td><td><b>${{p.IMD_Decile}}</b> <small>(${{p.norm_imd}})</small></td></tr>
                            <tr><td>Crimes</td><td><b>${{p.Crime_Count}}</b> <small>(${{p.norm_crime}})</small></td></tr>
                            <tr><td>Stop & Search</td><td><b>${{p.StopSearch_Count}}</b> <small>(${{p.norm_stop}})</small></td></tr>
                        </table>
                    </div>
                `;
                layer.bindPopup(content);
                layer.bindTooltip(p.LSOA21NM, {{sticky: true}});
            }}
        }}).addTo(mapInstance);

        // Add Brownfields
        brownfieldData.forEach(function(site) {{
            var circle = L.circleMarker([site.lat, site.lng], {{
                pane: 'brownfield_pane', 
                radius: 4 + (site.ha * 0.4),
                fillColor: '#333',
                color: '#000',
                weight: 1,
                opacity: 1,
                fillOpacity: 0.9
            }});
            
            var bfContent = `
                <div style="font-family: sans-serif; font-size: 12px; min-width: 200px;">
                    <h4 class="popup-header">${{site.name}}</h4>
                    <div style="margin-bottom:8px; color:#555;">${{site.address}}</div>
                    <table class="popup-table" style="width:100%;">
                        <tr><td><strong>Size:</strong></td><td>${{site.ha}} ha</td></tr>
                        <tr><td><strong>Dwellings:</strong></td><td>${{site.dwellings}}</td></tr>
                        <tr><td><strong>Status:</strong></td><td>${{site.status}}</td></tr>
                    </table>
                </div>
            `;
            
            circle.bindPopup(bfContent);
            circle.lsoaCode = site.lsoa; 
            bfLayerGroup.addLayer(circle);
        }});

        updateWeights();
    }}

    // 5. UPDATE FUNCTION
    function updateWeights() {{
        if (!geojsonLayer || !bfLayerGroup) return;

        weights.fuel = parseFloat(document.getElementById('w_fuel').value);
        weights.imd = parseFloat(document.getElementById('w_imd').value);
        weights.crime = parseFloat(document.getElementById('w_crime').value);
        weights.stop = parseFloat(document.getElementById('w_stop').value);

        document.getElementById('val_fuel').innerText = weights.fuel.toFixed(1);
        document.getElementById('val_imd').innerText = weights.imd.toFixed(1);
        document.getElementById('val_crime').innerText = weights.crime.toFixed(1);
        document.getElementById('val_stop').innerText = weights.stop.toFixed(1);

        var totalWeight = weights.fuel + weights.imd + weights.crime + weights.stop;
        if (totalWeight <= 0) totalWeight = 1; 

        // Update Polygons
        geojsonLayer.eachLayer(function(layer) {{
            var p = layer.feature.properties;
            var score = (
                ((p.norm_fuel || 0) * weights.fuel) + 
                ((p.norm_imd || 0) * weights.imd) + 
                ((p.norm_crime || 0) * weights.crime) + 
                ((p.norm_stop || 0) * weights.stop)
            ) / totalWeight;

            lsoaScores[p.LSOA21CD] = score;

            layer.setStyle({{
                fillColor: getColor(score)
            }});
        }});

        // Update Brownfields
        bfLayerGroup.eachLayer(function(layer) {{
            var code = layer.lsoaCode;
            var areaScore = lsoaScores[code];
            if (areaScore === undefined) areaScore = 0; 
            layer.setStyle({{
                fillColor: getColor(areaScore)
            }});
        }});
    }}

    initMap();
</script>
"""

m.get_root().html.add_child(Element(custom_ui))

# IMPORTANT: Add Layer Control at the end so it sees all layers
folium.LayerControl(collapsed=False).add_to(m)

# ================= 5. Output =================
outfile = "LCR_Integrated_Full_Map.html"
m.save(outfile)
print(f">>> Step 6/6: Done! Saved to {outfile}")
