import json
from pathlib import Path

import folium
import geopandas as gpd
import pandas as pd
from branca.element import Element

# ================= 0. Configuration =================
bound_path = Path("CAUTH_MAY_2025_EN_BSC_1271049543973882786.geojson")
brownfield_path = Path("brownfield-land.geojson")
lsoa_path = Path(
	"Lower_layer_Super_Output_Areas_December_2021_Boundaries_EW_BSC_V4_-4299016806856585929.geojson"
)
fuel_excel_path = Path("Sub-regional_fuel_poverty_statistics_2023.xlsx")
imd_excel_path = Path("Deprivation_2025.xlsx")
police_dir = Path("police_data")

# ================= 1. Data Processing (Backend) =================
print(">>> Step 1/5: Loading and Processing Data...")

# 1.1 Boundaries
print("   -> Loading LCR Boundary...")
bound_gdf = gpd.read_file(bound_path)
if bound_gdf.crs is not None and bound_gdf.crs.to_epsg() != 4326:
    bound_gdf = bound_gdf.to_crs(epsg=4326)

lcr = bound_gdf[bound_gdf["CAUTH25NM"] == "Liverpool City Region"].copy()
if lcr.empty:
	raise ValueError("ERROR: Could not find 'Liverpool City Region' in boundary file!")

# Fix deprecation warning: use union_all() for newer geopandas
try:
	centre = lcr.geometry.union_all().centroid
except AttributeError:
	# Fallback for older geopandas
	centre = lcr.geometry.unary_union.centroid
centre_lat, centre_lon = centre.y, centre.x

# 1.2 LSOA Polygons
print("   -> Loading LSOA Polygons...")
lsoa_gdf = gpd.read_file(lsoa_path)
if lsoa_gdf.crs is not None and lsoa_gdf.crs.to_epsg() != 4326:
    lsoa_gdf = lsoa_gdf.to_crs(epsg=4326)

# Spatial Join
lsoa_lcr = gpd.sjoin(
	lsoa_gdf, lcr[["geometry"]], how="inner", predicate="intersects"
).drop(columns=["index_right"])
print(f"      LSOAs inside Liverpool City Region: {len(lsoa_lcr)}")

valid_lsoa_codes = set(lsoa_lcr["LSOA21CD"].unique())

# 1.3 Load Metrics
print("   -> Loading Metrics...")
# Fuel
try:
    fuel_data = pd.read_excel(fuel_excel_path, sheet_name="Table 4", header=2)
    fuel_data = fuel_data.rename(
	    columns={
		    "LSOA Code"                             : "LSOA21CD",
		    "Proportion of households fuel poor (%)": "Fuel_Poverty_Rate",
	    }
    )
    fuel_data["Fuel_Poverty_Rate"] = pd.to_numeric(
	    fuel_data["Fuel_Poverty_Rate"], errors="coerce"
    )
except:
    fuel_data = pd.DataFrame(columns=["LSOA21CD", "Fuel_Poverty_Rate"])

# IMD
try:
    imd_data = pd.read_excel(imd_excel_path, sheet_name="IMD25", header=0)
    imd_col = next(
	    (c for c in imd_data.columns if "Decile" in str(c) and "IMD" in str(c)), None
    )
    if imd_col:
	    imd_data = imd_data.rename(
		    columns={"LSOA code (2021)": "LSOA21CD", imd_col: "IMD_Decile"}
	    )
        imd_data["IMD_Decile"] = pd.to_numeric(imd_data["IMD_Decile"], errors="coerce")
    else:
        imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])
except:
    imd_data = pd.DataFrame(columns=["LSOA21CD", "IMD_Decile"])

# Police
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
            mask = (df["Latitude"].between(53.0, 54.0)) & (
	            df["Longitude"].between(-3.5, -2.0)
            )
            df_trimmed = df[mask]
            if not df_trimmed.empty:
                stop_search_points.append(df_trimmed)
    except:
        continue

df_street = pd.DataFrame(
	list(street_counts.items()), columns=["LSOA21CD", "Crime_Count"]
)

df_stop_agg = pd.DataFrame(columns=["LSOA21CD", "StopSearch_Count"])
if stop_search_points:
    all_stops = pd.concat(stop_search_points, ignore_index=True)
    gdf_stops = gpd.GeoDataFrame(
	    all_stops,
	    geometry=gpd.points_from_xy(all_stops.Longitude, all_stops.Latitude),
	    crs="EPSG:4326",
    )
    stops_with_lsoa = gpd.sjoin(
	    gdf_stops, lsoa_lcr[["LSOA21CD", "geometry"]], predicate="within"
    )
    stop_counts = stops_with_lsoa["LSOA21CD"].value_counts().reset_index()
    stop_counts.columns = ["LSOA21CD", "StopSearch_Count"]
    df_stop_agg = stop_counts

# 1.4 Merge
master_gdf = lsoa_lcr.merge(fuel_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(imd_data, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_street, on="LSOA21CD", how="left")
master_gdf = master_gdf.merge(df_stop_agg, on="LSOA21CD", how="left")

master_gdf["Fuel_Poverty_Rate"] = master_gdf["Fuel_Poverty_Rate"].fillna(0)
master_gdf["IMD_Decile"] = master_gdf["IMD_Decile"].fillna(10)
master_gdf["Crime_Count"] = master_gdf["Crime_Count"].fillna(0)
master_gdf["StopSearch_Count"] = master_gdf["StopSearch_Count"].fillna(0)

# 1.5 Normalization
print(">>> Step 2/5: Normalizing Metrics...")
f_min, f_max = (
	master_gdf["Fuel_Poverty_Rate"].min(),
	master_gdf["Fuel_Poverty_Rate"].max(),
)
master_gdf["norm_fuel"] = (master_gdf["Fuel_Poverty_Rate"] - f_min) / (
		f_max - f_min + 1e-9
)
master_gdf["norm_imd"] = (11 - master_gdf["IMD_Decile"]) / 10.0
c_max = master_gdf["Crime_Count"].quantile(0.95)
master_gdf["norm_crime"] = (master_gdf["Crime_Count"] / c_max).clip(upper=1.0)
s_max = master_gdf["StopSearch_Count"].quantile(0.95)
master_gdf["norm_stop"] = (master_gdf["StopSearch_Count"] / s_max).clip(upper=1.0)

cols_to_round = [
	"norm_fuel",
	"norm_imd",
	"norm_crime",
	"norm_stop",
	"Fuel_Poverty_Rate",
]
master_gdf[cols_to_round] = master_gdf[cols_to_round].round(3)

print("   -> Simplifying Geometries...")
master_gdf.geometry = master_gdf.geometry.simplify(
	tolerance=0.0001, preserve_topology=True
)

# 1.6 Brownfield Data
print(">>> Step 3/5: Preparing Brownfield Data...")
brownfield = gpd.read_file(brownfield_path).to_crs(epsg=4326)
bf_joined = gpd.sjoin(
	brownfield, master_gdf[["geometry", "LSOA21CD"]], predicate="within"
)
bf_joined["hectares"] = pd.to_numeric(bf_joined["hectares"], errors="coerce").fillna(0)

bf_data_list = []
for _, row in bf_joined.iterrows():
	if row.geometry.is_empty:
		continue
	bf_data_list.append(
		{
			"lat" : row.geometry.centroid.y,
			"lng" : row.geometry.centroid.x,
			"ha"  : round(row["hectares"], 2),
			"name": row.get("SiteName", "Unknown"),
			"lsoa": row["LSOA21CD"],
		}
	)

# ================= 2. Generating Map with Custom JS =================
print(">>> Step 4/5: Building Interactive Map...")

m = folium.Map(
	location=[centre_lat, centre_lon], zoom_start=11, tiles="CartoDB positron"
)

# 2.1 Add Custom Panes
folium.map.CustomPane("poly_pane", z_index=400).add_to(m)
folium.map.CustomPane("brownfield_pane", z_index=650).add_to(m)

# 2.2 Inject Raw GeoJSON
geojson_str = master_gdf.to_json()

# 2.3 The HTML/CSS/JS Interface
custom_ui = f"""
<style>
    #controls {{
        position: absolute; bottom: 30px; right: 10px; width: 300px;
        background: rgba(255, 255, 255, 0.95); backdrop-filter: blur(5px);
        padding: 15px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.3);
        z-index: 1000; font-family: sans-serif;
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
</style>

<div id="controls">
    <h4>Layer Weights (权重)</h4>
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
    var mapInstance = null; // Store the Found map
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

    // 4. FIND FOLIUM MAP (The Fix)
    function findFoliumMap() {{
        for (var key in window) {{
            // Folium maps are instances of L.Map and usually stored in global variables starting with 'map_'
            if (window[key] instanceof L.Map) {{
                console.log("Found Folium Map Instance: " + key);
                return window[key];
            }}
        }}
        return null;
    }}

    // 5. INIT LAYERS
    function initMap() {{
        // Try to find the map
        mapInstance = findFoliumMap();
        
        if (!mapInstance) {{
            console.log("Map not ready yet, retrying...");
            setTimeout(initMap, 200); // Retry loop
            return;
        }}

        console.log("Map ready. initializing layers...");

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
                var content = "<b>" + p.LSOA21NM + "</b><br>" +
                              "Fuel Norm: " + (p.norm_fuel || 0) + "<br>" +
                              "IMD Norm: " + (p.norm_imd || 0);
                layer.bindTooltip(content);
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
            circle.bindPopup("<b>Site:</b> " + site.name + "<br>Ha: " + site.ha);
            circle.lsoaCode = site.lsoa; 
            bfLayerGroup.addLayer(circle);
        }});

        // Trigger first calculation
        updateWeights();
    }}

    // 6. UPDATE FUNCTION
    function updateWeights() {{
        // Ensure layers are initialized before trying to update them
        if (!geojsonLayer || !bfLayerGroup) return;

        weights.fuel = parseFloat(document.getElementById('w_fuel').value);
        weights.imd = parseFloat(document.getElementById('w_imd').value);
        weights.crime = parseFloat(document.getElementById('w_crime').value);
        weights.stop = parseFloat(document.getElementById('w_stop').value);

        // Update UI Text
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

    // Start Initialization Loop
    initMap();
</script>
"""

m.get_root().html.add_child(Element(custom_ui))

# ================= 3. Output =================
outfile = "LCR_Interactive_Weighted_Map.html"
m.save(outfile)
print(f">>> Step 5/5: Done! Saved to {outfile}")
