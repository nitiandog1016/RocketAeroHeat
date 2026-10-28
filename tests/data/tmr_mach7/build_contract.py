"""Build the bounded Mach-7/20-degree data contract from unchanged NASA/TMR files.

This script performs parsing and unit conversion only; it does not run CFD.
"""
import hashlib
import json
import math
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent


def tecplot_numeric(name, columns):
    rows = []
    for line in (ROOT / "raw" / name).read_text().splitlines():
        if re.match(r"^\s*[-+0-9.]", line):
            values = [float(value) for value in line.split()]
            if len(values) != columns:
                raise ValueError("Unexpected numeric column count: " + name)
            rows.append(values)
    return rows


def zones(name):
    result = {}
    current = None
    for line in (ROOT / "raw" / name).read_text().splitlines():
        match = re.search(r'zone\s+t="([^"]+)"', line, re.I)
        if match:
            current = match.group(1)
            result[current] = []
        elif re.match(r"^\s*[-+0-9.]", line) and current is not None:
            result[current].append([float(value) for value in line.split()])
    return result


def wall_point(surface_cm):
    surface_m = surface_cm * 0.01
    angle = math.radians(20)
    x_source = surface_m if surface_m <= 0 else surface_m * math.cos(angle)
    radius = 0.1015 + max(surface_m, 0) * math.sin(angle)
    return surface_m, x_source, -x_source, radius


def main():
    pressure = tecplot_numeric("pw_exp.dat", 2)
    heat = tecplot_numeric("qw_exp.dat", 2)
    upstream = tecplot_numeric("upstream_exp.dat", 3)
    originals = zones("KussoyCompCorner.dat")
    inlet_original = zones("KussoyBLProfile.dat")["KussoyBL20_IncomingBL"]
    unique_p = []
    seen = set()
    for index, row in enumerate(pressure, 1):
        key = tuple(row)
        selected = key not in seen
        seen.add(key)
        s, xs, xp, radius = wall_point(row[0])
        unique_p.append(dict(source_row=index, S_m=s, source_x_m=xs,
                             project_X_m=xp, radius_m=radius,
                             ratio=row[1], pressure_Pa=row[1]*576,
                             independent_comparison_point=selected))
    q_data = []
    for index, row in enumerate(heat, 1):
        s, xs, xp, radius = wall_point(row[0])
        q_data.append(dict(source_row=index, S_m=s, source_x_m=xs,
                           project_X_m=xp, radius_m=radius, ratio=row[1],
                           reported_heating_magnitude_W_m2=row[1]*9300))
    profile = [dict(source_row=i, y_from_wall_m=row[0]*0.01,
                    velocity_ratio=row[1], total_temperature_ratio=row[2],
                    measured_M=inlet_original[i-1][1],
                    pressure_ratio=inlet_original[i-1][2])
               for i, row in enumerate(upstream, 1)]
    # Cross-check the two official distributions, retaining their actual precision.
    selected = [pressure[row["source_row"]-1] for row in unique_p
                if row["independent_comparison_point"]]
    max_p = max(max(abs(a[0]-b[0]), abs(a[1]-b[1]))
                for a, b in zip(selected, originals["KussoyP20"]))
    max_q = max(max(abs(a[0]-b[0]), abs(a[1]-b[2]))
                for a, b in zip(heat, originals["KussoyQ20"]))
    max_u = max(max(abs(a[0]-b[0]),abs(a[1]-b[3]),abs(a[2]-b[4]))
                for a,b in zip(upstream,inlet_original))
    assert len(pressure)==27 and len(originals["KussoyP20"])==24
    assert len(heat)==27 and len(profile)==21
    assert max(max_p,max_q,max_u)<1e-6
    manifest_hash = hashlib.sha256((ROOT/"source_manifest.json").read_bytes()).hexdigest()
    contract = {
        "schema": "RocketAeroHeat.tmr-mach7-acc.v1",
        "case_id": "NASA-TMR-ACC-SBLI-M7-20deg",
        "status": "source_data_and_unit_contract_only",
        "source_manifest_sha256": manifest_hash,
        "scope": {"flare_angle_deg":20,"angle_of_attack_deg":0,
                  "axisymmetric":True,"planar_2d_equivalent":False,
                  "solver_validation_complete":False,"new_physical_area_m2":0},
        "sources": {
            "case_url":"https://tmbwg.github.io/turbmodels/axiswblim7_val.html",
            "original_report_url":"https://ntrs.nasa.gov/api/citations/19890010729/downloads/19890010729.pdf",
            "original_report_local_sha256":None,
            "original_report_retrieval":"NTRS PDF read via web tool; local HTTP download returned 403; no PDF-byte hash claimed",
            "experiment_archive_url":"https://www.nasa.gov/wp-content/uploads/2025/10/kussoycompcorner.zip",
            "experiment_archive_identity":"NASA/TM-2013-216604 database, same Kussoy/Horstman experiment, original NASA-TM-101075",
            "source_locator_conditions":"raw/KussoyConditions.txt; raw/A5.0.ReadmeKussoyCompCorner.txt; NASA-TM-101075 Table I",
            "source_locator_uncertainty":"NASA-TM-101075 printed p.5, Experimental Uncertainties",
            "source_locator_measurements":"raw/*_exp.dat; raw/KussoyCompCorner.dat zones KussoyP20/KussoyQ20; NASA-TM-101075 Table IV(a)"},
        "geometry": {
            "cylinder_radius_m":0.1015,"cylinder_diameter_m":0.203,
            "flare_start_distance_from_nose_m":1.39,"nose_cone_half_angle_deg":10,
            "nose_arc_radius_m":None,"flare_end_S_m":None,
            "source_coordinates":"x increases downstream; x=S=0 at cylinder/flare junction; y is radius from axis",
            "source_surface_coordinate":"S<0: source x=S,r=R; S>=0: source x=S*cos(20deg),r=R+S*sin(20deg)",
            "project_coordinates":"X=-source_x, junction X=0; nose is at +X, free-stream velocity -X; radius unchanged",
            "profile_coordinate":"y_prime is distance from cylinder wall; source radius=y_prime+R",
            "wall_condition":{"type":"prescribed_no_slip_isothermal","temperature_K":311},
            "geometry_truth_limit":"Boundary grid gives the computational truncation, not full forebody nose-arc dimensions"},
        "experimental_local_reference": {
            "location_source_x_m":-0.06,"distance_from_nose_m":1.33,
            "definition":"Outside the boundary layer immediately upstream of interaction; inside the forebody/inflow shock layer",
            "Mach":7.05,"temperature_K":81.2,"pressure_Pa":576,
            "density_kg_m3":0.0252,"velocity_m_s":1274,
            "wall_temperature_K":311,"reported_wall_heating_magnitude_W_m2":9300,
            "wall_shear_Pa":25,"unit_Re_m_inv":5.8e6,
            "boundary_layer_thickness_m":0.025,"displacement_thickness_m":0.0074,
            "momentum_thickness_m":0.00065,"Cf":0.00122,
            "total_temperature_K":None,
            "independent_reported_rounded_quantities":True},
        "computational_truncated_inflow": {
            "identity":"TMR prescribed calculation boundary, distinct from measured profile station",
            "source_x_m":-0.7839,"project_X_m":0.7839,"Mach":7.11,
            "temperature_K_figure":80,"unit_Re_m_inv_figure":5706000,
            "total_pressure_psi_WIND_input":364.8,
            "total_pressure_Pa_WIND_input":364.8*6894.757293168,
            "total_temperature_rankine_WIND_input":1600,
            "total_temperature_K_WIND_input":1600*5/9,
            "wall_temperature_rankine_WIND_input":559.8,
            "angle_of_attack_deg":0,"sideslip_deg":0,
            "top_boundary":"supersonic_inflow on both blocks (NMF)",
            "outflow":"extrapolate from interior, NMF; WIND sample first order extrapolation",
            "azimuthal_boundary":"rotated periodic planes separated by 1 degree; alternatively consistent axisymmetric equations",
            "incoming_turbulence_profiles":None,
            "WIND_unit_source":"raw/wind_freestream.html: total pressure psi and temperature Rankine",
            "not_a_measured_complete_profile":True},
        "normalization": {
            "pressure_ratio":"TMR Pw/Pw-inf; author archive Pw/Pinf, with local Pinf=576 Pa",
            "pressure_Pa":"pressure_ratio * 576 Pa (same local reference)",
            "heat_ratio":"Qw/Qw-inf, experimental heating magnitude reference 9300 W/m^2",
            "heating_magnitude_W_m2":"heat_ratio * 9300 W/m^2",
            "temperature_profile":"TT/TTinf is total-temperature ratio, not static-temperature ratio",
            "signed_Fourier_warning":"TMR defines Qw=-k*dT/dy at wall; raw positive ratios do not determine the solver normal/sign. The dimensional field here is reported positive heating magnitude, not a signed Fourier conversion.",
            "forbidden_reference_substitution":"Do not use truncated M7.11 inlet state as the denominator of the x=-6 cm local experimental ratios"},
        "measurement_uncertainties": {
            "surface_pressure_relative_estimate":0.10,
            "surface_heat_relative_estimate":0.10,
            "velocity_relative_estimate":0.03,
            "total_temperature_relative_estimate":0.015,
            "profile_y_absolute_m":0.0002,
            "confidence_level":None,"per_point_covariance":None,
            "do_not_infer_ratio_independent_uncertainty":True},
        "measurements": {"upstream_profile":profile,"pressure_rows":unique_p,
                         "heating_rows":q_data},
        "data_quality": {
            "pressure_raw_rows":27,"pressure_independent_points":24,
            "heating_raw_rows":27,"upstream_rows":21,
            "pressure_duplicate_source_rows_1based":[25,26,27],
            "pressure_duplicate_policy":"Keep raw rows and their identity; compare only the first occurrence of an exactly repeated (S,ratio) pair. Do not give repeated tail rows extra statistical weight.",
            "author_first_heat_position_typo":"raw/KussoyAxiCompCorner.txt prints +12.06 cm; TMR qw_exp.dat, NASA database KussoyQ20 and original Table IV(a) give -12.06 cm. Contract uses the latter and retains both sources unchanged.",
            "tecplot_placeholder_warning":"Zero columns QQInf/Tau_TauInf in KussoyP20 and PPInf/Tau_TauInf in KussoyQ20 are unavailable-field placeholders, not measurements.",
            "separation_label_disagreement":"2013 database README calls 20deg attached; original report/Conditions give separation at S=0, reattachment undetermined. Do not treat either as a measured finite separation bubble or reattachment coordinate.",
            "cross_source_max_abs_native_units":{"pressure":max_p,"heat":max_q,"upstream":max_u}},
        "unresolved_before_physical_validation": [
            "Choose and freeze caloric/transport/Prandtl and turbulence-model variants; do not silently substitute current NASA7 for source gas/transport closure",
            "Generate the boundary layer from TMR M7.11 inflow and demonstrate agreement at x=-6 cm, or separately label measured-profile local-domain validation",
            "Measured full turbulence profiles, Reynolds stresses and a complete transverse velocity profile are not supplied by these files",
            "Experimental total-temperature reference is not supplied numerically in the selected source files; retain normalized TT, never substitute calculation inlet TT",
            "Convergence, near-wall resolution, grid and outer/inlet/outlet sensitivity still require an actual internal viscous solution",
            "Match wall normal and heat sign explicitly; compare source locations and normalized values without treating duplicate rows as independent observations"]}
    (ROOT/"contract.json").write_text(json.dumps(contract,indent=2,ensure_ascii=True)+"\n",encoding="utf-8")
    print("Contract built: 21 profile rows, 24 independent pressure points, 27 heating points.")


if __name__ == "__main__":
    main()
