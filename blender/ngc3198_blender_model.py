"""Build an interactive Blender model of NGC 3198 from the local SPARC data.

Run this file from Blender's Scripting workspace. The generated scene has
200 test stars, a timeline animation, and an N-panel dark-matter control.
One Blender unit represents one kiloparsec.

This file also works as a Blender add-on. Install it from Preferences > Add-ons
to restore the controls automatically whenever the saved .blend is reopened.
"""

from __future__ import annotations

import bisect
import csv
import math
from pathlib import Path
from typing import Any, TypedDict

try:
    import bpy as _bpy  # type: ignore[reportMissingImports]
    from mathutils import Euler as _Euler, Vector as _Vector  # type: ignore[reportMissingImports]
except ImportError:
    _bpy = None
    _Euler = None
    _Vector = None

bpy: Any = _bpy
Euler: Any = _Euler
Vector: Any = _Vector

bl_info = {
    "name": "NGC 3198 Dark Matter Model",
    "author": "OpenAI",
    "version": (1, 0, 0),
    "blender": (4, 0, 0),
    "location": "View3D > Sidebar > Galaxy",
    "description": "Build and interact with a SPARC-based NGC 3198 orbit model",
    "category": "3D View",
}


KM_S_TO_KPC_MYR = 0.0010227121650537077
UPSILON_DISK = 0.50
UPSILON_BULGE = 0.70
MAX_DARK_TO_BARYON = 8.0
COSMIC_DARK_TO_BARYON = 5.4
FRAME_START = 1
FRAME_END = 301
STAR_COUNT = 200
SPARC_RELATIVE_PATH = Path("data") / "sparc" / "table2.dat"
COLLECTION_NAME = "NGC 3198 Interactive Model"
STAR_OBJECT_PREFIX = "NGC3198_Star_"
TRAIL_OBJECT_PREFIX = "NGC3198_Orbit_"

class RotationProfile(TypedDict):
    radii: list[float]
    baryonic_v: list[float]
    dark_v2: list[float]
    outer_ratio: float


State = tuple[float, float, float, float]
_MODEL: dict[str, Any] | None = None
_UPDATING_PROPERTIES = False


def _project_dir() -> Path:
    candidates = []
    if bpy is not None:
        scene = bpy.context.scene
        stored_directory = scene.get("ngc3198_source_directory", "")
        configured_directory = getattr(scene, "ngc3198_data_directory", "")
        if stored_directory:
            candidates.append(Path(stored_directory).expanduser())
        if configured_directory:
            candidates.append(Path(configured_directory).expanduser())
        if bpy.data.filepath:
            candidates.append(Path(bpy.data.filepath).resolve().parent)
    if "__file__" in globals():
        candidates.append(Path(__file__).resolve().parent)
    candidates.append(Path.cwd())
    for candidate in candidates:
        # The script and .blend live in blender/, so also try the folder above.
        for folder in (candidate, candidate.parent):
            if (folder / SPARC_RELATIVE_PATH).is_file():
                return folder.resolve()
    raise FileNotFoundError(
        "Select the EPQ Rotation Curves project folder (it contains data/sparc/table2.dat and results/)."
    )


def _load_sparc_profile(path: Path) -> RotationProfile:
    """Read NGC 3198's SPARC table and form baryonic and residual speeds."""
    radii = []
    baryonic_speeds = []
    observed_speeds = []
    with path.open("r", encoding="utf-8", errors="replace") as source:
        for line in source:
            columns = line.split()
            if len(columns) < 10 or columns[0].upper() != "NGC3198":
                continue
            try:
                radius, observed, gas, disk, bulge = map(
                    float, (columns[2], columns[3], columns[5], columns[6], columns[7])
                )
            except ValueError:
                continue
            if radius <= 0.0:
                continue
            baryonic_squared = (
                gas * abs(gas) + UPSILON_DISK * disk * disk + UPSILON_BULGE * bulge * bulge
            )
            radii.append(radius)
            baryonic_speeds.append(math.sqrt(max(baryonic_squared, 0.0)))
            observed_speeds.append(abs(observed))

    if len(radii) < 3:
        raise ValueError(f"Could not read NGC 3198 SPARC rows from {path}")
    rows = sorted(zip(radii, baryonic_speeds, observed_speeds))
    radii, baryonic_speeds, observed_speeds = map(list, zip(*rows))
    if any(right <= left for left, right in zip(radii, radii[1:])):
        raise ValueError("SPARC radii must be strictly increasing for NGC 3198")

    dark_squared = [
        max(observed * observed - baryonic * baryonic, 0.0)
        for observed, baryonic in zip(observed_speeds, baryonic_speeds)
    ]
    outer_ratio = dark_squared[-1] / max(baryonic_speeds[-1] ** 2, 1e-12)
    outer_ratio = max(outer_ratio, 1e-6)
    return {
        "radii": radii,
        "baryonic_v": baryonic_speeds,
        "dark_v2": dark_squared,
        "outer_ratio": outer_ratio,
    }


def _load_baryonic_trajectories(path: Path) -> dict[int, list[tuple[float, float, float]]]:
    trajectories: dict[int, list[tuple[float, float, float]]] = {}
    with path.open("r", newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            star_id = int(row["star_id"])
            trajectories.setdefault(star_id, []).append(
                (float(row["x_kpc"]), float(row["y_kpc"]), float(row["z_kpc"]))
            )
    if len(trajectories) < STAR_COUNT:
        raise ValueError(f"Expected at least {STAR_COUNT} stars in {path}")
    for positions in trajectories.values():
        if len(positions) != FRAME_END:
            raise ValueError(f"Each trajectory must contain {FRAME_END} frames: {path}")
    return trajectories


def _interpolate(radii: list[float], values: list[float], radius: float) -> float:
    """Piecewise-linear interpolation with physical inner/outer extensions."""
    if radius <= radii[0]:
        return values[0] * radius / radii[0]
    if radius >= radii[-1]:
        return values[-1]
    index = bisect.bisect_right(radii, radius) - 1
    fraction = (radius - radii[index]) / (radii[index + 1] - radii[index])
    return values[index] + fraction * (values[index + 1] - values[index])


def _baryonic_speed(profile: RotationProfile, radius: float) -> float:
    if radius > profile["radii"][-1]:
        return profile["baryonic_v"][-1] * math.sqrt(profile["radii"][-1] / radius)
    return _interpolate(profile["radii"], profile["baryonic_v"], radius)


def _dark_speed_squared(profile: RotationProfile, radius: float) -> float:
    radii = profile["radii"]
    values = profile["dark_v2"]
    if radius <= radii[0]:
        return values[0] * (radius / radii[0]) ** 2
    return _interpolate(radii, values, radius)


def _accelerations(
    x: float, y: float, profile: RotationProfile, dark_ratio: float
) -> tuple[float, float]:
    radius = max(math.hypot(x, y), 1e-5)
    baryonic_v = _baryonic_speed(profile, radius)
    dark_v2 = _dark_speed_squared(profile, radius)
    dark_scale = dark_ratio / profile["outer_ratio"]
    v2_kpc_myr2 = (baryonic_v**2 + dark_scale * dark_v2) * KM_S_TO_KPC_MYR**2
    acceleration = v2_kpc_myr2 / radius
    return -acceleration * x / radius, -acceleration * y / radius


def _derivative(
    state: State,
    profile: RotationProfile,
    dark_ratio: float,
) -> State:
    x, y, vx, vy = state
    ax, ay = _accelerations(x, y, profile, dark_ratio)
    return vx, vy, ax, ay


def _rk4_step(
    state: State,
    profile: RotationProfile,
    dark_ratio: float,
    step_myr: float,
) -> State:
    k1 = _derivative(state, profile, dark_ratio)
    s2: State = (
        state[0] + step_myr * k1[0] / 2.0,
        state[1] + step_myr * k1[1] / 2.0,
        state[2] + step_myr * k1[2] / 2.0,
        state[3] + step_myr * k1[3] / 2.0,
    )
    k2 = _derivative(s2, profile, dark_ratio)
    s3: State = (
        state[0] + step_myr * k2[0] / 2.0,
        state[1] + step_myr * k2[1] / 2.0,
        state[2] + step_myr * k2[2] / 2.0,
        state[3] + step_myr * k2[3] / 2.0,
    )
    k3 = _derivative(s3, profile, dark_ratio)
    s4: State = (
        state[0] + step_myr * k3[0],
        state[1] + step_myr * k3[1],
        state[2] + step_myr * k3[2],
        state[3] + step_myr * k3[3],
    )
    k4 = _derivative(s4, profile, dark_ratio)
    return (
        state[0] + step_myr * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0,
        state[1] + step_myr * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0,
        state[2] + step_myr * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6.0,
        state[3] + step_myr * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) / 6.0,
    )


def _integrate_stars(
    profile: RotationProfile,
    initial_conditions: list[tuple[float, float, float]],
    dark_ratio: float,
) -> dict[int, list[tuple[float, float, float]]]:
    """Integrate massless disk tracers in the fixed, axisymmetric mid-plane."""
    trajectories = {}
    step_myr = 1.0
    for star_id, (x, y, initial_speed_kms) in enumerate(initial_conditions):
        radius = max(math.hypot(x, y), 1e-8)
        speed = initial_speed_kms * KM_S_TO_KPC_MYR
        state = (x, y, -speed * y / radius, speed * x / radius)
        positions = [(x, y, 0.0)]
        for _ in range(FRAME_END - 1):
            state = _rk4_step(state, profile, dark_ratio, step_myr)
            positions.append((state[0], state[1], 0.0))
        trajectories[star_id] = positions
    return trajectories


def _initial_conditions(path: Path) -> list[tuple[float, float, float]]:
    by_id = {}
    with path.open("r", newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            if int(row["frame"]) == 1:
                by_id[int(row["star_id"])] = (
                    float(row["x_kpc"]),
                    float(row["y_kpc"]),
                    float(row["v_obs_kms"]),
                )
    if len(by_id) < STAR_COUNT:
        raise ValueError(f"Expected initial conditions for {STAR_COUNT} stars in {path}")
    return [by_id[index] for index in sorted(by_id)[:STAR_COUNT]]


def _material(name: str, color: tuple[float, float, float, float], emission: float = 0.0):
    material = bpy.data.materials.new(name)
    material.diffuse_color = color
    material.use_nodes = True
    shader = material.node_tree.nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = color
    emission_color = shader.inputs.get("Emission Color") or shader.inputs.get("Emission")
    if emission_color:
        emission_color.default_value = color
    shader.inputs["Emission Strength"].default_value = emission
    alpha_socket = shader.inputs.get("Alpha")
    if alpha_socket and color[3] < 1.0:
        alpha_socket.default_value = color[3]
    return material


def _new_collection(scene):
    old_collection = bpy.data.collections.get(COLLECTION_NAME)
    if old_collection:
        for obj in list(old_collection.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(old_collection)
    collection = bpy.data.collections.new(COLLECTION_NAME)
    scene.collection.children.link(collection)
    return collection


def _velocity_color_group():
    group_name = "NGC3198_Velocity_Color_Ramp"
    group = bpy.data.node_groups.get(group_name)
    if group is None:
        group = bpy.data.node_groups.new(group_name, "ShaderNodeTree")
    for socket in list(group.interface.items_tree):
        if getattr(socket, "item_type", None) == "SOCKET":
            group.interface.remove(socket)
    group.interface.new_socket(name="Vrot", in_out="INPUT", socket_type="NodeSocketFloat")
    group.interface.new_socket(
        name="Emission Strength", in_out="INPUT", socket_type="NodeSocketFloat"
    )
    group.interface.new_socket(name="Shader", in_out="OUTPUT", socket_type="NodeSocketShader")
    nodes = group.nodes
    nodes.clear()
    input_node = nodes.new("NodeGroupInput")
    ramp = nodes.new("ShaderNodeValToRGB")
    ramp.color_ramp.elements[0].position = 0.0
    ramp.color_ramp.elements[0].color = (0.02, 0.38, 1.0, 1.0)
    ramp.color_ramp.elements[1].position = 1.0
    ramp.color_ramp.elements[1].color = (1.0, 0.025, 0.01, 1.0)
    middle = ramp.color_ramp.elements.new(0.5)
    middle.color = (1.0, 0.82, 0.06, 1.0)
    emission = nodes.new("ShaderNodeEmission")
    output_node = nodes.new("NodeGroupOutput")
    group.links.new(input_node.outputs["Vrot"], ramp.inputs["Fac"])
    group.links.new(input_node.outputs["Emission Strength"], emission.inputs["Strength"])
    group.links.new(ramp.outputs["Color"], emission.inputs["Color"])
    group.links.new(emission.outputs["Emission"], output_node.inputs["Shader"])
    return group


def _velocity_color(value: float) -> tuple[float, float, float, float]:
    stops = (
        (0.0, (0.02, 0.38, 1.0)),
        (0.5, (1.0, 0.82, 0.06)),
        (1.0, (1.0, 0.025, 0.01)),
    )
    value = max(0.0, min(1.0, value))
    for (start, start_color), (end, end_color) in zip(stops, stops[1:]):
        if value <= end:
            factor = (value - start) / (end - start)
            red, green, blue = (
                start_color[index] + factor * (end_color[index] - start_color[index])
                for index in range(3)
            )
            return red, green, blue, 1.0
    return (*stops[-1][1], 1.0)


def _velocity_material(
    star_id: int,
    normalized_velocity: float,
    emission_strength: float,
    color_group,
):
    material = bpy.data.materials.new(f"NGC3198_StarMaterial_{star_id:03d}")
    material.diffuse_color = _velocity_color(normalized_velocity)
    material.use_nodes = True
    nodes = material.node_tree.nodes
    nodes.clear()
    group_node = nodes.new("ShaderNodeGroup")
    group_node.node_tree = color_group
    group_node.inputs["Vrot"].default_value = normalized_velocity
    group_node.inputs["Emission Strength"].default_value = emission_strength
    output_node = nodes.new("ShaderNodeOutputMaterial")
    material.node_tree.links.new(group_node.outputs["Shader"], output_node.inputs["Surface"])
    return material


def _make_star_mesh():
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=1, radius=0.065)
    mesh = bpy.context.object.data
    mesh.name = "NGC3198_Shared_Star_Mesh"
    bpy.data.objects.remove(bpy.context.object, do_unlink=True)
    return mesh


def _make_halo_material():
    material = bpy.data.materials.new("NGC3198_Dark_Matter_Halo")
    halo_color = (0.28, 0.2, 0.72, 1.0)
    material.diffuse_color = (*halo_color[:3], 0.08)
    material.use_nodes = True
    nodes = material.node_tree.nodes
    nodes.clear()
    output_node = nodes.new("ShaderNodeOutputMaterial")
    transparent = nodes.new("ShaderNodeBsdfTransparent")
    emission = nodes.new("ShaderNodeEmission")
    emission.inputs["Color"].default_value = halo_color
    emission.inputs["Strength"].default_value = 0.45
    mix = nodes.new("ShaderNodeMixShader")
    mix.inputs["Fac"].default_value = 0.08
    material.node_tree.links.new(transparent.outputs["BSDF"], mix.inputs[1])
    material.node_tree.links.new(emission.outputs["Emission"], mix.inputs[2])
    material.node_tree.links.new(mix.outputs["Shader"], output_node.inputs["Surface"])
    if hasattr(material, "surface_render_method"):
        material.surface_render_method = "BLENDED"
    elif hasattr(material, "blend_method"):
        material.blend_method = "BLEND"
    return material


def _configure_world_camera_and_glow(scene, collection):
    scene.render.engine = "BLENDER_EEVEE"
    scene.render.resolution_x = 1440
    scene.render.resolution_y = 1080
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"

    world = scene.world or bpy.data.worlds.new("NGC3198_Deep_Space")
    scene.world = world
    world.use_nodes = True
    background = world.node_tree.nodes.get("Background")
    background.inputs["Color"].default_value = (0.0, 0.0, 0.0, 1.0)
    background.inputs["Strength"].default_value = 1.0

    scene.use_nodes = True
    if hasattr(scene, "compositing_node_group"):
        compositor = scene.compositing_node_group
        if compositor is None:
            compositor = bpy.data.node_groups.new(
                "NGC3198_Deep_Space_Compositor", "CompositorNodeTree"
            )
        for socket in list(compositor.interface.items_tree):
            if getattr(socket, "item_type", None) == "SOCKET":
                compositor.interface.remove(socket)
        compositor.interface.new_socket(
            name="Image", in_out="INPUT", socket_type="NodeSocketColor"
        )
        compositor.interface.new_socket(
            name="Image", in_out="OUTPUT", socket_type="NodeSocketColor"
        )
        compositor.nodes.clear()
        render_layers = compositor.nodes.new("CompositorNodeRLayers")
        glow = compositor.nodes.new("CompositorNodeGlare")
        glow.inputs["Type"].default_value = "Fog Glow"
        glow.inputs["Threshold"].default_value = 1.0
        glow.inputs["Strength"].default_value = 0.8
        glow.inputs["Size"].default_value = 0.85
        group_output = compositor.nodes.new("NodeGroupOutput")
        compositor.links.new(render_layers.outputs["Image"], glow.inputs["Image"])
        compositor.links.new(
            glow.outputs["Image"],
            group_output.inputs["Image"],
        )
        scene.compositing_node_group = compositor
    else:
        compositor = scene.node_tree
        compositor.nodes.clear()
        render_layers = compositor.nodes.new("CompositorNodeRLayers")
        glow = compositor.nodes.new("CompositorNodeGlare")
        glow.glare_type = "FOG_GLOW"
        glow.threshold = 1.0
        glow.size = 8
        composite = compositor.nodes.new("CompositorNodeComposite")
        compositor.links.new(render_layers.outputs["Image"], glow.inputs["Image"])
        compositor.links.new(glow.outputs["Image"], composite.inputs["Image"])

    camera_data = bpy.data.cameras.new("NGC3198_Render_Camera")
    camera = bpy.data.objects.new("NGC3198_Render_Camera", camera_data)
    collection.objects.link(camera)
    camera.location = (0.0, -82.0, 57.4)
    target_direction = Vector((0.0, 0.0, 0.0)) - camera.location
    camera.rotation_euler = target_direction.to_track_quat("-Z", "Y").to_euler()
    camera_data.lens = 35.0
    camera_data.clip_start = 0.1
    camera_data.clip_end = 500.0
    scene.camera = camera

    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type != "VIEW_3D":
                continue
            space = area.spaces.active
            space.show_region_ui = True
            space.shading.type = "MATERIAL"
            space.shading.use_scene_world = True
            space.shading.use_scene_lights = True
            if hasattr(space.shading, "use_compositor"):
                space.shading.use_compositor = "ALWAYS"
            space.overlay.show_floor = False
            space.overlay.show_axis_x = False
            space.overlay.show_axis_y = False
            space.region_3d.view_location = (0.0, 0.0, 0.0)
            space.region_3d.view_rotation = camera.rotation_euler.to_quaternion()
            space.region_3d.view_distance = 92.0


def _make_orbit_curve(collection, star_id: int, positions, material):
    curve = bpy.data.curves.new(f"NGC3198_OrbitData_{star_id:03d}", "CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 1
    curve.bevel_depth = 0.008
    curve.bevel_resolution = 0
    spline = curve.splines.new("POLY")
    spline.points.add(len(positions) - 1)
    for point, position in zip(spline.points, positions):
        point.co = (*position, 1.0)
    orbit = bpy.data.objects.new(f"{TRAIL_OBJECT_PREFIX}{star_id:03d}", curve)
    collection.objects.link(orbit)
    curve.materials.append(material)
    return orbit


def _create_scene(trajectories, profile, initial_conditions):
    scene = bpy.context.scene
    default_cube = bpy.data.objects.get("Cube")
    if default_cube is not None:
        bpy.data.objects.remove(default_cube, do_unlink=True)
    collection = _new_collection(scene)
    orbit_material = _material("Orbit traces", (0.12, 0.2, 0.28, 0.08), 0.008)
    if hasattr(orbit_material, "surface_render_method"):
        orbit_material.surface_render_method = "BLENDED"
    elif hasattr(orbit_material, "blend_method"):
        orbit_material.blend_method = "BLEND"
    core_material = _material("Galactic Core Glow", (1.0, 0.86, 0.52, 1.0), 1.4)
    guide_material = _material("Scale guides", (0.16, 0.23, 0.28, 1.0), 0.0)
    halo_material = _make_halo_material()
    star_mesh = _make_star_mesh()
    color_group = _velocity_color_group()
    velocity_values = [condition[2] for condition in initial_conditions]
    minimum_velocity = min(velocity_values)
    velocity_span = max(max(velocity_values) - minimum_velocity, 1e-8)

    bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=20, radius=0.62)
    core = bpy.context.object
    core.name = "Galactic Core"
    core.data.materials.append(core_material)
    for polygon in core.data.polygons:
        polygon.use_smooth = True
    for linked_collection in list(core.users_collection):
        linked_collection.objects.unlink(core)
    collection.objects.link(core)

    bpy.ops.mesh.primitive_uv_sphere_add(segments=20, ring_count=10, radius=40.0)
    halo = bpy.context.object
    halo.name = "NGC3198_Dark_Matter_Halo_Shell"
    halo.data.materials.append(halo_material)
    wireframe = halo.modifiers.new("Translucent Halo Wireframe", "WIREFRAME")
    wireframe.thickness = 0.045
    wireframe.use_even_offset = True
    halo.show_in_front = False
    for linked_collection in list(halo.users_collection):
        linked_collection.objects.unlink(halo)
    collection.objects.link(halo)

    star_objects = []
    trail_objects = []
    for star_id in range(STAR_COUNT):
        star = bpy.data.objects.new(f"{STAR_OBJECT_PREFIX}{star_id:03d}", star_mesh)
        star.location = trajectories[star_id][0]
        speed_kms = velocity_values[star_id]
        normalized_velocity = (speed_kms - minimum_velocity) / velocity_span
        initial_radius = math.hypot(initial_conditions[star_id][0], initial_conditions[star_id][1])
        normalized_radius = min(initial_radius / profile["radii"][-1], 1.0)
        emission_strength = 2.0 + 3.0 * (1.0 - normalized_radius) ** 2
        star_material = _velocity_material(
            star_id, normalized_velocity, emission_strength, color_group
        )
        if star_id == 0:
            star_mesh.materials.append(star_material)
        star.material_slots[0].link = "OBJECT"
        star.material_slots[0].material = star_material
        star["v_rot_kms"] = speed_kms
        collection.objects.link(star)
        star_objects.append(star)
        trail_objects.append(
            _make_orbit_curve(collection, star_id, trajectories[star_id], orbit_material)
        )

    outer_radius = profile["radii"][-1]
    for guide_radius in (10.0, 20.0, 30.0, 40.0, outer_radius):
        if guide_radius > outer_radius:
            continue
        bpy.ops.curve.primitive_bezier_circle_add(radius=guide_radius)
        guide = bpy.context.object
        guide.name = f"Radius guide {guide_radius:.0f} kpc"
        guide.data.dimensions = "3D"
        guide.data.resolution_u = 1
        guide.data.bevel_depth = 0.006
        guide.data.bevel_resolution = 0
        guide.data.materials.append(guide_material)
        for linked_collection in list(guide.users_collection):
            linked_collection.objects.unlink(guide)
        collection.objects.link(guide)

    scene.frame_start = FRAME_START
    scene.frame_end = FRAME_END
    scene.frame_set(FRAME_START)
    scene.render.fps = 30
    scene.unit_settings.system = "NONE"
    _configure_world_camera_and_glow(scene, collection)
    scene["galaxy"] = "NGC 3198"
    scene["model_note"] = (
        "Fixed Newtonian axisymmetric mid-plane potential; 200 non-interacting "
        "test particles; no mutual stellar gravity or vertical dynamics. SPARC "
        "2016 rotation curve; M/L_disk=0.50 and M/L_bulge=0.70. The slider is "
        "a spherical-equivalent effective dark-to-baryonic ratio at the outer "
        "SPARC radius, not a unique mass decomposition for a disk galaxy."
    )
    return scene, star_objects, trail_objects


def _write_positions(trajectories, star_objects, trail_objects, frame_index: int) -> None:
    for star_id, (star, orbit) in enumerate(zip(star_objects, trail_objects)):
        positions = trajectories[star_id]
        star.location = positions[frame_index]
        spline = orbit.data.splines[0]
        for point, position in zip(spline.points, positions):
            point.co = (*position, 1.0)


def _ratio_update(self, context):
    global _UPDATING_PROPERTIES
    if _UPDATING_PROPERTIES or _MODEL is None:
        return
    ratio = max(0.0, min(float(self.dark_to_baryon), MAX_DARK_TO_BARYON))
    _UPDATING_PROPERTIES = True
    self.use_cosmic_mean = abs(ratio - COSMIC_DARK_TO_BARYON) < 0.05
    self["ngc3198_dark_to_baryon"] = ratio
    self["ngc3198_use_cosmic_mean"] = self.use_cosmic_mean
    _UPDATING_PROPERTIES = False
    trajectories = (
        _MODEL["baryonic_trajectories"]
        if ratio < 1e-6
        else _integrate_stars(_MODEL["profile"], _MODEL["initial_conditions"], ratio)
    )
    _MODEL["trajectories"] = trajectories
    _write_positions(
        trajectories,
        _MODEL["stars"],
        _MODEL["trails"],
        max(0, min(context.scene.frame_current - FRAME_START, FRAME_END - 1)),
    )


def _preset_update(self, context):
    global _UPDATING_PROPERTIES
    if _UPDATING_PROPERTIES or _MODEL is None:
        return
    if not self.use_cosmic_mean and abs(self.dark_to_baryon - COSMIC_DARK_TO_BARYON) < 0.05:
        self["ngc3198_use_cosmic_mean"] = False
        return
    if self.use_cosmic_mean:
        _UPDATING_PROPERTIES = True
        self.dark_to_baryon = COSMIC_DARK_TO_BARYON
        _UPDATING_PROPERTIES = False
    _ratio_update(self, context)


def _frame_change(scene, depsgraph=None):
    if _MODEL is None:
        return
    frame_index = max(0, min(scene.frame_current - FRAME_START, FRAME_END - 1))
    for star_id, star in enumerate(_MODEL["stars"]):
        star.location = _MODEL["trajectories"][star_id][frame_index]
        star.update_tag()


if bpy is not None:
    _frame_change = bpy.app.handlers.persistent(_frame_change)


def _draw_panel(self, context):
    layout = self.layout
    scene = context.scene
    if _MODEL is None:
        column = layout.column(align=True)
        column.prop(scene, "ngc3198_data_directory", text="Data folder")
        column.operator("ngc3198.build_model", text="Build NGC 3198 model", icon="FILE_NEW")
        return
    column = layout.column(align=True)
    column.label(text="Effective dark / baryonic ratio")
    column.prop(scene, "dark_to_baryon", slider=True, text="DM : baryons")
    column.prop(scene, "use_cosmic_mean", text="Cosmic mean (about 5.4:1)")
    column.operator("ngc3198.set_sparc_ratio", text="Set NGC 3198 SPARC value")
    column.label(text=f"SPARC outer-radius effective fit: {_MODEL['recommended_ratio']:.2f}:1")
    column.separator()
    row = column.row(align=True)
    row.operator("screen.animation_play", text="Play", icon="PLAY")
    row.operator("screen.animation_play", text="Reverse", icon="PLAY_REVERSE").reverse = True
    column.label(text="Timeline: 1 Myr per frame")
    column.label(text="Space: kpc | 2D mid-plane dynamics")


if bpy is not None:
    class NGC3198_PT_dark_matter(bpy.types.Panel):
        bl_label = "NGC 3198 Dynamics"
        bl_idname = "NGC3198_PT_dark_matter"
        bl_space_type = "VIEW_3D"
        bl_region_type = "UI"
        bl_category = "Galaxy"

        draw = _draw_panel


    class NGC3198_OT_sparc_ratio(bpy.types.Operator):
        bl_idname = "ngc3198.set_sparc_ratio"
        bl_label = "Set NGC 3198 SPARC Ratio"
        bl_description = "Set the dark-to-baryonic ratio inferred from NGC 3198's outer SPARC rotation point"

        def execute(self, context):
            if _MODEL is None:
                self.report({"ERROR"}, "Run the NGC 3198 model script first")
                return {"CANCELLED"}
            context.scene.dark_to_baryon = _MODEL["recommended_ratio"]
            return {"FINISHED"}


    class NGC3198_OT_build_model(bpy.types.Operator):
        bl_idname = "ngc3198.build_model"
        bl_label = "Build NGC 3198 Model"
        bl_description = "Read the local SPARC and trajectory files and build the interactive galaxy scene"

        def execute(self, context):
            try:
                build_model()
            except Exception as error:
                self.report({"ERROR"}, str(error))
                return {"CANCELLED"}
            return {"FINISHED"}


if bpy is not None:
    @bpy.app.handlers.persistent
    def _restore_model_after_load(_unused):
        global _MODEL, _UPDATING_PROPERTIES
        scene = bpy.context.scene
        if scene.get("galaxy") != "NGC 3198":
            _MODEL = None
            return
        try:
            base_dir = Path(scene.get("ngc3198_source_directory", ""))
            profile = _load_sparc_profile(base_dir / SPARC_RELATIVE_PATH)
            trajectory_path = base_dir / "results" / "star_trajectories_baryonic.csv"
            baryonic_trajectories = _load_baryonic_trajectories(trajectory_path)
            initial_conditions = _initial_conditions(trajectory_path)
            collection = bpy.data.collections.get(COLLECTION_NAME)
            if collection is None:
                raise ValueError("The NGC 3198 model collection is missing from this scene")
            stars_by_id = {
                int(obj.name.removeprefix(STAR_OBJECT_PREFIX)): obj
                for obj in collection.objects
                if obj.name.startswith(STAR_OBJECT_PREFIX)
            }
            trails_by_id = {
                int(obj.name.removeprefix(TRAIL_OBJECT_PREFIX)): obj
                for obj in collection.objects
                if obj.name.startswith(TRAIL_OBJECT_PREFIX)
            }
            if len(stars_by_id) < STAR_COUNT or len(trails_by_id) < STAR_COUNT:
                raise ValueError("The saved scene does not contain all 200 stars and orbit trails")
            _MODEL = {
                "profile": profile,
                "initial_conditions": initial_conditions,
                "baryonic_trajectories": baryonic_trajectories,
                "trajectories": baryonic_trajectories,
                "stars": [stars_by_id[index] for index in range(STAR_COUNT)],
                "trails": [trails_by_id[index] for index in range(STAR_COUNT)],
                "recommended_ratio": profile["outer_ratio"],
            }
            ratio = float(scene.get("ngc3198_dark_to_baryon", 0.0))
            cosmic_mean = bool(scene.get("ngc3198_use_cosmic_mean", False))
            _UPDATING_PROPERTIES = True
            scene.dark_to_baryon = ratio
            scene.use_cosmic_mean = cosmic_mean
            _UPDATING_PROPERTIES = False
            trajectories = (
                baryonic_trajectories
                if ratio < 1e-6
                else _integrate_stars(profile, initial_conditions, ratio)
            )
            _MODEL["trajectories"] = trajectories
            _write_positions(
                trajectories,
                _MODEL["stars"],
                _MODEL["trails"],
                max(0, min(scene.frame_current - FRAME_START, FRAME_END - 1)),
            )
        except Exception as error:
            _MODEL = None
            print(f"NGC 3198 add-on could not restore model controls: {error}")


def _register_runtime() -> None:
    if bpy is None:
        raise RuntimeError("Run this add-on from Blender.")
    scene_type = bpy.types.Scene
    if not hasattr(scene_type, "ngc3198_data_directory"):
        scene_type.ngc3198_data_directory = bpy.props.StringProperty(
            name="NGC 3198 data folder",
            description="Project folder containing data/sparc/table2.dat and results/star_trajectories_baryonic.csv",
            subtype="DIR_PATH",
            default=str(Path(__file__).resolve().parent),
        )
    if not hasattr(scene_type, "dark_to_baryon"):
        scene_type.dark_to_baryon = bpy.props.FloatProperty(
            name="Dark-to-baryon ratio",
            description=(
                "Spherical-equivalent effective dark/baryonic ratio at the last "
                "SPARC radius; scales the halo circular-speed-squared profile"
            ),
            min=0.0,
            max=MAX_DARK_TO_BARYON,
            default=0.0,
            precision=2,
            update=_ratio_update,
        )
    if not hasattr(scene_type, "use_cosmic_mean"):
        scene_type.use_cosmic_mean = bpy.props.BoolProperty(
            name="Cosmic mean dark-to-baryon ratio",
            description="Set the ratio to the cosmological mean of about 5.4 dark-matter to 1 baryonic mass",
            default=False,
            update=_preset_update,
        )
    for operator_class in (
        NGC3198_PT_dark_matter,
        NGC3198_OT_sparc_ratio,
        NGC3198_OT_build_model,
    ):
        if not hasattr(bpy.types, operator_class.__name__):
            bpy.utils.register_class(operator_class)
    for handlers, callback in (
        (bpy.app.handlers.frame_change_post, _frame_change),
        (bpy.app.handlers.load_post, _restore_model_after_load),
    ):
        if callback not in handlers:
            handlers.append(callback)


def _restore_model_deferred():
    _restore_model_after_load(None)
    return None


def register() -> None:
    _register_runtime()
    if not bpy.app.timers.is_registered(_restore_model_deferred):
        bpy.app.timers.register(_restore_model_deferred, first_interval=0.1)


def unregister() -> None:
    if bpy is None:
        return
    for handlers, callback in (
        (bpy.app.handlers.frame_change_post, _frame_change),
        (bpy.app.handlers.load_post, _restore_model_after_load),
    ):
        if callback in handlers:
            handlers.remove(callback)
    for class_name in (
        "NGC3198_OT_build_model",
        "NGC3198_OT_sparc_ratio",
        "NGC3198_PT_dark_matter",
    ):
        registered_class = getattr(bpy.types, class_name, None)
        if registered_class:
            bpy.utils.unregister_class(registered_class)
    for property_name in ("use_cosmic_mean", "dark_to_baryon", "ngc3198_data_directory"):
        if hasattr(bpy.types.Scene, property_name):
            delattr(bpy.types.Scene, property_name)


def build_model() -> None:
    global _MODEL
    _register_runtime()
    base_dir = _project_dir()
    sparc_path = base_dir / SPARC_RELATIVE_PATH
    trajectory_path = base_dir / "results" / "star_trajectories_baryonic.csv"
    profile = _load_sparc_profile(sparc_path)
    baryonic_trajectories = _load_baryonic_trajectories(trajectory_path)
    initial_conditions = _initial_conditions(trajectory_path)
    scene, stars, trails = _create_scene(baryonic_trajectories, profile, initial_conditions)
    scene["ngc3198_source_directory"] = str(base_dir)
    scene["ngc3198_dark_to_baryon"] = 0.0
    scene["ngc3198_use_cosmic_mean"] = False
    _MODEL = {
        "profile": profile,
        "initial_conditions": initial_conditions,
        "baryonic_trajectories": baryonic_trajectories,
        "trajectories": baryonic_trajectories,
        "stars": stars,
        "trails": trails,
        "recommended_ratio": profile["outer_ratio"],
    }

    scene.dark_to_baryon = 0.0
    scene.use_cosmic_mean = False
    scene.frame_set(FRAME_START)
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == "VIEW_3D":
                area.spaces.active.show_region_ui = True
                area.spaces.active.region_3d.view_location = (0.0, 0.0, 0.0)
                area.spaces.active.region_3d.view_rotation = Euler(
                    (math.radians(58.0), 0.0, math.radians(32.0)), "XYZ"
                ).to_quaternion()
                area.spaces.active.region_3d.view_distance = 105.0

    output_path = base_dir / "blender" / "ngc3198_dark_matter_interactive.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(output_path))
    print(f"Built NGC 3198 scene: {output_path}")
    print(
        "SPARC-derived enclosed dark-to-baryonic mass ratio at "
        f"{profile['radii'][-1]:.2f} kpc: {profile['outer_ratio']:.2f}:1"
    )


if __name__ == "__main__":
    if bpy is None:
        raise RuntimeError("Run this script from Blender's Scripting workspace.")
    build_model()