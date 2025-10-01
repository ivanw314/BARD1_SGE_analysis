import pandas as pd
from pymol import cmd

# === Parameters ===
pdb_file = "8G4L_subset.cif"
mybpc3_chain = "ao"
c10_start = 1164
output_pml = "c10_visualization.pml"
output_pse = "c10_visualization.pse"
score_threshold = 0.4259  # Threshold for counting variants
base_sphere_size = 0.5  # Starting size for sphere
size_increment = 0.1  # Size increase per variant below threshold

# Load your score data
score_df = pd.read_csv('df_plot_SASA.csv')

classification_df = pd.read_csv('position_classifications.csv')  # Replace with your actual file


# Calculate sphere sizes based on variants below threshold at each position
variants_below_threshold = {}
if 'average_score_norm' in score_df.columns and 'aa_pos' in score_df.columns:
    below_threshold = score_df[score_df['average_score_norm'] < score_threshold].copy()
    below_threshold['pos'] = below_threshold['aa_pos'].astype(str)
    
    for pos, group in below_threshold.groupby('pos'):
        variants_below_threshold[str(int(float(pos)))] = len(group)
    
    print(f"Found {len(variants_below_threshold)} positions with variants below threshold {score_threshold}")

# Add sphere size to classification DataFrame
classification_df['sphere_size'] = classification_df['resi'].apply(
    lambda pos: base_sphere_size + (variants_below_threshold.get(str(pos), 1) - 1) * size_increment
)

# Only keep positions that have variants below threshold
classification_df = classification_df[classification_df['resi'].astype(str).isin(variants_below_threshold.keys())]

print(f"Visualizing {len(classification_df)} positions with spheres")

# Initialize PyMOL
cmd.reinitialize()
cmd.load(pdb_file, "structure")

# Create objects
cmd.create("chain_ao_only", f"structure and chain {mybpc3_chain}")
cmd.create("myh7_surface", f"structure and not chain {mybpc3_chain}")
cmd.hide("everything")

# Show MYBPC3 C10 domain cartoon
c10_sel = f"chain_ao_only and resi {c10_start}-"
cmd.show("cartoon", c10_sel)
cmd.set("cartoon_fancy_helices", 1)
cmd.color("neon", f"{c10_sel} and name N+C+O+CA")

# Show MYH7 surface
cmd.show("surface", "myh7_surface")
cmd.color("gray90", "myh7_surface")
cmd.set("transparency", 0.3, "myh7_surface")

# Add spheres for each position in DataFrame
for idx, row in classification_df.iterrows():
    pos = int(row['resi'])
    color = row['color']
    size = row['sphere_size']
    
    # Check if glycine (use CA) or other (use CB)
    stored.resn_list = []
    cmd.iterate(f"chain_ao_only and resi {pos} and name CA", "stored.resn_list.append(resn)")
    
    if stored.resn_list and stored.resn_list[0] == 'GLY':
        atom_sel = f"chain_ao_only and resi {pos} and name CA"
    else:
        atom_sel = f"chain_ao_only and resi {pos} and name CB"
    
    # Show sphere with specified color and size
    cmd.show("spheres", atom_sel)
    cmd.color(color, atom_sel)
    cmd.set("sphere_scale", size, atom_sel)
    cmd.set("sphere_transparency", 0.0, atom_sel)

# Zoom and styling
cmd.zoom(c10_sel)
cmd.bg_color("white")
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.3)

# Generate PML script
pml_commands = f"""# PyMOL visualization script
# Sphere size = {base_sphere_size} + {size_increment} * (variants below {score_threshold} - 1)
load {pdb_file}
create chain_ao_only, chain {mybpc3_chain}
create myh7_surface, not chain {mybpc3_chain}
hide everything

# Show MYBPC3 C10 domain cartoon
show cartoon, chain_ao_only and resi {c10_start}-
color neon, chain_ao_only and name N+C+O+CA
set cartoon_fancy_helices, 1

# Show MYH7 surface
show surface, myh7_surface
color gray90, myh7_surface
set transparency, 0.3, myh7_surface

"""

# Add sphere commands from DataFrame
for idx, row in classification_df.iterrows():
    pos = int(row['resi'])
    color = row['color']
    size = row['sphere_size']
    classification = row.get('classification', 'unknown')
    n_variants = variants_below_threshold.get(str(pos), 1)
    
    # Determine atom name
    stored.resn_list = []
    cmd.iterate(f"chain_ao_only and resi {pos} and name CA", "stored.resn_list.append(resn)")
    atom_name = 'CA' if (stored.resn_list and stored.resn_list[0] == 'GLY') else 'CB'
    
    pml_commands += f"# Position {pos}: {classification} ({n_variants} variants)\n"
    pml_commands += f"show spheres, chain_ao_only and resi {pos} and name {atom_name}\n"
    pml_commands += f"color {color}, chain_ao_only and resi {pos} and name {atom_name}\n"
    pml_commands += f"set sphere_scale, {size}, chain_ao_only and resi {pos} and name {atom_name}\n"
    pml_commands += f"set sphere_transparency, 0.0, chain_ao_only and resi {pos} and name {atom_name}\n\n"

pml_commands += f"""
zoom chain_ao_only and resi {c10_start}-
bg_color white
set ray_shadows, 0
set ambient, 0.3
"""

# Save files
with open(output_pml, 'w') as f:
    f.write(pml_commands)
print(f"Saved: {output_pml}")

cmd.save(output_pse)
print(f"Saved: {output_pse}")

print(f"\nGenerated visualization for {len(classification_df)} positions")
print(f"Sphere sizes range from {classification_df['sphere_size'].min():.1f} to {classification_df['sphere_size'].max():.1f}")