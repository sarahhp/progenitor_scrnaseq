
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import yaml

# -------------------------------
# 1. Interaction helpers
# -------------------------------
def ask(question, options=None):
    print(question)
    if options:
        print(f"Options: {', '.join(options)}")
    return input("> ").strip()

def info(msg):
    print(f"[INFO] {msg}")

def confirm(msg):
    return input(f"{msg} (Y/N): ").strip().lower() == "y"

# -------------------------------
# 2. Templates
# -------------------------------
def load_templates(path="templates.yaml"):
    with open(path, "r") as f:
        return yaml.safe_load(f)

# -------------------------------
# 3. Scan source project
# -------------------------------
def find_rmds(source_root):
    # If you’re porting out of workflowr, Rmds are under source_root/analysis
    analysis_dir = os.path.join(source_root, "analysis")
    rmds = []
    for root, _, files in os.walk(analysis_dir):
        for fn in files:
            if fn.endswith(".Rmd"):
                rmds.append(os.path.join(root, fn))
    return rmds

# -------------------------------
# 4. Extraction utilities
# -------------------------------
def extract_libraries(rmd_path):
    libs = set()
    with open(rmd_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith("library("):
                pkg = s.split("(")[1].split(")")[0]
                pkg = pkg.replace('"', '').replace("'", "")
                libs.add(pkg)
    return sorted(libs)

def extract_palettes(rmd_path):
    pal, clusters = None, None
    with open(rmd_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith("pal=") or s.startswith("pal ="):
                pal = s.split("=", 1)[1].strip()
            if s.startswith("clusters=") or s.startswith("clusters ="):
                clusters = s.split("=", 1)[1].strip()
    return pal, clusters

def extract_params_chunk(rmd_path):
    """Return dict of key=value from {r params} chunk; values are raw strings without quotes."""
    meta = {}
    in_params = False
    with open(rmd_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith("```{r") and "params" in s:
                in_params = True
                continue
            if in_params:
                if s.startswith("```"):
                    break
                if "=" in s and not s.startswith("#"):
                    key, val = s.split("=", 1)
                    key = key.strip()
                    val = val.strip().replace('"', '').replace("'", "")
                    meta[key] = val
    return meta

def detect_module(rmd_path, categories):
    meta = extract_params_chunk(rmd_path)
    module = meta.get("module")
    if module and module in categories:
        return module
    # Fallback to filename hint
    base = os.path.basename(rmd_path).lower()
    for cat in categories:
        if cat in base:
            return cat
    return "uncategorized"

# -------------------------------
# 5. Rendering blocks
# -------------------------------
def build_libs_chunk(libs, chunk_opts):
    libs_code = "\n".join(["library({})".format(x) for x in libs])
    return "```{r " + chunk_opts + "}\n" + libs_code + "\n```"

def render_params_block(tpls, metadata, optional_params):
    tmpl = tpls["params_block"]["content"]
    block = tmpl.replace("{module}", metadata["module"]) \
                .replace("{dataname}", metadata["dataname"]) \
                .replace("{load_data}", metadata["load_data"]) \
                .replace("{subset_data}", metadata["subset_data"]) \
                .replace("{optional_params}", optional_params)
    return block

def render_setup_block(tpls, seed, rmd_filename):
    tmpl = tpls["setup_block"]["content"]
    return tmpl.replace("{seed}", str(seed)).replace("{rmd_filename}", rmd_filename)

def render_themes_block(tpls, pal, clusters):
    tmpl = tpls["themes_block"]["content"]
    return tmpl.replace("{pal}", pal).replace("{clusters}", clusters)

# -------------------------------
# 6. Filename generation & conflicts
# -------------------------------
def infer_analysis_type(original_fn):
    base = os.path.basename(original_fn).lower().replace(".rmd", "")
    for hint in ["initial", "qc", "integration", "plots", "tss", "beige", "adipogenesis", "bulk", "progenitors"]:
        if hint in base:
            return "{}_analysis".format(hint) if not hint.endswith("_analysis") else hint
    return "analysis"

def generate_new_filename(module, dataname, original_fn):
    return "{}_{}_{}.Rmd".format(module, dataname, infer_analysis_type(original_fn))

def resolve_conflict(original_fn, new_fn, target_dir):
    target_path = os.path.join(target_dir, new_fn)
    if not os.path.exists(target_path):
        return new_fn
    print("[WARNING] Conflict detected!")
    print("Original file:", os.path.basename(original_fn))
    print("Proposed new name:", new_fn)
    print("Options:\n  [1] Enter a new suffix\n  [2] Auto-append a number\n  [3] Skip this file")
    choice = input("Choose an option (1/2/3): ").strip()
    if choice == "1":
        suffix = input("Enter new suffix: ").strip()
        return new_fn.replace(".Rmd", "_{}.Rmd".format(suffix))
    elif choice == "2":
        i = 2
        candidate = new_fn.replace(".Rmd", "_{}.Rmd".format(i))
        while os.path.exists(os.path.join(target_dir, candidate)):
            i += 1
            candidate = new_fn.replace(".Rmd", "_{}.Rmd".format(i))
        return candidate
    elif choice == "3":
        return None
    else:
        print("[INFO] Invalid choice; skipping file.")
        return None

# -------------------------------
# 7. Migration core
# -------------------------------
def migrate_file(src_rmd, dest_root, tpls, categories, default_seed=1234, dry_run=True, report=None):
    report = report if report is not None else []
    cats = list(categories.keys())
    module = detect_module(src_rmd, cats)
    params = extract_params_chunk(src_rmd)
    dataname = params.get("dataname", "dataset")
    load_data = params.get("load_data", "T")
    subset_data = params.get("subset_data", "F")

    # Optional params beyond core keys (e.g., resolution, etc.)
    optional_keys = [k for k in params.keys() if k not in ["module", "dataname", "load_data", "subset_data"]]
    optional_lines = ["{} = {}".format(k, params[k]) for k in optional_keys]
    optional_params_rendered = "\n".join(optional_lines)

    libs = extract_libraries(src_rmd)
    pal, clusters = extract_palettes(src_rmd)
    if not pal or not clusters:
        pal = categories.get(module, {}).get("pal", 'c("#1F78B4")')
        clusters = categories.get(module, {}).get("clusters", "ggsci::pal_simpsons()(12)")

    new_fn = generate_new_filename(module, dataname, src_rmd)

    # Destination directories (no extra Rmd dir)
    rmd_dir = os.path.join(dest_root, "analysis", module, dataname)
    figs_dir = os.path.join(dest_root, "analysis", module, dataname, "figures")
    os.makedirs(rmd_dir, exist_ok=True)
    os.makedirs(figs_dir, exist_ok=True)

    new_fn_resolved = resolve_conflict(src_rmd, new_fn, rmd_dir)
    if new_fn_resolved is None:
        info("Skipped {}".format(src_rmd))
        report.append({"original": os.path.basename(src_rmd), "new": "-", "action": "Skipped"})
        return report

    if not confirm("Migrate {} as {}?".format(os.path.basename(src_rmd), new_fn_resolved)):
        info("Skipped {}".format(src_rmd))
        report.append({"original": os.path.basename(src_rmd), "new": new_fn_resolved, "action": "Skipped"})
        return report

    # Read original content
    with open(src_rmd, "r", encoding="utf-8") as f:
        content = f.read()

    # Split YAML; remove workflowr output if present
    yaml_header = ""
    body = content
    if "---" in content:
        parts = content.split("---")
        yaml_header = parts[1] if len(parts) > 1 else ""
        body = "---".join(parts[2:]) if len(parts) > 2 else content
        yaml_header = yaml_header.replace("workflowr::wflow_html", "html_document")

    # Ensure output: html_document
    new_header = yaml_header
    if "output:" not in new_header:
        new_header += "\noutput: html_document\n"

    # Render blocks
    libs_block = ""
    for block in tpls.get("global", []):
        if block["name"] == "load_libs" and block["type"] == "dynamic":
            libs_block = build_libs_chunk(libs, block["chunk_options"])

    meta_block = render_params_block(
        tpls,
        {
            "module": module,
            "dataname": dataname,
            "load_data": load_data,
            "subset_data": subset_data
        },
        optional_params_rendered
    )

    rmd_filename_var = os.path.splitext(os.path.basename(new_fn_resolved))[0]
    setup_block = render_setup_block(tpls, default_seed, rmd_filename_var)
    themes_block = render_themes_block(tpls, pal, clusters)

    joined_blocks = "\n\n".join([libs_block, meta_block, setup_block, themes_block])
    new_content = "---" + new_header + "---\n\n" + joined_blocks + "\n\n" + body

    dest_rmd = os.path.join(rmd_dir, new_fn_resolved)
    if dry_run:
        info("Dry-run; preview first 500 chars:")
        print("-" * 40)
        print(new_content[:500] + "...")
        print("-" * 40)
    else:
        with open(dest_rmd, "w", encoding="utf-8") as f:
            f.write(new_content)
        info("Migrated → {}".format(dest_rmd))

    report.append({"original": os.path.basename(src_rmd), "new": new_fn_resolved, "action": "Migrated"})
    return report

# -------------------------------
# 8. Report
# -------------------------------
def print_report(entries):
    print("\nMigration Report")
    for e in entries:
        print("- {} → {} [{}]".format(e["original"], e["new"], e["action"]))

# -------------------------------
# 9. Main
# -------------------------------
def main():
    info("Port‑out migration: workflowr → flexible Rmd")
    root = ask("Enter path to source repository (current workflowr or similar):")
    dest = ask("Enter path to destination repository (target):")
    dry = confirm("Enable dry-run mode?")
    seed = ask("Default RNG seed (integer)?")
    seed = int(seed) if seed.isdigit() else 1234

    tpls = load_templates("templates.yaml")
    rmds = find_rmds(root)
    cats = tpls["categories"]
    report = []

    info("Found {} Rmd files".format(len(rmds)))
    for rmd in rmds:
        report = migrate_file(rmd, dest, tpls, cats, default_seed=seed, dry_run=dry, report=report)

    print_report(report)
    info("Done.")

if __name__ == "__main__":
    main()

