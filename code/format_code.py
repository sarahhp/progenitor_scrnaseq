
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
def find_rmds(source_root, port_dir="analysis", pattern=""):
    print("[INFO] seraching source dir",source_root)
    # If you’re porting out of workflowr, Rmds are under source_root/analysis
    subdirs = [d for d in os.listdir(source_root) if os.path.isdir(d)]
    if "analysis" in subdirs:
        rmd_dir = os.path.join(source_root, port_dir)
    else: rmd_dir = source_root
    rmds = []
    for root, _, files in os.walk(rmd_dir):
        for fn in files:
            if fn.endswith(".Rmd"):
                if pattern == "":
                    rmds.append(os.path.join(root, fn))
                else:
                    if pattern.strip().lower() in fn.lower():
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
    pal, clusters,spectral = None, None, False
    with open(rmd_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith("pal=") or s.startswith("pal ="):
                pal = s.split("=", 1)[1].strip()
            if s.startswith("clusters=") or s.startswith("clusters ="):
                clusters = s.split("=", 1)[1].strip()
                #if "(" in clusters and ")" not in clusters:
                    #grab next line
                    
            if s.startswith("spectral=") or s.startswith("spectral ="):
                spectral = True
    return pal, clusters, spectral

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
            if s.startswith("module=") or s.startswith("module ="):
                in_params = True
            if in_params:
                if s.startswith("```"):
                    break
                if "=" in s and not s.startswith("#") and not "{" in s and not "i_am" in s:
                    key, val = s.split("=", 1)
                    key = key.strip()
                    val = val.strip().replace('"', '\"').replace("'", "\'")
                    meta[key] = val
            else:
                #Fallback to filename hint
                fi = rmd_path.split("_")
                meta["module"] = fi[0]
                meta["dataname"] = fi[1]
    #print("Extracted params {}".format(meta))
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
        elif "white_only" in base:
            return "adipogenesis"
        elif "tss" in base:
            return "beige_vs_white"
    return "uncategorized"

# -------------------------------
# 5. Rendering blocks
# -------------------------------
def build_libs_chunk(libs, chunk_opts):
    libs_code = "\n".join(["library({})".format(x) for x in libs])
    return "```{r libs, " + chunk_opts + "}\n" + libs_code + "\n```"

def render_params_block(tpls, metadata, optional_params):
    tmpl = tpls["params_block"]["content"]
    block = tmpl.replace("{rmd_filename}", metadata["rmd_filename"]) \
                .replace("{module}", metadata["module"]) \
                .replace("{dataname}", metadata["dataname"]) \
                .replace("{load_data}", metadata["load_data"]) \
                .replace("{subset_data}", metadata["subset_data"]) \
                .replace("{odir}", metadata["odir"]) \
                .replace("{optional_params}", optional_params)
    #block_fix = block.replace("\"file.path","file.path") \
     #                 .replace("dataname)\"","dataname)")
    return block

def render_setup_block(tpls, seed):
    tmpl = tpls["setup_block"]["content"]
    return tmpl.replace("{seed}", str(seed))

def render_themes_block(tpls, pal, clusters, spectral):
    tmpl = tpls["themes_block"]["content"]
    if spectral:
        spectral = "spectral = rev(RColorBrewer::brewer.pal(n = 11, name = 'Spectral'))"
        block = tmpl.replace("\n```\n", "\n" + spectral +"\n```\n") 
        return block.replace("{pal}", pal).replace("{clusters}", clusters)
    else:
        return tmpl.replace("{pal}", pal).replace("{clusters}", clusters)


# -------------------------------
# 6. Remove duplicate chunks
# -------------------------------

def remove_duplicate_chunks(body):
    rmd_list = body.split("\n")
    
    print("Creating indexes")
    chunks_to_remove = []    
    chunk_start_index = [i for i,line in enumerate(rmd_list) if line.startswith("```{r")]
    chunk_end_index = [i for i,line in enumerate(rmd_list) if line.endswith("```")]    
    library_call_index = [i for i,line in enumerate(rmd_list) if line.startswith("library(")]
    assert len(chunk_start_index) == len(chunk_end_index), "Different number of chunk start to end strings. Check if there are any additional ``` in Rmd document"
    
    print("Checking for chunks to remove...")
    for i,cs in enumerate(chunk_start_index):
        ce = chunk_end_index[i]
        
        # Check for names
        header = rmd_list[cs].split(" ")
        if len(header) > 1:
            name=header[1]
            if name in ["libs","themes","load_libs","params","setup"]:
                chunks_to_remove.append((cs,ce))
                
        # Check for unnamed library chunks
        else: 
            lib_calls_in_chunk = 0
            for lib_call in library_call_index:
                if lib_call > cs & lib_call < ce: lib_calls_in_chunk += 1
            if lib_calls_in_chunk > 2:
               chunks_to_remove.append((cs,ce))
    
    print("Converting chunks to lines")
    # Convert chunks to lines
    lines_to_remove = []
    for cs,ce in chunks_to_remove:
        i = cs
        while i <= ce:
            lines_to_remove.append(i)
            i+=1
            
    print("Removing lines")
    new_body = [line for j,line in enumerate(rmd_list) if j not in lines_to_remove]
    new_body = "\n".join(new_body)
    return(new_body)

# -------------------------------
# 7. Filename generation & conflicts
# -------------------------------
def infer_analysis_type(original_fn, module, dataname):
    base = os.path.basename(original_fn).lower().replace(".rmd", "")
    #Check for common analysis types
    suffix = base.replace(module +"_", "") \
                .replace(dataname+"_","") \
                .replace(dataname,"")
    print(suffix)
    if suffix.strip() == "":
      return None
    
    for hint in ["qc", "integration", "plots", "GO"]:
        if hint in base:
            return hint 
    # If not use end of filename
    fns = suffix.split("_")
    return "_".join(fns[2:])
    

# def generate_new_filename(module, dataname, original_fn):
#     analysis_type = infer_analysis_type(original_fn, module, dataname)
#     if analysis_type == None:
#         return "{}_{}.Rmd".format(module, dataname)
#     else:
#         return "{}_{}_{}.Rmd".format(module, dataname, analysis_type)

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
# 8. Migration core
# -------------------------------
def migrate_file(src_rmd, dest_root, tpls, categories, default_seed=1234, dry_run=True, report=None):
    report = report if report is not None else []
    cats = list(categories.keys())
    module = detect_module(src_rmd, cats)
    params = extract_params_chunk(src_rmd)
    dataname = params.get("dataname", "dataset")
    odir = params.get("odir","here('output',module, dataname)")
    load_data = params.get("load_data", "T")
    subset_data = params.get("subset_data", "F")

    # Optional params beyond core keys (e.g., resolution, etc.)
    optional_keys = [k for k in params.keys() if k not in ["module", "dataname", "load_data", "subset_data","odir"]]
    optional_lines = ["{} = {}".format(k, params[k]) for k in optional_keys]
    optional_params_rendered = "\n".join(optional_lines)

    libs = extract_libraries(src_rmd)
    pal, clusters,spectral = extract_palettes(src_rmd)
    if not pal or not clusters:
        pal = categories.get(module, {}).get("pal", 'c("#1F78B4")')
        clusters = categories.get(module, {}).get("clusters", "ggsci::pal_simpsons()(12)")

    new_fn = os.path.basename(src_rmd)#generate_new_filename(module, dataname, src_rmd)

    # Destination directories
    rmd_dir = os.path.join(dest_root, "analysis", module)
    figs_dir = os.path.join(dest_root, "analysis", module, dataname, "figures")

    new_fn_resolved = resolve_conflict(src_rmd, new_fn, rmd_dir)
    if new_fn_resolved is None:
        info("Skipped {}".format(src_rmd))
        report.append({"original": os.path.basename(src_rmd), "new": "-", "action": "Skipped"})
        return report

    if not confirm("Migrate {} as {}?".format(os.path.basename(src_rmd), new_fn_resolved)):
        add_new_name = ask("Enter a new file name or press Enter to skip file migration.")
        if add_new_name != "":
          new_fn_resolved = resolve_conflict(src_rmd, add_new_name, rmd_dir)
        else:
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
        
    # Extract preamble for document top
    preamble = ""
    parts = body.split("```{r")
    if "```" not in parts[0]:
        preamble = parts[0]
        body = "```{r" + "```{r".join(parts[1:])
        print("Adding preamble to the top of document")
        

    # Render blocks
    libs_block = ""
    for block in tpls.get("global", []):
        if block["name"] == "load_libs" and block["type"] == "dynamic":
            libs_block = build_libs_chunk(libs, block["chunk_options"])
    
    rmd_filename_var = os.path.splitext(os.path.basename(new_fn_resolved))[0]
    meta_block = render_params_block(
        tpls,
        {
            "rmd_filename": rmd_filename_var,
            "module": module,
            "dataname": dataname,
            "load_data": load_data,
            "subset_data": subset_data,
            "odir":odir
            
        },
        optional_params_rendered
    )

    setup_block = render_setup_block(tpls, default_seed)
    themes_block = render_themes_block(tpls, pal, clusters, spectral)
    
    new_body = body
    #new_body = remove_duplicate_chunks(body)
    
    ## Find and replace edits
    search_strs = tpls["strings_to_replace"]
    for s in search_strs.keys():
        if s in new_body:
            new_body = new_body.replace(s, search_strs[s])
            print("Replacing \"{}\" {} times".format(s,new_body.count(s)))
            
    ##Build Rmd  
    joined_blocks = "\n\n".join([libs_block, meta_block, setup_block, themes_block])
    new_content = "---" + new_header + "---\n" + preamble + \
                    joined_blocks + "\n\n" + new_body

    dest_rmd = os.path.join(rmd_dir, new_fn_resolved)
    if dry_run:
        info("Dry-run; preview first 2000 chars:")
        print("-" * 40)
        print(new_content[:2000] + "...")
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
    is_cwd_root = confirm("Is the current working dir the root dir for the workflow?")
    if is_cwd_root :
        root = "./"
        dest = "./"
    else:
        root = ask("Enter path to source repository (current workflowr or similar):")
        dest = ask("Enter path to destination repository (target):")
        
    pattern = ask("Would you like to migrate all files? Press enter for yes or enter a pattern to migrate file names containing that pattern.  ")
    dry = confirm("Enable dry-run mode?")
    seed = 8

    tpls = load_templates("code/r_chunk_templates.yaml")
    rmds = find_rmds(root, pattern=pattern)
    cats = tpls["categories"]
    report = []

    info("Found {} Rmd files".format(len(rmds)))
    for rmd in rmds:
        report = migrate_file(rmd, dest, tpls, cats, default_seed=seed, dry_run=dry, report=report)

    print_report(report)
    info("Done.")

if __name__ == "__main__":
    main()

