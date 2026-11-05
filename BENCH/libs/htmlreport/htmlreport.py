##########################################################
#  CROCO under CeCILL-C
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
import json
import shutil
import base64
import glob

# internal
from ..helpers import Messaging


def format_status(value):
    if value is True:
        format_td = '<td class="status-true">Done</td>'
    elif value is False:
        format_td = '<td class="status-false">Error</td>'
    else:
        format_td = '<td class="status-null">N/A</td>'
    return format_td


def generate_status_html(data):
    html_rows = ""
    for key, entry in data.items():
        case = entry["case"]
        if entry["restarted"]:
            variant = entry["variant"] + "-rst"
        else:
            variant = entry["variant"]
        status = entry["status"]
        html_rows += f"""
            <tr>
               <td>{case}</td>
               <td>{variant}</td>
               {format_status(status["build"]["status"])}
               {format_status(status["run"]["status"])}
               {format_status(status["check"]["status"])}
               {format_status(status["plotphy"]["status"])}
            </tr>
            """
    html_table = f"""
        <table>
            <thead>
                  <tr>
                     <th>Case</th>
                     <th>Variant</th>
                     <th>Build Status</th>
                     <th>Run Status</th>
                     <th>Check Status</th>
                     <th>Plotphy Status</th>
                  </tr>
            </thead>
            <tbody>
                  {html_rows}
            </tbody>
        </table>"""

    return html_table


def _encode_image_b64(img_path):
    """Return a data-URI for a PNG so the gallery is self-contained."""
    with open(img_path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f"data:image/png;base64,{data}"


def _plotphy_gallery_html(result_dir, root_path, ref_result_dir=None):
    """
    Build an inline plotphy gallery for one case+variant leaf directory.

    result_dir      – directory that may contain PNG files (the current variant)
    root_path       – root of the whole results tree (used to build relative URLs
                      for the iframe links)
    ref_result_dir  – if provided, its PNGs are shown side-by-side with the
                      current ones (current left, reference right)

    Returns an HTML string, or "" when no PNG files are found in result_dir.
    """
    png_paths = sorted(glob.glob(os.path.join(result_dir, "*.png")))
    if not png_paths:
        return ""

    def rel(p):
        return os.path.relpath(p, root_path).replace(os.sep, "/")

    if ref_result_dir and os.path.isdir(ref_result_dir):
        # ---- side-by-side mode ------------------------------------------------
        items_html = ""
        for png in png_paths:
            fname = os.path.basename(png)
            ref_png = os.path.join(ref_result_dir, fname)

            cur_src = _encode_image_b64(png)
            cur_link = rel(png)

            if os.path.isfile(ref_png):
                ref_src = _encode_image_b64(ref_png)
                ref_link = rel(ref_png)
                ref_cell = f"""
                    <div class="plot-cell">
                        <div class="plot-label plot-label-ref">ref</div>
                        <a href="{ref_link}" target="main-frame">
                            <img src="{ref_src}" alt="{fname} (ref)" loading="lazy">
                        </a>
                    </div>"""
            else:
                ref_cell = '<div class="plot-cell plot-missing">ref not found</div>'

            items_html += f"""
            <div class="plot-pair">
                <div class="plot-pair-title">{fname}</div>
                <div class="plot-pair-cols">
                    <div class="plot-cell">
                        <div class="plot-label plot-label-cur">current</div>
                        <a href="{cur_link}" target="main-frame">
                            <img src="{cur_src}" alt="{fname}" loading="lazy">
                        </a>
                    </div>
                    {ref_cell}
                </div>
            </div>"""

        return f'<div class="plotphy-gallery plotphy-compare">{items_html}</div>'

    else:
        # ---- single-run mode --------------------------------------------------
        items_html = ""
        for png in png_paths:
            fname = os.path.basename(png)
            src = _encode_image_b64(png)
            link = rel(png)
            items_html += f"""
            <div class="plot-item">
                <a href="{link}" target="main-frame">
                    <img src="{src}" alt="{fname}" loading="lazy">
                </a>
                <div class="plot-caption">{fname}</div>
            </div>"""

        return f'<div class="plotphy-gallery">{items_html}</div>'


def _find_ref_variant_dir(variant_dir, ref_variant_name, use_ref_path=None):
    """
    Resolve the reference directory for a given variant leaf directory.

    - use_ref_path set  → look for <use_ref_path>/<case>/<variant_name>
                          (mirrors the same layout in the external ref tree)
    - use_ref_path None → sibling variant: <case_dir>/<ref_variant_name>

    Returns None when the candidate does not exist, when no ref is configured,
    or when variant_dir is already the reference variant itself.
    """
    current_variant = os.path.basename(variant_dir)
    case_dir = os.path.dirname(variant_dir)

    if use_ref_path:
        # External reference tree: same case/variant layout rooted at use_ref_path
        case_name = os.path.basename(case_dir)
        candidate = os.path.join(use_ref_path, case_name)
        return candidate if os.path.isdir(candidate) else None
    else:
        if ref_variant_name is None:
            return None
        if current_variant == ref_variant_name:
            return None  # don't compare the ref against itself
        candidate = os.path.join(case_dir, ref_variant_name)
        return candidate if os.path.isdir(candidate) else None


def _generate_all_plots_section(root_path, ref_variant_name=None, use_ref_path=None):
    """
    Build a flat <details> block listing every case+variant that has PNG files.
    Structure:
        <details id="plots-section">
          <summary>Plots</summary>
          <button onclick="window.print()">Export PDF</button>
    root_path = os.path.abspath(root_path)
    if use_ref_path:
        use_ref_path = os.path.abspath(use_ref_path)
          <h2>cas : BASIN - variant : sequential</h2>
          <div class="plotphy-gallery">...</div>
          <h2>cas : BASIN - variant : mpi-4</h2>
          ...
        </details>

    Returns "" if no PNG is found anywhere under root_path.
    """
    blocks = []

    for case_name in sorted(os.listdir(root_path), key=str.lower):
        case_dir = os.path.join(root_path, case_name)
        if not os.path.isdir(case_dir) or os.path.islink(case_dir):
            continue

        for variant_name in sorted(os.listdir(case_dir), key=str.lower):
            variant_dir = os.path.join(case_dir, variant_name)
            if not os.path.isdir(variant_dir) or os.path.islink(variant_dir):
                continue
            png_files = glob.glob(os.path.join(variant_dir, "*.png"))
            if not png_files:
                continue

            ref_dir = _find_ref_variant_dir(variant_dir, ref_variant_name, use_ref_path)
            gallery_html = _plotphy_gallery_html(variant_dir, root_path, ref_dir)
            if not gallery_html:
                continue

            blocks.append(
                f'<h2 class="plots-section-title">cas : {case_name} — variant : {variant_name}</h2>\n'
                f"{gallery_html}"
            )

    if not blocks:
        return ""

    inner = "\n".join(blocks)
    return (
        '<details id="plots-section">\n'
        "  <summary>Plots</summary>\n"
        '  <div class="plots-section-toolbar">'
        '<button class="pdf-btn" onclick="window.print()">⬇ Export PDF</button>'
        "</div>\n"
        f'  <div class="plots-section-body">{inner}</div>\n'
        "</details>\n"
    )


def generate_treeview_json(root_path, ref_variant_name=None, use_ref_path=None):
    """
    Walk root_path and build the treeview JSON structure.

    ref_variant_name – variant name used as reference for side-by-side comparison
                       (e.g. "sequential").
    use_ref_path     – path to an external reference results tree (--use-ref).

    # Normalize to absolute so all os.path.join / os.path.relpath calls are consistent
    root_path = os.path.abspath(root_path)
    if use_ref_path:
        use_ref_path = os.path.abspath(use_ref_path)
                       When set, comparisons are made against that tree instead
                       of a sibling variant in the current tree.

    Layout produced:
      <details> case_name
          [all plots of all variants, immediately visible]
          <details> variant_name   ← still there for other files / sub-nodes
              ...non-PNG files...
          </details>
      </details>
    """

    def traverse_directory(path, depth=0):
        nodes = []
        files = []
        directories = []
        status_report_node = None

        for name in sorted(os.listdir(path), key=str.lower):
            full_path = os.path.join(path, name)
            relative_path = os.path.relpath(full_path, root_path)
            url_path = relative_path.replace(os.sep, "/")

            if os.path.isdir(full_path):
                if not os.path.islink(full_path):
                    # Recurse – but strip PNG gallery from sub_nodes when we are
                    # at case level (the gallery is rendered as a grouped block).
                    sub_nodes = traverse_directory(full_path, depth + 1)
                    directories.append({"text": name, "nodes": sub_nodes})
            else:
                if name == "status_report.json":
                    status_report_node = {
                        "text": generate_status_html(json.load(open(full_path, "r")))
                    }
                elif ("treeview.html" not in url_path) and (
                    "styles.css" not in url_path
                ):
                    files.append(
                        {"text": f'<a href="{url_path}" target="main-frame">{name}</a>'}
                    )

        files = sorted(files, key=lambda x: x["text"].lower())
        directories = sorted(directories, key=lambda x: x["text"].lower())

        if status_report_node:
            nodes.append(status_report_node)

        nodes.extend(files)
        nodes.extend(directories)
        return nodes

    return traverse_directory(root_path)


# Function to convert JSON to HTML
def json_to_html(data):
    html_content = ""
    for item in data:
        if isinstance(item, dict):
            text = item.get("text", "")
            if "nodes" in item:
                html_content += "<details>\n"
                html_content += f"  <summary>{text}</summary>\n"
                html_content += '  <div class="tree">\n'
                html_content += json_to_html(item["nodes"])
                html_content += "  </div>\n"
                html_content += "</details>\n"
            else:
                html_content += text + "\n"
    return html_content


# Function to convert JSON tree to HTML
def generate_html(
    data_path,
    output_file="treeview.html",
    navigation_content="",
    style_file="libs/htmlreport/styles.css",
    ref_variant_name=None,
    use_ref_path=None,
):
    json_data = generate_treeview_json(
        data_path,
        ref_variant_name=ref_variant_name,
        use_ref_path=use_ref_path,
    )

    all_plots_section = _generate_all_plots_section(
        data_path,
        ref_variant_name=ref_variant_name,
        use_ref_path=use_ref_path,
    )

    html_output = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CROCO test report</title>
    <link rel="stylesheet" href="styles.css">
    <style>
    @media print {
        /* Force color rendering — browsers strip colors by default to save ink */
        * {
            -webkit-print-color-adjust: exact !important;
            print-color-adjust: exact !important;
            color-adjust: exact !important;
        }
        /* Hide everything except the plots section content */
        nav, .tree, iframe, h1, .plots-section-toolbar { display: none !important; }
        body { display: block; margin: 0; }
        main { padding: 0; }
        #plots-section { display: block !important; }
        #plots-section > summary { display: none; }
        .plots-section-body { display: block; }
        .plots-section-title {
            font-size: 14pt;
            margin: 16pt 0 6pt 0;
            page-break-after: avoid;
        }
        /* Gallery: wrap images, constrain size so they don't overflow */
        .plotphy-gallery, .plotphy-compare {
            display: flex !important;
            flex-wrap: wrap !important;
            gap: 6pt;
            page-break-inside: avoid;
            width: 100%%;
            box-sizing: border-box;
        }
        .plot-item, .plot-cell {
            flex: 0 1 auto;
            max-width: 45%%;
            box-sizing: border-box;
        }
        .plot-item img, .plot-cell img {
            width: 100%% !important;
            max-width: 100%% !important;
            height: auto !important;
            border: 1px solid #ccc;
            display: block;
        }
        .plot-pair {
            page-break-inside: avoid;
            width: 100%%;
            box-sizing: border-box;
        }
        .plot-pair-cols {
            display: flex !important;
            flex-wrap: wrap !important;
            gap: 6pt;
        }
    }
    </style>
</head>
<body>
    %s
    <main>
    <h1>CROCO test report - %s</h1>

    %s

    <div class="tree">
""" % (
        navigation_content,
        os.path.split(data_path)[-1],
        all_plots_section,
    )

    html_output += json_to_html(json_data)

    html_output += """
    </div>

    <iframe id="main-frame" name="main-frame"></iframe>
    </main>
</body>
</html>
"""

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_output)
    shutil.copyfile(
        style_file,
        os.path.join(data_path, "styles.css"),
    )

    Messaging.step(f"HTML file generated: {output_file}")


# Function to create index html file for all results
def generate_global_html(
    base_dir,
    style_file="libs/htmlreport/styles.css",
    ref_name=None,
    use_ref_path=None,
):
    """
    Generate one treeview.html per run-folder found in base_dir, plus an
    index.html at the root.

    ref_name     – the variant name used as reference (e.g. "sequential").
    use_ref_path – path to an external reference results tree (from --use-ref).
    """
    output_file = os.path.join(base_dir, "index.html")

    folders = [
        folder
        for folder in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, folder))
        and not os.path.islink(os.path.join(base_dir, folder))
    ]

    nav_in_folder = f"""  <nav>
      <h2><a href="../index.html">Results</a></h2>
      <ul>
         {
        "".join(
            f'<li><a href="../{folder}/treeview.html">View {folder}</a></li>'
            for folder in folders
        )
    }
      </ul>
      </nav>"""

    for folder in folders:
        folder_path = os.path.join(base_dir, folder)
        if not os.path.islink(folder_path):
            generate_html(
                folder_path,
                output_file=os.path.join(folder_path, "treeview.html"),
                navigation_content=nav_in_folder,
                style_file=style_file,
                ref_variant_name=ref_name,
                use_ref_path=use_ref_path,
            )

    html_content = f"""<!DOCTYPE html>
   <html lang="en">
   <head>
   <meta charset="UTF-8">
   <meta name="viewport" content="width=device-width, initial-scale=1.0">
   <title>CROCO test report</title>
   <link rel="stylesheet" href="styles.css">
   </head>
   <body>
   <nav>
      <h2><a href="index.html">Results</a></h2>
      <ul>
         {
        "".join(
            f'<li><a href="{folder}/treeview.html">View {folder}</a></li>'
            for folder in folders
        )
    }
      </ul>
   </nav>
   <main>
      <h1>CROCO test report</h1>
      <p>Use the navigation menu on the left to explore the treeviews in each folder.</p>
   </main>
   </body>
   </html>
   """

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)
    shutil.copyfile(
        style_file,
        os.path.join(base_dir, "styles.css"),
    )

    Messaging.step(f"Index file generated: {os.path.abspath(output_file)}")
