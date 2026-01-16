#!/usr/bin/env python3
"""
Patches MetaPop_Colab.ipynb to add a cell that copies BAM files 
from Google Drive to local Colab storage.

Run this script once to update the notebook.
"""

import json
import os

# Path to the notebook
notebook_path = os.path.join(os.path.dirname(__file__), "MetaPop_Colab.ipynb")

# New cells to insert after "Mount Google Drive" section
new_markdown_cell = {
    "cell_type": "markdown",
    "metadata": {},
    "source": [
        "## 2.5 Copy BAM Files to Local Storage (Recommended)\n",
        "\n",
        "**Important:** BAM files accessed directly from Google Drive can cause `Exec format error` with samtools.\n",
        "Run this cell to copy your BAM files to local Colab storage for reliable processing."
    ]
}

new_code_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [
        "import shutil\n",
        "import os\n",
        "\n",
        "# === CONFIGURE YOUR PATHS HERE ===\n",
        "drive_bam_dir = '/content/drive/MyDrive/metapop_data/bams'  # Your Google Drive BAM folder\n",
        "drive_ref_dir = '/content/drive/MyDrive/metapop_data/references'  # Your Google Drive reference folder\n",
        "local_bam_dir = '/content/local_bams'  # Local Colab storage for BAMs\n",
        "local_ref_dir = '/content/local_refs'  # Local Colab storage for references\n",
        "\n",
        "def copy_files_locally(src_dir, dst_dir, extensions):\n",
        "    \"\"\"Copy files with specified extensions from src to dst.\"\"\"\n",
        "    if not os.path.exists(src_dir):\n",
        "        print(f\"Warning: Source directory not found: {src_dir}\")\n",
        "        return 0\n",
        "    \n",
        "    os.makedirs(dst_dir, exist_ok=True)\n",
        "    copied = 0\n",
        "    for f in os.listdir(src_dir):\n",
        "        if any(f.lower().endswith(ext) for ext in extensions):\n",
        "            src = os.path.join(src_dir, f)\n",
        "            dst = os.path.join(dst_dir, f)\n",
        "            if not os.path.exists(dst):\n",
        "                print(f\"  Copying {f}...\")\n",
        "                shutil.copy2(src, dst)\n",
        "                copied += 1\n",
        "            else:\n",
        "                print(f\"  {f} already exists locally\")\n",
        "    return copied\n",
        "\n",
        "print(\"Copying BAM files to local storage...\")\n",
        "bam_count = copy_files_locally(drive_bam_dir, local_bam_dir, ['.bam', '.bai'])\n",
        "print(f\"Copied {bam_count} BAM/BAI files\\n\")\n",
        "\n",
        "print(\"Copying reference files to local storage...\")\n",
        "ref_count = copy_files_locally(drive_ref_dir, local_ref_dir, ['.fa', '.fasta', '.fna', '.fai'])\n",
        "print(f\"Copied {ref_count} reference files\\n\")\n",
        "\n",
        "print(\"=\"*50)\n",
        "print(\"Done! Update the paths in Section 4 to:\")\n",
        "print(f\"  BAM Directory: {local_bam_dir}\")\n",
        "print(f\"  Reference Directory: {local_ref_dir}\")\n",
        "print(\"=\"*50)\n",
        "\n",
        "# Verify samtools can read the files\n",
        "if os.path.exists(local_bam_dir) and os.listdir(local_bam_dir):\n",
        "    bam_files = [f for f in os.listdir(local_bam_dir) if f.endswith('.bam')]\n",
        "    if bam_files:\n",
        "        first_bam = bam_files[0]\n",
        "        print(f\"\\nVerifying samtools can read {first_bam}...\")\n",
        "        !samtools quickcheck {local_bam_dir}/{first_bam} && echo \"✅ BAM file is valid\" || echo \"❌ BAM file check failed\""
    ]
}

def patch_notebook():
    # Load the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook = json.load(f)
    
    # Find the index of the cell after "Mount Google Drive" code cell
    insert_index = None
    for i, cell in enumerate(notebook['cells']):
        if cell['cell_type'] == 'code':
            source = ''.join(cell['source']) if isinstance(cell['source'], list) else cell['source']
            if 'drive.mount' in source:
                insert_index = i + 1
                break
    
    if insert_index is None:
        print("ERROR: Could not find the Google Drive mount cell.")
        return False
    
    # Check if already patched
    for cell in notebook['cells']:
        source = ''.join(cell['source']) if isinstance(cell['source'], list) else cell['source']
        if 'Copy BAM Files to Local Storage' in source:
            print("Notebook already patched!")
            return True
    
    # Insert the new cells
    notebook['cells'].insert(insert_index, new_markdown_cell)
    notebook['cells'].insert(insert_index + 1, new_code_cell)
    
    # Save the notebook
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=1)
    
    print(f"Successfully patched {notebook_path}")
    print("Added new section 2.5: Copy BAM Files to Local Storage")
    return True

if __name__ == "__main__":
    patch_notebook()
