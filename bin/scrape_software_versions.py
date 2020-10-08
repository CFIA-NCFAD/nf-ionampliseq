#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    "peterk87/nf-ionampliseq": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "FastQC": ["v_fastqc.txt", r"FastQC v(\S+)"],
    "MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
    "samtools": ["v_samtools.txt", r"samtools (\S+)"],
    "bcftools": ["v_bcftools.txt", r"bcftools (\S+)"],
    "Torrent Mapper (TMAP)": ["v_tmap.txt", r"Version\:\s+(\S+)"],
    "Torrent Variant Caller (TVC)": ["v_tvc.txt", r"tvc (\S+)"],
    "Mash": ["v_mash.txt", r"(\S+)"],
    "pigz": ["v_pigz.txt", r"pigz (\S+)"],
}
results = OrderedDict()
results["peterk87/nf-ionampliseq"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["Torrent Mapper (TMAP)"] = '<span style="color:#999999;">N/A</span>'
results["Torrent Variant Caller (TVC)"] = '<span style="color:#999999;">N/A</span>'
results["samtools"] = '<span style="color:#999999;">N/A</span>'
results["bcftools"] = '<span style="color:#999999;">N/A</span>'
results["Mash"] = '<span style="color:#999999;">N/A</span>'
results["FastQC"] = '<span style="color:#999999;">N/A</span>'
results["MultiQC"] = '<span style="color:#999999;">N/A</span>'
results["pigz"] = '<span style="color:#999999;">N/A</span>'

# regex to remove ANSI escape sequences from strings
reaesc = re.compile(r'\x1b[^m]*m')

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            versions = reaesc.sub('', versions)
            match = re.search(v[1], versions)
            if match:
                results[k] = f"v{match.group(1)}"
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'peterk87/nf-ionampliseq Software Versions'
section_href: 'https://github.com/peterk87/nf-ionampliseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
# adjust width and left margin of dt and dd elements to fit text content without ellipsis
for k, v in results.items():
    print(f"""        <dt style=\"width:200px !important;\">
        {k}
        </dt>
        <dd style=\"margin-left:220px !important;\"><samp>
        {v}
        </samp></dd>""")
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write(f"{k}\t{v}\n")
