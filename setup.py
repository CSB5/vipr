from distutils.core import setup
#, Extension
#import os
import sys
import subprocess

# see also http://docs.python.org/distutils/setupscript.html

def which(prog):
    """make sure prog is in PATH
    """
        
    try:
        subprocess.call([prog], 
                        stderr=subprocess.PIPE, 
                        stdout=subprocess.PIPE,)
        return True
    except OSError:
        return False

DEPS = ['samtools', 'bwa', 'java', 'make']# there are more    
missing = []
for prog in DEPS:
    if not which(prog):
        missing.append(prog)
if len(missing):
    sys.stderr.write("#\nWARNING: couldn't find required program/s '%s'" % ', '.join(missing))
    raw_input("Press Enter to continue anyway.")


scripts = []
scripts.append("./src/bam2cons.py")
scripts.append("./src/bam2cons_iter.sh")
scripts.append("./src/base_qual_calib_wrapper.sh")
scripts.append("./src/bwa_unique.sh")
scripts.append("./src/coverage_plot.py")
scripts.append("./src/mapping_success.sh")
scripts.append("./src/mark_primer.py")
scripts.append("./src/mask_primer.py")
scripts.append("./src/primer_pos_from_seq.sh")
scripts.append("./src/primer_pos_to_bed.py")
scripts.append("./src/vipr.py")
scripts.append("./src/snp_coldspots.py")
scripts.append("./src/snp_hotspots.py")

# To get an approximate list of Python dependencies:
# grep ^src/import | cut -f 2 -d : | sort -u
modules = []
modules.append("./src/primer_pos.py")
modules = [m.replace(".py", "") for m in modules]

setup(name = 'Viral Analysis Pipeline',
      version = '0.1.1',
      description="GIS internal viral analysis pipeline",
      author="Andreas Wilm",
      author_email='wilma@gis.a-star.edu.sg',
      url='https://github.com/CSB5/vipr/',
      scripts = scripts,
      py_modules = modules,
      keywords='bioinformatics'
      )
