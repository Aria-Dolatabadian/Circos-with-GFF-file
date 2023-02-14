#install pycirclize
import matplotlib.pyplot as plt
from pycirclize import Circos
from pycirclize.utils import load_prokaryote_example_file
from pycirclize.parser import Gff

# Load GFF file
gff_file = ("escherichia_coli.gff.gz")
gff = Gff(gff_file)
# Initialize circos sector by genome size
circos = Circos(sectors={gff.name: gff.range_size})
circos.text("Escherichia coli\n(NC_000913)", size=15)
sector = circos.sectors[0]
# Outer track
outer_track = sector.add_track((98, 100))
outer_track.axis(fc="lightgrey")
outer_track.xticks_by_interval(500000, label_formatter=lambda v: f"{v / 1000000:.1f} Mb")
outer_track.xticks_by_interval(100000, tick_length=1, show_label=False)
# Forward CDS genomic features
f_cds_track = sector.add_track((88, 95), r_pad_ratio=0.1)
f_cds_track.genomic_features(gff.extract_features("CDS", 1),  fc="red")
# Reverse CDS genomic features
r_cds_track = sector.add_track((81, 88), r_pad_ratio=0.1)
r_cds_track.genomic_features(gff.extract_features("CDS", -1), fc="blue")
# rRNA genomic features
rrna_track = sector.add_track((74, 81), r_pad_ratio=0.1)
rrna_track.genomic_features(gff.extract_features("rRNA"), fc="limegreen")
# tRNA genomic features
trna_track = sector.add_track((67, 74), r_pad_ratio=0.1)
trna_track.genomic_features(gff.extract_features("tRNA"), color="magenta", lw=0.1)

fig = circos.plotfig()


plt.show()

