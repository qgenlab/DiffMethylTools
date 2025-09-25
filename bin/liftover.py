import sys
from pyliftover import LiftOver


input_bed = sys.argv[1]
chain_file = sys.argv[2]
base_output = sys.argv[3]

output_bed = f"{base_output}.bed"
unmapped_bed = f"{base_output}.unmapped.bed"

lo = LiftOver(chain_file)

with open(input_bed, "r") as infile, open(output_bed, "w") as out_mapped, open(unmapped_bed, "w") as out_unmapped:
    for line in infile:
        line = line.rstrip()
        if line.startswith("#") or line == "":
            continue
        fields = line.split("\t")
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        lifted_start = lo.convert_coordinate(chrom, start)
        lifted_end = lo.convert_coordinate(chrom, end)
        if lifted_start and lifted_end:
            new_chrom, new_start, strand, _ = lifted_start[0]
            _, new_end, _, _ = lifted_end[0]
            fields[0] = new_chrom
            fields[1] = str(new_start)
            fields[2] = str(new_end)
            out_mapped.write("\t".join(fields) + "\n")
        else:
            out_unmapped.write(line + "\n")


print(f"Liftover complete:")
print(f"Mapped BED file: {output_bed}")
print(f"Unmapped lines: {unmapped_bed}")


