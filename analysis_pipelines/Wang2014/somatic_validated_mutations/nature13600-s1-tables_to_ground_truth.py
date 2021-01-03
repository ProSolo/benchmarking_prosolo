import pandas
import re

from pyliftover import LiftOver
liftover = LiftOver('hg18', 'hg19')

clonal = pandas.read_csv("nature13600-s1-table-s6-clonal.tsv", sep = "\t")
subclonal = pandas.read_csv("nature13600-s1-table-s7-subclonal.tsv", sep = "\t")

all = pandas.concat( [clonal, subclonal] )

# keep only validated variants, i.e. those with a Duplex_P_val smaller 0.01 and not NA
filtered = all[ all['Duplex_P_val'].str.contains('0,0[123456789]|NA|0,[123456789]') == False ]

with open("Wang2014_ground_truth_non_synonymous_variants.hg18_to_hg19.tsv", mode='w') as out:
    print( '\t'.join( filtered.columns.values.tolist() ), file=out, end='\n')
    for index, row in filtered.iterrows():
        print( row['chrom'], row['pos'], row['REF'], row['VAR'] )
        lo = liftover.convert_coordinate(row['chrom'],row['pos'] - 1)[0]
        row['chrom'] = lo[0]
        row['pos'] = lo[1] + 1
        print( '\t'.join(map(str, row)) )
        print( '\t'.join(map(str, row)), file=out, end='\n' )

