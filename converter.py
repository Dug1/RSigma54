import sys
import subprocess
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

def meme(filenames, outdir):
    return subprocess.run(["meme-meme"] + filenames + ["-dna","-oc", outdir])

def rpstblastn(input_file_name, output_file_name, database_path):
    return subprocess.run(["rpstblastn", "-query", input_file_name, "-out", output_file_name, "-db", database_path, "-outfmt", "5"])

def blastp(queryfile, db, outfile):
    return subprocess.run(["blastp", "-query", queryfile, "-db", db, "-out", outfile, "-outfmt", "5"])

def makeblastdb(f, db_type):
    return subprocess.run(["makeblastdb", "-in", f, "-dbtype", db_type])


def filter_matches(filename, e):
    """ Parses the the matches results of rpstblastn and save every match 
    above a certain evalue"""
    transcription_factors = []
    with open(filename) as f:
        for record in NCBIXML.parse(f):
            best = e
            if record.alignments:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < best:
                            best = hsp.expect

            if best < e:
                iden = record.query.split(" ")[0]
                locus = iden.split("|")[1]
                transcription_factors.append(locus)

    return transcription_factors

def find_best_match(filename, e):
    """ Find the best match for every query protein above a certain e-value
    above a certain evalue"""
    best_matches = {}
    with open(filename) as f:
        for record in NCBIXML.parse(f):
            best = e
            best_name = ""
            if record.alignments:
                for alignment in record.alignments:
                    print(alignment.hit_def)
                    locus = alignment.hit_def.split(" ")[0].split("|")[1]
                    for hsp in alignment.hsps:
                        if hsp.expect < best:
                            best_name = locus
                            best = hsp.expect

            if best < e:
                iden = record.query.split(" ")[0]
                best_matches[iden] = best_name

    return best_matches 

def extract_genes(filename, fmt, outname, outfmt, locuses=None):
    """ Extra all proteins coded in the file """

    genome = SeqIO.parse(filename, fmt)
    proteins = []
    for record in genome:
        for feature in filter(lambda f: f.type == "CDS", record.features):
            qualifiers = feature.qualifiers
            if locuses is None or qualifiers["locus_tag"][0] in locuses:
                name = qualifiers["gene"][0] if "gene" in qualifiers else ""
                desc = qualifiers["product"][0] if "product" in  qualifiers else ""
                protein = SeqRecord(Seq(qualifiers["translation"][0], IUPAC.protein), 
                    id="{0}|{1}".format(record.id, qualifiers["locus_tag"][0]),
                    name=name,
                    dbxrefs=qualifiers["db_xref"],
                    description=desc)

                proteins.append(protein)
    
    SeqIO.write(proteins, outname, outfmt)


def find_transcription_factors(query_id, filename, output_file, profiledb, e):
    """ Finds all transcription factors which match the [profiledb] 
    and saves them to a fasta file"""

    temp_file = "{0}-tf.fsa".format(query_id)
    extract_genes(filename + ".gbk", "genbank", temp_file, "fasta")
    outfile = "{0}-tf.xml".format(query_id)
    print(rpstblastn(temp_file, outfile, profiledb))
    locuses = filter_matches(outfile, e)
    print(locuses)
    extract_genes(filename + ".gbk", "genbank", output_file + ".fsa", "fasta", locuses)
    
def find_orthologs(job_id, test_against_matrix):
    """ Finds the orthologs of in a given set of each query set of genes against a list
    of others."""

    files = []
    for query, data in test_against_matrix.items():
        original_file, test_against = data
        makeblastdb(original_file, "prot")
        orthologs = dict()
        for filename in test_against:
            makeblastdb(filename, "prot")
            out_file = "{0}-ortho".format(job_id)
            blastp(query,filename, out_file)
            best_locuses = find_best_match(out_file, 1)
            print(best_locuses)

            query_file = "{0}-query.fsa".format(job_id)            
            print(filename)
            print(best_locuses.values())
            extract_genes(filename, "fasta", query_file, "fasta",  set(best_locuses.values())) 

            blastp(query_file, original_file, out_file)
            reverse_matricies = find_best_match(out_file)
            
            if locus in best_locuses:
                other_locus = best_locuses[locus]
                if reverse_matricies[other_locus] == locus:
                    if not locus in orthologs:
                        orthologs[locus] = []
                    orthologs[locus].append(other_locus)

        for original_locus, orthologs in orthologs.items():
            request = [{"file":original_file, "locuses":[original_locus]}]
            for ortho_locus, other_file in zip(orthologs, test_against):
                request.append({"file": other_file, "locuses":[ortho_locus]})
            
            result_name = "{0}-memereq.json".format(len(files))
            files.append(result_name)
            with open(result_name, "w") as f:
                f.write(json.dumps(request))

    return files

def extract_upstream(indicies, genome, amount, overlap, min_length=8):
    """ For each index in [indicies], extract [amount] codons before it or
        up until the previous genes is overlap is false """

    records = []
    prev_end = -1
    index = 0
    for feature in filter(lambda f: f.type == "CDS", genome.features):
        if index in indicies:
            end = int(feature.location.start)
            start = max(end - amount, 0)
            if not overlap:
                start = max(start, prev_end)

            if (end - start) > min_length:
                upstream = genome[start:end]
                upstream.id =  "{0}|{1}".format(genome.id, feature.qualifiers["locus_tag"][0])
                records.append(upstream)

        index += 1
        prev_end = int(feature.location.end)

    return records


def extract_upstream_for_meme(genomes, locuses, upstream, radius, overlap):
    """ Given genomes and locuses, returns all [upstream] number
        of codons before the gene and its [radius] number of
        surrounding genes. [overlap] determines whether or not to allow the upstream to
        bleed into the previous CDS."""

    records = []
    for genome in genomes:
        feature_len = len(genome.features)

        index = 0
        locations = set()
        for feature in filter(lambda f: f.type == "CDS", genome.features):
            locus = feature.qualifiers["locus_tag"][0] 
            if locus in locuses:
                locations.add(index)
                for i in range(index - radius, index + radius):
                    locations.add(i)
            
            index += 1

        print(locations)
        records += extract_upstream(locations, genome, upstream, overlap)

    return records

    
def find_tfbs(job_id, filename, outdir, upstream, radius, overlap):
    """ Finds the the best motif using using the query specified in [filename]
        and outputs the files into the directory specified in [outdir]"""

    records = []
    with open(filename) as f:
        jobs = []
        for req in json.loads(f.read()):
            print(req["file"])
            genomes = list(SeqIO.parse(req["file"], "genbank"))
            records += extract_upstream_for_meme(genomes, req["locuses"], upstream, radius, overlap)

    temp_file = "{0}-meme.fsa".format(job_id)
    SeqIO.write(records, temp_file, "fasta")
    meme([temp_file], outdir)
        
        
#find_transcription_factors("1", "511145", "results", "domain-matching/db/Sigma.pin", 1)
#find_tfbs(1, "meme/request.json", "output", 150, 0, False)
jobs = frozenset(["211586", "225849", "318167", "319224", "323850", "325240", "326297", "351745", "392500", "398579", "425104", "458817", "60480", "60481", "94122"])
matrix = dict()
for job in jobs:
    test_against = ["meme/{0}".format(other) for other in jobs - set([job])]
    matrix["query/" + job + ".fsa"] = ("meme/{0}.fsa".format(job), test_against)

find_orthologs(1, matrix)
